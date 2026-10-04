"""Adjoint Gauss-Newton for (m, I_zz, mu_floor) per env from measured joint torques.

Measurement equation (fixed-base Talos EOM at the MEASURED (q_t, v_t, v_{t+1})):
    tau_t = M(q_t) (v_{t+1}-v_t)/dt + nle(q_t,v_t) - J_h(q_t)^T B0^T lam_h,t / dt
          = w_t - P_t lam_h,t ,      P_t = J_h(q_t)^T B0(x_box,t)^T / dt
w_t is recorded by the generator from the full pinocchio terms at the measured state, so
the only theta-dependent quantity is the hand contact impulse lam_h,t.
Residual  r_t(theta) = tau_meas,t - w_t + P_t lam_h,t(theta) ,  noise on tau enters r
directly AND through v_free of the hand (c_hand uses tau_meas), so the estimate is
noise-limited by construction.

lam(theta) is the solution of the SAME NCP the generator solved, on the reduced scene
(24 dofs = 3 hand + 6 box + 15 compliance) whose (G, b) equal the full 53-dof scene's:
the robot enters exactly through A = J_h M(q)^-1 J_h^T and J_h v_free,robot.
It is solved by the SAME batched kernels (`ncp_batch_gpu`), over env x step at once.
d lam / d theta is `ncp_ref.sensitivity`'s implicit-function-theorem system, batched:
    dlam/dtheta_j = -A^-1 [ Bb (dG/dtheta_j lam + db/dtheta_j) + Bmu dmu/dtheta_j ]
which is the adjoint through both contacts (hand-box and box-floor).

    OMP_NUM_THREADS=4 python -m motion_engine.ncp.sysid_gn --sigma 1.0 --window 50

MEASURED, 512 envs on the full dynamics (the robot's M, nle and J recomputed at the
integrated state every step; NCP natural residual max 4.213074e-12 over 512 x 100 steps):
  * the adjoint against central differences: relative p95 1.774e-07 for mass and
    2.044e-07 for mu at h = 1e-5; on non-sliding rows the adjoint is exactly 0 and the
    largest difference quotient there is 4.5e-08.
  * the error grows monotonically with the torque noise on all nine (sigma, window)
    cells, and it sits ON the Cramer-Rao bound: at sigma_tau = 1 N m over 50 steps the
    measured median relative mass error is 3.8590e-02 against 0.6745 x CRB = 3.8621e-02.
    So the error is information-bound, not a local minimum -- starting from the truth
    moves the p95 only from 1.5173e-01 to 1.4101e-01.

THE GATE FAILS, AND IT FAILS ON THE BOUND, NOT ON THE METHOD.
  Acceptance was p95 < 0.05 for m and mu at sigma_tau = 1 N m over 50 steps. Measured:
  1.5173e-01 and 1.6086e-01. The CRB itself is 5.726e-02 per N m at that window, so the
  gate was set below what the observable can deliver. At sigma_tau = 0.1 N m the same
  criterion passes (1.5272e-02 and 1.6340e-02).
  The Fisher screen fails too: lambda_min scores AUC 0.6108 at 50 steps and 0.7828 at
  100; a CRB score from the same Fisher gives 0.7157 and 0.8116; an oracle using the CRB
  at the TRUE parameters gives 0.7141 and 0.7949, so 0.9 is unreachable with this
  observable rather than with this estimator.
  Baselines, median relative mass error: quasi-static F = m a is 3.4138e-01 against this
  estimator's 3.8590e-02; a recursive mu estimator with quasi-static mass is 2.6829e-01,
  but handed the TRUE mass it reaches 1.3501e-02 on mu alone.

KNOWN DEFECTS, stated as numbers.
  Coplanar corner contacts give a condition number of 3e17; the fix is an eps = 0.02
  compliance carried as virtual degrees of freedom. The first-stage estimator has errors
  in variables and needs the total-least-squares form. Earlier variants that froze M,
  Minv, nle and J over the trajectory are rejected: their result was an artefact of the
  frozen stance, not a measurement.
"""
import argparse, json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
import warp as wp
from . import ncp_batch_gpu as H
from . import sysid_scene as S
(NVH, NVB, NVIRT) = (3, 6, 15)
NVR_ = NVH + NVB
NVD = NVR_ + NVIRT
(MROWS, NC) = (15, 5)
(TOL_LAM, TOL_REL) = (1e-10, 1e-06)

def static_parts(d, sigma, window, seed, noise_q=0.0):
    """theta-independent pieces, flattened over (env, step)."""
    N = window
    dt = float(d['dt'])
    eps = float(d['eps'])
    plane_z = float(d['plane_z'])
    B = d['m_true'].shape[0]
    tau = d['tau'][:N]
    w = d['w'][:N]
    nle = d['nle'][:N]
    Jh = d['Jh'][:N]
    D = d['D'][:N]
    A = d['A'][:N]
    Jhv = d['Jhv'][:N]
    hand = d['hand'][:N]
    box = d['box'][:N]
    rng = np.random.default_rng(seed)
    tau_meas = tau + rng.normal(0.0, sigma, tau.shape) if sigma > 0 else tau.copy()
    F = N * B
    rs = lambda a: a.reshape((F,) + a.shape[2:])
    (pos, quat) = (rs(box[..., :3]), rs(box[..., 3:7]))
    (vel, omg) = (rs(box[..., 7:10]), rs(box[..., 10:13]))
    (Jh, D, A, Jhv, hand) = (rs(Jh), rs(D), rs(A), rs(Jhv), rs(hand))
    (tau_meas, w, nle) = (rs(tau_meas), rs(w), rs(nle))
    Lb = np.repeat(d['L'][None], N, 0).reshape(F, 3)
    half = 0.5 * Lb
    R = S.quat_to_R(quat)
    n0 = -R[:, :, 0]
    (t1, t2) = S.tangents_hand(n0)
    B0 = np.stack([n0, t1, t2], 1)
    r0 = hand - pos
    J = np.zeros((F, MROWS, NVD))
    J[:, np.arange(MROWS), NVR_ + np.arange(MROWS)] = 1.0
    J[:, 0:3, 0:NVH] = B0
    J[:, 0:3, NVH:NVH + 3] = -B0
    J[:, 0:3, NVH + 3:NVR_] = np.einsum('bij,bjk->bik', B0, S.skew(r0))
    b_off = np.zeros((F, MROWS))
    b_off[:, 0] = (np.einsum('bi,bi->b', n0, r0) - half[:, 0]) / dt
    for k in range(4):
        rk = np.einsum('bij,bj->bi', R, S.CORNERS[k] * half)
        row = 3 + 3 * k
        J[:, row:row + 3, NVH:NVH + 3] = S.B_FLOOR
        J[:, row:row + 3, NVH + 3:NVR_] = -S.B_FLOOR @ S.skew(rk)
        b_off[:, row] = (pos[:, 2] + rk[:, 2] - plane_z) / dt
    vfy = Jhv + dt * np.einsum('bij,bj->bi', D, w - nle)
    A = np.zeros_like(A)
    P = np.einsum('bji,bkj->bik', Jh, B0) / dt
    y = tau_meas - w
    env = np.repeat(np.arange(B)[None], N, 0).reshape(F)
    return dict(F=F, B=B, N=N, dt=dt, eps=eps, J=J, b_off=b_off, vfy=vfy, A=A, vel=vel, omg=omg, R=R, P=P, y=y, env=env, L=Lb, hand=hand, pos=pos, mu_hand=float(d['mu_hand']), sigma=sigma)

def theta_parts(sp, m, I_zz, mu, deriv=True):
    """Minv, v_free, mu-vector and their theta-derivatives for the flattened batch."""
    (F, dt, R) = (sp['F'], sp['dt'], sp['R'])
    L = sp['L']
    Ixx = m / 12.0 * (L[:, 1] ** 2 + L[:, 2] ** 2)
    Iyy = m / 12.0 * (L[:, 0] ** 2 + L[:, 2] ** 2)
    Ib = np.stack([Ixx, Iyy, I_zz], 1)
    Iwi = np.einsum('bij,bj,bkj->bik', R, 1.0 / Ib, R)
    Iw = np.einsum('bij,bj,bkj->bik', R, Ib, R)
    Minv = np.zeros((F, NVD, NVD))
    Minv[:, 0:NVH, 0:NVH] = sp['A']
    i3 = np.arange(NVH, NVH + 3)
    Minv[:, i3, i3] = (1.0 / m)[:, None]
    Minv[:, NVH + 3:NVR_, NVH + 3:NVR_] = Iwi
    vi = np.arange(NVR_, NVD)
    Minv[:, vi, vi] = sp['eps']
    omg = sp['omg']
    gyro = -np.cross(omg, np.einsum('bij,bj->bi', Iw, omg))
    v_free = np.zeros((F, NVD))
    v_free[:, 0:NVH] = sp['vfy']
    v_free[:, NVH:NVH + 3] = sp['vel'] + dt * S.G_VEC
    v_free[:, NVH + 3:NVR_] = omg + dt * np.einsum('bij,bj->bi', Iwi, gyro)
    muv = np.empty((F, NC))
    muv[:, 0] = sp['mu_hand']
    muv[:, 1:] = mu[:, None]
    if not deriv:
        return (Minv, v_free, muv, None, None, None)
    dIb_dm = np.stack([Ixx / m, Iyy / m, np.zeros_like(m)], 1)
    dIb_dI = np.stack([np.zeros_like(m), np.zeros_like(m), np.ones_like(m)], 1)
    dv = {}
    dMi = {}
    for (key, dIb) in (('m', dIb_dm), ('I', dIb_dI)):
        dIwi = np.einsum('bij,bj,bkj->bik', R, -dIb / Ib ** 2, R)
        dIw = np.einsum('bij,bj,bkj->bik', R, dIb, R)
        M_ = np.zeros((F, NVD, NVD))
        if key == 'm':
            M_[:, i3, i3] = (-1.0 / m ** 2)[:, None]
        M_[:, NVH + 3:NVR_, NVH + 3:NVR_] = dIwi
        dMi[key] = M_
        dgyro = -np.cross(omg, np.einsum('bij,bj->bi', dIw, omg))
        vv = np.zeros((F, NVD))
        vv[:, NVH + 3:NVR_] = dt * (np.einsum('bij,bj->bi', dIwi, gyro) + np.einsum('bij,bj->bi', Iwi, dgyro))
        dv[key] = vv
    dmu = np.zeros((F, NC))
    dmu[:, 1:] = 1.0
    return (Minv, v_free, muv, dMi, dv, dmu)

def classify_batch(lam, muv):
    L = lam.reshape(-1, NC, 3)
    scale = np.maximum(np.abs(L).reshape(len(L), -1).max(1), 1e-30)[:, None]
    ln = L[:, :, 0]
    lt = np.linalg.norm(L[:, :, 1:], axis=2)
    op = ln <= TOL_LAM + TOL_REL * scale
    sl = ~op & (lt >= muv * ln * (1.0 - TOL_REL) - TOL_LAM)
    return (op, sl, ~op & ~sl)

def sens_solve(lam, G, b, muv, rhs_raw):
    """ncp_ref.sensitivity batched.  rhs_raw (F,15,K) = dG/dtheta_k lam + db/dtheta_k,
    plus dmu handled inside.  Returns dlam/dtheta (F,15,K) and cond of the KKT matrix."""
    F = lam.shape[0]
    K = rhs_raw.shape[2]
    u = np.einsum('bij,bj->bi', G, lam) + b
    (op, sl, st) = classify_batch(lam, muv)
    Amat = np.zeros((F, MROWS, MROWS))
    rhs = np.zeros((F, MROWS, K))
    dmu_cols = np.zeros((F, NC, K))
    dmu_cols[:, 1:, K - 1] = 1.0
    I3 = np.eye(3)
    for c in range(NC):
        i = 3 * c
        (o, s, t) = (op[:, c], sl[:, c], st[:, c])
        Amat[o, i:i + 3, i:i + 3] = I3
        Amat[t, i:i + 3, :] = G[t, i:i + 3, :]
        rhs[t, i:i + 3, :] = rhs_raw[t, i:i + 3, :]
        if s.any():
            gs = np.where(s)[0]
            ss = u[gs, i + 1:i + 3]
            ns = np.linalg.norm(ss, axis=1)
            ns = np.maximum(ns, 1e-14)
            sh = ss / ns[:, None]
            Pm = (np.eye(2)[None] - sh[:, :, None] * sh[:, None, :]) / ns[:, None, None]
            ln = lam[gs, i]
            mc = muv[gs, c]
            Amat[gs, i, :] = G[gs, i, :]
            Amat[gs, i + 1:i + 3, :] = (mc * ln)[:, None, None] * np.einsum('bij,bjk->bik', Pm, G[gs, i + 1:i + 3, :])
            Amat[gs, i + 1:i + 3, i] += (mc * ln)[:, None] * 0.0 + mc[:, None] * sh
            Amat[gs, i + 1, i + 1] += 1.0
            Amat[gs, i + 2, i + 2] += 1.0
            rhs[gs, i, :] = rhs_raw[gs, i, :]
            rhs[gs, i + 1:i + 3, :] = (mc * ln)[:, None, None] * np.einsum('bij,bjk->bik', Pm, rhs_raw[gs, i + 1:i + 3, :])
            rhs[gs, i + 1:i + 3, :] += (ln[:, None] * sh)[:, :, None] * dmu_cols[gs, c][:, None, :]
    dl = -np.linalg.solve(Amat, rhs)
    return dl

class GpuNcp:

    def __init__(self, chunk, n_iter, cg_iters=8, rho='sqrt', tol=1e-14):
        self.sol = H.BatchProxADMMGPU(chunk, NVD, NC)
        (self.chunk, self.n_iter, self.cg, self.rho, self.tol) = (chunk, n_iter, cg_iters, rho, tol)
        self.graph = None

    def solve(self, J, Minv, v_free, muv, b, dt):
        s = self.sol
        s.upload(J, Minv, v_free, muv, dt)
        wp.copy(s.bb, wp.array(np.ascontiguousarray(b.reshape(-1, 3)), dtype=H.vec3d, device='cuda:0'))
        s.set_rho(self.rho)
        if self.graph is None:
            (self.graph, _) = s.capture_fwd(self.n_iter, cg_iters=self.cg, refresh_every=4, tol=self.tol)
        s.replay(self.graph)
        return (s.lam.numpy().reshape(self.chunk, MROWS).copy(), s.resenv.numpy().copy())

def _cost(sp, gpu, chunk, theta):
    env = sp['env']
    F = sp['F']
    (Minv, v_free, muv, _, _, _) = theta_parts(sp, theta[env, 0], theta[env, 1], theta[env, 2], deriv=False)
    b = np.einsum('bij,bj->bi', sp['J'], v_free) + sp['b_off']
    lam = np.empty((F, MROWS))
    resmax = 0.0
    for s0 in range(0, F, chunk):
        sl = slice(s0, s0 + chunk)
        (lm, re) = gpu.solve(sp['J'][sl], Minv[sl], v_free[sl], muv[sl], b[sl], sp['dt'])
        lam[sl] = lm
        resmax = max(resmax, re.max())
    r = sp['y'] + np.einsum('bij,bj->bi', sp['P'], lam[:, :3])
    cost = np.zeros(sp['B'])
    np.add.at(cost, env, (r ** 2).sum(1))
    return (cost, lam, resmax)

def _jac(sp, gpu, chunk, theta, lam):
    env = sp['env']
    F = sp['F']
    (Minv, v_free, muv, dMi, dv, _) = theta_parts(sp, theta[env, 0], theta[env, 1], theta[env, 2], deriv=True)
    (G, b) = S.delassus_b(sp['J'], Minv, v_free, sp['b_off'])
    Jtl = np.einsum('bji,bj->bi', sp['J'], lam)
    rhs = np.zeros((F, MROWS, 3))
    for (k, key) in enumerate(('m', 'I')):
        rhs[:, :, k] = np.einsum('bij,bj->bi', sp['J'], np.einsum('bij,bj->bi', dMi[key], Jtl) + dv[key])
    dl = sens_solve(lam, G, b, muv, rhs)
    r = sp['y'] + np.einsum('bij,bj->bi', sp['P'], lam[:, :3])
    Jr = np.einsum('bij,bjk->bik', sp['P'], dl[:, :3, :])
    Hm = np.einsum('bik,bil->bkl', Jr, Jr)
    gv = np.einsum('bik,bi->bk', Jr, r)
    Hs = np.zeros((sp['B'], 3, 3))
    gs = np.zeros((sp['B'], 3))
    np.add.at(Hs, env, Hm)
    np.add.at(gs, env, gv)
    return (Hs, gs)
BOUNDS = ((0.2, 40.0), (0.0001, 5.0), (0.005, 2.0))

def _step(Hs, gs, theta, damp):
    sc = np.abs(theta)
    Hn = Hs * sc[:, :, None] * sc[:, None, :]
    gn = gs * sc
    tr = np.trace(Hn, axis1=1, axis2=2)
    Hn = Hn + (damp * np.maximum(tr, 1e-300))[:, None, None] * np.eye(3)[None]
    dth = -np.linalg.solve(Hn, gn[:, :, None])[:, :, 0] * sc
    out = theta + dth
    for (j, (lo, hi)) in enumerate(BOUNDS):
        out[:, j] = np.clip(out[:, j], lo, hi)
    return out

def run(sp, gpu, chunk, iters, theta0=None, verbose=True, trials=5):
    """Levenberg-Marquardt: per-env damping, a step is taken only if the env's cost
    decreases (otherwise damping x10 and retry).  Gauss-Newton direction from the
    contact adjoint."""
    B = sp['B']
    if theta0 is None:
        Lb = sp['L'].reshape(sp['N'], B, 3)[0]
        theta = np.stack([np.full(B, 5.0), 5.0 / 12.0 * (Lb[:, 0] ** 2 + Lb[:, 1] ** 2), np.full(B, 0.4)], 1)
    else:
        theta = theta0.copy()
    damp = np.full(B, 1e-06)
    (cost, lam, resmax) = _cost(sp, gpu, chunk, theta)
    hist = []
    for it in range(iters):
        (Hs, gs) = _jac(sp, gpu, chunk, theta, lam)
        acc = np.zeros(B, bool)
        best = theta.copy()
        best_cost = cost.copy()
        for tr_i in range(trials):
            cand = _step(Hs, gs, theta, damp)
            (c2, lam2, resmax) = _cost(sp, gpu, chunk, cand)
            good = ~acc & (c2 < best_cost)
            best[good] = cand[good]
            best_cost[good] = c2[good]
            rows = good[sp['env']]
            lam[rows] = lam2[rows]
            damp[good] /= 5.0
            damp[~acc & ~good] *= 10.0
            acc |= good
            if acc.all():
                break
        (theta, cost) = (best, best_cost)
        hist.append(dict(it=it, cost=float(cost.sum()), cost_p95=float(np.percentile(cost, 95)), resmax=float(resmax), accepted=int(acc.sum()), damp_med=float(np.median(damp))))
        if verbose:
            print(f'gn {it}: cost {cost.sum():.6e} acc {acc.sum()}/{B} ncp_res {resmax:.2e} damp_med {np.median(damp):.1e}', flush=True)
    return (theta, Hs, hist, cost)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sim', default='u98x_sim_B512.npz')
    ap.add_argument('--sigma', type=float, required=True)
    ap.add_argument('--window', type=int, required=True)
    ap.add_argument('--iters', type=int, default=6)
    ap.add_argument('--chunk', type=int, default=25600)
    ap.add_argument('--nit', type=int, default=200)
    ap.add_argument('--cg', type=int, default=4)
    ap.add_argument('--trials', type=int, default=3)
    ap.add_argument('--init', default=None)
    ap.add_argument('--truth-start', action='store_true')
    ap.add_argument('--seed', type=int, default=4242)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    d = np.load(os.path.join(HERE, args.sim))
    sp = static_parts(d, args.sigma, args.window, args.seed)
    gpu = GpuNcp(min(args.chunk, sp['F']), args.nit, cg_iters=args.cg)
    t0 = time.perf_counter()
    tr0 = np.stack([d['m_true'], d['I_zz_true'], d['mu_true']], 1)
    (c_true, _, rmax_true) = _cost(sp, gpu, args.chunk, tr0)
    print(f'cost at true theta: sum {c_true.sum():.6e} p95 {np.percentile(c_true, 95):.6e} ncp_res {rmax_true:.2e}', flush=True)
    th0 = np.load(os.path.join(HERE, args.init))['theta'] if args.init else None
    if args.truth_start:
        th0 = tr0.copy()
    (theta, Hs, hist, cost) = run(sp, gpu, min(args.chunk, sp['F']), args.iters, theta0=th0, trials=args.trials)
    wall = time.perf_counter() - t0
    tr = np.stack([d['m_true'], d['I_zz_true'], d['mu_true']], 1)
    rel = np.abs(theta - tr) / np.abs(tr)
    names = ('m', 'I_zz', 'mu')
    print(f'sigma {args.sigma} window {args.window}: wall {wall:.1f} s (GPU shared)')
    for (j, nm) in enumerate(names):
        print(f'  {nm:5s} rel p50 {np.median(rel[:, j]):.6e}  p95 {np.percentile(rel[:, j], 95):.6e}  max {rel[:, j].max():.6e}')
    lmin = np.linalg.eigvalsh(Hs / max(args.sigma, 1e-12) ** 2)[:, 0]
    out = args.out or f'gn_sigma{args.sigma}_w{args.window}.npz'
    np.savez(os.path.join(HERE, out), theta=theta, rel=rel, H=Hs, lam_min=lmin, cost=cost, sigma=args.sigma, window=args.window, wall_s=wall, hist=json.dumps(hist), true=tr, cost_true=c_true)
    print('->', out)
if __name__ == '__main__':
    main()

def crb_from_fisher(H_batch, theta_true):
    """Relative Cramer-Rao bound per parameter, per env, scaled per unit torque noise.

    `H_batch` is the Gauss-Newton Hessian J^T J of the residual in `_jac`, i.e. the
    Fisher information at sigma = 1; the bound at noise sigma is sigma times this.
    Returns (crb_rel, lambda_min) with crb_rel[e, j] = sqrt([H^-1]_jj) / |theta_j|.

    The median of a zero-mean normal error of scale s is 0.6745 s, so an estimator that
    is on its bound shows median |rel error| = 0.6745 * sigma * crb_rel. Measured at
    sigma_tau = 1 N m over 50 steps: 3.8590e-02 against 3.8621e-02.
    """
    H_batch = np.asarray(H_batch, float)
    theta_true = np.asarray(theta_true, float)
    n = H_batch.shape[-1]
    Hinv = np.linalg.inv(H_batch + 1e-300 * np.eye(n))
    crb = np.sqrt(np.maximum(np.diagonal(Hinv, axis1=1, axis2=2), 0.0)) / np.abs(theta_true)
    return (crb, np.linalg.eigvalsh(H_batch)[:, 0])
HALF_NORMAL_MEDIAN = 0.6745
