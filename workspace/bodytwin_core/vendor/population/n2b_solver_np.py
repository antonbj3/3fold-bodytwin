"""N2b copy of results/N2/n2_np.py (solver part only: _amap, _J, _mv, _mtv, solve_poly, dual_lower_bound, solve_minmax),
sha256 of source file at copy time in n2b_sources.json. Unchanged algorithm and constants."""
import numpy as np
CAP = 1e-12

def _amap(s, p, U):
    sp = np.maximum(s, 0.0)
    a0 = sp / 2.0 if p == 2 else (sp / p) ** (1.0 / (p - 1.0))
    a = np.minimum(a0, U[:, None])
    inner = (a0 > 0) & (a0 < U[:, None])
    if p == 2:
        w = np.where(inner, 0.5, 0.0)
    else:
        # curvature da/ds = a/((p-1)s) is unbounded as s -> 0+; evaluate it at s >= CAP max(s) (bounded, still
        # a descent direction; exact Hessian away from the kink)
        sfl = np.maximum(s, CAP * sp.max(-1, keepdims=True))
        w = np.where(inner, (np.minimum(sfl, p * U[:, None] ** (p - 1)) / p) ** (1.0 / (p - 1.0)) / ((p - 1.0) * sfl), 0.0)
    return a, w


def _J(lam, s, a, b, p):
    return -(lam * b).sum(-1) - (a ** p - s * a).sum(-1)


def _mv(A, x):      # (B,m,n) @ (B,n)
    return (A @ x[..., None])[..., 0]   # N2b: matmul


def _mtv(A, y):     # (B,m,n)^T @ (B,m)
    return (y[:, None, :] @ A)[:, 0, :]   # N2b: matmul


def solve_poly(Ap, b, p, U=None, tol=1e-14, maxit=200, lam0=None, stop_infeasible=False):
    B, m, n = Ap.shape
    U = np.full(B, np.inf) if U is None else U
    sc = np.maximum(np.linalg.norm(b, axis=-1), 1e-300)
    bh = b / sc[:, None]; Ah = Ap / sc[:, None, None]
    G = Ah @ Ah.transpose(0, 2, 1)   # N2b: matmul (BLAS) instead of einsum, same arithmetic
    reg0 = 1e-12 * np.trace(G, axis1=1, axis2=2) / m
    eye = np.eye(m)
    if lam0 is None:
        lam = 2.0 * np.linalg.solve(G + reg0[:, None, None] * eye, bh[..., None])[..., 0]   # unconstrained p=2 optimum
        if p != 2:   # rescale so the start has the activity scale of the p=2 optimum: s = p a^(p-1)
            s2 = _mtv(Ah, lam); atyp = np.maximum(s2, 0).sum(-1) / 2.0 / np.maximum((s2 > 0).sum(-1), 1)
            lam = lam * (p * np.maximum(atyp, 1e-12) ** (p - 1.0) / (2.0 * np.maximum(atyp, 1e-12)))[:, None]
    else:
        lam = lam0.copy()
    mu = np.zeros(B)                       # Levenberg damping, raised when the line search fails
    s = _mtv(Ah, lam)
    a, w = _amap(s, p, U)
    grad = _mv(Ah, a) - bh
    J = _J(lam, s, a, bh, p)
    active = np.linalg.norm(grad, axis=-1) > tol
    iters = np.zeros(B, dtype=np.int64)
    for it in range(maxit):
        if not active.any():
            break
        idx = np.nonzero(active)[0]
        Ai, bi, li, gi, wi, Ji, Ui = Ah[idx], bh[idx], lam[idx], grad[idx], w[idx], J[idx], U[idx]
        H = (Ai * wi[:, None, :]) @ Ai.transpose(0, 2, 1)   # N2b: matmul
        trH = np.trace(H, axis1=1, axis2=2)
        reg = reg0[idx] + 1e-12 * trH + mu[idx] * (trH / m + reg0[idx])
        d = -np.linalg.solve(H + reg[:, None, None] * eye, gi[..., None])[..., 0]
        slope = (gi * d).sum(-1)
        t = np.ones(len(idx)); done = np.zeros(len(idx), bool)
        ln, Jn = li.copy(), Ji.copy()
        sn = np.empty_like(s[idx]); an = np.empty_like(sn); wn = np.empty_like(sn)
        gnorm = np.linalg.norm(gi, axis=-1)
        for _ in range(40):
            lt = li + t[:, None] * d
            st = _mtv(Ai, lt)
            at, wt = _amap(st, p, Ui)
            Jt = _J(lt, st, at, bi, p)
            gt = np.linalg.norm(_mv(Ai, at) - bi, axis=-1)
            ok = (~done) & ((Jt <= Ji + 1e-4 * t * slope) | ((gnorm < 1e-6) & (gt < (1 - 1e-4 * t) * gnorm)))
            ln[ok], sn[ok], an[ok], wn[ok], Jn[ok] = lt[ok], st[ok], at[ok], wt[ok], Jt[ok]
            done |= ok
            if done.all():
                break
            t = np.where(done, t, t * 0.5)
        if not done.all():
            nd = ~done
            sn[nd] = _mtv(Ai[nd], li[nd])
            an[nd], wn[nd] = _amap(sn[nd], p, Ui[nd])
        gn = _mv(Ai, an) - bi
        lam[idx], s[idx], a[idx], w[idx], J[idx], grad[idx] = ln, sn, an, wn, Jn, gn
        iters[idx] += 1
        mu[idx] = np.where(done, mu[idx] * 0.1, np.maximum(mu[idx] * 10.0, 1e-6))
        still = (np.linalg.norm(gn, axis=-1) > tol) & (done | (mu[idx] < 1e8))
        if stop_infeasible:
            cert = (ln * bi).sum(-1) > Ui * np.maximum(sn, 0).sum(-1) * (1 + 1e-12)
            still = still & ~cert
        active[idx] = still
    resid = np.linalg.norm(grad, axis=-1)
    return dict(a=a, lam=lam, s=s, resid=resid, iters=iters, conv=resid <= 1e-10, sc=sc)


def dual_lower_bound(Ap, b, lam_n):
    sc = np.maximum(np.linalg.norm(b, axis=-1), 1e-300)
    s = _mtv(Ap / sc[:, None, None], lam_n)
    den = np.maximum(s, 0).sum(-1); num = (lam_n * b / sc[:, None]).sum(-1)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(den > 0, num / np.where(den > 0, den, 1), np.where(num > 0, np.inf, 0.0))


def solve_minmax(Ap, b, rel_gap=1e-9, max_bisect=60, qp_maxit=80):
    base = solve_poly(Ap, b, 2)
    cone_ok = base['conv']
    hi = np.where(cone_ok, base['a'].max(-1), np.inf)
    lo = np.maximum(dual_lower_bound(Ap, b, base['lam']), 0.0)
    lo = np.where(cone_ok, np.minimum(lo, hi), lo)
    lam_ws = base['lam'].copy()
    for k in range(max_bisect):
        act = cone_ok & ((hi - lo) > rel_gap * hi)
        if not act.any():
            break
        idx = np.nonzero(act)[0]
        mid = 0.5 * (lo[idx] + hi[idx])
        r = solve_poly(Ap[idx], b[idx], 2, U=mid, maxit=qp_maxit, lam0=lam_ws[idx], stop_infeasible=True)
        lb = dual_lower_bound(Ap[idx], b[idx], r['lam'])
        feas = r['conv']
        hi[idx] = np.where(feas, np.minimum(hi[idx], np.minimum(r['a'].max(-1), mid)), hi[idx])
        lo[idx] = np.where(feas, np.maximum(lo[idx], np.minimum(lb, mid)), np.maximum(lo[idx], mid))
        lam_ws[idx] = np.where(feas[:, None], r['lam'], lam_ws[idx])
    Ut = np.where(cone_ok, hi * (1 + 1e-6), 1.0)
    r = solve_poly(Ap, b, 2, U=Ut, lam0=lam_ws)
    return dict(a=r['a'], lam=r['lam'], s=r['s'], resid=r['resid'], conv=r['conv'] & cone_ok, U_ub=hi, U_lb=lo,
                cone_ok=cone_ok, U=Ut)


# ---------------- problem ----------------
class Data:
    def __init__(self):
        z = np.load(_find('k1_hip.npz')); nd = np.load(_find('n2_nodes.npz'))
        self.names = list(z['names']); self.t = z['t']
        self.gW, self.rW, self.kW = z['g'], z['r'], z['k']
        self.FJ, self.RJ, self.N = z['Ft'], z['R'], z['N'][0]
        self.rect = np.array([m.startswith('RectusFemoris') for m in self.names]); self.free = ~self.rect
        self.c, self.kc, self.kax, self.isb = z['c'], z['kc'], z['kax'], z['isb']
        self.P0, self.P1 = nd['P'][:, :, 0], nd['P'][:, :, 1]
        self.fsd = nd['first_seg_end_distal']; self.fsk = nd['knee'][:, 1]
        self.gidx = nd['gidx']; self.groups = list(nd['groups'])
        self.bfull = np.concatenate([np.einsum('tnj,tn->tj', self.rW, self.FJ), (self.kW * self.FJ).sum(-1, keepdims=True)], -1)
        self.GJ = np.einsum('tnj,tn->tj', self.gW, self.FJ)


