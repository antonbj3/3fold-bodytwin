"""CPU float64 reference: exact-cone NCP contact step with analytical sensitivities.

CPU side of the exact-cone contact pair; the GPU side is `ncp_gpu.py`.

Formulation (velocity level, one rigid scene, restitution 0), per contact c with
components ordered (n, t1, t2), rows contact-major:

    u = G lam + b,      G = J Minv J^T (Delassus),   b = J v_free (+ offset)

    NCP (exact Coulomb cone K_mu = {lam : ||lam_t|| <= mu lam_n}):
        K_mu  ∋  lam_c   ⊥   (u_c + Gamma(u_c))  ∈  K_mu*        (de Saxce)
        Gamma(u, mu) = (mu ||u_t||, 0, 0)                        [2304.06372 eq. 14]
        K_mu* = K_{1/mu}                                         [2304.06372 §III-A]

    CCP / Anitescu convex relaxation: same with Gamma := 0       [2304.06372 eq. 16]
    Pyramid: K replaced by the 4-facet pyramid.

Exact local (3x3) case analysis used by the Gauss-Seidel sweeps, derived from the
complementarity condition (derivation recorded in tests/test_ncp_reference.py):
    take-off :  lam_c = 0                     admissible iff w_n >= 0    (NCP)
                                              admissible iff w ∈ K*      (CCP)
    stick    :  lam_c = -G_cc^{-1} w          admissible iff lam_c ∈ K   (u_c = 0)
    slip     :  lam_c = lam_n (1, -mu d),  d unit 2-vector,  u_t ∥ +d,
                NCP:  u_n = 0                 (no hover)
                CCP:  u_n = mu (u_t·d)        (hover  =  the gliding artifact)
    The slip branch reduces to a scalar root find in the slip angle theta, because
    lam_n is linear in w once d is fixed.  Solved by grid-scan + bisection.

References: Le Lidec & Carpentier et al., arXiv:2304.06372 (Alg. 1/2/3/5/6, eq. 14-19);
            "Simple", arXiv:2405.17020 / RSS 2024 (ADMM + proximal on the Delassus).

Self-contained: numpy only, except `solve_pyramid_qp(..., method="scipy")` which
uses scipy.optimize.minimize (SLSQP) as an independent cross-check. ADMM preserves
the legacy Cholesky/LU arithmetic by default. Explicit xstep="triangular" factors
G + rho I with scipy.linalg.cholesky once per rho, keeps the factor in Fortran order
and solves each step by two BLAS triangular substitutions (dtrsv), with the legacy
numpy path as fallback when scipy is missing. This opt-in path is numerically close
but is not bit-identical to the legacy iteration.

Audited numbers this module reproduces (tests/test_ncp_reference.py):
  * complementarity of the ADMM solution on 6 scenes: eps_p = 0, eps_d <= 4.9e-13.
  * sum lam_n = M g dt to 3.8e-11; pushed-cube a_x against (F - mu m g)/m to 1.4e-15.
  * gliding gap: NCP 1.1e-14 m; Anitescu convex relaxation 5.0e-4 m mean / 4.3e-3 m peak,
    ratio to the analytical dt mu ||u_t|| 1.044; the relaxation takes 7.9x more tangential
    momentum at touchdown.
  * dlam/dmu and dlam/db by the implicit function theorem against central FD: 3.1e-7
    (the FD floor); two independent derivations agree to 3e-16.
  * ADMM (Simple-style) converges on 6/6 scenes in 64-299 iterations; Cholesky
    is rebuilt whenever rho changes (up to 13 factors on those scenes).
  * PGS fails on 4/6: stack-1000 stalls at 1.5e-3, the other three are only slow. Those
    four rows are xfail in the test file with the measured residuals, not removed.

Audit reservations carried as numbers, not adjectives:
  * with four coplanar closed box contacts, the four normal branch rows have
    rank at most three, including all-slip states. The branch Jacobian is singular.
    `sensitivity` withholds ambiguous per-contact derivatives by default; an
    explicit minimum_norm selection returns only compatible linear responses. A readout C is
    unique on a compatible linearized fiber only if ker(dF/dlam) is in ker(C).
    Net body impulse is not automatically unique for the full Coulomb NCP.
  * `scene_tripod`'s docstring claim "determinate" is wrong: rank 6/9.
  * boxed-LCP-PGS (the ODE/Bullet fixed point) is not the pyramid-QP optimum: gap 4-51 %.
"""
from __future__ import annotations
import math
from dataclasses import dataclass, field
import numpy as np
from .contact_rank import rank_factor, guarded_response
__all__ = ['Scene', 'proj_cone', 'proj_dual_cone', 'desaxce', 'natural_residual', 'complementarity', 'classify', 'solve_ncp_pgs', 'solve_ncp_admm', 'solve_convex_pgs', 'solve_pyramid_qp', 'sensitivity', 'sensitivity_fd', 'contact_diagnostics', 'ContactResult', 'Box', 'build_scene', 'simulate', 'scene_pushed_cube', 'scene_sliding_box', 'scene_massratio_stack', 'scene_tower']

@dataclass
class Scene:
    """J (3 n_c x n_v), Minv (n_v x n_v or diagonal (n_v,)), v_free (n_v),
    mu (n_c), body_pairs (n_c x 2 int, -1 = static), dt."""
    J: np.ndarray
    Minv: np.ndarray
    v_free: np.ndarray
    mu: np.ndarray
    body_pairs: np.ndarray
    dt: float
    b_offset: np.ndarray | None = None

    @property
    def n_c(self) -> int:
        return self.J.shape[0] // 3

    @property
    def n_v(self) -> int:
        return self.J.shape[1]

    def minv_mat(self) -> np.ndarray:
        M = np.asarray(self.Minv, dtype=np.float64)
        return np.diag(M) if M.ndim == 1 else M

    def delassus(self) -> np.ndarray:
        """G = J Minv J^T, formed explicitly (allowed for n_c <= 4096)."""
        J = np.asarray(self.J, dtype=np.float64)
        M = np.asarray(self.Minv, dtype=np.float64)
        MJt = M[:, None] * J.T if M.ndim == 1 else M @ J.T
        return J @ MJt

    def b(self) -> np.ndarray:
        bb = np.asarray(self.J, dtype=np.float64) @ np.asarray(self.v_free, dtype=np.float64)
        if self.b_offset is not None:
            bb = bb + np.asarray(self.b_offset, dtype=np.float64)
        return bb

    def apply(self, lam: np.ndarray) -> np.ndarray:
        """v = v_free + Minv J^T lam."""
        J = np.asarray(self.J, dtype=np.float64)
        M = np.asarray(self.Minv, dtype=np.float64)
        imp = J.T @ np.asarray(lam, dtype=np.float64)
        dv = M * imp if M.ndim == 1 else M @ imp
        return np.asarray(self.v_free, dtype=np.float64) + dv

def proj_cone(x: np.ndarray, mu) -> np.ndarray:
    """Euclidean projection onto the second-order Coulomb cone K_mu."""
    x = np.asarray(x, dtype=np.float64)
    single = x.ndim == 1
    X = np.atleast_2d(x).copy()
    m = np.broadcast_to(np.atleast_1d(np.asarray(mu, dtype=np.float64)), (X.shape[0],))
    xn = X[:, 0]
    xt = X[:, 1:3]
    nt = np.linalg.norm(xt, axis=1)
    out = np.zeros_like(X)
    zero_mu = m <= 0.0
    if np.any(zero_mu):
        out[zero_mu, 0] = np.maximum(0.0, xn[zero_mu])
    pos = ~zero_mu
    if np.any(pos):
        inside = pos & (nt <= m * xn)
        polar = pos & (m * nt <= -xn)
        mid = pos & ~inside & ~polar
        out[inside] = X[inside]
        if np.any(mid):
            a = (m[mid] * nt[mid] + xn[mid]) / (1.0 + m[mid] ** 2)
            out[mid, 0] = a
            out[mid, 1:3] = (m[mid] * a / nt[mid])[:, None] * xt[mid]
    return out[0] if single else out

def proj_dual_cone(x: np.ndarray, mu) -> np.ndarray:
    """Projection onto the dual cone K_mu* = K_{1/mu}."""
    x = np.asarray(x, dtype=np.float64)
    single = x.ndim == 1
    X = np.atleast_2d(x).copy()
    m = np.broadcast_to(np.atleast_1d(np.asarray(mu, dtype=np.float64)), (X.shape[0],))
    out = np.empty_like(X)
    zero_mu = m <= 0.0
    if np.any(zero_mu):
        out[zero_mu] = X[zero_mu]
        out[zero_mu, 0] = np.maximum(0.0, X[zero_mu, 0])
    pos = ~zero_mu
    if np.any(pos):
        out[pos] = proj_cone(X[pos], 1.0 / m[pos])
    return out[0] if single else out

def desaxce(u: np.ndarray, mu) -> np.ndarray:
    """de Saxce correction Gamma(u, mu) = (mu ||u_t||, 0, 0)  [2304.06372 eq. 14]."""
    u = np.asarray(u, dtype=np.float64)
    single = u.ndim == 1
    U = np.atleast_2d(u)
    m = np.broadcast_to(np.atleast_1d(np.asarray(mu, dtype=np.float64)), (U.shape[0],))
    g = np.zeros_like(U)
    g[:, 0] = m * np.linalg.norm(U[:, 1:3], axis=1)
    return g[0] if single else g

def _rho_of(G: np.ndarray, n_c: int) -> np.ndarray:
    """Per-contact step rho_c = 1/||G_cc||_2 for the natural-map residual."""
    rho = np.empty(n_c)
    for c in range(n_c):
        Gcc = G[3 * c:3 * c + 3, 3 * c:3 * c + 3]
        s = np.linalg.norm(Gcc, ord=2)
        rho[c] = 1.0 / s if s > 0 else 1.0
    return rho

def natural_residual(lam, G, b, mu, rho=None, desax: bool=True) -> float:
    """||lam - proj_K(lam - rho (u + Gamma(u)))||_inf   (SPEC A.1)."""
    lam = np.asarray(lam, dtype=np.float64)
    n_c = len(mu)
    if rho is None:
        rho = _rho_of(G, n_c)
    u = (G @ lam + b).reshape(n_c, 3)
    L = lam.reshape(n_c, 3)
    s = desaxce(u, mu) if desax else 0.0
    P = proj_cone(L - np.asarray(rho).reshape(-1, 1) * (u + s), mu)
    return float(np.max(np.abs(L - P))) if n_c else 0.0

def complementarity(lam, G, b, mu, desax: bool=True):
    """Per-contact (eps_p, eps_d, eps_c) of 2304.06372 §III-A:
    eps_p = dist_K(lam), eps_d = dist_K*(u + Gamma(u)), eps_c = |<lam, u + Gamma(u)>|."""
    lam = np.asarray(lam, dtype=np.float64)
    n_c = len(mu)
    L = lam.reshape(n_c, 3)
    u = (G @ lam + b).reshape(n_c, 3)
    w = u + (desaxce(u, mu) if desax else 0.0)
    eps_p = np.linalg.norm(L - proj_cone(L, mu), axis=1)
    eps_d = np.linalg.norm(w - proj_dual_cone(w, mu), axis=1)
    eps_c = np.abs(np.sum(L * w, axis=1))
    return (eps_p, eps_d, eps_c)

def classify(lam, G, b, mu, tol_lam=1e-10, tol_rel=1e-06):
    """Active-set labels per contact: 'open' / 'stick' / 'slip'."""
    lam = np.asarray(lam, dtype=np.float64)
    n_c = len(mu)
    L = lam.reshape(n_c, 3)
    scale = max(np.max(np.abs(L)), 1e-30)
    labels = np.empty(n_c, dtype='<U5')
    for c in range(n_c):
        (ln, lt) = (L[c, 0], np.linalg.norm(L[c, 1:3]))
        if ln <= tol_lam + tol_rel * scale:
            labels[c] = 'open'
        elif lt >= mu[c] * ln * (1.0 - tol_rel) - tol_lam:
            labels[c] = 'slip'
        else:
            labels[c] = 'stick'
    return labels
_THETA_GRID = 64
_BISECT_IT = 60

def _local_residual(l3, A, w, mu, desax, rho):
    u = A @ l3 + w
    s = desaxce(u, mu) if desax else np.zeros(3)
    return float(np.max(np.abs(l3 - proj_cone(l3 - rho * (u + s), mu))))

def _slip_eval(A, w, mu, th, desax):
    """Vectorised slip branch at angles `th`.  Returns (cross, dot, lam_n, lam (N,3)).

    lam = lam_n (1, -mu d),  d = (cos th, sin th).  lam_n is fixed by the branch
    equation (NCP: u_n = 0;  CCP: u_n = mu (u_t . d)), which is LINEAR in lam_n
    once d is fixed.  The remaining scalar equation is cross(d, u_t) = 0 with
    u_t . d > 0 (slip direction aligned with the tangential contact velocity)."""
    th = np.atleast_1d(np.asarray(th, dtype=np.float64))
    d = np.stack([np.cos(th), np.sin(th)], axis=1)
    a_nn = A[0, 0]
    a_nt = A[0, 1:3]
    a_tn = A[1:3, 0]
    A_tt = A[1:3, 1:3]
    k_n = a_nn - mu * (d @ a_nt)
    if desax:
        (den, num) = (k_n, -w[0] * np.ones_like(k_n))
    else:
        k_t = d @ a_tn - mu * np.einsum('ij,jk,ik->i', d, A_tt, d)
        den = k_n - mu * k_t
        num = mu * (d @ w[1:3]) - w[0]
    with np.errstate(divide='ignore', invalid='ignore'):
        ln = np.where(np.abs(den) > 1e-300, num / den, np.nan)
    lam_t = -mu * ln[:, None] * d
    u_t = ln[:, None] * a_tn[None, :] + lam_t @ A_tt.T + w[1:3][None, :]
    cross = d[:, 0] * u_t[:, 1] - d[:, 1] * u_t[:, 0]
    dot = np.sum(d * u_t, axis=1)
    lam = np.concatenate([ln[:, None], lam_t], axis=1)
    return (cross, dot, ln, lam)

def _slip_scalar(P, mu, th, desax):
    """Scalar slip-branch evaluation (no numpy overhead; hot inner loop).

    P = (a_nn, a_nt0, a_nt1, a_tn0, a_tn1, A00, A01, A10, A11, w0, w1, w2).
    Returns (cross, dot, lam_n, lam_t0, lam_t1) or None."""
    (a_nn, a_nt0, a_nt1, a_tn0, a_tn1, A00, A01, A10, A11, w0, w1, w2) = P
    cd = math.cos(th)
    sd = math.sin(th)
    k_n = a_nn - mu * (a_nt0 * cd + a_nt1 * sd)
    if desax:
        den = k_n
        num = -w0
    else:
        dAd = A00 * cd * cd + (A01 + A10) * cd * sd + A11 * sd * sd
        k_t = a_tn0 * cd + a_tn1 * sd - mu * dAd
        den = k_n - mu * k_t
        num = mu * (w1 * cd + w2 * sd) - w0
    if abs(den) < 1e-300:
        return None
    ln = num / den
    lt0 = -mu * ln * cd
    lt1 = -mu * ln * sd
    ut0 = ln * a_tn0 + A00 * lt0 + A01 * lt1 + w1
    ut1 = ln * a_tn1 + A10 * lt0 + A11 * lt1 + w2
    return (cd * ut1 - sd * ut0, cd * ut0 + sd * ut1, ln, lt0, lt1)

def _local_solve(A, w, mu, desax, lam_init, rho, tol=1e-14):
    """Exact local NCP/CCP solve for one contact (others frozen).  Returns lam (3,)."""
    if desax:
        ok = w[0] >= 0.0
    else:
        wt = math.hypot(w[1], w[2])
        ok = mu * wt <= w[0] if mu > 0 else w[0] >= 0.0
    if ok:
        return np.zeros(3)
    lam_st = None
    try:
        lam_st = -np.linalg.solve(A, w)
        if lam_st[0] >= 0.0 and math.hypot(lam_st[1], lam_st[2]) <= mu * lam_st[0]:
            return lam_st
    except np.linalg.LinAlgError:
        lam_st = None
    cands = []
    if lam_st is not None and np.all(np.isfinite(lam_st)):
        cands.append(lam_st)
    if mu > 0.0:
        P = (A[0, 0], A[0, 1], A[0, 2], A[1, 0], A[2, 0], A[1, 1], A[1, 2], A[2, 1], A[2, 2], w[0], w[1], w[2])
        step = 2.0 * math.pi / _THETA_GRID
        prev = _slip_scalar(P, mu, 0.0, desax)
        for i in range(1, _THETA_GRID + 1):
            th = i * step
            cur = _slip_scalar(P, mu, th, desax)
            if prev is not None and cur is not None and (prev[0] * cur[0] <= 0.0):
                (lo, hi, flo) = (th - step, th, prev[0])
                for _ in range(_BISECT_IT):
                    mid = 0.5 * (lo + hi)
                    r = _slip_scalar(P, mu, mid, desax)
                    if r is None:
                        break
                    if flo * r[0] <= 0.0:
                        hi = mid
                    else:
                        (lo, flo) = (mid, r[0])
                r = _slip_scalar(P, mu, 0.5 * (lo + hi), desax)
                if r is not None and r[2] > 0.0 and (r[1] > 0.0):
                    cands.append(np.array([r[2], r[3], r[4]]))
            prev = cur
    (best, best_r) = (None, np.inf)
    for cd_ in cands:
        if not np.all(np.isfinite(cd_)):
            continue
        r = _local_residual(cd_, A, w, mu, desax, rho)
        if r < best_r:
            (best, best_r) = (cd_, r)
    if best is not None and best_r <= tol:
        return best
    l = np.array(lam_init, dtype=np.float64)
    for _ in range(200):
        u = A @ l + w
        s = desaxce(u, mu) if desax else np.zeros(3)
        l2 = proj_cone(l - rho * (u + s), mu)
        if np.max(np.abs(l2 - l)) < 1e-16:
            l = l2
            break
        l = l2
    for cd_ in (l, np.array(lam_init, dtype=np.float64)):
        if not np.all(np.isfinite(cd_)):
            continue
        r = _local_residual(cd_, A, w, mu, desax, rho)
        if r < best_r:
            (best, best_r) = (cd_, r)
    return best if best is not None else np.zeros(3)

def _pgs(G, b, mu, iters, tol, warm, desax):
    G = np.asarray(G, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    n_c = len(mu)
    lam = np.zeros(3 * n_c) if warm is None else np.array(warm, dtype=np.float64).copy()
    if n_c == 0:
        return (lam, np.zeros(0), np.empty(0, dtype='<U5'))
    rho = _rho_of(G, n_c)
    blocks = [G[3 * c:3 * c + 3, 3 * c:3 * c + 3] for c in range(n_c)]
    hist = []
    for _ in range(int(iters)):
        for c in range(n_c):
            sl = slice(3 * c, 3 * c + 3)
            w = b[sl] + G[sl, :] @ lam - blocks[c] @ lam[sl]
            lam[sl] = _local_solve(blocks[c], w, mu[c], desax, lam[sl], rho[c])
        r = natural_residual(lam, G, b, mu, rho, desax)
        hist.append(r)
        if r < tol:
            break
    return (lam, np.array(hist), classify(lam, G, b, mu))

def solve_ncp_pgs(G, b, mu, iters=2000, tol=1e-12, warm=None):
    """A.1 PGS with de Saxce; tuple-compatible ContactResult has .info and .force_selection."""
    return ContactResult(_pgs(G, b, mu, iters, tol, warm, desax=True), G, b, mu, 'PGS_ORDER_WARM_START')

def solve_convex_pgs(G, b, mu, iters=2000, tol=1e-12, warm=None):
    """A.3  Anitescu / MuJoCo-style convex relaxation (no Gamma)."""
    return ContactResult(_pgs(G, b, mu, iters, tol, warm, desax=False), G, b, mu, 'PGS_ORDER_WARM_START', model='CCP')

def _admm_core(G, b, mu, iters, tol, warm, rho, adapt_every, desax, xstep='legacy', diagnostics=None, physical=None):
    """One ADMM run; see solve_ncp_admm for the algorithm and the warm-start caveat."""
    G = np.asarray(G, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    n_c = len(mu)
    n = 3 * n_c
    diagnostics = {} if diagnostics is None else diagnostics
    diagnostics.update(status='MAX_ITERATIONS', iterations=0, n_factor=0)
    if n_c == 0:
        diagnostics['status'] = 'CONVERGED'
        return (np.zeros(0), np.zeros(0), np.empty(0, dtype='<U5'))
    z = np.zeros(n) if warm is None else np.array(warm, dtype=np.float64).copy()
    lam = z.copy()
    if warm is None:
        gamma = np.zeros(n)
    else:
        u0 = (G @ z + b).reshape(n_c, 3)
        gamma = -(G @ z + b + (desaxce(u0, mu).reshape(-1) if desax else 0.0))
    z_prev = z.copy()
    if rho is None:
        rho = float(max(np.trace(G) / (3.0 * n_c), 1e-12))
    rho_floor = 64.0 * np.finfo(float).eps * float(np.linalg.norm(G, ord=np.inf))
    if physical is not None:
        rho = max(rho, rho_floor)
    rho_nat = _rho_of(G, n_c)
    if xstep not in ('legacy', 'triangular'):
        raise ValueError("xstep must be 'legacy' or 'triangular'")
    try:
        if xstep == 'legacy':
            raise ImportError('legacy arithmetic requested')
        from scipy.linalg import cholesky as _chol
        from scipy.linalg.blas import dtrsv as _dtrsv

        def _factor(A):
            return np.asfortranarray(_chol(A, lower=True, check_finite=False))

        def _solve_with(c, rhs):
            return _dtrsv(c, _dtrsv(c, rhs, lower=1), lower=1, trans=1)
    except ImportError:

        def _factor(A):
            return np.linalg.cholesky(A)

        def _solve_with(c, rhs):
            return np.linalg.solve(c.T, np.linalg.solve(c, rhs))
    chol = _factor(G + rho * np.eye(n))
    diagnostics['n_factor'] += 1
    if not np.isfinite(chol).all():
        raise FloatingPointError('non-finite ADMM Cholesky factor')
    hist = []
    best_r = np.inf
    stale = 0

    def _solve_sys(rhs):
        return _solve_with(chol, rhs)
    for k in range(int(iters)):
        u = (G @ z + b).reshape(n_c, 3)
        g = b + (desaxce(u, mu).reshape(-1) if desax else 0.0)
        lam = -_solve_sys(g + gamma - rho * z)
        if not all((np.isfinite(a).all() for a in (u, g, lam))):
            diagnostics['status'] = 'NUMERICAL_FAILURE'
            break
        z_prev = z
        z_next = proj_cone((lam + gamma / rho).reshape(n_c, 3), mu).reshape(-1)
        gamma_next = gamma + rho * (lam - z_next)
        if not np.isfinite(z_next).all() or not np.isfinite(gamma_next).all():
            diagnostics['status'] = 'NUMERICAL_FAILURE'
            break
        if physical is None:
            r = natural_residual(z_next, G, b, mu, rho_nat, desax)
        else:
            (Gp, bp, dp, rp) = physical
            r = natural_residual(z_next * dp, Gp, bp, mu, rp, desax)
        if not np.isfinite(r):
            diagnostics['status'] = 'NUMERICAL_FAILURE'
            break
        (z, gamma) = (z_next, gamma_next)
        hist.append(r)
        diagnostics['iterations'] = k + 1
        if r < tol:
            diagnostics['status'] = 'CONVERGED'
            break
        stale = 0 if r < best_r else stale + 1
        best_r = min(best_r, r)
        impulse_scale = np.max(np.abs(z if physical is None else z * physical[2]))
        if stale >= 100 and best_r <= 64.0 * np.finfo(float).eps * max(1.0, impulse_scale):
            diagnostics['status'] = 'PRECISION_LIMIT'
            break
        if adapt_every and (k + 1) % adapt_every == 0:
            pr = float(np.linalg.norm(lam - z))
            dr = max(rho * float(np.linalg.norm(z - z_prev)), 1e-300)
            new_rho = rho
            if physical is not None:
                pr *= rho
            if pr > 10.0 * dr:
                new_rho = rho * 2.0
            elif dr > 10.0 * pr:
                new_rho = rho / 2.0
            if new_rho < rho_floor:
                if physical is None:
                    diagnostics['status'] = 'PRECISION_LIMIT'
                    break
                new_rho = rho_floor
            if new_rho != rho:
                if physical is None:
                    gamma = gamma * (new_rho / rho)
                rho = new_rho
                chol = _factor(G + rho * np.eye(n))
                diagnostics['n_factor'] += 1
                if not np.isfinite(chol).all():
                    diagnostics['status'] = 'NUMERICAL_FAILURE'
                    break
    diagnostics.update(rho=float(rho), rho_floor=float(rho_floor))
    return (z, np.array(hist), classify(z, G, b, mu))

def solve_ncp_admm(G, b, mu, iters=5000, tol=1e-12, warm=None, rho=None, adapt_every=25, desax=True, cold_retry=True, *, xstep='legacy'):
    """A.2  ADMM on   min 1/2 lam^T G lam + (b + Gamma(u))^T lam   s.t. lam in K.

    lam <- -(G + rho I)^{-1} (g + gamma - rho z);  z <- proj_K(lam + gamma/rho);
    gamma <- gamma + rho (lam - z).                       [2304.06372 Alg. 3]
    The de Saxce shift Gamma(u) is refreshed from u = G z + b each iteration; its
    fixed point is exactly the NCP natural map (see module docstring).

    MEASURED NEGATIVE: a warm lam from the previous time step can
    drive this one-loop ADMM into a non-converging regime on stiff scenes (the
    mass-ratio-100 stack diverges within 6 steps), even with a consistent dual
    init.  With cold_retry=True a failed warm solve is redone from zero, so the
    returned residual history includes both attempts. A failed legacy solve
    is retried in contact-equilibrated coordinates, with a scale-matched
    penalty and a roundoff floor. One positive scalar scales all three rows
    of each contact, preserving the Coulomb cone and de Saxce correction.
    Every stopping residual is evaluated in the original impulse coordinates.

    Returns tuple-compatible ContactResult(lam, residual_history, labels).
    xstep='legacy' preserves the original arithmetic and bytes. Opt-in
    xstep='triangular' accelerates the x-step when scipy is available; it can
    change rounding, residual histories and stopping iterations. The numpy
    fallback remains available when scipy cannot be imported.
    .solver_status reports CONVERGED, PRECISION_LIMIT, MAX_ITERATIONS or
    NUMERICAL_FAILURE; a finite impulse does not imply convergence.
    .solver_diagnostics records all attempts and the returned residual.
    Invalid/non-finite input raises ValueError; failed finite arithmetic never
    returns NaN impulses. .info additionally diagnoses local rank on request."""
    G = np.asarray(G, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    n = 3 * len(mu)
    if G.shape != (n, n) or b.shape != (n,) or mu.ndim != 1:
        raise ValueError('contact arrays have incompatible shapes')
    if not all((np.isfinite(a).all() for a in (G, b, mu))) or np.any(mu < 0):
        raise ValueError('ADMM inputs must be finite with nonnegative friction')
    if warm is not None:
        warm = np.asarray(warm, dtype=np.float64)
        if warm.shape != (n,) or not np.isfinite(warm).all():
            raise ValueError('ADMM warm start must be a finite impulse vector')
    if not np.isfinite(tol) or tol <= 0 or int(iters) < 0:
        raise ValueError('ADMM tolerance must be positive and iteration budget nonnegative')
    if rho is not None and (not np.isfinite(rho) or rho <= 0):
        raise ValueError('ADMM rho must be positive and finite')
    if xstep not in ('legacy', 'triangular'):
        raise ValueError("xstep must be 'legacy' or 'triangular'")
    attempts = []
    outputs = []

    def attempt(Ag, ab, aw, arho, selection, physical=None):
        diag = {'selection': selection}
        with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
            try:
                vals = _admm_core(Ag, ab, mu, iters, tol, aw, arho, adapt_every, desax, xstep, diag, physical)
            except (np.linalg.LinAlgError, FloatingPointError) as exc:
                z = np.zeros(n) if aw is None else aw.copy()
                vals = (z, np.zeros(0), np.full(len(mu), 'open', dtype='<U5'))
                diag.update(status='NUMERICAL_FAILURE', reason=str(exc))
            (lam, hist, labels) = vals
            if physical is not None:
                lam = lam * physical[2]
                labels = classify(lam, G, b, mu) if n else labels
            residual = natural_residual(lam, G, b, mu, desax=desax) if n else 0.0
        if not np.isfinite(lam).all() or not np.isfinite(residual):
            raise FloatingPointError('ADMM could not produce a finite residual and impulse')
        diag['residual'] = residual
        attempts.append(diag)
        outputs.append((lam, hist, labels))
        return residual < tol
    ok = attempt(G, b, warm, rho, 'ADMM_WARM_START')
    if not ok and cold_retry and (warm is not None):
        ok = attempt(G, b, None, rho, 'ADMM_COLD_RETRY')
    if not ok and n and (int(iters) > 0):
        norms = np.array([np.linalg.norm(G[i:i + 3, i:i + 3], 2) for i in range(0, n, 3)])
        if np.all(norms > 0) and np.isfinite(norms).all():
            d = np.repeat(1.0 / np.sqrt(norms), 3)
            Gh = G * d[:, None] * d[None, :]
            bh = b * d
            physical = (G, b, d, _rho_of(G, len(mu)))
            ok = attempt(Gh, bh, None, None, 'ADMM_EQUILIBRATED_RETRY', physical)
    pick = min(range(len(attempts)), key=lambda i: attempts[i]['residual'])
    (lam, _, labels) = outputs[pick]
    hist = np.concatenate([o[1] for o in outputs])
    residual = attempts[pick]['residual']
    if len(hist) and hist[-1] != residual:
        hist = np.append(hist, residual)
    status = 'CONVERGED' if residual < tol else attempts[pick]['status']
    diag = dict(status=status, residual=residual, selected_attempt=pick, attempts=attempts, iterations=sum((a['iterations'] for a in attempts)), n_factor=sum((a['n_factor'] for a in attempts)))
    return ContactResult((lam, hist, labels), G, b, mu, attempts[pick]['selection'], model='NCP' if desax else 'CCP', solver_info=diag)

def pyramid_bound(mu, scale):
    return scale * mu

def solve_pyramid_qp(G, b, mu, iters=5000, tol=1e-12, warm=None, scale=1.0 / math.sqrt(2.0), method='scipy'):
    """A.4  4-facet pyramid  P = {lam_n >= 0, |lam_t1| <= s mu lam_n, |lam_t2| <= s mu lam_n}.

    scale = 1/sqrt(2)  -> pyramid INSCRIBED in K_mu (never over-estimates friction);
    scale = 1          -> circumscribed box (ODE/Bullet convention).
    method='pgs'   : coordinate PGS on the QP  min 1/2 lam^T G lam + b^T lam, lam in P.
    method='scipy' : scipy.optimize.minimize(SLSQP), independent cross-check.
    Returns (lam, residual_history, labels)."""
    G = np.asarray(G, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    n_c = len(mu)
    n = 3 * n_c
    if n_c == 0:
        return (np.zeros(0), np.zeros(0), np.empty(0, dtype='<U5'))
    if method == 'scipy':
        from scipy.optimize import minimize
        A = np.zeros((5 * n_c, n))
        for c in range(n_c):
            (i, s) = (3 * c, scale * mu[c])
            A[5 * c + 0, i] = 1.0
            (A[5 * c + 1, i], A[5 * c + 1, i + 1]) = (s, -1.0)
            (A[5 * c + 2, i], A[5 * c + 2, i + 1]) = (s, 1.0)
            (A[5 * c + 3, i], A[5 * c + 3, i + 2]) = (s, -1.0)
            (A[5 * c + 4, i], A[5 * c + 4, i + 2]) = (s, 1.0)
        x0 = np.zeros(n) if warm is None else np.array(warm, dtype=np.float64)
        res = minimize(lambda x: 0.5 * x @ (G @ x) + b @ x, x0, jac=lambda x: G @ x + b, constraints=[{'type': 'ineq', 'fun': lambda x: A @ x, 'jac': lambda x: A}], method='SLSQP', options={'maxiter': int(iters), 'ftol': 1e-14})
        lam = res.x
        return (lam, np.array([_pyr_residual(lam, G, b, mu, scale)]), classify(lam, G, b, mu))
    lam = np.zeros(n) if warm is None else np.array(warm, dtype=np.float64).copy()
    d = np.diag(G).copy()
    d[d <= 0] = 1.0
    hist = []
    for _ in range(int(iters)):
        prev = lam.copy()
        for c in range(n_c):
            i = 3 * c
            r = G[i, :] @ lam + b[i]
            lam[i] = max(0.0, lam[i] - r / d[i])
            bnd = scale * mu[c] * lam[i]
            for k in (1, 2):
                r = G[i + k, :] @ lam + b[i + k]
                lam[i + k] = min(bnd, max(-bnd, lam[i + k] - r / d[i + k]))
        hist.append(float(np.max(np.abs(lam - prev))))
        if hist[-1] < tol:
            break
    return (lam, np.array(hist), classify(lam, G, b, mu))

def _pyr_residual(lam, G, b, mu, scale):
    """PGS fixed-point residual for the pyramid QP (stationarity measure)."""
    n_c = len(mu)
    G = np.asarray(G)
    x = np.array(lam, dtype=np.float64)
    d = np.diag(G).copy()
    d[d <= 0] = 1.0
    prev = x.copy()
    for c in range(n_c):
        i = 3 * c
        x[i] = max(0.0, x[i] - (G[i, :] @ x + b[i]) / d[i])
        bnd = scale * mu[c] * x[i]
        for k in (1, 2):
            x[i + k] = min(bnd, max(-bnd, x[i + k] - (G[i + k, :] @ x + b[i + k]) / d[i + k]))
    return float(np.max(np.abs(x - prev)))

def _active_system(lam, G, b, mu, labels):
    """Assemble fixed-branch residual Jacobian and parameter maps."""
    lam = np.asarray(lam, dtype=np.float64)
    G = np.asarray(G, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    n_c = len(mu)
    n = 3 * n_c
    if G.shape != (n, n) or b.shape != (n,) or lam.shape != (n,) or (len(labels) != n_c):
        raise ValueError('contact arrays and labels have incompatible shapes')
    if any((lab not in ('open', 'stick', 'slip') for lab in labels)):
        raise ValueError('labels must be open, stick or slip')
    if not all((np.isfinite(x).all() for x in (lam, G, b, mu))):
        raise ValueError('contact inputs must be finite')
    u = G @ lam + b
    A = np.zeros((n, n))
    Bmu = np.zeros((n, n_c))
    Bb = np.zeros((n, n))
    I2 = np.eye(2)
    boundary = []
    for c in range(n_c):
        i = 3 * c
        lab = labels[c]
        if lab == 'open':
            A[i:i + 3, i:i + 3] = np.eye(3)
            if abs(u[i]) < 1e-09:
                boundary.append((c, 'open_with_zero_normal_velocity'))
        elif lab == 'stick':
            A[i:i + 3, :] = G[i:i + 3, :]
            Bb[i:i + 3, i:i + 3] = np.eye(3)
        else:
            A[i, :] = G[i, :]
            Bb[i, i] = 1.0
            s = u[i + 1:i + 3]
            ns = float(np.linalg.norm(s))
            if ns < 1e-14:
                boundary.append((c, 'slip_with_zero_slip_speed'))
                ns = 1e-14
            sh = s / ns
            P = (I2 - np.outer(sh, sh)) / ns
            rows = slice(i + 1, i + 3)
            A[rows, :] = mu[c] * lam[i] * (P @ G[rows, :])
            A[rows, i] += mu[c] * sh
            A[i + 1, i + 1] += 1.0
            A[i + 2, i + 2] += 1.0
            Bb[rows, rows] = mu[c] * lam[i] * P
            Bmu[rows, c] = lam[i] * sh
        (ln, lt) = (lam[i], float(np.linalg.norm(lam[i + 1:i + 3])))
        scale = max(abs(ln), 1e-12)
        if lab != 'open' and abs(lt - mu[c] * ln) < 1e-07 * mu[c] * scale and (lab == 'stick'):
            boundary.append((c, 'stick_on_cone_boundary'))
        if lab != 'open' and ln < 1e-09:
            boundary.append((c, 'near_zero_normal_impulse'))
    return (A, Bmu, Bb, boundary)

def contact_diagnostics(lam, G, b, mu, labels, rcond=1e-12, observable=None):
    """Numerical local rank, nullspace and optional invariant observable rows.

    Labels alone do not establish feasibility or differentiability. Boundary
    diagnostics invalidate differentiability even when the matrix has full rank.
    """
    (A, Bmu, Bb, boundary) = _active_system(lam, G, b, mu, labels)
    factors = rank_factor(A, rcond)
    info = factors[3]
    info.update(boundary=boundary, labels=list(labels), residual_differentiable=not boundary, response_scope='NONSMOOTH_BRANCH' if boundary else 'LINEARIZED_RESPONSE_ONLY' if info['singular'] else 'IMPLICIT_RESPONSE_UNDER_VALID_BRANCH_ASSUMPTIONS')
    if observable is not None:
        (_, check) = guarded_response(A, np.zeros((A.shape[0], 0)), factors, observable=observable)
        info.update(observable_invariant=check['observable_invariant'], observable_nullspace_residual=check['observable_nullspace_residual'])
    return info

class ContactResult(tuple):
    """Three-item legacy tuple with an explicit algorithm-selected representative.

    .info performs the dense local rank diagnostic only on request. Iteration
    order, warm start and stopping tolerance select lam; this is not a proof
    that a per-contact force distribution is unique. Large callers can consume
    body motion without paying for a dense sensitivity diagnostic. ADMM also
    exposes .solver_status and .solver_diagnostics, distinct from the local
    rank status in .info. Query .info before mutating input or result arrays.
    """

    def __new__(cls, values, G, b, mu, selection, model='NCP', solver_info=None):
        out = super().__new__(cls, values)
        out.force_selection = selection
        out.contact_model = model
        out.solver_diagnostics = solver_info
        out.solver_status = None if solver_info is None else solver_info['status']
        out._problem = (G, b, mu)
        out._info = None
        return out

    def __reduce__(self):
        (G, b, mu) = self._problem
        return (type(self), (tuple(self), G, b, mu, self.force_selection, self.contact_model, self.solver_diagnostics))

    @property
    def info(self):
        if self._info is None:
            if self.contact_model != 'NCP':
                self._info = {'status': 'UNASSESSED_CONTACT_MODEL', 'scope': self.contact_model}
            else:
                (G, b, mu) = self._problem
                self._info = contact_diagnostics(self[0], G, b, mu, self[2])
            self._info['force_selection'] = self.force_selection
            if self.solver_diagnostics is not None:
                self._info['solver_status'] = self.solver_status
                self._info['solver_diagnostics'] = self.solver_diagnostics
        return self._info

def sensitivity(lam, G, b, mu, labels, rcond=1e-12, *, selection=None, observable=None, db_dparam=None):
    """Fixed-branch impulse sensitivity with a numerical uniqueness guard.

    open: lam_c=0; stick: u_c=0; slip: u_n=0 and
    lam_t + mu lam_n u_t/||u_t||=0, with u=G lam+b.
    Returns (dlam_dmu, dlam_db, info). Ambiguous force derivatives are NaN
    unless selection='minimum_norm' is explicitly requested. Inconsistent
    perturbations and nonsmooth branch boundaries remain NaN in either mode.
    The selection is a linear-response convention, not a physical compliance.

    observable=Q exposes derivatives of Q lam in info only for invariant rows
    Q N=0 and compatible RHS. db_dparam maps physical parameters to b, e.g.
    Scene.J maps free body velocity; arbitrary independent corner b changes
    need not be compatible with a closed redundant contact branch. These are
    local responses at the supplied impulse representative: d(Q lam)/dmu can
    depend on that representative even when Q lam itself is invariant.
    Compatibility alone does not establish a differentiable nonlinear solution
    branch. Singular observable outputs are formal linearized responses, as
    reported by info["response_scope"].
    """
    (A, Bmu, Bb, boundary) = _active_system(lam, G, b, mu, labels)
    if db_dparam is not None:
        P = np.asarray(db_dparam, dtype=np.float64)
        if P.ndim != 2 or P.shape[0] != Bb.shape[1] or (not np.isfinite(P).all()):
            raise ValueError('db_dparam must be a finite matrix with 3*n_c rows')
        Bb = Bb @ P
    factors = rank_factor(A, rcond)
    (dm, im) = guarded_response(A, Bmu, factors, selection, observable)
    (db, ib) = guarded_response(A, Bb, factors, selection, observable)
    info = factors[3]
    info.update(boundary=boundary, labels=list(labels), residual_differentiable=not boundary, response_scope='NONSMOOTH_BRANCH' if boundary else 'LINEARIZED_RESPONSE_ONLY' if info['singular'] else 'IMPLICIT_RESPONSE_UNDER_VALID_BRANCH_ASSUMPTIONS', selection=im['selection'], mu_compatible=im['compatible'], b_compatible=ib['compatible'], mu_response_finite=im['response_finite'], b_response_finite=ib['response_finite'], mu_compatibility_residual=im['compatibility_residual'], b_compatibility_residual=ib['compatibility_residual'])
    if observable is not None:
        info.update(observable_dmu=im['observable'], observable_db=ib['observable'], observable_mu_valid=im['observable_valid'], observable_b_valid=ib['observable_valid'], observable_invariant=im['observable_invariant'], observable_nullspace_residual=im['observable_nullspace_residual'])
    if boundary:
        dm[:] = np.nan
        db[:] = np.nan
        if observable is not None:
            info['observable_dmu'][:] = np.nan
            info['observable_db'][:] = np.nan
            info['observable_mu_valid'][:] = False
            info['observable_b_valid'][:] = False
    return (dm, db, info)

def sensitivity_fd(G, b, mu, lam_star, h_mu=1e-06, h_b=1e-06, solver=None, iters=20000, tol=1e-14):
    """Central finite differences of the NCP solution w.r.t. mu and b."""
    solver = solver or solve_ncp_pgs
    G = np.asarray(G, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    n_c = len(mu)
    n = 3 * n_c
    dmu = np.zeros((n, n_c))
    db = np.zeros((n, n))
    for c in range(n_c):
        (mp, mm) = (mu.copy(), mu.copy())
        mp[c] += h_mu
        mm[c] -= h_mu
        lp = solver(G, b, mp, iters, tol, warm=lam_star)[0]
        lm = solver(G, b, mm, iters, tol, warm=lam_star)[0]
        dmu[:, c] = (lp - lm) / (2.0 * h_mu)
    for j in range(n):
        (bp, bm) = (b.copy(), b.copy())
        bp[j] += h_b
        bm[j] -= h_b
        lp = solver(G, bp, mu, iters, tol, warm=lam_star)[0]
        lm = solver(G, bm, mu, iters, tol, warm=lam_star)[0]
        db[:, j] = (lp - lm) / (2.0 * h_b)
    return (dmu, db)
GRAVITY = np.array([0.0, 0.0, -9.81])

def _skew(r):
    return np.array([[0.0, -r[2], r[1]], [r[2], 0.0, -r[0]], [-r[1], r[0], 0.0]])

def _quat_to_R(q):
    (w, x, y, z) = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)], [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)], [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])

def _quat_integrate(q, w, dt):
    ww = np.array([0.0, w[0], w[1], w[2]])
    (qw, qx, qy, qz) = q
    (vw, vx, vy, vz) = ww
    dq = 0.5 * np.array([qw * vw - qx * vx - qy * vy - qz * vz, qw * vx + qx * vw + qy * vz - qz * vy, qw * vy - qx * vz + qy * vw + qz * vx, qw * vz + qx * vy - qy * vx + qz * vw])
    q = q + dt * dq
    return q / np.linalg.norm(q)

@dataclass
class Box:
    half: np.ndarray
    mass: float
    pos: np.ndarray
    quat: np.ndarray = field(default_factory=lambda : np.array([1.0, 0.0, 0.0, 0.0]))
    vel: np.ndarray = field(default_factory=lambda : np.zeros(3))
    omega: np.ndarray = field(default_factory=lambda : np.zeros(3))

    def I_body(self):
        (hx, hy, hz) = self.half
        m = self.mass
        return np.diag([m / 3.0 * (hy ** 2 + hz ** 2), m / 3.0 * (hx ** 2 + hz ** 2), m / 3.0 * (hx ** 2 + hy ** 2)])

    def R(self):
        return _quat_to_R(self.quat)

    def I_world_inv(self):
        R = self.R()
        return R @ np.linalg.inv(self.I_body()) @ R.T

    def vertices(self):
        (R, h) = (self.R(), self.half)
        out = []
        for sx in (-1, 1):
            for sy in (-1, 1):
                for sz in (-1, 1):
                    out.append(self.pos + R @ (np.array([sx, sy, sz]) * h))
        return np.array(out)

def _tangents(n):
    a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    t1 = np.cross(n, a)
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(n, t1)
    return (t1, t2)

def contacts_box_plane(box, idx, margin=0.01, plane_z=0.0):
    """Vertex-plane contacts.  normal = +z (ground -> box).  Returns list of dicts."""
    out = []
    for v in box.vertices():
        gap = v[2] - plane_z
        if gap < margin:
            out.append({'p': v, 'n': np.array([0.0, 0.0, 1.0]), 'gap': gap, 'a': idx, 'b': -1})
    return out
_BOX_AXES = [np.eye(3)[i] for i in range(3)]

def _face_polygon(box, axis_local, sign):
    """4 world vertices of the face whose outward local normal is sign*e_axis."""
    (R, h, c) = (box.R(), box.half, box.pos)
    others = [i for i in range(3) if i != axis_local]
    pts = []
    for s0 in (-1, 1):
        for s1 in (-1, 1):
            l = np.zeros(3)
            l[axis_local] = sign * h[axis_local]
            l[others[0]] = s0 * h[others[0]]
            l[others[1]] = s1 * h[others[1]]
            pts.append(c + R @ l)
    return [pts[0], pts[1], pts[3], pts[2]]

def _clip_poly(poly, plane_n, plane_d, eps=1e-09):
    """Keep the part of poly with  plane_n . x <= plane_d  (Sutherland-Hodgman).

    `eps` treats a vertex lying on the plane as inside and suppresses the spurious
    intersection vertex it would otherwise generate.  Without it, face-face contact
    between exactly aligned boxes emits a degenerate manifold (a duplicated vertex
    displacing a real corner) as soon as positions carry ~1e-13 m of rounding; that
    manifold is rank-deficient in a way the solvers cannot resolve (measured: the
    mass-ratio-100 stack lost a corner at step 5)."""
    out = []
    m = len(poly)
    for i in range(m):
        (a, bpt) = (poly[i], poly[(i + 1) % m])
        da = float(plane_n @ a) - plane_d
        db = float(plane_n @ bpt) - plane_d
        if da <= eps:
            out.append(a)
        if da > eps and db < -eps:
            out.append(a + da / (da - db) * (bpt - a))
        elif da < -eps and db > eps:
            out.append(a + da / (da - db) * (bpt - a))
    return out

def _dedup_points(pts, tol=1e-07):
    """Drop points that coincide within `tol` (clipping emits duplicates on
    exactly-aligned faces)."""
    out = []
    for p in pts:
        if not any((float(np.linalg.norm(p - q)) <= tol for q in out)):
            out.append(p)
    return out

def contacts_box_box(A, ia, B, ib, margin=0.01):
    """SAT over the 6 face axes + reference-face clipping, <= 4 points.
    normal n points from B to A.  Edge-edge contacts are NOT generated
    (documented limitation: adequate for the stacked-box scenes of SPEC A.6)."""
    (RA, RB) = (A.R(), B.R())
    axesA = [RA[:, i] for i in range(3)]
    axesB = [RB[:, i] for i in range(3)]
    d = A.pos - B.pos

    def overlap(ax):
        ra = float(np.sum(np.abs(np.array([float(ax @ axesA[i]) for i in range(3)])) * A.half))
        rb = float(np.sum(np.abs(np.array([float(ax @ axesB[i]) for i in range(3)])) * B.half))
        return ra + rb - abs(float(ax @ d))
    (best, best_ax, best_from) = (np.inf, None, None)
    for (i, ax) in enumerate(axesA):
        o = overlap(ax)
        if o < best:
            (best, best_ax, best_from) = (o, ax, ('A', i))
    for (i, ax) in enumerate(axesB):
        o = overlap(ax)
        if o < best:
            (best, best_ax, best_from) = (o, ax, ('B', i))
    if best_ax is None or best < -margin:
        return []
    n = best_ax.copy()
    if float(n @ d) < 0:
        n = -n
    (which, ai) = best_from
    if which == 'A':
        (ref_box, inc_box) = (A, B)
        ref_n = -n
        ref_sign = -1.0 if float(axesA[ai] @ n) > 0 else 1.0
    else:
        (ref_box, inc_box) = (B, A)
        ref_n = n
        ref_sign = 1.0 if float(axesB[ai] @ n) > 0 else -1.0
    inc_axes = [inc_box.R()[:, i] for i in range(3)]
    (bi, bs, bv) = (0, 1.0, np.inf)
    for i in range(3):
        for s in (-1.0, 1.0):
            val = float(s * inc_axes[i] @ ref_n)
            if val < bv:
                (bv, bi, bs) = (val, i, s)
    poly = _face_polygon(inc_box, bi, bs)
    ref_poly = _face_polygon(ref_box, ai, ref_sign)
    m = len(ref_poly)
    for i in range(m):
        e = ref_poly[(i + 1) % m] - ref_poly[i]
        sn = np.cross(e, ref_n)
        nl = np.linalg.norm(sn)
        if nl < 1e-12:
            continue
        sn = sn / nl
        if float(sn @ (ref_box.pos - ref_poly[i])) > 0:
            sn = -sn
        poly = _clip_poly(poly, sn, float(sn @ ref_poly[i]))
        if not poly:
            return []
    poly = _dedup_points(poly)
    p0 = ref_poly[0]
    out = []
    for p in poly:
        gap = float(ref_n @ (p - p0))
        if gap > margin:
            continue
        out.append({'p': p - 0.5 * gap * ref_n, 'n': n, 'gap': gap, 'a': ia, 'b': ib})
    if len(out) > 4:
        out.sort(key=lambda cc: cc['gap'])
        keep = [out[0]]
        rest = out[1:]
        while len(keep) < 4 and rest:
            j = max(range(len(rest)), key=lambda i: min((float(np.linalg.norm(rest[i]['p'] - k['p'])) for k in keep)))
            keep.append(rest.pop(j))
        out = keep
    return out

def build_scene(bodies, contacts, mu, dt, v_free=None, speculative=True):
    """Assemble Scene (J, Minv diagonal-block, v_free, mu, body_pairs, dt)."""
    nb = len(bodies)
    n_v = 6 * nb
    n_c = len(contacts)
    J = np.zeros((3 * n_c, n_v))
    pairs = np.zeros((n_c, 2), dtype=np.int64)
    off = np.zeros(3 * n_c)
    Minv = np.zeros((n_v, n_v))
    for (k, bd) in enumerate(bodies):
        Minv[6 * k:6 * k + 3, 6 * k:6 * k + 3] = np.eye(3) / bd.mass
        Minv[6 * k + 3:6 * k + 6, 6 * k + 3:6 * k + 6] = bd.I_world_inv()
    for (c, ct) in enumerate(contacts):
        n = ct['n']
        (t1, t2) = _tangents(n)
        Bm = np.vstack([n, t1, t2])
        (ia, ib) = (ct['a'], ct['b'])
        pairs[c] = (ia, ib)
        if ia >= 0:
            r = ct['p'] - bodies[ia].pos
            J[3 * c:3 * c + 3, 6 * ia:6 * ia + 3] = Bm
            J[3 * c:3 * c + 3, 6 * ia + 3:6 * ia + 6] = -Bm @ _skew(r)
        if ib >= 0:
            r = ct['p'] - bodies[ib].pos
            J[3 * c:3 * c + 3, 6 * ib:6 * ib + 3] = -Bm
            J[3 * c:3 * c + 3, 6 * ib + 3:6 * ib + 6] = Bm @ _skew(r)
        if speculative and ct['gap'] > 0:
            off[3 * c] = ct['gap'] / dt
    if v_free is None:
        v_free = np.zeros(n_v)
    mu_arr = np.full(n_c, mu, dtype=np.float64) if np.isscalar(mu) else np.asarray(mu, float)
    return Scene(J=J, Minv=Minv, v_free=np.asarray(v_free, float), mu=mu_arr, body_pairs=pairs, dt=dt, b_offset=off)

def collect_contacts(bodies, ground=True, margin=0.01):
    ct = []
    for (i, bd) in enumerate(bodies):
        if ground:
            ct += contacts_box_plane(bd, i, margin=margin)
    for i in range(len(bodies)):
        for j in range(i + 1, len(bodies)):
            ct += contacts_box_box(bodies[j], j, bodies[i], i, margin=margin)
    return ct

def simulate(bodies, steps, dt, mu, solver=solve_ncp_pgs, iters=500, tol=1e-12, force_fn=None, margin=0.01, ground=True, record=None, warm=True):
    """Semi-implicit step loop.  Returns a dict of recorded diagnostics."""
    hist = {'t': [], 'min_gap': [], 'max_pen': [], 'res': [], 'iters': [], 'n_c': [], 'pos': [], 'vel': [], 'lam_sum': []}
    lam_prev = None
    for s in range(int(steps)):
        t = s * dt
        v_free = np.zeros(6 * len(bodies))
        for (k, bd) in enumerate(bodies):
            f = bd.mass * GRAVITY.copy()
            tau = np.zeros(3)
            if force_fn is not None:
                (ff, tt) = force_fn(t, k, bd)
                f = f + ff
                tau = tau + tt
            Iw = np.linalg.inv(bd.I_world_inv())
            v_free[6 * k:6 * k + 3] = bd.vel + dt * f / bd.mass
            v_free[6 * k + 3:6 * k + 6] = bd.omega + dt * bd.I_world_inv() @ (tau - np.cross(bd.omega, Iw @ bd.omega))
        ct = collect_contacts(bodies, ground=ground, margin=margin)
        gaps = np.array([c['gap'] for c in ct]) if ct else np.array([np.inf])
        if ct:
            sc = build_scene(bodies, ct, mu, dt, v_free=v_free)
            (G, b) = (sc.delassus(), sc.b())
            w = lam_prev if warm and lam_prev is not None and (len(lam_prev) == 3 * len(ct)) else None
            (lam, res, labels) = solver(G, b, sc.mu, iters, tol, w)
            v = sc.apply(lam)
            lam_prev = lam
            hist['res'].append(float(res[-1]) if len(res) else 0.0)
            hist['iters'].append(int(len(res)))
            hist['lam_sum'].append(float(np.sum(lam.reshape(-1, 3)[:, 0])))
        else:
            v = v_free
            lam_prev = None
            hist['res'].append(0.0)
            hist['iters'].append(0)
            hist['lam_sum'].append(0.0)
        for (k, bd) in enumerate(bodies):
            bd.vel = v[6 * k:6 * k + 3]
            bd.omega = v[6 * k + 3:6 * k + 6]
            bd.pos = bd.pos + dt * bd.vel
            bd.quat = _quat_integrate(bd.quat, bd.omega, dt)
        hist['t'].append(t)
        hist['n_c'].append(len(ct))
        hist['min_gap'].append(float(np.min(gaps)))
        hist['max_pen'].append(float(max(0.0, -np.min(gaps))) if ct else 0.0)
        hist['pos'].append(np.array([b.pos.copy() for b in bodies]))
        hist['vel'].append(np.array([b.vel.copy() for b in bodies]))
        if record is not None:
            record(s, bodies, hist)
    for k in ('t', 'min_gap', 'max_pen', 'res', 'iters', 'n_c', 'lam_sum'):
        hist[k] = np.array(hist[k])
    hist['pos'] = np.array(hist['pos'])
    hist['vel'] = np.array(hist['vel'])
    return hist

def scene_pushed_cube(dt=0.004, mu=0.5, size=0.2, mass=15.0, force=147.15, period=0.5, ramp=0.3):
    """(i) ODYNSim pushed cube: 0.20 m, 15 kg, mu=0.5, dt=0.004,
    horizontal force alternating x/y every `period` s, ramped over `ramp` s."""
    h = size / 2.0
    box = Box(half=np.array([h, h, h]), mass=mass, pos=np.array([0.0, 0.0, h]))

    def force_fn(t, k, bd):
        idx = int(t // period)
        d = np.array([1.0, 0.0, 0.0]) if idx % 2 == 0 else np.array([0.0, 1.0, 0.0])
        tau_in = t - idx * period
        s = min(1.0, tau_in / ramp) if ramp > 0 else 1.0
        return (force * s * d, np.zeros(3))
    return ([box], force_fn, dict(dt=dt, mu=mu))

def scene_sliding_box(dt=0.004, mu=0.3, v0=1.0, size=0.2, mass=1.0):
    """(ii) box sliding on a plane with initial velocity v0 -> gliding artifact."""
    h = size / 2.0
    box = Box(half=np.array([h, h, h]), mass=mass, pos=np.array([0.0, 0.0, h]), vel=np.array([v0, 0.0, 0.0]))
    return ([box], None, dict(dt=dt, mu=mu))

def scene_massratio_stack(ratio=1000.0, H=0.2, dt=1.0 / 240.0, mu=0.5, n=3):
    """(iii) 1000:1 mass-ratio stack of 3 boxes (geometry of
    an independent mass-ratio study script: H=0.20,
    ratios 1/100/1000, top box heavy, resting exactly at z_k = (k+0.5) H)."""
    h = H / 2.0
    bodies = []
    for k in range(n):
        m = float(ratio) if k == n - 1 else 1.0
        bodies.append(Box(half=np.array([h, h, h]), mass=m, pos=np.array([0.0, 0.0, (k + 0.5) * H])))
    return (bodies, None, dict(dt=dt, mu=mu))

def scene_tower(n=8, H=0.2, dt=0.004, mu=0.5, mass=1.0):
    """(iv) 8-box tower K8."""
    h = H / 2.0
    bodies = [Box(half=np.array([h, h, h]), mass=mass, pos=np.array([0.0, 0.0, (k + 0.5) * H])) for k in range(n)]
    return (bodies, None, dict(dt=dt, mu=mu))

def scene_tripod(mu=0.4, mass=2.0, dt=0.004, vx=0.8, vy=0.3, seed=0):
    """3 point contacts under a single body -> statically determinate (G invertible,
    active-set Jacobian non-singular) -> clean sensitivity/FD comparison."""
    box = Box(half=np.array([0.1, 0.1, 0.1]), mass=mass, pos=np.array([0.0, 0.0, 0.1]), vel=np.array([vx, vy, 0.0]), omega=np.array([0.0, 0.0, 0.2]))
    pts = [np.array([0.1, 0.0, 0.0]), np.array([-0.05, 0.0866, 0.0]), np.array([-0.05, -0.0866, 0.0])]
    ct = [{'p': p, 'n': np.array([0.0, 0.0, 1.0]), 'gap': 0.0, 'a': 0, 'b': -1} for p in pts]
    v_free = np.concatenate([box.vel + dt * GRAVITY, box.omega])
    sc = build_scene([box], ct, np.full(3, mu), dt, v_free=v_free)
    return sc

def scene_random(n_c=6, n_b=3, mu=0.35, dt=0.004, seed=1):
    """Random well-conditioned synthetic scene for solver cross-checks."""
    rng = np.random.default_rng(seed)
    n_v = 6 * n_b
    J = rng.standard_normal((3 * n_c, n_v))
    Minv = np.diag(np.exp(rng.uniform(-1.0, 1.0, n_v)))
    v_free = rng.standard_normal(n_v) * 0.3
    v_free[2::6] -= 0.5
    pairs = np.zeros((n_c, 2), dtype=np.int64)
    pairs[:, 0] = rng.integers(0, n_b, n_c)
    pairs[:, 1] = -1
    return Scene(J=J, Minv=Minv, v_free=v_free, mu=np.full(n_c, mu), body_pairs=pairs, dt=dt)
