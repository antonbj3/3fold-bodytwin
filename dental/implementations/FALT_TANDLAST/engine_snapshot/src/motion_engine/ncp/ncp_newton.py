"""Body-space Newton solvers for the exact-cone (de Saxce) Coulomb NCP, next to ncp_ref's PGS/ADMM.

Same problem and residual as `ncp_ref`:  u = G lam + b,  G = J Minv J^T,
    K_mu ∋ lam_c ⊥ u_c + Gamma(u_c) ∈ K_mu*,   Gamma(u) = (mu |u_t|, 0, 0),
stopping test the natural-map residual ||lam - P_K(lam - rho (u + Gamma(u)))||_inf, rho_c = 1/||G_cc||_2.

Pipeline `solve_ncp_newton` (measured against Siconos 4.3.1 FC3D NSGS):
  (a) primal-dual interior point (Mehrotra, Nesterov-Todd scaling) on the conic QP
      min 1/2|v - v_f|_M^2  s.t.  J v + off + (s(v),0,0) ∈ K_mu*, with the de Saxce term s(v) = mu|u_t(v)|
      inside the Newton system; each step solves one 6 n_b x 6 n_b body-space system
      M + J^T W^-2 (J + dS J)  (the vertex pencil with the NT scaling as contact compliance);
  (b) semismooth Newton on the natural map, Levenberg-Marquardt shift for hyperstatic (singular)
      directions, each step by Woodbury in body space: one factor of M + J^T Bd^-1 C J;
  (c) fold escape: when Newton stalls, the contacts carrying the near-null Newton mode are tried on
      the take-off branch; an answer is accepted only on the FULL problem's residual;
  (d) friction continuation (mu ramp with branch jumps) as last resort.
`uniqueness_report` diagnoses force fibres and motion linearization at a converged
solution. Floating-point rank is numerical evidence, not a uniqueness certificate.
A singular motion linearization does not prove multiple nonlinear solutions. The
second-start probe reports numerical witnesses after checking both residuals;
agreement cannot establish global uniqueness. Supported input: positive friction
and independent rigid-body mass blocks with six degrees of freedom per body.

Linear algebra: MKL PARDISO through pypardiso when it is installed, otherwise scipy SuperLU.
"""
from __future__ import annotations
import time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

def ldet(x):
    return x[:, 0] ** 2 - np.einsum('ij,ij->i', x[:, 1:], x[:, 1:])

def jprod(x, y):
    """Jordan product x∘y = (x·y, x0 y1 + y0 x1)."""
    out = np.empty_like(x)
    out[:, 0] = np.einsum('ij,ij->i', x, y)
    out[:, 1:] = x[:, :1] * y[:, 1:] + y[:, :1] * x[:, 1:]
    return out

def jsolve(x, d):
    """u with x∘u = d."""
    det = ldet(x)
    u = np.empty_like(d)
    u[:, 0] = (x[:, 0] * d[:, 0] - np.einsum('ij,ij->i', x[:, 1:], d[:, 1:])) / det
    u[:, 1:] = (d[:, 1:] - u[:, :1] * x[:, 1:]) / x[:, :1]
    return u

def max_step(x, d):
    """largest alpha >= 0 with x + alpha d in L (x interior), per cone, then min."""
    a = ldet(d)
    b = x[:, 0] * d[:, 0] - np.einsum('ij,ij->i', x[:, 1:], d[:, 1:])
    c = ldet(x)
    disc = np.maximum(b * b - a * c, 0.0)
    sq = np.sqrt(disc)
    alpha = np.full(len(x), np.inf)
    with np.errstate(divide='ignore', invalid='ignore'):
        q = -(b + np.sign(b + (b == 0)) * sq)
        r1 = np.where(a != 0, q / a, np.inf)
        r2 = np.where(q != 0, c / q, np.inf)
    for r in (r1, r2):
        ok = np.isfinite(r) & (r > 0)
        alpha = np.where(ok, np.minimum(alpha, r), alpha)
    with np.errstate(divide='ignore', invalid='ignore'):
        t0 = np.where(d[:, 0] < 0, -x[:, 0] / d[:, 0], np.inf)
    alpha = np.minimum(alpha, t0)
    return float(alpha.min()) if len(alpha) else np.inf

def nt_scaling(s, z):
    """Nesterov-Todd scaling per cone: returns W, Winv (m,3,3) and lam = W z = W^-1 s."""
    (ds, dz) = (np.sqrt(ldet(s)), np.sqrt(ldet(z)))
    (sb, zb) = (s / ds[:, None], z / dz[:, None])
    gam = np.sqrt(0.5 * (1.0 + np.einsum('ij,ij->i', sb, zb)))
    wb = np.empty_like(s)
    wb[:, 0] = (sb[:, 0] + zb[:, 0]) / (2 * gam)
    wb[:, 1:] = (sb[:, 1:] - zb[:, 1:]) / (2 * gam[:, None])
    eta = np.sqrt(ds / dz)
    m = len(s)
    W = np.empty((m, 3, 3))
    W[:, 0, 0] = wb[:, 0]
    W[:, 0, 1:] = wb[:, 1:]
    W[:, 1:, 0] = wb[:, 1:]
    W[:, 1:, 1:] = np.eye(2)[None] + np.einsum('ij,ik->ijk', wb[:, 1:], wb[:, 1:]) / (1.0 + wb[:, 0])[:, None, None]
    Winv = W.copy()
    Winv[:, 0, 1:] *= -1
    Winv[:, 1:, 0] *= -1
    W *= eta[:, None, None]
    Winv /= eta[:, None, None]
    return (W, Winv)

def bmv(A, x):
    return np.einsum('ijk,ik->ij', A, x)

def _get_ps(mtype):
    """Each factor owner needs its own mutable PARDISO handle."""
    import pypardiso
    return pypardiso.PyPardisoSolver(mtype=mtype)

class BodySolver:
    """Sparse direct solves of the 6 n_b x 6 n_b body-space matrices.  kind = "auto" uses MKL PARDISO
    through pypardiso when installed (reusing the symbolic analysis while the pattern is fixed),
    otherwise scipy SuperLU.  symmetric=False for the de Saxce-coupled / semismooth matrices."""

    def __init__(self, n, kind='auto', symmetric=True):
        self.n = n
        self.symmetric = symmetric
        self.t_fact = self.t_solve = 0.0
        self.n_fact = self.n_solve = 0
        if kind in ('auto', 'pardiso'):
            try:
                import pypardiso
                kind = 'pardiso'
            except Exception:
                kind = 'superlu'
        self.kind = kind
        self._ps = None
        self._A = None

    def factor(self, H):
        t = time.perf_counter()
        if self.kind == 'pardiso':
            A = sp.triu(H, format='csr') if self.symmetric else H.tocsr()
            A.sort_indices()
            same = self._ps is not None and self._A is not None and (A.nnz == self._A.nnz) and np.array_equal(A.indptr, self._A.indptr) and np.array_equal(A.indices, self._A.indices)
            if self._ps is None:
                self._ps = _get_ps(2 if self.symmetric else 11)
                self._ps.free_memory(everything=True)
                same = False
            if not same:
                self._ps.free_memory(everything=True)
                self._ps.set_phase(11)
                self._ps._call_pardiso(A, np.zeros(self.n))
            self._ps.set_phase(22)
            self._ps._call_pardiso(A, np.zeros(self.n))
            self._A = A
        else:
            self._lu = splu(H.tocsc(), permc_spec='MMD_AT_PLUS_A' if self.symmetric else 'COLAMD')
        self.t_fact += time.perf_counter() - t
        self.n_fact += 1

    def solve(self, r):
        t = time.perf_counter()
        if self.kind == 'pardiso':
            self._ps.set_phase(33)
            x = self._ps._call_pardiso(self._A, np.ascontiguousarray(r, float))
        else:
            x = self._lu.solve(np.asarray(r, float))
        self.t_solve += time.perf_counter() - t
        self.n_solve += 1
        return x

    def free(self):
        if self._ps is not None:
            try:
                self._ps.free_memory(everything=True)
            except Exception:
                pass
            self._ps = None

class Problem:
    """J (3nc x 6nb csr), M, Minv (block diag csr), v_free, off, mu."""

    def __init__(self, J, M, Minv, v_free, off, mu):
        self.J = J.tocsr()
        self.M = M.tocsr()
        self.Minv = Minv.tocsr()
        self.v_free = np.asarray(v_free, float)
        self.off = np.asarray(off, float)
        self.mu = np.asarray(mu, float)
        if not np.all(np.isfinite(self.mu)) or not np.all(self.mu > 0):
            raise ValueError('Newton solver requires finite positive friction coefficients')
        self.nc = len(self.mu)
        self.nv = self.J.shape[1]
        sc = np.ones(3 * self.nc)
        sc[0::3] = 1.0 / self.mu
        self.Jh = sp.diags(sc) @ self.J
        self.Jh = self.Jh.tocsr()
        self._G = None
        self._rho = None

    @classmethod
    def from_scene(cls, sc):
        """from a sparse scene with attributes J, M, Minv, v_free, b, mu."""
        return cls(sc.J, sc.M, sc.Minv, sc.v_free, sc.b - sc.J @ sc.v_free, sc.mu)

    @classmethod
    def from_engine_scene(cls, scene):
        """from motion_engine.ncp.ncp_ref.Scene (dense J, Minv full or diagonal, b_offset)."""
        J = sp.csr_matrix(np.asarray(scene.J, float))
        Mi = np.asarray(scene.Minv, float)
        Mi = np.diag(Mi) if Mi.ndim == 1 else Mi
        nv = J.shape[1]
        if nv % 6 or Mi.shape != (nv, nv):
            raise ValueError('Newton solver requires six degrees of freedom per rigid body')
        if not np.isfinite(Mi).all() or not np.allclose(Mi, Mi.T, rtol=0, atol=1e-14):
            raise ValueError('Inverse mass matrix must be finite and symmetric')
        block_ids = np.arange(nv) // 6
        if np.any(Mi[block_ids[:, None] != block_ids[None, :]] != 0):
            raise ValueError('Newton scene conversion requires independent rigid-body mass blocks')
        if not np.isfinite(J.data).all() or not np.isfinite(scene.v_free).all():
            raise ValueError('Scene data must be finite')
        if np.asarray(scene.v_free).shape != (nv,) or np.asarray(scene.mu).shape != (J.shape[0] // 3,) or J.shape[0] % 3:
            raise ValueError('Incompatible scene dimensions')
        (blocks_i, blocks) = ([], [])
        for k in range(nv // 6):
            B = Mi[6 * k:6 * k + 6, 6 * k:6 * k + 6]
            try:
                np.linalg.cholesky(B)
            except np.linalg.LinAlgError as ex:
                raise ValueError('Rigid-body inverse mass blocks must be positive definite') from ex
            blocks_i.append(B)
            blocks.append(np.linalg.inv(B))
        off = np.zeros(J.shape[0]) if scene.b_offset is None else np.asarray(scene.b_offset, float)
        if off.shape != (J.shape[0],) or not np.isfinite(off).all():
            raise ValueError('Contact offset must be finite and match contact rows')
        mass = sp.block_diag(blocks, format='csr') if blocks else sp.csr_matrix((nv, nv))
        inverse_mass = sp.block_diag(blocks_i, format='csr') if blocks_i else sp.csr_matrix((nv, nv))
        return cls(J, mass, inverse_mass, np.asarray(scene.v_free, float), off, np.asarray(scene.mu, float))

    @property
    def G(self):
        if self._G is None:
            self._G = (self.J @ self.Minv @ self.J.T).tocsr()
        return self._G

    @property
    def b(self):
        return self.J @ self.v_free + self.off

    def rho_nat(self):
        """rho_c = 1/||G_cc||_2; large inputs currently assemble sparse G."""
        if self._rho is None:
            J = self.J.tocoo()
            nc = self.nc
            Gcc = np.zeros((nc, 3, 3))
            Jd = self.J.toarray() if self.nv * 3 * nc < 4000000.0 else None
            if Jd is not None:
                Mi = self.Minv.toarray()
                for c in range(nc):
                    Jc = Jd[3 * c:3 * c + 3]
                    nz = np.where(np.any(Jc != 0, axis=0))[0]
                    Gcc[c] = Jc[:, nz] @ Mi[np.ix_(nz, nz)] @ Jc[:, nz].T
            else:
                G = self.G.tocoo()
                m = G.row // 3 == G.col // 3
                Gcc[G.row[m] // 3, G.row[m] % 3, G.col[m] % 3] = G.data[m]
            s = np.linalg.norm(Gcc, ord=2, axis=(1, 2))
            self._rho = np.where(s > 0, 1.0 / np.where(s > 0, s, 1.0), 1.0)
        return self._rho

    def lam_from_z(self, z):
        lam = z.copy()
        lam[:, 0] = z[:, 0] / self.mu
        return lam.reshape(-1)

    def natres(self, lam, v=None):
        """engine natural residual on lam (u from lam, not from the primal iterate)."""
        if v is None:
            v = self.v_free + self.Minv @ (self.J.T @ lam)
        u = (self.J @ v + self.off).reshape(self.nc, 3)
        L = lam.reshape(self.nc, 3)
        s = np.zeros_like(u)
        s[:, 0] = self.mu * np.linalg.norm(u[:, 1:], axis=1)
        X = L - self.rho_nat()[:, None] * (u + s)
        return float(np.max(np.abs(L - proj_K(X, self.mu)))) if self.nc else 0.0

def proj_K(X, mu):
    (xn, xt) = (X[:, 0], X[:, 1:])
    nt = np.linalg.norm(xt, axis=1)
    P = np.zeros_like(X)
    inside = nt <= mu * xn
    polar = mu * nt <= -xn
    mid = ~inside & ~polar
    P[inside] = X[inside]
    a = (mu[mid] * nt[mid] + xn[mid]) / (1 + mu[mid] ** 2)
    P[mid, 0] = a
    P[mid, 1:] = (mu[mid] * a / np.maximum(nt[mid], 1e-300))[:, None] * xt[mid]
    return P

def _blockdiag3(B):
    nc = len(B)
    rows = (3 * np.arange(nc))[:, None, None] + np.arange(3)[None, :, None] + np.zeros((1, 1, 3), int)
    cols = (3 * np.arange(nc))[:, None, None] + np.zeros((1, 3, 1), int) + np.arange(3)[None, None, :]
    return sp.csr_matrix((B.ravel(), (rows.ravel(), cols.ravel())), shape=(3 * nc, 3 * nc))

def _dsax_rows(pb, x, eps):
    """de Saxce term in scaled coordinates and its Jacobian.
    c_hat_n = (off_n + mu |u_t|_eps)/mu,  d c_hat_n / dv = that^T J_t  (that = u_t/|u_t|_eps).
    |u_t|_eps = sqrt(|u_t|^2 + eps^2) - eps  (smoothing, eps -> 0 with the barrier)."""
    nc = pb.nc
    u = (pb.J @ x + pb.off).reshape(nc, 3)
    nt = np.sqrt(np.einsum('ij,ij->i', u[:, 1:], u[:, 1:]) + eps * eps)
    svec = pb.mu * (nt - eps)
    that = u[:, 1:] / nt[:, None]
    rows = np.repeat(3 * np.arange(nc), 2)
    cols = (3 * np.arange(nc))[:, None] + np.array([1, 2])[None]
    Dt = sp.csr_matrix((that.ravel(), (rows, cols.ravel())), shape=(3 * nc, 3 * nc))
    return (svec, Dt)

def solve(pb: Problem, tol=1e-10, max_iter=200, desaxce='newton', linsolver='auto', verbose=False, eps_fac=1.0, anitescu_phase=True, time_limit=None, hook=None, s_fixed=None, conv_tol=None, init=None):
    """Primal-dual IPM (Mehrotra, NT scaling) in body space for the de Saxce NCP.

    desaxce = "newton": the de Saxce term s(v) = mu |u_t(v)|_eps is part of the Newton system
              (constraint Jacobian Jh + Dt J, reduced matrix M + Jh^T W^-2 (Jh + Dt J),
              nonsymmetric 6 n_b x 6 n_b), smoothing eps ~ sqrt(barrier) -> 0;
              "fixed":  s(v) refreshed after each step but not differentiated (symmetric H,
              linear rate of the Cadoux fixed point);
              "none":   Anitescu convex relaxation (s = 0) -- NOT the NCP, for reference.
    anitescu_phase: start with s = 0 until the relative gap < 1e-2 (global phase on the convex
              problem), then switch the de Saxce coupling on.
    Stop: engine natural residual (lam from the dual iterate) < tol."""
    t0 = time.perf_counter()
    (nc, nv) = (pb.nc, pb.nv)
    mu = pb.mu
    Mvf = pb.M @ pb.v_free
    sym = desaxce != 'newton'
    lin = BodySolver(nv, linsolver, symmetric=True)
    lin_ns = BodySolver(nv, linsolver, symmetric=False)
    e = np.zeros((nc, 3))
    e[:, 0] = 1.0
    zero_Dt = sp.csr_matrix((3 * nc, 3 * nc))

    def chat(svec):
        c = pb.off.copy()
        c[0::3] += svec
        c[0::3] /= mu
        return c.reshape(nc, 3)
    ch = chat(np.zeros(nc))
    lin.factor((pb.M + pb.Jh.T @ pb.Jh).tocsr())
    x = lin.solve(Mvf - pb.Jh.T @ ch.reshape(-1))
    s = (pb.Jh @ x).reshape(nc, 3) + ch

    def shift_in(w):
        a = np.max(np.linalg.norm(w[:, 1:], axis=1) - w[:, 0])
        if a >= 0:
            w = w + (1.0 + a) * e
        return w
    s = shift_in(s)
    z = s.copy()
    hist = []
    coupled = desaxce != 'none' and (not anitescu_phase)
    zscale = None
    it = 0
    resn = np.inf
    t_natres = 0.0
    n_ns = 0
    status = 'max_iter'
    for it in range(1, max_iter + 1):
        gap = float(np.sum(s * z))
        mug = gap / nc
        if zscale is None:
            zscale = mug
        if coupled:
            eps = eps_fac * mug / max(float(np.mean(z[:, 0])), 1e-300)
            (svec, Dt) = _dsax_rows(pb, x, eps)
            ch = chat(svec)
            if desaxce != 'newton':
                Dt = zero_Dt
        else:
            ch = chat(np.zeros(nc) if s_fixed is None else s_fixed)
            Dt = zero_Dt
        rx = pb.M @ x - Mvf - pb.Jh.T @ z.reshape(-1)
        rz = -(pb.Jh @ x).reshape(nc, 3) + s - ch
        (W, Winv) = nt_scaling(s, z)
        lam = bmv(W, z)
        D = _blockdiag3(np.einsum('ijk,ikl->ijl', Winv, Winv))
        Jeff = pb.Jh + Dt @ pb.J if coupled and desaxce == 'newton' else pb.Jh
        H = (pb.M + pb.Jh.T @ D @ Jeff).tocsr()
        L = lin_ns if coupled and desaxce == 'newton' else lin
        try:
            L.factor(H)
        except Exception as ex:
            status = f'factor_failed: {ex}'
            break
        if L is lin_ns:
            n_ns += 1

        def newton(ds):
            ldd = jsolve(lam, ds)
            t1 = rz + bmv(W, ldd)
            Wi2t = bmv(Winv, bmv(Winv, t1))
            dx = L.solve(-rx + pb.Jh.T @ Wi2t.reshape(-1))
            dz = bmv(Winv, bmv(Winv, -(Jeff @ dx).reshape(nc, 3) + t1))
            dsv = bmv(W, ldd - bmv(W, dz))
            return (dx, dsv, dz)
        (dx_a, ds_a, dz_a) = newton(-jprod(lam, lam))
        a_a = min(1.0, max_step(s, ds_a), max_step(z, dz_a))
        gap_a = float(np.sum((s + a_a * ds_a) * (z + a_a * dz_a)))
        sigma = min(1.0, max(0.0, gap_a / gap)) ** 3
        corr = jprod(bmv(Winv, ds_a), bmv(W, dz_a))
        (dx, dsv, dz) = newton(-jprod(lam, lam) - corr + sigma * mug * e)
        alpha = min(1.0, 0.99 * min(max_step(s, dsv), max_step(z, dz)))
        x = x + alpha * dx
        s = s + alpha * dsv
        z = z + alpha * dz
        tt = time.perf_counter()
        lamv = pb.lam_from_z(z)
        resn = pb.natres(lamv)
        t_natres += time.perf_counter() - tt
        rec = dict(it=it, alpha=alpha, sigma=sigma, gap=gap, rx=float(np.linalg.norm(rx)), rz=float(np.linalg.norm(rz)), natres=resn, coupled=coupled)
        hist.append(rec)
        if verbose:
            print(rec, flush=True)
        if resn < tol:
            status = 'converged'
            break
        if conv_tol is not None and gap / nc < conv_tol and (float(np.linalg.norm(rz)) < conv_tol * 1000.0):
            status = 'convex_converged'
            break
        if hook is not None and hook(it, z, x, resn):
            status = 'hook_done'
            break
        if not coupled and desaxce != 'none' and (gap / nc < 0.01 * zscale):
            coupled = True
        if time_limit is not None and time.perf_counter() - t0 > time_limit:
            status = 'time_limit'
            break
    lin.free()
    lin_ns.free()
    lamv = pb.lam_from_z(z)
    v = pb.v_free + pb.Minv @ (pb.J.T @ lamv)
    T = time.perf_counter() - t0
    fin = pb.natres(lamv)
    return dict(lam=lamv, v=v, x=x, s=s, z=z, natres=fin, iters=it, hist=hist, time=T, t_fact=lin.t_fact + lin_ns.t_fact, t_solve=lin.t_solve + lin_ns.t_solve, n_fact=lin.n_fact + lin_ns.n_fact, n_nonsym=n_ns, t_natres=t_natres, converged=bool(fin < tol), status=status, linsolver=lin.kind)

def _branch_maps(pb, lab, d):
    """U = E J (rows), Vt = T^T J (rows) so V = J^T T = Vt^T; ny; index maps."""
    nc = pb.nc
    S = np.where(lab == 1)[0]
    P = np.where(lab == 2)[0]
    erow = np.concatenate([(3 * S[:, None] + np.arange(3)).ravel(), 3 * P]).astype(int)
    ny = len(erow)
    E = sp.csr_matrix((np.ones(ny), (np.arange(ny), erow)), shape=(ny, 3 * nc))
    r_s = np.repeat(np.arange(3 * len(S)), 1)
    c_s = (3 * S[:, None] + np.arange(3)).ravel()
    v_s = np.ones(len(c_s))
    k0 = 3 * len(S)
    r_p = np.repeat(k0 + np.arange(len(P)), 3)
    c_p = (3 * P[:, None] + np.arange(3)).ravel()
    v_p = np.stack([np.ones(len(P)), -pb.mu[P] * d[:, 0], -pb.mu[P] * d[:, 1]], 1).ravel()
    Tt = sp.csr_matrix((np.concatenate([v_s, v_p]), (np.concatenate([r_s, r_p]), np.concatenate([c_s, c_p]))), shape=(ny, 3 * nc))
    return (S, P, E, Tt, ny)

def classify_natmap(pb, lam):
    """labels from the natural map X = lam - rho (u + Gamma(u)):  X in int K -> stick,
    X in polar cone -> open, else slip with direction d = -X_t/|X_t| (Alart-Curnier active set)."""
    nc = pb.nc
    L = lam.reshape(nc, 3)
    v = pb.v_free + pb.Minv @ (pb.J.T @ lam)
    u = (pb.J @ v + pb.off).reshape(nc, 3)
    g = np.zeros_like(u)
    g[:, 0] = pb.mu * np.linalg.norm(u[:, 1:], axis=1)
    X = L - pb.rho_nat()[:, None] * (u + g)
    (xn, xt) = (X[:, 0], X[:, 1:])
    nt = np.linalg.norm(xt, axis=1)
    lab = np.full(nc, 2, int)
    lab[nt <= pb.mu * xn] = 1
    lab[pb.mu * nt <= -xn] = 0
    d = -xt / np.maximum(nt, 1e-300)[:, None]
    return (lab, d, u)

def polish(pb, lam0, v0=None, tol=1e-10, eps_rel=1e-09, max_rounds=8, max_refine=6, linsolver='auto', log=None):
    """Crossover from an IPM iterate (semismooth active-set Newton with proximal regularisation).
    Each round: labels and slip directions from the natural map at the current lam; branch
    equations  E u(lam) = 0,  lam = T y  solved by proximal refinement
        y <- y + (eps I + A)^-1 (-E u),   A = E J M^-1 J^T T,
    every prox solve in BODY space by Woodbury with ONE factor of  eps M + J^T T E J  (the
    vertex pencil restricted to the active set).  On hyperstatic sets A is singular; the
    prox iteration leaves ker A untouched, so lam stays at the IPM's analytic centre instead
    of jumping to an arbitrary vertex of the force fibre.  Globalised by monotonicity: a round
    is accepted only if the engine natural residual decreases; otherwise the polish returns
    'fail' and the caller continues the interior-point phase."""
    nc = pb.nc
    t0 = time.perf_counter()
    lam = lam0.copy()
    best = pb.natres(lam)
    out = dict(rounds=[], status='fail')
    nfact = 0
    gd = float(np.median(np.asarray(pb.J.multiply(pb.J) @ pb.Minv.diagonal())))
    eps = eps_rel * gd
    for rnd in range(max_rounds):
        (lab, dall, u) = classify_natmap(pb, lam)
        L = lam.reshape(nc, 3)
        P = np.where(lab == 2)[0]
        (S, P, E, Tt, ny) = _branch_maps(pb, lab, dall[P])
        y = np.concatenate([L[S].ravel(), np.maximum(L[P, 0], 0.0)])
        U = (E @ pb.J).tocsr()
        Vt = (Tt @ pb.J).tocsr()
        K = (eps * pb.M + Vt.T @ U).tocsr()
        sym = len(P) == 0
        ls = BodySolver(pb.nv, linsolver, symmetric=sym)
        try:
            ls.factor(K)
        except Exception as ex:
            out['status'] = f'factor_failed: {ex}'
            break
        nfact += 1
        res_hist = []
        lam_try = Tt.T @ y
        for k in range(max_refine):
            vv = pb.v_free + pb.Minv @ (pb.J.T @ lam_try)
            r = -(U @ vv + E @ pb.off)
            dl = (r - U @ ls.solve(Vt.T @ r)) / eps
            y = y + dl
            lam_try = Tt.T @ y
            rn = pb.natres(lam_try)
            res_hist.append(rn)
            if rn < tol or (k > 0 and rn > 0.5 * res_hist[-2]):
                break
        ls.free()
        rn = pb.natres(lam_try)
        out['rounds'].append(dict(round=rnd, n_stick=int(len(S)), n_slip=int(len(P)), n_open=int((lab == 0).sum()), res=res_hist, eps=eps, accepted=bool(rn < best)))
        if log:
            log(out['rounds'][-1])
        if rn < best:
            (lam, best) = (lam_try, rn)
        else:
            break
        if best < tol:
            out['status'] = 'converged'
            break
    out['lam'] = lam
    out['v'] = pb.v_free + pb.Minv @ (pb.J.T @ lam)
    out['natres'] = best
    out['labels'] = classify_natmap(pb, lam)[0]
    out['time'] = time.perf_counter() - t0
    out['n_fact'] = nfact
    return out

def _natmap_parts(pb, lam):
    nc = pb.nc
    L = lam.reshape(nc, 3)
    v = pb.v_free + pb.Minv @ (pb.J.T @ lam)
    u = (pb.J @ v + pb.off).reshape(nc, 3)
    ut = u[:, 1:]
    nut = np.linalg.norm(ut, axis=1)
    g = np.zeros_like(u)
    g[:, 0] = pb.mu * nut
    rho = pb.rho_nat()
    X = L - rho[:, None] * (u + g)
    P = proj_K(X, pb.mu)
    F = L - P
    return (F, X, u, nut, v)

def _dproj(X, mu):
    """Clarke-Jacobian element of P_K at X, (m,3,3)."""
    m = len(X)
    (xn, xt) = (X[:, 0], X[:, 1:])
    r = np.linalg.norm(xt, axis=1)
    D = np.zeros((m, 3, 3))
    inside = r <= mu * xn
    polar = mu * r <= -xn
    mid = ~inside & ~polar
    D[inside] = np.eye(3)
    if np.any(mid):
        (mm, rr) = (mu[mid], r[mid])
        xb = xt[mid] / rr[:, None]
        a = (mm * rr + xn[mid]) / (1 + mm ** 2)
        k = 1.0 / (1 + mm ** 2)
        Dm = np.zeros((mid.sum(), 3, 3))
        Dm[:, 0, 0] = k
        Dm[:, 0, 1:] = (mm * k)[:, None] * xb
        Dm[:, 1:, 0] = (mm * k)[:, None] * xb
        xx = np.einsum('ij,ik->ijk', xb, xb)
        Dm[:, 1:, 1:] = (mm ** 2 * k)[:, None, None] * xx + (mm * a / rr)[:, None, None] * (np.eye(2)[None] - xx)
        D[mid] = Dm
    return (D, inside, polar, mid)

def _jac_blocks(pb, X, u, nut, dir_from_X='small', stick_g0=True, weak_open=None):
    nc = pb.nc
    mu = pb.mu
    rho = pb.rho_nat()
    (Dp, inside, polar, mid) = _dproj(X, mu)
    if weak_open is not None and np.any(weak_open):
        Dp[weak_open] = 0.0
        mid = mid & ~weak_open
        polar = polar | weak_open
    DG = np.zeros((nc, 3, 3))
    tt = np.where(nut[:, None] > 1e-300, u[:, 1:] / np.maximum(nut, 1e-300)[:, None], 0.0)
    DG[:, 0, 1:] = mu[:, None] * tt
    if stick_g0:
        DG[inside] = 0.0
    sel = mid & (rho * nut < 0.001 * np.linalg.norm(X[:, 1:], axis=1))
    if dir_from_X == 'always':
        sel = mid
    if dir_from_X and np.any(sel):
        xt = X[sel, 1:]
        DG[sel, 0, 1:] = mu[sel, None] * (-xt / np.maximum(np.linalg.norm(xt, axis=1), 1e-300)[:, None])
    C = np.einsum('ijk,ikl->ijl', Dp, np.eye(3)[None] + DG) * rho[:, None, None]
    return (Dp, C, (inside, polar, mid))

def ssn(pb, lam0, tol=1e-10, max_iter=60, linsolver='auto', eps0=0.001, verbose=False, time_limit=None, stick_g0=True, polish_on_fail=True, dir_from_X='small', nonmono=1, retries=2, polish_below=None, weak_frac=0.05):
    """Semismooth Newton on F(lam) = lam - P_K(lam - rho (u + Gamma(u))), u = G lam + b.
    Generalised Jacobian  dF = (I - D_P) + D_P rho (I + D_Gamma) G  (block-diagonal parts
    times G = J M^-1 J^T), Levenberg-Marquardt shift eps_k = min(eps0, |F|_inf) for the
    singular hyperstatic directions, and the step solved in BODY space by Woodbury:
        (Bd + C J M^-1 J^T)^-1 = Bd^-1 - Bd^-1 C J (M + J^T Bd^-1 C J)^-1 J^T Bd^-1,
        Bd = (1+eps) I - D_P,  C = D_P rho (I + D_Gamma).
    Globalisation: Armijo backtracking on 1/2 |F|_2^2 (optionally non-monotone); when it
    fails, the Jacobian element is re-chosen at the full-step trial point (branch switch
    across a kink), up to `retries` times; then the proximal active-set polish is tried."""
    t0 = time.perf_counter()
    nc = pb.nc
    lam = lam0.copy()
    last_dir = None
    (F, X, u, nut, v) = _natmap_parts(pb, lam)
    phi = 0.5 * float(F.ravel() @ F.ravel())
    hist = []
    lin = BodySolver(pb.nv, linsolver, symmetric=False)
    status = 'max_iter'
    it = 0
    phis = []
    n_polish = 0
    for it in range(1, max_iter + 1):
        res = float(np.max(np.abs(F)))
        phis.append(phi)
        if res < tol:
            status = 'converged'
            it -= 1
            break
        if polish_below is not None and res < polish_below:
            po = polish(pb, lam, tol=tol, linsolver=linsolver)
            n_polish += 1
            hist.append(dict(it=it, polish=po['status'], polish_natres=po['natres'], polish_rounds=len(po['rounds'])))
            if po['natres'] < tol:
                lam = po['lam']
                (F, X, u, nut, v) = _natmap_parts(pb, lam)
                break
        eps = min(eps0, max(res, 1e-14))
        (Xj, uj, nutj) = (X, u, nut)
        acc = False
        Lm = lam.reshape(nc, 3)
        act = Lm[:, 0] > 0
        lnmed = float(np.median(Lm[act, 0])) if np.any(act) else 0.0
        (_, _, _, mid0) = _dproj(X, pb.mu)
        weak = mid0 & (Lm[:, 0] < weak_frac * lnmed)
        plan = ['std'] + (['weak'] if np.any(weak) and weak_frac > 0 else []) + ['trial'] * retries
        plan = plan[:1 + max(retries, 0) + (1 if np.any(weak) and weak_frac > 0 else 0)]
        for (attempt, kind) in enumerate(plan):
            (Dp, C, (inside, polar, mid)) = _jac_blocks(pb, Xj, uj, nutj, dir_from_X, stick_g0, weak_open=weak if kind == 'weak' else None)
            Bd = (1.0 + eps) * np.eye(3)[None] - Dp
            Bi = np.linalg.inv(Bd)
            BiC = np.einsum('ijk,ikl->ijl', Bi, C)
            H = (pb.M + pb.J.T @ _blockdiag3(BiC) @ pb.J).tocsr()
            try:
                lin.factor(H)
            except Exception as ex:
                status = f'factor_failed: {ex}'
                break
            y = bmv(Bi, -F).reshape(-1)
            w = lin.solve(pb.J.T @ y)
            dlam = y - bmv(BiC, (pb.J @ w).reshape(nc, 3)).reshape(-1)
            alpha = 1.0
            ref = max([phi] + phis[-nonmono:])
            full = None
            for ls in range(30):
                lt = lam + alpha * dlam
                (F2, X2, u2, nut2, v2) = _natmap_parts(pb, lt)
                if full is None:
                    full = (X2, u2, nut2)
                phi2 = 0.5 * float(F2.ravel() @ F2.ravel())
                if phi2 <= (1 - 0.0001 * alpha) * ref:
                    acc = True
                    break
                alpha *= 0.5
                if alpha < 1e-06:
                    break
            if acc:
                break
            if kind != 'weak':
                (Xj, uj, nutj) = full
        last_dir = dlam
        hist.append(dict(it=it, res=res, eps=eps, alpha=alpha, n_stick=int(inside.sum()), n_slip=int(mid.sum()), n_open=int(polar.sum()), ls=ls, attempt=attempt, kind=kind))
        if verbose:
            print(hist[-1], flush=True)
        if status.startswith('factor_failed'):
            break
        if not acc:
            status = 'line_search_failed'
            if polish_on_fail:
                po = polish(pb, lam, tol=tol, linsolver=linsolver)
                n_polish += 1
                hist.append(dict(it=it, polish=po['status'], polish_natres=po['natres'], polish_rounds=len(po['rounds'])))
                if po['natres'] < res:
                    lam = po['lam']
                    (F, X, u, nut, v) = _natmap_parts(pb, lam)
                    phi = 0.5 * float(F.ravel() @ F.ravel())
                    status = 'polished'
                    if po['natres'] < tol:
                        break
                    continue
            break
        (lam, F, X, u, nut, v, phi) = (lt, F2, X2, u2, nut2, v2, phi2)
        if time_limit is not None and time.perf_counter() - t0 > time_limit:
            status = 'time_limit'
            break
    lin.free()
    res = pb.natres(lam)
    if res < tol:
        status = 'converged'
    return dict(lam=lam, v=pb.v_free + pb.Minv @ (pb.J.T @ lam), natres=res, iters=it, status=status, hist=hist, time=time.perf_counter() - t0, n_fact=lin.n_fact, t_fact=lin.t_fact, n_polish=n_polish, last_dir=last_dir)

def subproblem(pb, keep):
    rows = (3 * keep[:, None] + np.arange(3)).ravel()
    p2 = Problem(pb.J[rows], pb.M, pb.Minv, pb.v_free, pb.off[rows], pb.mu[keep])
    p2._rho = pb.rho_nat()[keep]
    return (p2, rows)

def ssn_escape(pb, lam0, tol=1e-10, max_iter=40, n_cand=6, max_escapes=8, log=None, time_limit=None, **kw):
    """Semismooth Newton with fold escape.  When Newton stalls, its last direction is dominated
    by the near-null mode of the active-set Jacobian (a fold of the stick-slip branch, typically
    a degenerate slip contact whose slip element is rank deficient).  The contacts carrying the
    mode are candidates; each is put on the take-off branch (lambda_c = 0, contact removed) in
    order of increasing lambda_n and Newton is rerun warm on the reduced problem.
      - reduced converges and the FULL natural residual < tol            -> done;
      - reduced converges but the removed contact wants to close (u_n<0) -> candidate rejected;
      - reduced stalls at a lower residual                               -> keep, look again
        (nested folds).
    Removed contacts that later violate u_n >= 0 are put back.  Every accepted answer is checked
    on the FULL problem; nothing is accepted on the reduced residual alone."""
    t0 = time.perf_counter()
    nc = pb.nc
    active = np.ones(nc, bool)
    lam = lam0.copy()
    escapes = []
    tot = dict(it=0, f=0)

    def run(mask, lam_full):
        idx = np.where(mask)[0]
        (p2, rows) = subproblem(pb, idx) if len(idx) < nc else (pb, np.arange(3 * nc))
        r = ssn(p2, lam_full[rows], tol=tol, max_iter=max_iter, **kw)
        tot['it'] += r['iters']
        tot['f'] += r['n_fact']
        lf = np.zeros(3 * nc)
        lf[rows] = r['lam']
        return (r, lf, idx)
    (r, lam, idx) = run(active, lam)
    for k in range(max_escapes):
        full = pb.natres(lam)
        if full < tol:
            break
        if time_limit is not None and time.perf_counter() - t0 > time_limit:
            break
        if r['natres'] < tol:
            u = (pb.J @ (pb.v_free + pb.Minv @ (pb.J.T @ lam)) + pb.off).reshape(nc, 3)
            back = ~active & (u[:, 0] < 0)
            if not np.any(back):
                break
            active = active | back
            escapes.append(dict(kind='readd', contacts=np.where(back)[0].tolist()))
            (r, lam, idx) = run(active, lam)
            continue
        if r.get('last_dir') is None:
            break
        d = r['last_dir'].reshape(-1, 3)
        Lr = r['lam'].reshape(-1, 3)
        top = np.argsort(-np.abs(d).max(1))[:n_cand]
        top = top[np.argsort(Lr[top, 0])]
        best = None
        for ci in top:
            if Lr[ci, 0] <= 0:
                continue
            c = int(idx[ci])
            m2 = active.copy()
            m2[c] = False
            (r3, lam3, idx3) = run(m2, lam)
            full3 = pb.natres(lam3)
            rec = dict(kind='drop', contact=c, lam_n=float(Lr[ci, 0]), sub_status=r3['status'], sub_natres=r3['natres'], natres_full=full3, iters=r3['iters'])
            escapes.append(rec)
            if log:
                log(rec)
            if full3 < tol:
                best = (m2, r3, lam3, idx3)
                break
            if r3['natres'] < tol:
                continue
            if r3['natres'] < 0.5 * r['natres'] and (best is None or r3['natres'] < best[1]['natres']):
                best = (m2, r3, lam3, idx3)
        if best is None:
            break
        (active, r, lam, idx) = best
    full = pb.natres(lam)
    return dict(lam=lam, v=pb.v_free + pb.Minv @ (pb.J.T @ lam), natres=full, converged=bool(full < tol), iters=tot['it'], n_fact=tot['f'], escapes=escapes, time=time.perf_counter() - t0, removed=np.where(~active)[0].tolist(), status='converged' if full < tol else 'not_converged')

def with_mu(pb, m):
    p = Problem(pb.J, pb.M, pb.Minv, pb.v_free, pb.off, np.full(pb.nc, m) if np.isscalar(m) else m)
    p._rho = pb.rho_nat()
    return p

def ssn_continuation(pb, lam0=None, tol=1e-10, mu0=0.1, dmu0=0.1, dmu_min=0.005, max_iter=25, escape=True, time_limit=None, log=None, **kw):
    """Continuation in the friction coefficient: mu_k = theta_k mu_target, theta from mu0/mu_max to 1,
    each step a warm semismooth Newton solve; step halved on failure.  At the minimum step a
    failure means the followed stick-slip branch ends in a fold: then the step is taken anyway
    from the last good solution (branch jump, as the true solution set is followed by its other
    branch) and, at the target, the fold escape is used.  Returns path log incl. fold events."""
    t0 = time.perf_counter()
    mu_t = pb.mu
    mmax = float(np.max(mu_t))
    th = min(1.0, mu0 / mmax)
    dth = dmu0 / mmax
    dmin = dmu_min / mmax
    lam = np.zeros(3 * pb.nc) if lam0 is None else lam0.copy()
    path = []
    tot = dict(it=0, f=0)
    folds = []
    while True:
        p = with_mu(pb, th * mu_t)
        if th >= 1.0 and escape:
            r = ssn_escape(p, lam, tol=tol, max_iter=max_iter, **kw)
        else:
            r = ssn(p, lam, tol=tol, max_iter=max_iter, **kw)
        tot['it'] += r['iters']
        tot['f'] += r['n_fact']
        ok = r['natres'] < tol
        path.append(dict(theta=th, mu=th * mmax, ok=ok, iters=r['iters'], natres=r['natres'], status=r['status']))
        if log:
            log(path[-1])
        if ok:
            lam = r['lam']
            if th >= 1.0:
                break
            th_new = min(1.0, th + dth)
            dth = min(dth * 1.5, dmu0 / mmax)
            th = th_new
            continue
        if th >= 1.0 and dth <= dmin:
            break
        if time_limit is not None and time.perf_counter() - t0 > time_limit:
            break
        th_prev = th - dth
        if dth > dmin:
            dth = max(dth / 2.0, dmin)
            th = min(1.0, th_prev + dth)
        else:
            folds.append(th * mmax)
            th = min(1.0, th + dth)
    full = pb.natres(lam)
    return dict(lam=lam, v=pb.v_free + pb.Minv @ (pb.J.T @ lam), natres=full, converged=bool(full < tol), iters=tot['it'], n_fact=tot['f'], path=path, folds=folds, time=time.perf_counter() - t0, status='converged' if full < tol else 'not_converged')

def solve_newton(pb, tol=1e-10, time_limit=200, ipm_iters=30, log=None):
    """The pipeline: (a) body-space IPM (global phase, de Saxce in the Newton system) until
    natres < 1e-4, then semismooth Newton with fold escape; (b) if not converged: semismooth
    Newton with fold escape from lam = 0; (c) if not converged: friction continuation.  Every
    stage is timed and counted; the result says which stage produced the answer."""
    t0 = time.perf_counter()
    stages = []
    lam = np.zeros(3 * pb.nc)
    best = (np.inf, lam)
    r = solve(pb, tol=0.0001, max_iter=ipm_iters, desaxce='newton', eps_fac=0.01, time_limit=time_limit)
    stages.append(dict(stage='ipm', iters=r['iters'], n_fact=r['n_fact'], natres=r['natres'], time=r['time']))
    if r['natres'] < 0.001:
        e = ssn_escape(pb, r['lam'], tol=tol, max_iter=40, eps0=0.01, dir_from_X='', retries=2, time_limit=time_limit)
        stages.append(dict(stage='ssn_after_ipm', iters=e['iters'], n_fact=e['n_fact'], natres=e['natres'], time=e['time'], escapes=len(e['escapes']), removed=e['removed']))
        best = min(best, (e['natres'], e['lam']), key=lambda q: q[0])
    if best[0] >= tol and time.perf_counter() - t0 < time_limit:
        e = ssn_escape(pb, lam, tol=tol, max_iter=40, eps0=0.01, dir_from_X='', retries=2, time_limit=time_limit)
        stages.append(dict(stage='ssn_cold', iters=e['iters'], n_fact=e['n_fact'], natres=e['natres'], time=e['time'], escapes=len(e['escapes']), removed=e['removed']))
        best = min(best, (e['natres'], e['lam']), key=lambda q: q[0])
    if best[0] >= tol and time.perf_counter() - t0 < time_limit:
        c = ssn_continuation(pb, tol=tol, eps0=0.01, dir_from_X='', retries=2, time_limit=max(1.0, time_limit - (time.perf_counter() - t0)))
        stages.append(dict(stage='mu_continuation', iters=c['iters'], n_fact=c['n_fact'], natres=c['natres'], time=c['time'], folds=c['folds']))
        best = min(best, (c['natres'], c['lam']), key=lambda q: q[0])
    lamb = best[1]
    full = pb.natres(lamb)
    return dict(lam=lamb, v=pb.v_free + pb.Minv @ (pb.J.T @ lamb), natres=full, converged=bool(full < tol), stages=stages, iters=int(sum((s['iters'] for s in stages))), n_fact=int(sum((s['n_fact'] for s in stages))), time=time.perf_counter() - t0, status='converged' if full < tol else 'not_converged', solved_by=stages[-1]['stage'] if full < tol else None)

def labels_at(pb, lam, tol_rel=1e-07):
    nc = pb.nc
    L = lam.reshape(nc, 3)
    v = pb.v_free + pb.Minv @ (pb.J.T @ lam)
    u = (pb.J @ v + pb.off).reshape(nc, 3)
    scale = max(float(np.abs(L[:, 0]).max()), 1e-300)
    ln = L[:, 0]
    lt = np.linalg.norm(L[:, 1:], axis=1)
    ut = np.linalg.norm(u[:, 1:], axis=1)
    lab = np.full(nc, 'stick', dtype='<U5')
    lab[ln <= tol_rel * scale] = 'open'
    lab[(lab != 'open') & (lt >= pb.mu * ln * (1 - 1e-06))] = 'slip'
    return (lab, u, scale, ln, lt, ut)

def uniqueness_report(pb, lam, dense_max=4000, null_tol=1e-09, near_tol=1e-06, tol=1e-10):
    """Numerical local diagnostics, conditional on a finite residual below tol.

    A null vector of the nonlinear linearization is first-order evidence only.
    The dimension lower bound is structural; SVD ranks are not exact certificates.
    """
    _check_tol(tol)
    lam = np.asarray(lam, float)
    if lam.shape != (3 * pb.nc,) or not np.isfinite(lam).all():
        return dict(verdict='INVALID_SOLUTION', certificate=False, motion_checked=False)
    residual = pb.natres(lam)
    if not np.isfinite(residual) or residual >= tol:
        return dict(verdict='NOT_CONVERGED_INCONCLUSIVE', natres=residual, certificate=False, motion_checked=False)
    nc = pb.nc
    if not nc:
        return dict(verdict='EMPTY_CONTACT_SET', natres=0.0, certificate=False, force_null_dim=0, force_null_exact=False, motion_checked=False)
    (lab, u, scale, ln, lt, ut) = labels_at(pb, lam)
    uscale = max(float(np.abs(u).max()), float(np.abs(pb.b).max()), 1e-300)
    S = np.where(lab == 'stick')[0]
    P = np.where(lab == 'slip')[0]
    O = np.where(lab == 'open')[0]
    rep = dict(n_stick=int(len(S)), n_slip=int(len(P)), n_open=int(len(O)), natres=residual, certificate=False, evidence='floating_point_linearization')
    deg = []
    for c in S:
        if abs(lt[c] - pb.mu[c] * ln[c]) < 1e-07 * pb.mu[c] * ln[c]:
            deg.append((int(c), 'stick_on_cone_boundary'))
        if ln[c] < 1e-09 * scale:
            deg.append((int(c), 'stick_near_zero_normal'))
    for c in P:
        if ut[c] < 1e-10 * uscale:
            deg.append((int(c), 'slip_zero_speed'))
        if ln[c] < 1e-09 * scale:
            deg.append((int(c), 'slip_near_zero_normal'))
    for c in O:
        if abs(u[c, 0]) < 1e-10 * uscale:
            deg.append((int(c), 'open_zero_gap_velocity'))
    rep['degenerate'] = deg
    d = np.zeros((len(P), 2))
    if len(P):
        d = u[P, 1:] / np.maximum(ut[P], 1e-300)[:, None]
    cols = []
    J = pb.J.tocsr()
    ny = 3 * len(S) + len(P)
    rep['ny'] = int(ny)
    touched = set()
    Jc = J.tocoo()
    act = np.zeros(nc, bool)
    act[S] = True
    act[P] = True
    m = act[Jc.row // 3]
    touched = np.unique(Jc.col[m] // 6)
    rep['force_null_lower_bound'] = int(max(0, ny - 6 * len(touched)))
    if ny and max(ny, pb.nv, 3 * nc) <= dense_max:
        Jd = J.toarray()
        F = np.zeros((pb.nv, ny))
        k = 0
        for c in S:
            F[:, k:k + 3] = Jd[3 * c:3 * c + 3].T
            k += 3
        for (i, c) in enumerate(P):
            F[:, k] = Jd[3 * c] - pb.mu[c] * (d[i, 0] * Jd[3 * c + 1] + d[i, 1] * Jd[3 * c + 2])
            k += 1
        sv = np.linalg.svd(F, compute_uv=False)
        r = int(np.sum(sv > null_tol * sv[0])) if len(sv) else 0
        rep['force_null_dim'] = int(ny - r)
        rep['force_null_exact'] = False
        rep['force_rank_computed'] = True
    else:
        rep['force_null_dim'] = rep['force_null_lower_bound']
        rep['force_null_exact'] = False
    if 3 * nc <= dense_max:
        G = pb.G.toarray()
        n = 3 * nc
        A = np.zeros((n, n))
        for c in O:
            A[3 * c:3 * c + 3, 3 * c:3 * c + 3] = np.eye(3)
        for c in S:
            A[3 * c:3 * c + 3, :] = G[3 * c:3 * c + 3, :]
        for (i, c) in enumerate(P):
            A[3 * c, :] = G[3 * c, :]
            sh = d[i]
            ns = max(ut[c], 1e-300)
            Pm = (np.eye(2) - np.outer(sh, sh)) / ns
            A[3 * c + 1:3 * c + 3, :] = pb.mu[c] * ln[c] * (Pm @ G[3 * c + 1:3 * c + 3, :])
            A[3 * c + 1:3 * c + 3, 3 * c] += pb.mu[c] * sh
            A[3 * c + 1, 3 * c + 1] += 1.0
            A[3 * c + 2, 3 * c + 2] += 1.0
        rs = np.maximum(np.abs(A).max(1), 1e-300)
        As = A / rs[:, None]
        (U_, sv, Vt) = np.linalg.svd(As)
        smax = sv[0]
        null = sv < null_tol * smax
        near = (sv < near_tol * smax) & ~null
        JT = J.T.toarray() if pb.nv * n < 40000000.0 else None
        Jn = np.linalg.norm(JT, 2) if JT is not None else 1.0

        def motion_frac(idx):
            if not np.any(idx) or JT is None:
                return 0.0
            Z = Vt[idx].T
            return float(np.linalg.norm(JT @ Z, 2) / Jn)
        rep['A_sigma_min_rel'] = float(sv[-1] / smax)
        rep['A_null_dim'] = int(null.sum())
        rep['A_near_null_dim'] = int(near.sum())
        rep['motion_frac_null'] = motion_frac(null)
        rep['motion_frac_near_null'] = motion_frac(near)
        motion_null = rep['motion_frac_null'] > 1e-06
        motion_near = rep['motion_frac_near_null'] > 1e-06
        rep['motion_checked'] = JT is not None
    else:
        motion_null = motion_near = False
        rep['motion_checked'] = False
    if any((k == 'slip_zero_speed' for (_, k) in deg)):
        v = 'DEGENERATE_INCONCLUSIVE'
    elif rep.get('motion_checked') and motion_null:
        v = 'MOTION_LINEARIZATION_SINGULAR'
    elif deg:
        v = 'DEGENERATE_INCONCLUSIVE'
    elif rep.get('motion_checked') and motion_near:
        v = 'MOTION_ILL_CONDITIONED'
    elif rep['force_null_dim'] > 0:
        v = 'FORCE_AMBIGUOUS_MOTION_REGULAR' if rep.get('motion_checked') else 'FORCE_AMBIGUOUS_MOTION_UNCHECKED'
    elif ny and (not rep.get('force_rank_computed')):
        v = 'RANK_UNCHECKED_INCONCLUSIVE'
    else:
        v = 'LOCALLY_REGULAR' if rep.get('motion_checked') else 'FORCE_RANK_FULL_MOTION_UNCHECKED'
    rep['verdict'] = v
    return rep

def multiplicity_probe(pb, lam, first_stage, tol=1e-10, time_limit=200, absolute_motion_tol=1e-06):
    """Second-start numerical witness; both answers must pass the same residual gate.

    Floating-point residuals do not enclose exact roots. Agreement proves nothing.
    Resolved differences also need an absolute motion gap: by default 1e-6 m/s
    in translation or 1e-6 rad/s in rotation, for the engine's SI scene units.
    """
    _check_tol(tol)
    if not np.isfinite(absolute_motion_tol) or absolute_motion_tol <= 0:
        raise ValueError('absolute_motion_tol must be finite and strictly positive')
    lam = np.asarray(lam, float)
    if lam.shape != (3 * pb.nc,) or not np.isfinite(lam).all():
        return dict(verdict='FIRST_START_INVALID', converged=False, certificate=False)
    first_residual = pb.natres(lam)
    if not np.isfinite(first_residual) or first_residual >= tol:
        return dict(verdict='FIRST_START_NOT_CONVERGED', converged=False, first_natres=float(first_residual), certificate=False)
    if not pb.nc:
        return dict(verdict='SECOND_START_AGREES', converged=True, natres=0.0, first_natres=0.0, dv_rel=0.0, certificate=False, lam_second=lam.copy())
    if first_stage == 'ssn_after_ipm':
        r = ssn_escape(pb, np.zeros(3 * pb.nc), tol=tol, max_iter=40, eps0=0.01, dir_from_X='', retries=2, time_limit=time_limit)
        start = 'cold'
    else:
        r0 = solve(pb, tol=0.0001, max_iter=30, desaxce='newton', eps_fac=0.01, time_limit=time_limit)
        r = ssn_escape(pb, r0['lam'], tol=tol, max_iter=40, eps0=0.01, dir_from_X='', retries=2, time_limit=time_limit)
        start = 'interior_point'
    v1 = pb.v_free + pb.Minv @ (pb.J.T @ lam)
    out = dict(start=start, converged=bool(np.isfinite(r['natres']) and r['natres'] < tol), natres=float(r['natres']), first_natres=float(first_residual), certificate=False)
    if out['converged']:
        v2 = pb.v_free + pb.Minv @ (pb.J.T @ r['lam'])
        delta = v2 - v1
        out['dv_rel'] = float(np.sqrt(max(float(delta @ (pb.M @ delta)), 0.0)) / max(np.sqrt(max(float(v1 @ (pb.M @ v1)), 0.0)), 1e-300))
        out['dv_rel_metric'] = 'mass_norm'
        per_body = delta.reshape(-1, 6)
        out['dv_max_translation'] = float(np.linalg.norm(per_body[:, :3], axis=1).max(initial=0.0))
        out['dv_max_angular'] = float(np.linalg.norm(per_body[:, 3:], axis=1).max(initial=0.0))
        (e1, e2) = (float(0.5 * v1 @ (pb.M @ v1)), float(0.5 * v2 @ (pb.M @ v2)))
        out['dE_rel'] = float(abs(e2 - e1) / max(abs(e1), 1e-300))
        out['dE_abs'] = abs(e2 - e1)
        out['motion_resolution'] = dict(translation_m_per_s=absolute_motion_tol, angular_rad_per_s=absolute_motion_tol)
        resolved = max(out['dv_max_translation'], out['dv_max_angular']) > absolute_motion_tol
        out['verdict'] = 'NUMERICAL_MULTIPLICITY_WITNESS' if resolved and out['dv_rel'] > 1e-08 else 'SECOND_START_AGREES'
        out['agreement_scope'] = 'declared_motion_budget'
        out['lam_second'] = r['lam']
    else:
        out['verdict'] = 'SECOND_START_NOT_CONVERGED'
    return out

def _check_tol(tol):
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError('tol must be finite and strictly positive')

def solve_ncp_newton(scene, tol=1e-10, time_limit=200, report=True, probe_multiplicity=False):
    """Engine entry point next to solve_ncp_pgs / solve_ncp_admm.  `scene` is an ncp_ref.Scene or a
    Problem.  Returns (lam, info); info carries the residual, iteration / factorisation counts, the
    stage that produced the answer, (report=True) the local uniqueness verdict and
    (probe_multiplicity=True) a second-start numerical witness, not an exact certificate."""
    _check_tol(tol)
    if not np.isfinite(time_limit) or time_limit <= 0:
        raise ValueError('time_limit must be finite and strictly positive')
    pb = scene if isinstance(scene, Problem) else Problem.from_engine_scene(scene)
    if not pb.nc:
        lam = np.empty(0)
        info = dict(natres=0.0, converged=True, iters=0, n_fact=0, time=0.0, status='converged', solved_by='free_motion', stages=[])
        if report:
            info['uniqueness'] = uniqueness_report(pb, lam, tol=tol)
        if probe_multiplicity:
            info['multiplicity'] = multiplicity_probe(pb, lam, 'free_motion', tol=tol)
        return (lam, info)
    r = solve_newton(pb, tol=tol, time_limit=time_limit)
    info = {k: r[k] for k in ('natres', 'converged', 'iters', 'n_fact', 'time', 'status', 'solved_by', 'stages')}
    if report:
        info['uniqueness'] = uniqueness_report(pb, r['lam'], tol=tol)
    if probe_multiplicity and r['converged']:
        info['multiplicity'] = multiplicity_probe(pb, r['lam'], r['solved_by'], tol=tol, time_limit=time_limit)
    return (r['lam'], info)
