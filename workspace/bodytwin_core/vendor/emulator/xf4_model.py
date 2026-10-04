"""BT-N29 core: exact p=2 hip-recruitment model on K1 data + parameter sweep.

Builds on inputs/K1/k1_tests.py (solve_poly, p=2), inputs/K1/k1_extract.py (k1_hip.npz).
The p=2 exact solve is replaced by a 4-D dual Newton active-set iteration whose
KKT residual is <=3.3e-15 and whose peak R matches K1's B0_p2 peak_R_N
(2583.268667071279 N) to ~1e-12 (see n29_bench.json / n29_core.check_k1).
No source code is imported across agents; functions reimplemented with citation.
"""
import numpy as np

K1 = 'inputs/BT-N29/inputs/K1/k1_hip.npz'


class Model:
    def __init__(self, k1=K1, n_dom=10):
        z = np.load(k1)
        self.z = z
        self.names = [str(x) for x in z['names']]
        self.T = len(z['t'])
        self.n = len(self.names)
        self.rect = np.array([m.startswith('RectusFemoris') for m in self.names])
        self.free = ~self.rect
        self.N0 = z['N'][0]
        self.g = z['g']; self.r = z['r']; self.k = z['k']
        self.FJ = z['Ft']; self.RJ = z['R']; self.isb = z['isb']
        # baseline p=2 solution at all frames -> dominant muscles by peak moment contribution
        self.base = self.solve_params(np.zeros(3), np.ones(self.n), np.ones(self.n))
        peak = int(np.argmax(np.linalg.norm(self.RJ, axis=1)))
        self.peak = peak
        contrib = np.linalg.norm(self.r[peak] * self.base['f'][peak][:, None], axis=1)
        order = np.argsort(-contrib)
        self.dom = [i for i in order if self.free[i]][:n_dom]
        self.dom_names = [self.names[i] for i in self.dom]
        self.d = 3 + 2 * n_dom
        self.hjc_mm = 10.0
        self.scale_frac = 0.15          # scales in [1-f, 1+f]

    # ---------------- exact model ----------------
    def _shifted(self, delta):
        r = self.r
        if np.any(delta):
            r = r - np.cross(delta, self.g)
        return r

    def _scaled(self, r, s_arm):
        r = r.copy()
        if np.any(s_arm != 1.0):
            r[:, self.dom, :] *= s_arm[self.dom][None, :, None]
        return r

    def solve_params(self, delta, s_arm, s_str):
        """Exact p=2 solve for one parameter setting.  delta in metres (3,),
        s_arm/s_str full-length per-muscle scale vectors (92,).  Returns dict.
        b uses the TRUE (unscaled) shifted geometry, exactly as K1 L2 kept the
        demand from the wrapped reference geometry; A uses the scaled arms."""
        r0 = self._shifted(delta)          # for demand b (unscaled arms)
        r = self._scaled(r0, s_arm)        # for constraint A
        N = self.N0 * s_str
        g, k = self.g, self.k
        FJ, RJ = self.FJ, self.RJ
        free = self.free
        T = self.T
        R = np.empty((T, 3)); Rmag = np.empty(T); nu = np.empty((T, 4))
        act = np.empty((T, self.n)); f = np.empty((T, self.n))
        for ti in range(T):
            rt = r[ti]
            Aall = np.vstack([rt.T, k[ti][None, :]])
            b = np.concatenate([r0[ti].T @ FJ[ti] - np.cross(delta, RJ[ti]), [k[ti] @ FJ[ti]]])
            b = b - Aall[:, ~free] @ FJ[ti][~free]
            C = Aall[:, free] * N[free][None, :]
            a, v = _dual_newton(C, b)
            fi = FJ[ti].copy(); fi[free] = a * N[free]
            R[ti] = RJ[ti] + g[ti].T @ FJ[ti] - g[ti].T @ fi
            Rmag[ti] = np.linalg.norm(R[ti]); nu[ti] = v
            act[ti][free] = (C.T @ v > 0) * (a * N[free])
            f[ti] = fi
        return dict(R=R, Rmag=Rmag, nu=nu, act=act, f=f)

    def targets(self, x):
        delta, s_arm, s_str = self.decode(x)
        return self.solve_params(delta[0], s_arm[0], s_str[0])

    def decode(self, x):
        x = np.atleast_2d(x)
        k = len(self.dom)
        n = x.shape[0]
        delta = x[:, :3] * (self.hjc_mm / 1000.0)                       # metres
        sarm = np.ones((n, self.n)); sarm[:, self.dom] = 1.0 + x[:, 3:3 + k] * self.scale_frac
        sstr = np.ones((n, self.n)); sstr[:, self.dom] = 1.0 + x[:, 3 + k:3 + 2 * k] * self.scale_frac
        return delta, sarm, sstr

    def _C(self, ti, r, N):
        Aall = np.vstack([r[ti].T, self.k[ti][None, :]])
        return Aall[:, self.free] * N[self.free][None, :]

    def active_masks(self, x, nu_pred):
        """Activity-sign masks from predicted dual nu_pred, per frame (T,n)."""
        delta, s_arm, s_str = self.decode(x)
        delta = delta[0]; s_arm = s_arm[0]; s_str = s_str[0]
        r = self._scaled(self._shifted(delta), s_arm); N = self.N0 * s_str
        masks = np.zeros((self.T, self.n), bool)
        for ti in range(self.T):
            C = self._C(ti, r, N)
            masks[ti, self.free] = C.T @ nu_pred[ti] > 0
        return masks

    def margin(self, x, nu_pred):
        """Intrinsic active-set margin min_i |C^T nu| / |b| per frame (T,)."""
        delta, s_arm, s_str = self.decode(x)
        delta = delta[0]; s_arm = s_arm[0]; s_str = s_str[0]
        r = self._scaled(self._shifted(delta), s_arm); N = self.N0 * s_str
        r0 = self._shifted(delta)
        out = np.empty(self.T)
        for ti in range(self.T):
            C = self._C(ti, r, N)
            b = np.concatenate([r0[ti].T @ self.FJ[ti] - np.cross(delta, self.RJ[ti]),
                                [self.k[ti] @ self.FJ[ti]]])
            b = b - np.vstack([r[ti].T, self.k[ti][None, :]])[:, ~self.free] @ self.FJ[ti][~self.free]
            s = C.T @ nu_pred[ti]
            out[ti] = np.min(np.abs(s)) / max(1e-12, np.linalg.norm(b))
        return out

    def reconstruct(self, x, nu_pred):
        """Physics reconstruction of R from predicted dual nu_pred (T,4)."""
        delta, s_arm, s_str = self.decode(x)
        delta = delta[0]; s_arm = s_arm[0]; s_str = s_str[0]
        r = self._scaled(self._shifted(delta), s_arm); N = self.N0 * s_str
        g = self.g; FJ = self.FJ; RJ = self.RJ; free = self.free
        R = np.empty((self.T, 3))
        for ti in range(self.T):
            C = self._C(ti, r, N)
            a = np.maximum(C.T @ nu_pred[ti], 0.0)
            fi = FJ[ti].copy(); fi[free] = a * N[free]
            R[ti] = RJ[ti] + g[ti].T @ FJ[ti] - g[ti].T @ fi
        return R


def _dual_newton(C, b, iters=40, tol=1e-12):
    """Exact solution of min 1/2||a||^2 s.t. C a = b, a>=0 (p=2, N folded into C).
    Dual Newton with active-set mask; returns (a, nu)."""
    CCt = C @ C.T + 1e-12 * np.eye(C.shape[0])
    nu = np.linalg.solve(CCt, b)
    for _ in range(iters):
        s = C.T @ nu
        mask = s > 0
        a = np.where(mask, s, 0.0)
        grad = C @ a - b
        if np.linalg.norm(grad) <= tol * max(1.0, np.linalg.norm(b)):
            break
        Cm = C[:, mask]
        d = np.linalg.solve(Cm @ Cm.T + 1e-10 * np.eye(C.shape[0]), grad)
        p0 = 0.5 * a @ a - nu @ b
        t = 1.0
        while t > 1e-6:
            nu2 = nu - t * d
            s2 = np.maximum(C.T @ nu2, 0.0)
            if 0.5 * s2 @ s2 - nu2 @ b <= p0 - 1e-4 * t * (grad @ d):
                break
            t *= 0.5
        nu = nu - t * d
    s = C.T @ nu
    return np.where(s > 0, s, 0.0), nu
