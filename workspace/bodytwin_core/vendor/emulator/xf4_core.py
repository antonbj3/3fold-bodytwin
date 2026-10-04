"""BT-XF4 core: region-wise emulator for the K1 p=2 hip reaction.

Vendors the exact model from BT-N29 (`n29_core.py`, model = K1 p=2 dual-Newton;
no cross-agent import, copied with citation).  Adds
  * a *region-frozen* p=2 solve (affine on a fixed active set S),
  * the analytic one-sided KKT derivative dR/dx on a frozen active set S
    (N13-style implicit differentiation, vectorised over frames),
  * region-classifier + local affine ("Taylor") and exact frozen-solve emulators.

Parameter vector x in [-1,1]^23 -> decode identical to BT-N29:
  delta_i = 0.01 * x_i (m),  s_arm_m = 1 + 0.15 x_m,  s_str_m = 1 + 0.15 x_m
  for the 10 dominant muscles.  Chain factors: 0.01 (delta), 0.15 (arms/strength).
"""
import numpy as np
from xf4_model import Model, _dual_newton  # noqa: F401  (vendored BT-N29)


class RegionEmulator:
    def __init__(self):
        self.m = Model()
        self.T = self.m.T
        self.n = self.m.n
        self.free = self.m.free
        self.NF = int(self.free.sum())
        self.freeidx = np.where(self.free)[0]
        self.dom = self.m.dom
        self.kd = len(self.dom)
        self.D = 3 + 2 * self.kd
        self.scale = np.empty(self.D)
        self.scale[:3] = 0.01          # d/dx = 0.01 d/ddelta
        self.scale[3:3 + self.kd] = 0.15
        self.scale[3 + self.kd:] = 0.15
        self.g = self.m.g
        self.k = self.m.k
        self.FJ = self.m.FJ
        self.RJ = self.m.RJ
        self.N0 = self.m.N0

    # -------- geometry --------
    def _geom(self, x):
        delta, sa, ss = self.m.decode(x)
        delta, sa, ss = delta[0], sa[0], ss[0]
        r0 = self.m._shifted(delta)
        r = self.m._scaled(r0, sa)
        N = self.N0 * ss
        return delta, sa, ss, r0, r, N

    def _C_b(self, r0, r, N, delta):
        Aall = np.empty((self.T, 4, self.n))
        Aall[:, :3, :] = r.transpose(0, 2, 1)
        Aall[:, 3, :] = self.k
        b = np.empty((self.T, 4))
        b[:, :3] = np.einsum('tma,tm->ta', r0, self.FJ) - np.cross(delta, self.RJ)
        b[:, 3] = np.einsum('tm,tm->t', self.k, self.FJ)
        b = b - np.einsum('tqm,tm->tq', Aall[:, :, ~self.free], self.FJ[:, ~self.free])
        Cfull = Aall[:, :, self.free] * N[self.free][None, None, :]
        return Aall, Cfull, b

    # -------- exact / frozen solves --------
    def masks_from_nu(self, Cfull, nu):
        return (np.einsum('tdf,td->tf', Cfull, nu) > 0)

    def masks_exact(self, x, nu=None):
        delta, sa, ss, r0, r, N = self._geom(x)
        _, Cfull, b = self._C_b(r0, r, N, delta)
        if nu is None:
            o = self.m.targets(x)
            nu = o['nu']
        return self.masks_from_nu(Cfull, nu)

    def frozen_R(self, x, masks):
        """Exact p=2 solution restricted to the given per-frame free-muscle mask."""
        delta, sa, ss, r0, r, N = self._geom(x)
        Aall, Cfull, b = self._C_b(r0, r, N, delta)
        Cact = Cfull * masks[:, None, :]
        M = Cact @ Cact.transpose(0, 2, 1) + 1e-14 * np.eye(4)
        nu = np.linalg.solve(M, b[:, :, None])[:, :, 0]
        a = np.einsum('tdf,td->tf', Cact, nu)
        fi = np.array(self.FJ, copy=True)
        fi[:, self.free] = 0.0
        fi[:, self.freeidx] = a * N[self.free]
        R = self.RJ + np.einsum('tia,ti->ta', self.g, self.FJ) - np.einsum('tia,ti->ta', self.g, fi)
        return R

    # -------- analytic one-sided gradient on frozen mask --------
    def gradient(self, x, masks, nu):
        delta, sa, ss, r0, r, N = self._geom(x)
        Aall, Cfull, b = self._C_b(r0, r, N, delta)
        Cact = Cfull * masks[:, None, :]
        M = Cact @ Cact.transpose(0, 2, 1)
        Minv = np.linalg.inv(M + 1e-12 * np.eye(4))
        gfree = self.g[:, self.free, :]
        Nfree = N[self.free]
        # activation a = C^T nu on the frozen set (f = a*N, so d f = da*N + a*dN)
        a_frozen = np.einsum('tdf,td->tf', Cact, nu)
        J = np.zeros((self.D, self.T, 3))
        zero4 = np.zeros((self.T, 4))

        def dR_from(dCfull, db, dNfree):
            dCact = dCfull * masks[:, None, :]
            dM = dCact @ Cact.transpose(0, 2, 1) + Cact @ dCact.transpose(0, 2, 1)
            dnu = (Minv @ (db[:, :, None] - dM @ nu[:, :, None]))[:, :, 0]
            da = np.einsum('tdf,td->tf', dCact, nu) + np.einsum('tdf,td->tf', Cact, dnu)
            df = da * Nfree[None, :] + a_frozen * dNfree[None, :]
            return -np.einsum('tfa,tf->ta', gfree, df)

        e = np.eye(3)
        for j in range(3):
            G = np.cross(e[j], self.g)                     # (T,n,3) e_j x g_i
            drd = -G                                       # d r0 / d delta_j
            dAall = np.zeros((self.T, 4, self.n))
            dAall[:, :3, :] = (drd * sa[None, :, None]).transpose(0, 2, 1)
            dC = dAall[:, :, self.free] * Nfree[None, None, :]
            db = np.zeros((self.T, 4))
            db[:, :3] = np.einsum('tma,tm->ta', drd, self.FJ) - np.cross(e[j], self.RJ)
            db = db - np.einsum('tqm,tm->tq', dAall[:, :, ~self.free], self.FJ[:, ~self.free])
            J[j] = dR_from(dC, db, np.zeros(self.NF))
        for q, mm in enumerate(self.dom):
            mf = int(np.where(self.freeidx == mm)[0][0])
            # moment-arm scale: dC = (d r / d s_arm) * N  (N unchanged)
            dAall = np.zeros((self.T, 4, self.n))
            dAall[:, :3, mm] = r0[:, mm, :]
            dC = dAall[:, :, self.free] * Nfree[None, None, :]
            J[3 + q] = dR_from(dC, zero4, np.zeros(self.NF))
            # strength scale: dN_m = N0_m, dC = Aall_m * N0_m (row 3 included)
            dAall2 = np.zeros((self.T, 4, self.n))
            dAall2[:, :, mm] = Aall[:, :, mm]
            dC2 = dAall2[:, :, self.free] * self.N0[mm]
            dNfree = np.zeros(self.NF)
            dNfree[mf] = self.N0[mm]
            J[3 + self.kd + q] = dR_from(dC2, zero4, dNfree)
        J *= self.scale[:, None, None]
        return J

    def state(self, x):
        """Exact targets + mask + analytic gradient at x (one call per centre)."""
        o = self.m.targets(x)
        delta, sa, ss, r0, r, N = self._geom(x)
        _, Cfull, _ = self._C_b(r0, r, N, delta)
        masks = self.masks_from_nu(Cfull, o['nu'])
        J = self.gradient(x, masks, o['nu'])
        return dict(R=o['R'], nu=o['nu'], masks=masks, J=J)

    # -------- emulator predictions --------
    def predict_affine(self, x, centres, centre_X, centre_R, centre_J, weights=None):
        """Equal (or weighted) average of region-frozen affine Taylor models."""
        R = np.zeros((self.T, 3))
        wsum = 0.0
        for c, i in enumerate(centres):
            w = 1.0 if weights is None else weights[c]
            R = R + w * (centre_R[i] + np.einsum('dta,d->ta', centre_J[i], x - centre_X[i]))
            wsum += w
        return R / wsum


if __name__ == '__main__':
    import time
    e = RegionEmulator()
    rng = np.random.default_rng(5)
    # control: analytic gradient vs frozen central FD on the same mask
    for xi in range(3):
        x = rng.uniform(-1, 1, e.D)
        st = e.state(x)
        h = 1e-5
        fd = np.zeros((e.D, e.T, 3))
        for j in range(e.D):
            xp = x.copy(); xp[j] += h
            xm = x.copy(); xm[j] -= h
            fd[j] = (e.frozen_R(xp, st['masks']) - e.frozen_R(xm, st['masks'])) / (2 * h)
        rel = np.linalg.norm(fd - st['J'], axis=2) / (np.linalg.norm(st['J'], axis=2) + 1e-30)
        print('grad control x%d median %.2e max %.2e' % (xi, np.median(rel), rel.max()))
