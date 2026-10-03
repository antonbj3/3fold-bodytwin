"""BT-XF4 region-wise emulator.

Mechanism (matches inputs/transfer.json):
  * region classifier on the KKT multiplier sign a_i = (C^T nu)_i  (equivalently
    the active-set membership sign)  -- a k-NN classifier in theta over the
    training targets' exact duals;
  * per-region one-sided KKT model: on a frozen active set S the p=2 system is
    affine,  (C_S C_S^T) nu = b ; a semismooth-Newton sweep re-selects
    S = {i : (C^T nu)_i > 0} (the one-sided/active-set derivative at the region
    boundary).  A fixed number K of sweeps polishes the classifier into the
    exact region without ever averaging a Taylor band across the kink.

Vendors Model / RegionEmulator (xf4_model.py / xf4_core.py, themselves vendored
from BT-N29 with citation).  No cross-agent import.

Fix vs the interrupted session: the analytic one-sided gradient's strength
branch omitted the a*dN term (d f = da*N + a*dN); corrected in xf4_core.gradient.
"""
import numpy as np
from xf4_core import RegionEmulator


class XF4Emulator:
    def __init__(self, K=6, knn=1):
        self.e = RegionEmulator()
        self.K = K
        self.knn = knn
        self.T = self.e.T

    def fit(self, Xtr, NUtr, Rtr=None):
        self.Xtr = np.asarray(Xtr)
        self.NUtr = np.asarray(NUtr)
        self.Rtr = None if Rtr is None else np.asarray(Rtr)

    def _classify(self, x):
        d2 = np.sum((self.Xtr - x) ** 2, axis=1)
        idx = np.argsort(d2)[:self.knn]
        return self.NUtr[idx].mean(0)

    def predict_mask(self, x, nu=None):
        if nu is None:
            nu = self._classify(x)
        delta, sa, ss, r0, r, N = self.e._geom(x)
        _, Cfull, _ = self.e._C_b(r0, r, N, delta)
        return np.einsum('tdf,td->tf', Cfull, nu) > 0, Cfull

    def predict(self, x, K=None):
        """Return (R_hat (T,3), info dict).  R_hat is the region-wise emulator."""
        K = self.K if K is None else K
        e = self.e
        delta, sa, ss, r0, r, N = e._geom(x)
        _, Cfull, b = e._C_b(r0, r, N, delta)
        nu = self._classify(x)
        for _ in range(K):
            masks = np.einsum('tdf,td->tf', Cfull, nu) > 0
            Cact = Cfull * masks[:, None, :]
            M = Cact @ Cact.transpose(0, 2, 1) + 1e-14 * np.eye(4)
            nu = np.linalg.solve(M, b[:, :, None])[:, :, 0]
        masks = np.einsum('tdf,td->tf', Cfull, nu) > 0
        Cact = Cfull * masks[:, None, :]
        M = Cact @ Cact.transpose(0, 2, 1) + 1e-14 * np.eye(4)
        nu = np.linalg.solve(M, b[:, :, None])[:, :, 0]
        a = np.einsum('tdf,td->tf', Cact, nu)
        fi = np.array(e.FJ, copy=True)
        fi[:, e.free] = 0.0
        fi[:, e.freeidx] = a * N[e.free]
        R = e.RJ + np.einsum('tia,ti->ta', e.g, e.FJ) - np.einsum('tia,ti->ta', e.g, fi)
        return R, dict(masks=masks, nu=nu)

    def predict_batch(self, X, K=None):
        R = np.empty((len(X), self.T, 3))
        for i, x in enumerate(X):
            R[i] = self.predict(x, K)[0]
        return R


def rel_error(Rhat, Rtrue):
    num = np.linalg.norm(Rhat - Rtrue, axis=2)
    den = np.linalg.norm(Rtrue, axis=2)
    return num / den
