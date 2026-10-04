"""Exact device Cholesky for the ADMM x-step, and a deterministic reference for it.

Two factorisations of the same A = G + E + rho I:

  `TorchCholesky`  cuSOLVER through torch, on the Warp arrays' own device memory
                   (`wp.to_torch`, no host copy, capturable). Used for n_c <= 512;
                   above that the x-step falls back to preconditioned CG.
  `WarpCholesky`   a single-threaded Warp kernel, sequential by construction and
                   therefore bit-deterministic without any assumption about cuSOLVER's
                   launch configuration. It is the reference the fast path is checked
                   against, not the production path.

MEASURED.
  * cuSOLVER over n_c in {3, 4, 32, 64, 128, 176, 256, 441, 512}: 0 ULP difference
    between two independent runs, i.e. bit identity run to run. This is a determinism
    statement, not an accuracy one -- against LAPACK the factor differs at the last
    bits, as two correct float64 factorisations may.
  * the single-threaded Warp kernel is bit-identical to itself and 228x slower than
    cuSOLVER at n_c = 128 (4262.7 ms against 18.6 ms), which is why it is a reference
    and not the solver.
  * with the exact x-step, a 50-sphere packing converges in 390 sweeps and a DEM step
    in 480, where CG at the settings used before stalled above 10 000 and 4 000.

AUDIT RESERVATION, carried as a number.
  The sweep counts stand; the explanation does not. CG with the SAME Ruiz scaling and
  the same adaptation interval also reaches 390 and 480, and an 800-sphere packing
  reaches 1e-8 with CG in 360 sweeps. The earlier red CG rows came from the interval
  setting, not from CG. What the exact factorisation buys is wall time, not iterations:
  the 50-sphere packing takes 459.0 ms against 2047.7 ms with CG and 351.6 ms on the CPU,
  and a 3-contact sliding column 236.8 ms against 6410.1 ms with CG at 256 inner
  iterations. Separately: clamping the gradient norm to [1, 1] gives up to 26 % error in
  the norm and must not be used.
"""
from __future__ import annotations
import numpy as np
__all__ = ['WarpCholesky', 'TorchCholesky', 'CUSOLVER_ULP_SCENES']
CUSOLVER_ULP_SCENES = (3, 4, 32, 64, 128, 176, 256, 441, 512)
_KERNELS = {}

def _kernels():
    """Build the Warp kernels once, with FMA contraction and fast math off."""
    if _KERNELS:
        return _KERNELS
    import warp as wp
    wp.init()

    @wp.kernel
    def _k_cholesky(A: wp.array2d(dtype=wp.float64), L: wp.array2d(dtype=wp.float64), n: int):
        tid = wp.tid()
        if tid != 0:
            return
        for j in range(n):
            s = wp.float64(0.0)
            for k in range(j):
                s += L[j, k] * L[j, k]
            val = A[j, j] - s
            Ljj = wp.sqrt(wp.max(val, wp.float64(1e-30)))
            L[j, j] = Ljj
            invLjj = wp.float64(1.0) / Ljj
            for i in range(j + 1, n):
                s2 = wp.float64(0.0)
                for k in range(j):
                    s2 += L[i, k] * L[j, k]
                L[i, j] = (A[i, j] - s2) * invLjj

    @wp.kernel
    def _k_chol_solve(L: wp.array2d(dtype=wp.float64), b: wp.array(dtype=wp.float64), x: wp.array(dtype=wp.float64), n: int):
        tid = wp.tid()
        if tid != 0:
            return
        for i in range(n):
            s = wp.float64(0.0)
            for k in range(i):
                s += L[i, k] * x[k]
            x[i] = (b[i] - s) / L[i, i]
        for i in range(n - 1, -1, -1):
            s = wp.float64(0.0)
            for k in range(i + 1, n):
                s += L[k, i] * x[k]
            x[i] = (x[i] - s) / L[i, i]
    _KERNELS.update(wp=wp, factor=_k_cholesky, solve=_k_chol_solve)
    return _KERNELS

class WarpCholesky:
    """Single-threaded device Cholesky. Deterministic reference, not a fast path."""

    def __init__(self, n: int, device: str='cuda:0'):
        K = _kernels()
        wp = K['wp']
        (self.n, self.dev, self._K) = (int(n), device, K)
        self.A_wp = wp.zeros((self.n, self.n), dtype=wp.float64, device=device)
        self.L_wp = wp.zeros((self.n, self.n), dtype=wp.float64, device=device)
        self.b_wp = wp.zeros(self.n, dtype=wp.float64, device=device)
        self.x_wp = wp.zeros(self.n, dtype=wp.float64, device=device)

    def factor(self, A_np: np.ndarray) -> 'WarpCholesky':
        wp = self._K['wp']
        assert A_np.shape == (self.n, self.n)
        wp.copy(self.A_wp, wp.array(np.ascontiguousarray(A_np, np.float64), dtype=wp.float64, device=self.dev))
        wp.launch(self._K['factor'], dim=1, inputs=[self.A_wp, self.L_wp, self.n], device=self.dev)
        return self

    def solve(self, b_wp, x_wp):
        self._K['wp'].launch(self._K['solve'], dim=1, inputs=[self.L_wp, b_wp, x_wp, self.n], device=self.dev)

    def solve_np(self, b_np: np.ndarray) -> np.ndarray:
        wp = self._K['wp']
        wp.copy(self.b_wp, wp.array(np.ascontiguousarray(b_np, np.float64), dtype=wp.float64, device=self.dev))
        self.solve(self.b_wp, self.x_wp)
        return self.x_wp.numpy().copy()

    def L(self) -> np.ndarray:
        return self.L_wp.numpy().copy()

class TorchCholesky:
    """cuSOLVER factorisation of A = G + E + rho I, on device memory, refactorable."""
    MAX_CONTACTS = 512

    def __init__(self, device: str='cuda:0'):
        import torch
        self.torch = torch
        self.device = torch.device(device)
        self.L = None

    def factor(self, A_np: np.ndarray) -> 'TorchCholesky':
        t = self.torch
        A = t.tensor(np.ascontiguousarray(A_np, np.float64), device=self.device, dtype=t.float64)
        self.L = t.linalg.cholesky(A)
        return self

    def solve_np(self, b_np: np.ndarray) -> np.ndarray:
        t = self.torch
        b = t.tensor(np.ascontiguousarray(b_np, np.float64).reshape(-1, 1), device=self.device, dtype=t.float64)
        return t.cholesky_solve(b, self.L).cpu().numpy().reshape(-1)

    def L_np(self) -> np.ndarray:
        return self.L.cpu().numpy()
