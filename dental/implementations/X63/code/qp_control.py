"""Independent convex QP epigraph control, identical physical information."""
import numpy as np
from scipy.optimize import minimize

def independent_qp(K, z, preload, load, kb=20000.0, kc=200000.0):
    n = len(z)
    nd = 2 * n
    e = preload * (1 / kb + 1 / kc)
    scale = np.r_[np.tile([0.01, 0.001], n), np.full(2 * n, 0.01)]
    H = np.zeros((4 * n, 4 * n))
    H[:nd, :nd] = K
    H[nd:nd + n, nd:nd + n] = np.eye(n) * kb
    H[nd + n:, nd + n:] = np.eye(n) * kc
    q = np.zeros(4 * n)
    q[:nd:2] = -load
    Hs = H * scale[:, None] * scale[None, :]
    qs = q * scale

    def fun(y):
        return (0.5 * y @ Hs @ y + qs @ y, Hs @ y + qs)
    C = np.zeros((2 * n, 4 * n))
    lower = np.r_[e - z, z]
    for j in range(n):
        C[j, nd + j] = 1
        C[j, 2 * j] = -1
        C[n + j, nd + n + j] = 1
        C[n + j, 2 * j] = 1
    Cs = C * scale[None, :]
    constraints = {'type': 'ineq', 'fun': lambda y: Cs @ y - lower, 'jac': lambda y: Cs}
    bounds = [(None, None)] * nd + [(0, None)] * (2 * n)
    initial = np.zeros(4 * n)
    initial[nd:nd + n] = np.maximum(e - z, 0) / 0.01
    initial[nd + n:] = np.maximum(z, 0) / 0.01
    r = minimize(fun, initial, jac=True, method='SLSQP', bounds=bounds, constraints=constraints, options={'ftol': 1e-13, 'maxiter': 4000})
    x = (r.x * scale)[:nd]
    g = x[::2] - z
    t = np.maximum(g + e, 0)
    c = np.maximum(-g, 0)
    grad = K @ x
    grad[::2] += kb * t - kc * c - load
    return {'x': x, 'bolt': kb * t, 'contact': kc * c, 'success': bool(r.success), 'message': str(r.message), 'constraint_violation_mm': float(max(0, -np.min(Cs @ r.x - lower))), 'gradient_N': float(np.max(abs(grad))), 'iterations': int(r.nit)}
