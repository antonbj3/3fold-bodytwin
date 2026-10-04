import numpy as np
from scipy.optimize import minimize

def control_box(center, entry, axis, L, r, sp):
    lo = center - sp / 2
    hi = center + sp / 2
    v = np.cross(axis, [1.0, 0, 0] if abs(axis[0]) < 0.9 else [0, 1.0, 0])
    v /= np.linalg.norm(v)
    w = np.cross(axis, v)

    def xyz(z):
        return entry + z[3] * axis + z[4] * v + z[5] * w

    def fun(z):
        dv = z[:3] - xyz(z)
        return dv @ dv

    def jac(z):
        dv = z[:3] - xyz(z)
        return np.r_[2 * dv, -2 * dv @ axis, -2 * dv @ v, -2 * dv @ w]
    c = lambda z: r * r - z[4] ** 2 - z[5] ** 2
    cj = lambda z: np.array([0.0, 0, 0, 0, -2 * z[4], -2 * z[5]])
    t = np.clip((center - entry) @ axis, 0, L)
    uv = np.array([(center - entry) @ v, (center - entry) @ w])
    uv *= min(1, r / max(np.linalg.norm(uv), 1e-30))
    res = minimize(fun, np.r_[center, t, uv], jac=jac, bounds=list(zip(lo, hi)) + [(0, L), (-r, r), (-r, r)], constraints=[{'type': 'ineq', 'fun': c, 'jac': cj}], method='SLSQP', options={'ftol': 1e-12, 'maxiter': 200})
    zz = res.x.copy()
    zz[:3] = np.clip(zz[:3], lo, hi)
    zz[3] = np.clip(zz[3], 0, L)
    zz[4:6] *= min(1, r / max(np.linalg.norm(zz[4:6]), 1e-30))
    return (float(np.sqrt(max(0, fun(zz)))), bool(res.success))
