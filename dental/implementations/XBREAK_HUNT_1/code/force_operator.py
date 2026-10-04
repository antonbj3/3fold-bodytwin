"""Force-balanced contact response; ordinary convex contact control included."""
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize, LinearConstraint

def balanced_basis(A):
    return null_space(A.T, rcond=1e-12)

def dual_contact(N, H, f0, h):
    G = np.linalg.inv(H)
    rhs = N.T @ h

    def fun(z):
        return 0.5 * z @ G @ z - rhs @ z

    def jac(z):
        return G @ z - rhs
    sol = minimize(fun, np.zeros(N.shape[1]), jac=jac, constraints=[LinearConstraint(N, -f0, np.full(len(f0), np.inf))], method='SLSQP', options={'ftol': 1e-13, 'maxiter': 1000})
    f = f0 + N @ sol.x
    active = f < 1e-06
    grad = jac(sol.x)
    if active.any():
        mu = np.linalg.lstsq(N[active].T, grad, rcond=1e-12)[0]
        stationarity = float(np.max(np.abs(grad - N[active].T @ mu)))
        dual_min = float(mu.min())
    else:
        stationarity = float(np.max(np.abs(grad)))
        dual_min = 0.0
    return (f, {'success': bool(sol.success), 'message': str(sol.message), 'iterations': int(sol.nit), 'min_force_N': float(f.min()), 'stationarity_mm': stationarity, 'dual_min_mm': dual_min})
