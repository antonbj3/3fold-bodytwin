"""Standard stationary active-set polish of independent QP sign discovery."""
import numpy as np
from scipy.linalg import solve
from qp_control import independent_qp

def independent_polished(K, z, preload, load, kb=20000.0, kc=200000.0):
    r = independent_qp(K, z, preload, load, kb, kc)
    x = r['x']
    n = len(z)
    e = preload * (1 / kb + 1 / kc)
    for it in range(8):
        g = x[::2] - z
        t = g + e
        bolt = t > 0
        seat = g < 0
        A = K.copy()
        rhs = np.zeros(2 * n)
        rhs[::2] = load + (kb * bolt + kc * seat) * z - kb * bolt * e
        A[np.arange(0, 2 * n, 2), np.arange(0, 2 * n, 2)] += kb * bolt + kc * seat
        x = solve(A, rhs, assume_a='sym')
        g = x[::2] - z
        t = g + e
        if np.all(np.where(bolt, t, -t) >= -1e-10) and np.all(np.where(seat, -g, g) >= -1e-10):
            break
    else:
        raise RuntimeError('Independent QP sign branch did not close')
    fb = kb * np.maximum(t, 0)
    fc = kc * np.maximum(-g, 0)
    grad = K @ x
    grad[::2] += fb - fc - load
    return {'x': x, 'bolt': fb, 'contact': fc, 'success': bool(np.max(abs(grad)) < 1e-07), 'message': 'QP branch followed by independently pivoted stationary solve', 'constraint_violation_mm': float(max(0, -min(np.min(np.where(bolt, t, -t)), np.min(np.where(seat, -g, g))))), 'gradient_N': float(np.max(abs(grad))), 'initial_QP_gradient_N': r['gradient_N'], 'polish_iterations': it + 1}
