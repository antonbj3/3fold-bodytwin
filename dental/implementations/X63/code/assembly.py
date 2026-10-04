"""Convex, unilateral four-support rail; units N, mm, MPa.

This is an unfolded arch Euler-Bernoulli surrogate, not a calibrated jaw FE.
Positive w opens a support. Bolts carry tension, seats carry compression.
"""
from fractions import Fraction as Q
from itertools import product
import numpy as np
from scipy.optimize import minimize

def beam_matrix(lengths, h=4, width=4, E=110000, exact=False):
    n = len(lengths) + 1
    K = [[Q(0) for _ in range(2 * n)] for _ in range(2 * n)] if exact else np.zeros((2 * n, 2 * n))
    q = Q if exact else float
    EI = q(E) * q(width) * q(h) ** 3 / q(12)
    for (j, value) in enumerate(lengths):
        L = q(str(value))
        rows = [[12, 6 * L, -12, 6 * L], [6 * L, 4 * L * L, -6 * L, 2 * L * L], [-12, -6 * L, 12, -6 * L], [6 * L, 2 * L * L, -6 * L, 4 * L * L]]
        ix = [2 * j, 2 * j + 1, 2 * j + 2, 2 * j + 3]
        for a in range(4):
            for b in range(4):
                K[ix[a]][ix[b]] += EI / L ** 3 * rows[a][b]
    return K

def rational_solve(A, b):
    n = len(b)
    M = [list(row) + [v] for (row, v) in zip(A, b)]
    for k in range(n):
        pivot = next((j for j in range(k, n) if M[j][k]), None)
        if pivot is None:
            raise ValueError('singular branch')
        (M[k], M[pivot]) = (M[pivot], M[k])
        p = M[k][k]
        M[k] = [v / p for v in M[k]]
        for j in range(n):
            if j != k and M[j][k]:
                a = M[j][k]
                M[j] = [x - a * y for (x, y) in zip(M[j], M[k])]
    return [r[-1] for r in M]

def state(K, z, preload, load, kb=20000.0, kc=200000.0):
    """Enumerate active spring branches; strictly convex relevant branches."""
    n = len(z)
    e = preload * (1 / kb + 1 / kc)
    best = None
    for ab in product([False, True], repeat=n):
        for ac in product([False, True], repeat=n):
            (b, c) = (np.array(ab), np.array(ac))
            A = K.copy()
            rhs = np.zeros(2 * n)
            rhs[::2] = load
            A[np.arange(0, 2 * n, 2), np.arange(0, 2 * n, 2)] += kb * b + kc * c
            rhs[::2] += (kb * b + kc * c) * z - kb * b * e
            try:
                x = np.linalg.solve(A, rhs)
            except np.linalg.LinAlgError:
                continue
            g = x[::2] - z
            t = g + e
            if np.any(np.where(b, t, -t) < -1e-10) or np.any(np.where(c, -g, g) < -1e-10):
                continue
            fb = kb * np.maximum(t, 0)
            fc = kc * np.maximum(-g, 0)
            grad = K @ x
            grad[::2] += fb - fc - load
            if np.max(np.abs(grad)) > 1e-07:
                continue
            energy = 0.5 * x @ K @ x + 0.5 * kb * np.maximum(t, 0) @ np.maximum(t, 0) + 0.5 * kc * np.minimum(g, 0) @ np.minimum(g, 0) - load @ x[::2]
            if best is None or energy < best['energy']:
                best = {'x': x, 'gap': np.maximum(g, 0), 'signed_gap': g, 'bolt': fb, 'contact': fc, 'active_bolts': b, 'active_contacts': c, 'energy': energy, 'residual_N': float(np.max(np.abs(grad)))}
    if best is None:
        raise RuntimeError('No verified branch')
    return best

def independent_state(K, z, preload, load, kb=20000.0, kc=200000.0):
    """Independent convex minimizer, no active-set solution as initial point."""
    n = len(z)
    e = preload * (1 / kb + 1 / kc)
    scale = np.ones(2 * n)
    scale[1::2] = 1 / 20

    def fun(y):
        x = y * scale
        g = x[::2] - z
        t = g + e
        val = 0.5 * x @ K @ x + 0.5 * kb * np.maximum(t, 0) @ np.maximum(t, 0) + 0.5 * kc * np.minimum(g, 0) @ np.minimum(g, 0) - load @ x[::2]
        grad = K @ x
        grad[::2] += kb * np.maximum(t, 0) + kc * np.minimum(g, 0) - load
        return (val, grad * scale)
    r = minimize(fun, np.zeros(2 * n), jac=True, method='BFGS', options={'gtol': 1e-10, 'maxiter': 3000})
    x = r.x * scale
    g = x[::2] - z
    return {'x': x, 'bolt': kb * np.maximum(g + e, 0), 'contact': kc * np.maximum(-g, 0), 'success': bool(r.success), 'message': str(r.message), 'gradient_N': float(np.max(np.abs(fun(r.x)[1])))}

def exact_first_event(lengths, z, preload, load_direction, h=4, kb=20000, kc=200000):
    """Exact rational affine response and validity interval on all-active branch.

    No floating enclosure is needed for the ideal *quantized* input model.
    Missing material/physical uncertainty is not enclosed by this calculation.
    """
    K = beam_matrix(lengths, h=h, exact=True)
    n = len(z)
    z = [Q(str(v)) for v in z]
    p = Q(str(preload))
    e = p * (Q(1, kb) + Q(1, kc))
    A = [r[:] for r in K]
    b = [Q(0)] * (2 * n)
    f = [Q(0)] * (2 * n)
    for j in range(n):
        A[2 * j][2 * j] += kb + kc
        b[2 * j] = (kb + kc) * z[j] - kb * e
        f[2 * j] = Q(str(load_direction[j]))
    x0 = rational_solve(A, b)
    v = rational_solve(A, f)
    g0 = [x0[2 * j] - z[j] for j in range(n)]
    t0 = [g + e for g in g0]
    if any((g > 0 for g in g0)) or any((t < 0 for t in t0)):
        return {'status': 'ASSEMBLY_BRANCH_INVALID', 'reason': 'An initially open support or slack bolt; use unilateral solver.'}
    roots = []
    for j in range(n):
        if v[2 * j] > 0:
            roots.append((-g0[j] / v[2 * j], j, 'contact_open'))
        if v[2 * j] < 0:
            roots.append((-t0[j] / v[2 * j], j, 'bolt_slack'))
    (lam, j, kind) = min(roots)
    residual = [sum((A[a][b] * x0[b] for b in range(2 * n))) - b0 for (a, b0) in enumerate(b)]
    return {'status': 'EXACT_RATIONAL_FIXED_MODEL', 'first_event_N': float(lam), 'first_event_fraction': str(lam), 'support_index': j, 'kind': kind, 'response_valid_interval_N': ['0', str(lam)], 'identity_residual_exact': str(max((abs(vv) for vv in residual))), 'initial_signed_gap_mm': [float(g) for g in g0], 'initial_contact_N': [float(-kc * g) for g in g0], 'x0_fractions': [str(x) for x in x0], 'response_fractions_per_N': [str(x) for x in v], 'physical_uncertainty_enclosure': 'MISSING; geometry quantization, beam reduction, E, kb, kc and load are conditional scenarios.'}
