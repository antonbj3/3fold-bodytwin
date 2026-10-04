"""Exact fixed-basis upper stress for a measured box intersecting a simplex."""
import itertools
import numpy as np
from scipy.optimize import linprog

def matrices(voigt):
    s = np.zeros(voigt.shape[:-1] + (3, 3))
    (s[..., 0, 0], s[..., 1, 1], s[..., 2, 2]) = (voigt[..., 0], voigt[..., 1], voigt[..., 2])
    s[..., 0, 1] = s[..., 1, 0] = voigt[..., 3]
    s[..., 1, 2] = s[..., 2, 1] = voigt[..., 4]
    s[..., 0, 2] = s[..., 2, 0] = voigt[..., 5]
    return s

def peak(tensor):
    return float(max(0, np.linalg.eigvalsh(matrices(tensor))[..., -1].max()))

def vertices(lower, upper, tol=1e-10):
    lower = np.asarray(lower, float)
    upper = np.asarray(upper, float)
    n = len(lower)
    if lower.shape != upper.shape or n < 1 or (not np.all(np.isfinite(lower))) or (not np.all(np.isfinite(upper))):
        raise ValueError('finite compatible bounds required')
    if np.any(lower < 0) or np.any(upper > 1) or np.any(lower > upper):
        raise ValueError('fraction bounds must be in [0,1]')
    if lower.sum() > 1 + tol or upper.sum() < 1 - tol:
        raise ValueError('empty pressure/simplex intersection')
    vals = {}
    for free in range(n):
        fixed = [j for j in range(n) if j != free]
        for choices in itertools.product((0, 1), repeat=n - 1):
            a = np.zeros(n)
            for (j, c) in zip(fixed, choices):
                a[j] = upper[j] if c else lower[j]
            a[free] = 1 - a.sum()
            if lower[free] - tol <= a[free] <= upper[free] + tol:
                a[free] = min(upper[free], max(lower[free], a[free]))
                vals[tuple(np.round(a, 12))] = a
    if not vals:
        raise ValueError('no feasible vertices')
    return np.stack(list(vals.values()))

def upper_stress(basis, center, epsilon):
    center = np.asarray(center, float)
    if len(center) != len(basis) or not np.all(np.isfinite(center)) or np.any(center < 0) or (abs(center.sum() - 1) > 1e-10):
        raise ValueError('registered fractions must be nonnegative and sum to one')
    if epsilon < 0 or not np.isfinite(epsilon):
        raise ValueError('finite positive fraction accuracy required')
    lower = np.maximum(0, center - epsilon)
    upper = np.minimum(1, center + epsilon)
    v = vertices(lower, upper)
    winner = None
    ub = -np.inf
    winner_element = None
    for a in v:
        tensor = np.einsum('j,jea->ea', a, basis, optimize=False)
        ev = np.linalg.eigvalsh(matrices(tensor))[..., -1]
        idx = int(np.argmax(ev))
        q = max(0, float(ev[idx]))
        if q > ub:
            ub = q
            winner = a.copy()
            winner_element = idx
    t = np.einsum('j,ja->a', winner, basis[:, winner_element])
    (e, vec) = np.linalg.eigh(matrices(t))
    direction = vec[:, -1]
    coeff = np.einsum('i,jik,k->j', direction, matrices(basis[:, winner_element]), direction)
    lp = linprog(-coeff, A_eq=np.ones((1, len(center))), b_eq=[1], bounds=list(zip(lower, upper)), method='highs')
    if not lp.success:
        raise RuntimeError('independent Rayleigh LP did not solve')
    lp_peak = max(0, float(-lp.fun))
    nominal = peak(np.einsum('j,jea->ea', center, basis))
    norms = np.max(np.abs(np.linalg.eigvalsh(matrices(basis))), axis=(-1, -2))
    lipschitz = nominal + float(np.sum(norms * np.maximum(center - lower, upper - center)))
    return {'epsilon_fraction': epsilon, 'lower': lower.tolist(), 'upper': upper.tolist(), 'vertex_count': len(v), 'upper_MPa': ub, 'witness_fraction': winner.tolist(), 'witness_element': winner_element, 'rayleigh_LP_MPa': lp_peak, 'rayleigh_LP_relative_error': abs(lp_peak - ub) / max(1, ub), 'conservative_spectral_upper_MPa': lipschitz, 'nominal_MPa': nominal, 'witness_feasible': bool(abs(winner.sum() - 1) < 1e-10 and np.all(winner >= lower - 1e-10) and np.all(winner <= upper + 1e-10)), 'physical_status': 'UNKNOWN'}
