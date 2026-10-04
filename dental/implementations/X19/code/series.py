"""Passive source-conditioned effective series compliance. Not a material law."""
import numpy as np
from scipy.optimize import lsq_linear

def fit(endpoints, control=False):
    s = np.array([r['sheet_mm'] for r in endpoints])
    K = np.array([r['K_N_per_mm'] for r in endpoints])
    x = 1 / s ** 3
    y = 1 / K
    if np.any(K <= 0):
        raise ValueError('positive measured stiffness required')
    if control:
        ab = lsq_linear(np.column_stack([x, np.ones(2)]), y, bounds=(0, np.inf), tol=1e-14).x
        (a, b) = ab
    else:
        a = (y[0] - y[1]) / (x[0] - x[1])
        b = y[0] - a * x[0]
    if a < 0 or b < 0:
        raise ValueError('nonpassive series compliance closure rejected')
    return {'a_mm4_per_N': float(a), 'b_mm_per_N': float(b), 'fit_scope': 'original-sheet thickness, source condition retained; no physical separation of bending and contact identified'}

def stiffness(f, s):
    return 1 / (f['a_mm4_per_N'] / s ** 3 + f['b_mm_per_N'])

def judge(pred, truth, tol):
    errors = [abs(a - b) / abs(b) for (a, b) in zip(pred, truth)]
    return {'relative_errors': errors, 'max_relative_error': max(errors), 'gate': 'PASS' if max(errors) <= tol else 'FAIL'}
