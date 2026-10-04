"""Three-height acquisition for one crown action under a linear full-contact closure."""
import numpy as np
from height_force import array, nonnegative_scalar

def validate(s):
    eps = nonnegative_scalar(s['raw_channel_error_N'], 'raw_channel_error_N')
    a = nonnegative_scalar(s['probe_height_mm'], 'probe_height_mm')
    if a <= 0:
        raise ValueError('Positive measured height step required')
    A = array(s['A'], 'A')
    b = array(s['baseline_force_N'], 'baseline_force_N')
    j = array(s['height_column_N_per_mm'], 'height_column_N_per_mm')
    w = array(s['w_N'], 'w_N')
    if A.shape != (len(b), 3) or j.shape != b.shape or w.shape != (3,):
        raise ValueError('Invalid corridor dimensions')
    if np.max(np.abs(A.T @ b - w)) > 1e-07 or np.max(np.abs(A.T @ j)) > 1e-07:
        raise ValueError('Measured column must preserve wrench balance')
    return (A, b, j, w, a, eps)

def error_contract(s, b, j, a, eps, h):
    ds = nonnegative_scalar(s.get('probe_step_error_mm', 0), 'probe_step_error_mm')
    dh = nonnegative_scalar(s.get('query_height_error_mm', 0), 'query_height_error_mm')
    if ds >= a:
        raise ValueError('Probe-height uncertainty must be smaller than its step')
    eta = np.sqrt(len(b)) * eps
    beta = (eta + ds * np.linalg.norm(j)) / (a - ds)
    radius = float(eta + abs(h) * beta + (np.linalg.norm(j) + beta) * dh)
    return (eta, beta, a - ds - dh, radius)

def predict(s, h, geometry_sha256=None, w_N=None):
    (A, b, j, w, a, eps) = validate(s)
    (eta, beta, corridor, radius) = error_contract(s, b, j, a, eps, h)
    if not np.isfinite(h):
        raise ValueError('Height must be finite')
    if geometry_sha256 is not None and geometry_sha256 != s['geometry_sha256']:
        raise ValueError('Geometry identity mismatch')
    if w_N is not None and (not np.array_equal(np.array(w_N), w)):
        return {'status': 'UNKNOWN', 'reason': 'Changed wrench needs a measured load column', 'force_interval_N': None}
    if not s.get('linear_full_contact_closure'):
        return {'status': 'UNKNOWN', 'reason': 'No validated linear/rate/history closure', 'force_interval_N': None}
    if corridor <= 0:
        return {'status': 'UNKNOWN', 'reason': 'Height metrology consumes the entire corridor', 'force_interval_N': None}
    Rmax = error_contract(s, b, j, a, eps, corridor)[3]
    if min(np.min(b - j * corridor - Rmax), np.min(b + j * corridor - Rmax)) <= 0:
        return {'status': 'UNKNOWN', 'reason': 'Full-contact corridor not certified; needs full energy or a new loaded measurement', 'force_interval_N': None}
    if abs(h) > corridor + 1e-12:
        return {'status': 'UNKNOWN', 'reason': 'Outside the measured positive contact corridor', 'force_interval_N': None}
    f = b + j * h
    return {'status': 'CONDITIONAL_CORRIDOR_BOUND', 'force_N': f.tolist(), 'force_interval_N': np.c_[f - radius, f + radius].tolist(), 'joint_l2_error_radius_N': radius, 'height_change_mm': float(h), 'fdi': s['action_fdi'], 'predicted_fdi': s['predicted_fdi'], 'claim_type': 'capability', 'resolution': 'PER_TOOTH', 'enclosure': 'Rigorous analytic L2 bound under linear full-contact closure; no formal float enclosure', 'physical_validity': 'UNKNOWN; simulated acquisition and prospective metrology'}

def inverse(s, target_N, tolerance_N=1):
    (A, b, j, w, a, eps) = validate(s)
    target = nonnegative_scalar(target_N, 'target_N')
    tol = nonnegative_scalar(tolerance_N, 'tolerance_N')
    i = s['predicted_fdi'].index(s['action_fdi'])
    (eta, beta, corridor, _) = error_contract(s, b, j, a, eps, 0)
    if j[i] - beta <= 0:
        return {'status': 'UNKNOWN', 'reason': 'Target slope is not identified as positive'}
    h = float((target - b[i]) / j[i])
    ans = predict(s, h)
    if ans.get('force_interval_N') is None:
        return ans
    endpoints = [(target - v) / slope for v in [b[i] - eta, b[i] + eta] for slope in [j[i] - beta, j[i] + beta]]
    hb = [float(min(endpoints)), float(max(endpoints))]
    (lo, hi) = ans['force_interval_N'][i]
    inside = bool(-corridor <= hb[0] <= hb[1] <= corridor and lo >= target - tol and (hi <= target + tol))
    return {'status': 'CONDITIONAL_TARGET_BAND' if inside else 'UNKNOWN', 'fdi': s['action_fdi'], 'target_N': target, 'nominal_height_change_mm': h, 'possible_exact_target_height_bracket_mm': hb, 'at_selected_height': ans, 'physical_validity': 'UNKNOWN', 'resolution': 'PER_TOOTH'}
