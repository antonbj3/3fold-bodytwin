"""Geometry-only decision guards. No clinical/physical calibration is implied.

Positive edge_error_mm decreases clearance (canal/root/pulp grows or support shrinks).
All inequalities accept equality. A flip distance is an infimum: crossing a passing
boundary requires a strictly greater adverse displacement.
"""
import numpy as np

def guard(distance_mm, threshold_mm, budget_mm=None, numerical_budget_mm=0.0):
    (d, t) = (float(distance_mm), float(threshold_mm))
    if not np.isfinite(d) or not np.isfinite(t):
        raise ValueError('Non-finite distance/threshold is not admissible evidence')
    if numerical_budget_mm < 0 or (budget_mm is not None and budget_mm < 0):
        raise ValueError('Error budgets must be nonnegative')
    m = d - t
    b = None if budget_mm is None else float(budget_mm) + float(numerical_budget_mm)
    status = 'UNCALIBRATED' if b is None else 'SCENARIO_STABLE_PASS' if m >= b else 'SCENARIO_STABLE_FAIL' if m < -b else 'ALREADY_UNCERTAIN_IN_SCENARIO'
    return {'distance_mm': d, 'threshold_mm': t, 'signed_margin_mm': m, 'nominal_pass': bool(m >= 0), 'flip_infimum_mm': abs(m), 'direction_to_flip': 'CLEARANCE_DECREASING' if m >= 0 else 'CLEARANCE_INCREASING', 'equality_accepted': True, 'scenario_budget_mm': b, 'scenario_status': status, 'physical_status': 'UNKNOWN_NO_MATCHED_LOCAL_TRUTH'}

def scalar_boundary_control(d, t):
    """Numerical predicate replay; independent bisection, no margin formula."""
    (lo, hi) = (-max(abs(d) + abs(t) + 1, 1), max(abs(d) + abs(t) + 1, 1))
    for _ in range(64):
        x = (lo + hi) / 2
        if d - x >= t:
            lo = x
        else:
            hi = x
    return (lo + hi) / 2

def fixed_design_record(distances, thresholds, budgets=None, numerical_budget_mm=0.0):
    missing = [k for (k, v) in distances.items() if v is None]
    rr = {k: guard(v, thresholds[k], None if budgets is None else budgets.get(k), numerical_budget_mm) for (k, v) in distances.items() if v is not None}
    if not rr:
        return {'guards': {}, 'nominal_pass': None, 'adverse_radius_mm': None, 'physical_status': 'UNKNOWN_MISSING_GEOMETRY'}
    b = min(rr, key=lambda k: rr[k]['signed_margin_mm'])
    nominal = False if any((not q['nominal_pass'] for q in rr.values())) else None if missing else True
    return {'guards': rr, 'missing_regions': missing, 'nominal_pass': nominal, 'adverse_radius_mm': None if missing else max(0.0, rr[b]['signed_margin_mm']), 'binding_region': b, 'physical_status': 'UNKNOWN_NO_MATCHED_LOCAL_TRUTH'}
