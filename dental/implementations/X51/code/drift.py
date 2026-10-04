"""Signed reference/coupon observation model. Physical transfer is a lab contract."""
import thread_budget
import numpy as np
from scipy.stats import norm

def limits(sigma=10, n=5, weeks=52, alpha=0.05):
    if sigma <= 0 or n < 1 or weeks < 1 or (not 0 < alpha < 1):
        raise ValueError('Invalid chart inputs')
    z = float(np.nextafter(norm.isf(alpha / (4 * weeks)), np.inf))
    while weeks * 4 * norm.sf(z) > alpha:
        z = float(np.nextafter(z, np.inf))
    se = np.array([sigma / np.sqrt(n), np.sqrt(2) * sigma / np.sqrt(n)])
    return {'z': z, 'reference_limit_um': float(z * se[0]), 'process_limit_um': float(z * se[1]), 'se_um': se.tolist(), 'annual_expected_channel_signals_bound': alpha, 'resolution': 'PER_SURFACE_REGION', 'scope': 'known Gaussian marginal variance, fixed baseline; alpha allocated to52 weeks x2 channels x2 tails'}

def paired(reference, coupon):
    a = np.asarray(reference, float)
    b = np.asarray(coupon, float)
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        raise ValueError('Finite signed metrology required')
    return np.stack([a, b - a], axis=-1)

def interval_hold(g0, dhat, measurement_error, ramp_bound, days_to_next, rejection=120):
    """Rigorous max enclosure for this additive gap model, NOT calibrated dental mechanics."""
    args = [g0, dhat, measurement_error, ramp_bound, days_to_next, rejection]
    if not all(np.isfinite(args)) or measurement_error < 0 or ramp_bound < 0 or (days_to_next < 0):
        raise ValueError('Invalid enclosure')
    upper = g0 + dhat + measurement_error + ramp_bound * days_to_next
    return dict(upper_gap_um=upper, verdict='HOLD' if upper > rejection else 'MODEL_WITHIN_LIMIT', scope='bounded-error, bounded-slope additive model; lab calibration required', resolution='PER_SURFACE_REGION')
