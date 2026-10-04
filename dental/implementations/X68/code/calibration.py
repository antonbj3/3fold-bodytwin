"""Research torque ports. This module never infers radial pressure from total torque."""
import math
import numpy as np

def monotone_inverse(points, target, measurement_bounds=None):
    """Outer set for a continuous nonincreasing torque response given calibration intervals.

    Exact means are not a physical uncertainty bound. A caller must supply independently
    established bounds to obtain a physical enclosure; otherwise status is conditional.
    """
    if measurement_bounds is not None:
        if len(measurement_bounds) != len(points):
            raise ValueError('One bound per input point required')
        zipped = sorted(zip(points, measurement_bounds), key=lambda v: v[0][0])
        points = [v[0] for v in zipped]
        measurement_bounds = [v[1] for v in zipped]
    else:
        points = sorted(points)
    ds = np.array([p[0] for p in points])
    ts = np.array([p[1] for p in points])
    if not np.isfinite(target) or not np.all(np.isfinite(ds)) or (not np.all(np.isfinite(ts))):
        raise ValueError('Finite measurements and target required')
    if len(ds) < 2 or len(set(ds)) != len(ds):
        raise ValueError('Need at least two distinct bore diameters')
    if measurement_bounds is None:
        lo = ts
        hi = ts
    else:
        lo = np.array([b[0] for b in measurement_bounds])
        hi = np.array([b[1] for b in measurement_bounds])
        if not (len(lo) == len(ds) and np.all(lo <= hi) and np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
            raise ValueError('Invalid independently specified bounds')
    if any((hi[i] < lo[j] for i in range(len(ds)) for j in range(i + 1, len(ds)))):
        return {'status': 'INCONSISTENT', 'Df_outer_interval_mm': None, 'radial_pressure_MPa': None}
    if target > hi[0] or target < lo[-1]:
        return {'status': 'OUTSIDE_CALIBRATION_RANGE', 'Df_outer_interval_mm': None, 'radial_pressure_MPa': None}
    left = float(ds[0])
    right = float(ds[-1])
    for (d, l, h) in zip(ds, lo, hi):
        if l > target:
            left = max(left, float(d))
        if h < target:
            right = min(right, float(d))
    if left > right:
        return {'status': 'INCONSISTENT', 'Df_outer_interval_mm': None, 'radial_pressure_MPa': None}
    return {'status': 'CONDITIONAL_ON_MEAN_CURVE_MONOTONICITY' if measurement_bounds is None else 'CONDITIONAL_ON_VALID_MEASUREMENT_BOUNDS_AND_MONOTONICITY', 'Df_outer_interval_mm': [left, right], 'width_mm': right - left, 'guaranteed_single_Df': None, 'radial_pressure_MPa': None, 'physical_uncertainty': 'UNKNOWN' if measurement_bounds is None else 'Caller-supplied bounds must cover true local response, not merely mean rounding or an SD'}

def lagrange(x, y, q):
    return sum((y[j] * math.prod(((q - x[k]) / (x[j] - x[k]) for k in range(len(x)) if k != j)) for j in range(len(x))))

def derivative_enclosure_area_poly(poly, D, dlo, dhi):
    """For quadratic T(u), u=D²-d²: T'(d)= -2d*(b+2cD²-2cd²).
    Cubic derivative reaches extrema at endpoints and at derivative' roots.
    """
    (c, b, a) = poly
    A = b + 2 * c * D * D
    candidates = [dlo, dhi]
    if c != 0:
        q = A / (6 * c)
        if q >= 0 and dlo <= math.sqrt(q) <= dhi:
            candidates.append(math.sqrt(q))
    values = [-2 * d * (A - 2 * c * d * d) for d in candidates]
    return {'dT_dDf_Ncm_per_mm': [min(values), max(values)], 'extremum_candidates_mm': candidates, 'scope': 'Exact extrema for fitted polynomial; physical approximation/remainder bound MISSING'}
