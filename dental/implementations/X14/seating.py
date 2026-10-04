"""Planar four-margin metrology → nonpenetrating rigid seating.

This is a small-angle 2D lower model. Full intaglio interference, cement,
elasticity and 3D registration remain outside its empirical validity.
"""
from itertools import combinations
import numpy as np
from scipy.optimize import linprog

def constraints(x_mm, z_mm):
    (x, z) = (np.asarray(x_mm, float), np.asarray(z_mm, float))
    if x.ndim != 1 or z.shape != x.shape or len(x) < 3 or (np.ptp(x) <= 0):
        raise ValueError('at least three finite points spanning x required')
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(z)):
        raise ValueError('finite observations required')
    a = np.vstack([np.column_stack([-np.ones(len(x)), -x, np.zeros(len(x))]), np.column_stack([np.ones(len(x)), x, -np.ones(len(x))]), [0.0, 0.0, -1.0]])
    b = np.concatenate([z, -z, [0.0]])
    return (x, z, a, b)

def seat_vertices(x_mm, z_mm):
    (x, z, a, b) = constraints(x_mm, z_mm)
    best = None
    for subset in combinations(range(len(b)), 3):
        mat = a[list(subset)]
        if abs(np.linalg.det(mat)) < 1e-12:
            continue
        q = np.linalg.solve(mat, b[list(subset)])
        if np.max(a @ q - b) <= 1e-10 and (best is None or q[2] < best[2]):
            best = q
    if best is None:
        raise ValueError('no feasible vertex found')
    gaps = z + best[0] + best[1] * x
    return dict(max_gap_mm=float(max(0.0, best[2])), translation_mm=float(best[0]), rigid_slope=float(best[1]), gaps_mm=gaps.tolist())

def seat_lp(x_mm, z_mm):
    (x, z, a, b) = constraints(x_mm, z_mm)
    result = linprog([0.0, 0.0, 1.0], A_ub=a, b_ub=b, bounds=[(None, None), (None, None), (0, None)], method='highs')
    if not result.success:
        raise ValueError(result.message)
    return dict(max_gap_mm=float(result.fun), translation_mm=float(result.x[0]), rigid_slope=float(result.x[1]), gaps_mm=(z + result.x[0] + result.x[1] * x).tolist())

def profile(alpha_deg, fraction, diameter_left_mm=7.0, diameter_right_mm=11.0, span_mm=24.0):
    """Declared synthetic profile, not a reproduction of source STL.

    Both retainer margin centers have z=0. The source relative angle is
    preserved while its unobserved allocation relative to the chord varies.
    """
    (left, right) = (-span_mm / 2, span_mm / 2)
    x = np.array([left - diameter_left_mm / 2, left + diameter_left_mm / 2, right - diameter_right_mm / 2, right + diameter_right_mm / 2])
    slopes = np.tan(np.radians([-fraction * alpha_deg, (1 - fraction) * alpha_deg]))
    z = np.array([-diameter_left_mm / 2 * slopes[0], diameter_left_mm / 2 * slopes[0], -diameter_right_mm / 2 * slopes[1], diameter_right_mm / 2 * slopes[1]])
    return (x, z)

def relative_angle_deg(x_mm, z_mm):
    (x, z) = (np.asarray(x_mm), np.asarray(z_mm))
    (slope1, slope2) = ((z[1] - z[0]) / (x[1] - x[0]), (z[3] - z[2]) / (x[3] - x[2]))
    return float(np.degrees(np.arctan(slope2) - np.arctan(slope1)))

def canonical_margins(x_mm, z_mm):
    """Two shape invariants modulo vertical translation and small rigid tilt.

    A measured relative angle alone supplies only one combination. Its
    missing chord-relative tilt must not be silently set to a symmetric value.
    """
    (x, z) = (np.asarray(x_mm, float), np.asarray(z_mm, float))
    if x.shape != (4,) or z.shape != (4,) or np.any(np.diff(x) <= 0):
        raise ValueError('four ordered margin coordinates required')
    centers = (x[:2].mean(), x[2:].mean())
    heights = (z[:2].mean(), z[2:].mean())
    chord = (heights[1] - heights[0]) / (centers[1] - centers[0])
    widths = (x[1] - x[0], x[3] - x[2])
    slopes = ((z[1] - z[0]) / widths[0] - chord, (z[3] - z[2]) / widths[1] - chord)
    normalized = np.array([-widths[0] * slopes[0] / 2, widths[0] * slopes[0] / 2, -widths[1] * slopes[1] / 2, widths[1] * slopes[1] / 2])
    return (normalized, dict(chord_relative_slopes=list(slopes), retainer_widths_mm=list(widths), center_span_mm=centers[1] - centers[0]))

def measured_seating_bound(x_mm, z_mm, height_error_mm=None):
    (canonical, invariants) = canonical_margins(x_mm, z_mm)
    answer = seat_vertices(x_mm, canonical)
    if height_error_mm is None:
        bound = None
        status = 'UNKNOWN_MEASUREMENT_ERROR'
    elif height_error_mm < 0:
        raise ValueError('nonnegative height error bound required')
    else:
        bound = [max(0, answer['max_gap_mm'] - 2 * height_error_mm), answer['max_gap_mm'] + 2 * height_error_mm]
        status = 'CONDITIONAL_PLANAR_BOUND'
    return dict(**answer, **invariants, measurement_error_bound_mm=height_error_mm, max_gap_bound_mm=bound, status=status, assumptions='fixed x; bounded z error; small rigid tilt; planar nonpenetration; no 3D intaglio contact')
