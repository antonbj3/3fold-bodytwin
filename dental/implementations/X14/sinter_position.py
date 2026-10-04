"""Observable-specific sintering model. No clinical or total-fit certification."""
import math
import numpy as np

def angular_gap_um(alpha_deg, diameter_mm):
    """Hirano section 4.5 equal-rotation assumption, with exact tangent."""
    if diameter_mm <= 0 or abs(alpha_deg) >= 90:
        raise ValueError('unsupported diameter or angle')
    return 1000 * diameter_mm * math.tan(math.radians(abs(alpha_deg)) / 2)

def distortion_decision(alpha_interval_deg, diameter_mm, other_gap_interval_um, geometry='Hirano_four_unit_FDP', threshold_um=120):
    if geometry != 'Hirano_four_unit_FDP':
        return dict(status='UNKNOWN', reason='geometry transfer unmeasured')
    if alpha_interval_deg is None or other_gap_interval_um is None:
        return dict(status='UNKNOWN', reason='angle uncertainty or measured other-process gap missing')
    (lo, hi) = alpha_interval_deg
    (blo, bhi) = other_gap_interval_um
    if lo > hi or blo < 0 or blo > bhi:
        raise ValueError('invalid interval')
    alo = 0 if lo <= 0 <= hi else min(abs(lo), abs(hi))
    ahi = max(abs(lo), abs(hi))
    (gl, gh) = (angular_gap_um(alo, diameter_mm) + blo, angular_gap_um(ahi, diameter_mm) + bhi)
    status = 'EXCEEDS_DISTORTION_BUDGET' if gl > threshold_um else 'WITHIN_DISTORTION_BUDGET' if gh <= threshold_um else 'UNKNOWN'
    return dict(status=status, conditional_gap_interval_um=[gl, gh], assumption='equal retainer rotation; additive externally supplied other-process bound')

def compensated_stretch(shrinkage_pct, enlargement_factor):
    values = np.asarray(shrinkage_pct, dtype=float)
    if np.any(values < 0) or np.any(values >= 100) or enlargement_factor <= 0:
        raise ValueError('nonphysical shrinkage/enlargement')
    return enlargement_factor * (1 - values / 100)

def deform(points_mm, stretch_xyz, kappa_per_mm=0):
    """Affine directional stretch + small-strain planar beam bending.

    Position-specific kappa needs same-batch calibration. It is not inferred
    from another study's shrinkage tensor or from a scalar material label.
    """
    p = np.asarray(points_mm, dtype=float)
    q = p * np.asarray(stretch_xyz)
    (x, z) = (q[:, 0].copy(), q[:, 2].copy())
    q[:, 0] = x - kappa_per_mm * x * z
    q[:, 2] = z + kappa_per_mm * x * x / 2
    return q

def calibration_design_matrix(span_mm=24, height_mm=5):
    """Four scalar features for log stretches XYZ and signed curvature.

    Three pre/post fiducial baseline changes identify log stretches. Signed
    relative retainer angle identifies kappa in the small-angle beam closure.
    Outputs are dimensionless and radians. Coordinates/units are explicit.
    """
    if span_mm <= 0 or height_mm <= 0:
        raise ValueError('invalid coupon geometry')
    return np.eye(4)

def infer_calibration(pre_lengths_mm, post_lengths_mm, angle_deg, span_mm=24):
    (pre, post) = (np.asarray(pre_lengths_mm), np.asarray(post_lengths_mm))
    if pre.shape != (3,) or post.shape != (3,) or np.any(pre <= 0) or np.any(post <= 0):
        raise ValueError('three positive independent pre/post baselines required')
    return dict(raw_stretch_xyz=(post / pre).tolist(), curvature_per_mm=math.radians(angle_deg) / span_mm, status='CALIBRATED_OBSERVABLES_ONLY', geometry_transfer='UNKNOWN until held-out crown/FDP scan and seating gap')
