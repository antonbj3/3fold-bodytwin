"""Installed tip displacement port. Exact rational interpolation bounds.

The curvature bound is REQUIRED. Two measurements alone do not justify affine
response. Calibrates only within the measured load interval and named pose.
"""
from fractions import Fraction as Q
import math

def q(x):
    if isinstance(x, Q):
        return x
    if isinstance(x, bool):
        raise ValueError('boolean is not a magnitude')
    if not math.isfinite(float(x)):
        raise ValueError('nonfinite magnitude')
    return Q(str(x))

def displacement_bound(port, F, pose_id):
    if port.get('pose_id') != pose_id:
        return dict(status='UNKNOWN_POSE_MISMATCH')
    M = port.get('curvature_bound_mm_N2')
    if M is None:
        return dict(status='UNKNOWN_NO_NONLINEAR_ENCLOSURE')
    M = q(M)
    if M < 0:
        raise ValueError('negative curvature bound')
    (a, b) = map(q, port['force_levels_N'])
    f = q(F)
    if a >= b:
        raise ValueError('force levels must increase')
    if f < a or f > b:
        return dict(status='UNKNOWN_EXTRAPOLATION')
    lo = port['tip_displacement_lower_mm']
    hi = port['tip_displacement_upper_mm']
    if any((q(l) < 0 or q(h) < q(l) for (l, h) in zip(lo, hi))):
        raise ValueError('bad measurement intervals')
    w = (f - a) / (b - a)
    rest = M * (f - a) * (b - f) / 2
    lower = max(Q(0), (1 - w) * q(lo[0]) + w * q(lo[1]) - rest)
    upper = (1 - w) * q(hi[0]) + w * q(hi[1]) + rest
    return dict(status='BOUNDED_WITHIN_DECLARED_CALIBRATION', lower_mm=float(lower), upper_mm=float(upper), lower_rational=str(lower), upper_rational=str(upper), remainder_rational=str(rest), rigorous_arithmetic_enclosure=True, curvature_source=port.get('curvature_source', 'UNKNOWN'), physical_status='UNKNOWN_UNQUALIFIED_CURVATURE_OR_MEASUREMENT' if port.get('data_kind') != 'independent_measurement' else 'CONDITIONAL_MEASUREMENT_BOUND', scope='Tip displacement magnitude in named installed pose and force interval; no actual force/CAM/whole-surface guarantee')

def calibration_controls():
    port = dict(pose_id='fixture_pose', force_levels_N=[0.5, 1.0], tip_displacement_lower_mm=[0.0015, 0.0035], tip_displacement_upper_mm=[0.0025, 0.0045], curvature_bound_mm_N2=0.0008, curvature_source='our_own_fixture declared derivative bound', data_kind='our_own_fixture')
    b = displacement_bound(port, 0.75, 'fixture_pose')
    exact = Q(1, 1000) * Q(3, 1) + Q(1, 40000)
    assert Q(b['upper_rational']) == Q(7, 2000) + Q(1, 40000)
    no = dict(port, curvature_bound_mm_N2=None)
    out = dict(baseline=b, injected_missing_curvature_rejected=displacement_bound(no, 0.75, 'fixture_pose')['status'].startswith('UNKNOWN'), injected_extrapolation_rejected=displacement_bound(port, 2, 'fixture_pose')['status'] == 'UNKNOWN_EXTRAPOLATION', injected_pose_rejected=displacement_bound(port, 0.75, 'other')['status'] == 'UNKNOWN_POSE_MISMATCH', injected50um_rejected=displacement_bound(dict(port, tip_displacement_lower_mm=[0.05, 0.05], tip_displacement_upper_mm=[0.05, 0.05]), 0.75, 'fixture_pose')['upper_mm'] > 0.02)

    def A(f):
        return f / 32

    def B(f):
        return A(f) + (f - Q(1, 2)) * (1 - f) / 32
    endpoints = [Q(1, 2), Q(1)]
    identity = max((abs(A(f) - B(f)) for f in endpoints))
    f = Q(3, 4)
    difference = B(f) - A(f)
    M = Q(1, 16)
    enclosure = M * (f - Q(1, 2)) * (1 - f) / 2
    out['sufficiency'] = dict(summary='two force-level displacement endpoints', summary_resolution='PER_POINT installed pose', identity_error_mm=float(identity), identity_rational=str(identity), midpoint_difference_mm=float(difference), midpoint_difference_rational=str(difference), rigorous_remainder_rational=str(enclosure), minimum_extension='Independent curvature bound, or measured full force-response envelope; no affine extrapolation', curve_A='u=F/32 mm', curve_B='u=F/32+(F-1/2)*(1-F)/32 mm', curvature_upper_mm_N2=float(M))
    assert identity == 0 and difference == enclosure and (difference > Q(1, 1000))
    out['injected_understated_enclosure_rejected'] = difference > enclosure / 2
    return out
