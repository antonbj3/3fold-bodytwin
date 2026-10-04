"""Two-mode mechanical measurement port; local bone reference is mandatory.

No dental constitutive calibration is implied by this linear operator.
The outward interval uses explicit IEEE-754 rounding and sign-preserving boxes.
"""
import math
import numpy as np

def project_intervals(qlo, qhi, z):
    """Enclose u+z*theta for independent axis-aligned q intervals."""
    vals = []
    for zz in z:
        if zz >= 0:
            (a, b) = (qlo[1], qhi[1])
        else:
            (a, b) = (qhi[1], qlo[1])
        lo = math.nextafter(qlo[0] + math.nextafter(zz * a, -math.inf), -math.inf)
        hi = math.nextafter(qhi[0] + math.nextafter(zz * b, math.inf), math.inf)
        vals.append([lo, hi])
    return np.array(vals)

def maximum_abs_interval(region_intervals):
    lo = max((0 if a <= 0 <= b else min(abs(a), abs(b)) for (a, b) in region_intervals))
    hi = max((max(abs(a), abs(b)) for (a, b) in region_intervals))
    return [math.nextafter(lo, -math.inf), math.nextafter(hi, math.inf)]

def local_threshold(interval_um, local_bone_reference, biological=False):
    if biological:
        return 'UNKNOWN_BIOLOGICAL_CLOSURE'
    if not local_bone_reference:
        return 'UNKNOWN_MISSING_LOCAL_BONE_REFERENCE'
    (a, b) = interval_um
    if b < 50:
        return 'BELOW_50_FOR_SAMPLED_REGIONS'
    if a > 150:
        return 'ABOVE_150_FOR_SAMPLED_REGIONS'
    return 'UNRESOLVED_50_150_BAND'

def identify(z_markers, wrenches, marker_displacements, bone_reference_displacements=None):
    if bone_reference_displacements is None:
        return {'status': 'UNKNOWN_MISSING_LOCAL_BONE_REFERENCE'}
    H = np.column_stack([np.ones(len(z_markers)), np.asarray(z_markers)])
    W = np.asarray(wrenches, dtype=float)
    Y = np.asarray(marker_displacements) - np.asarray(bone_reference_displacements)
    if np.linalg.matrix_rank(H) < 2 or np.linalg.matrix_rank(W) < 2:
        return {'status': 'UNKNOWN_RANK_DEFICIENT_MEASUREMENT'}
    Q = np.linalg.solve(H, Y)
    C = Q @ np.linalg.inv(W)
    sym_error = np.linalg.norm(C - C.T) / max(np.linalg.norm(C), 1e-300)
    if sym_error > 1e-10:
        return {'status': 'REJECT_NONRECIPROCAL_CALIBRATION', 'symmetry_error': float(sym_error)}
    if min(np.linalg.eigvalsh(C)) <= 0:
        return {'status': 'REJECT_NONPOSITIVE_COMPLIANCE'}
    return {'status': 'IDENTIFIED_LINEAR_PORT', 'C': C.tolist(), 'K': np.linalg.inv(C).tolist(), 'symmetry_error': float(sym_error)}

def stiffness(z, k):
    X = np.column_stack([np.ones(len(z)), np.asarray(z)])
    return X.T @ (np.asarray(k)[:, None] * X)
