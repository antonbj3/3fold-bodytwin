"""Exact fixed-design clearance for complete masks represented as closed voxels.

The API accepts any segmentation model's reference/prediction masks or voxel
coordinates on the same axis-aligned physical grid. It never infers anatomy
truth, label coverage or registration accuracy from the data array.
"""
import numpy as np
from scipy.spatial import cKDTree

def validate_points(a):
    a = np.asarray(a)
    if a.ndim != 2 or a.shape[1] != 3 or len(a) == 0:
        raise ValueError('Require nonempty N x 3 voxel coordinates; empty hazard is UNKNOWN')
    if not np.issubdtype(a.dtype, np.integer):
        raise ValueError('Voxel coordinates must be integers')
    return a.astype(np.int64, copy=False)

def validate_spacing(spacing):
    s = np.asarray(spacing, dtype=float)
    if s.shape != (3,) or not np.all(np.isfinite(s)) or (not np.all(s > 0)):
        raise ValueError('Spacing must be three finite positive mm values')
    return s

def local_cube_field(body, hazard, spacing=(0.3, 0.3, 0.3)):
    """Per-hazard-cube distance to all body cubes, with a nearest witness.

    Cube min distance <= min centre distance. A centre farther than
    candidate_upper + norm(spacing) cannot give a better cube pair. Every
    potential better pair is therefore enumerated; no nearest-k closure.
    """
    X = validate_points(body)
    Y = validate_points(hazard)
    s = validate_spacing(spacing)
    tree = cKDTree(X * s)
    (centre, nearest) = tree.query(Y * s, workers=1)
    base_delta = (X[nearest] - Y) * s
    base_sep = np.sign(base_delta) * np.maximum(abs(base_delta) - s, 0)
    upper = np.linalg.norm(base_sep, axis=1)
    dist = np.empty(len(Y))
    witness = np.empty_like(Y)
    direction = np.zeros((len(Y), 3))
    paircount = 0
    for lo in range(0, len(Y), 256):
        ids = tree.query_ball_point(Y[lo:lo + 256] * s, upper[lo:lo + 256] + np.linalg.norm(s) + 1e-10, workers=1)
        for (k, ii) in enumerate(ids, lo):
            ii = np.asarray(ii, dtype=int)
            paircount += len(ii)
            delta = (X[ii] - Y[k]) * s
            sep = np.sign(delta) * np.maximum(abs(delta) - s, 0)
            dd = np.einsum('ij,ij->i', sep, sep)
            j = int(np.argmin(dd))
            dist[k] = np.sqrt(dd[j])
            witness[k] = X[ii[j]]
            if dist[k] > 0:
                direction[k] = sep[j] / dist[k]
    return {'distance_mm': dist, 'body_witness_voxel': witness, 'approach_direction_zyx': direction, 'candidate_pairs': paircount, 'all_potential_better_pairs_covered': True}

def direct_cube_field(body, hazard, spacing=(0.3, 0.3, 0.3)):
    """Independent blockwise all-pairs control. No tree or covering radius."""
    X = validate_points(body)
    Y = validate_points(hazard)
    s = validate_spacing(spacing)
    out = np.full(len(Y), np.inf)
    for i in range(0, len(Y), 32):
        for j in range(0, len(X), 2048):
            dd = np.maximum(abs((Y[i:i + 32, None] - X[None, j:j + 2048]) * s) - s, 0)
            v = np.sum(dd * dd, axis=2).min(axis=1)
            out[i:i + 32] = np.minimum(out[i:i + 32], np.sqrt(v))
    return out

def decision_report(field, hazard, threshold_mm, budget_mm=None, spacing=(0.3, 0.3, 0.3), patch_mm=1.2):
    Y = validate_points(hazard)
    s = validate_spacing(spacing)
    d = np.asarray(field['distance_mm'])
    t = float(threshold_mm)
    if len(d) != len(Y) or not np.all(np.isfinite(d)) or (not np.isfinite(t)) or (t < 0):
        raise ValueError('Nonnegative finite threshold and complete finite field required')
    if patch_mm <= 0:
        raise ValueError('Patch scale must be positive')
    i = int(np.argmin(d))
    m = d - t
    out = {'gap_mm': float(d[i]), 'threshold_mm': t, 'nominal_pass': bool(d.min() >= t), 'adverse_flip_infimum_mm': max(0.0, float(m.min())), 'critical_hazard_voxel': Y[i].tolist(), 'critical_body_voxel': field['body_witness_voxel'][i].tolist(), 'adverse_direction_zyx': field['approach_direction_zyx'][i].tolist(), 'surface_cells': len(Y), 'geometry_coverage': 'COMPLETE_SUBMITTED_MASK', 'physical_status': 'UNKNOWN_NO_INDEPENDENT_ANATOMY', 'signed_distance_convention': 'distance between closed voxel cubes; zero for overlap/touch'}
    if budget_mm is not None:
        B = float(budget_mm)
        if not np.isfinite(B) or B < 0:
            raise ValueError('Nonnegative finite bound required')
        critical = m <= B
        pk = np.floor(Y * s / patch_mm).astype(int)
        patches = {tuple(v) for v in pk[critical]}
        measured = np.array([tuple(v) in patches for v in pk], bool)
        out.update(scenario_budget_mm=B, dangerous_cells=int(critical.sum()), dangerous_cell_fraction=float(critical.mean()), measurement_patches=len(patches), measurement_cells=int(measured.sum()), measurement_cell_fraction=float(measured.mean()), unmeasured_min_margin_mm=float(m[~measured].min()) if (~measured).any() else None, patch_indices_zyx=sorted([[int(x) for x in k] for k in patches]), patch_contract='All dangerous source cells measured independently; every unmeasured source cell has displacement <=B. New/missing hazards and design error excluded.')
    return out

def overlap_metrics(reference, prediction, spacing=(0.3, 0.3, 0.3)):
    """Dice and symmetric surface HD95 on supplied surface cell centres.

    Inputs are solid-voxel coordinates for Dice. Compute surface externally
    when reporting surface HD95 for nonsparse masks.
    """
    A = {tuple(v) for v in validate_points(reference)}
    B = {tuple(v) for v in validate_points(prediction)}
    return {'Dice': 2 * len(A & B) / (len(A) + len(B)), 'intersection_voxels': len(A & B), 'reference_voxels': len(A), 'prediction_voxels': len(B)}

def hd95_points(A, B, spacing=(0.3, 0.3, 0.3)):
    s = validate_spacing(spacing)
    x = cKDTree(A * s).query(B * s, workers=1)[0]
    y = cKDTree(B * s).query(A * s, workers=1)[0]
    return float(np.percentile(np.r_[x, y], 95))

def paired_evaluate(body, reference_voxel, prediction_voxel, spacing, threshold_mm):
    """Query both complete solid masks, including containment/overlap.

    Supply all solid hazard cells. Surface-only lists cannot distinguish a body
    contained in a solid hazard from a body in a cavity shell.
    """
    r = decision_report(local_cube_field(body, reference_voxel, spacing), reference_voxel, threshold_mm, spacing=spacing)
    p = decision_report(local_cube_field(body, prediction_voxel, spacing), prediction_voxel, threshold_mm, spacing=spacing)
    r['queried_solid_cells'] = r.pop('surface_cells')
    p['queried_solid_cells'] = p.pop('surface_cells')
    return {'reference': r, 'prediction': p, 'distance_error_mm': p['gap_mm'] - r['gap_mm'], 'decision_flip': bool(r['nominal_pass'] != p['nominal_pass']), 'compared_truth': 'USER_SUPPLIED_REFERENCE_GEOMETRY', 'physical_status': 'UNKNOWN_UNLESS_REFERENCE_INDEPENDENTLY_VALIDATED'}

def verify_patch_selection(distance_mm, threshold_mm, budget_mm, selected):
    """Can reject omission of any potentially decision-changing source cell."""
    d = np.asarray(distance_mm, float)
    selected = np.asarray(selected, bool)
    if d.shape != selected.shape or not np.all(np.isfinite(d)):
        return False
    required = d - float(threshold_mm) <= float(budget_mm)
    return bool(np.all(selected[required]))

def witness_subset_report(witness_distances_mm, threshold_mm):
    """A supported subset gives upper witnesses, no lower bound on full minimum."""
    v = np.asarray(witness_distances_mm, float)
    if len(v) == 0 or not np.all(np.isfinite(v)) or np.any(v < 0):
        return {'coverage': 'UNKNOWN', 'all_variants_close_witness': None, 'can_certify_full_mask_above_threshold': False}
    return {'coverage': 'SUPPORTED_SUBSET', 'all_variants_close_witness': bool(v.max() < threshold_mm), 'can_certify_full_mask_above_threshold': False, 'physical_status': 'UNKNOWN'}

def verify_witness_subset_claim(witness_distances_mm, threshold_mm, claim):
    expected = witness_subset_report(witness_distances_mm, threshold_mm)
    return all((claim.get(k) == v for (k, v) in expected.items()))
