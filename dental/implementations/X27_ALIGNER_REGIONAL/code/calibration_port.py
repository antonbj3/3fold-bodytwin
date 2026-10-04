"""Linear regional compliance identification with OBSERVED shell pose.

This port does not invent a shell pose, aligner gap, sensor zero, specimen ID,
or apex-to-centroid transformation when an experiment lacks them.
"""
import numpy as np
from scipy.optimize import nnls
REQUIRED = ['specimen_id', 'tooth_id', 'force_N', 'moment_Nmm', 'moment_origin', 'tooth_pose', 'shell_housing_pose', 'gap_mm', 'region_signatures', 'source_kind', 'baseline_seated', 'elapsed_s', 'temperature_C', 'optical_housing_rigidity_residual_mm']

def matrix(row):
    a = np.asarray(row['region_signatures'], float)
    pose = np.asarray(row['tooth_pose']) - np.asarray(row['shell_housing_pose'])
    compression = a @ pose - np.asarray(row['gap_mm'])
    return (-(a.T * np.maximum(compression, 0)), compression)

def fit(rows, length_scale_mm, tol=1e-07):
    for row in rows:
        missing = [k for k in REQUIRED if k not in row or row[k] is None]
        if missing:
            return dict(status='NEEDS_MATCHED_MEASUREMENT', missing=missing)
        if not row['baseline_seated']:
            return dict(status='NEEDS_SEATED_ZERO')
        if row['moment_origin'] != 'crown_centroid':
            return dict(status='NEEDS_EXACT_MOMENT_TRANSPORT', origin=row['moment_origin'])
        if row['optical_housing_rigidity_residual_mm'] > 0.05:
            return dict(status='RIGID_HOUSING_CLOSURE_REFUTED')
    if len({row['specimen_id'] for row in rows}) != 1:
        return dict(status='REJECTED_SPECIMEN_MIXING')
    if len({row['tooth_id'] for row in rows}) != 1:
        return dict(status='REJECTED_TOOTH_MIXING')
    scales = np.array([1, 1, 1, 1 / length_scale_mm, 1 / length_scale_mm, 1 / length_scale_mm])
    X = np.vstack([matrix(row)[0] * scales[:, None] for row in rows])
    y = np.concatenate([np.r_[row['force_N'], row['moment_Nmm']] * scales for row in rows])
    sv = np.linalg.svd(X, compute_uv=False)
    rank = int(np.sum(sv > sv[0] * 1e-09)) if sv[0] > 0 else 0
    if rank < X.shape[1]:
        null = np.linalg.svd(X, full_matrices=True)[2][rank:].T
        return dict(status='UNIDENTIFIABLE_REGIONAL_COMPLIANCE', rank=rank, parameters=X.shape[1], nullspace=null.tolist(), singular_values=sv.tolist())
    k = nnls(X, y, maxiter=1000)[0]
    residual = float(np.linalg.norm(X @ k - y, np.inf))
    return dict(status='CONDITIONAL_CALIBRATION' if residual <= tol else 'CONTACT_MODEL_REFUTED', stiffness_N_per_mm=k.tolist(), rank=rank, parameters=X.shape[1], singular_values=sv.tolist(), condition_number=float(sv[0] / sv[-1]), scaled_residual_N=residual, empirical_validation='UNKNOWN_REQUIRES_INDEPENDENT_SAME_SPECIMEN_HOLDOUT_AND_HOUSING_RIGIDITY_CHECK')

def predict(row, k):
    return matrix(row)[0] @ np.asarray(k)
