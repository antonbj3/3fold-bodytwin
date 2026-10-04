"""CLI: arbitrary model -> fixed-design decision error, not clinical advice.

NPZ inputs, same z/y/x physical grid:
  body_voxel, reference_hazard_voxel, prediction_hazard_voxel: integer N x 3
  spacing_zyx_mm: positive length-3 vector.
Alternatively supply body_mask, reference_hazard_mask, prediction_hazard_mask
as identically shaped 3-D binary masks and spacing_zyx_mm.
The reference may be annotated or independently measured; its provenance is
reported explicitly and never promoted to physical truth by this evaluator.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import argparse, json
import numpy as np
from scipy import ndimage as ndi
from voxel_decisions import paired_evaluate, overlap_metrics, hd95_points
from assessment_bounds import reference_interval

def mask_surface(v):
    return np.argwhere(v & ~ndi.binary_erosion(v))

def surface_of_coordinates(points):
    lo = points.min(0) - 1
    hi = points.max(0) + 2
    if int(np.prod(hi - lo)) > 100000000:
        raise ValueError('Coordinate bounding box >100M cells; supply a crop or use direct API')
    v = np.zeros(tuple(hi - lo), bool)
    v[tuple((points - lo).T)] = True
    return mask_surface(v) + lo

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', required=True)
    ap.add_argument('--threshold-mm', type=float, required=True)
    ap.add_argument('--reference-kind', choices=['annotation', 'independent_measurement'], default='annotation')
    ap.add_argument('--output')
    ap.add_argument('--bound-manifest', help='JSON complete-shape hard bounds; missing bounds remain UNKNOWN')
    a = ap.parse_args()
    f = np.load(a.input, allow_pickle=False)
    s = f['spacing_zyx_mm']
    if 'body_mask' in f:
        b = f['body_mask']
        r = f['reference_hazard_mask']
        p = f['prediction_hazard_mask']
        if b.ndim != 3 or r.shape != b.shape or p.shape != b.shape:
            raise ValueError('Mask grids differ')
        if any((not np.all((v == 0) | (v == 1)) for v in (b, r, p))):
            raise ValueError('Masks must be binary')
        X = np.argwhere(b)
        R = np.argwhere(r)
        Q = np.argwhere(p)
        yr = mask_surface(r.astype(bool))
        yq = mask_surface(p.astype(bool))
    else:
        X = f['body_voxel']
        R = f['reference_hazard_voxel']
        Q = f['prediction_hazard_voxel']
        for v in (X, R, Q):
            if v.ndim != 2 or v.shape[1] != 3 or len(v) == 0 or (not np.issubdtype(v.dtype, np.integer)):
                raise ValueError('Invalid or empty coordinates -> cannot certify')
            if len(np.unique(v, axis=0)) != len(v):
                raise ValueError('Duplicate voxel coordinates')
        yr = surface_of_coordinates(R)
        yq = surface_of_coordinates(Q)
    result = paired_evaluate(X, R, Q, s, a.threshold_mm)
    result.update(overlap_metrics(R, Q, s))
    result['surface_HD95_mm'] = hd95_points(yr, yq, s)
    result['reference_kind'] = a.reference_kind
    result['reference_validation'] = 'Caller-declared; modality, completeness and registration bounds must be independently verified'
    result['scope'] = 'Fixed geometry threshold, complete submitted masks; geometry safety only. Original stochastic guide risk, fabrication error and biological outcome not evaluated.'
    manifest = json.load(open(a.bound_manifest)) if a.bound_manifest else None
    result['bounded_reference_decision'] = reference_interval(result['reference']['gap_mm'], a.threshold_mm, manifest)
    result['physical_certified'] = False
    text = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if a.output:
        open(a.output, 'w').write(text)
    else:
        print(text, end='')
if __name__ == '__main__':
    main()
