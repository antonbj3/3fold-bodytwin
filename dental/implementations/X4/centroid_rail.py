"""A 3D centre rail of the mirrored missing surface, not an inferior border."""
import json
import numpy as np
from freeze import ROOT, write, freeze, sha, state

def rail(q, pre, anchor, direction=None):
    from planner import dist
    missing = q[dist(q, pre) > 3.0]
    if len(missing) < 200:
        raise ValueError('UNKNOWN_INSUFFICIENT_MIRROR_ADDITION')
    centre = np.median(missing, axis=0)
    if direction is None:
        (_, _, Vt) = np.linalg.svd(missing[::max(1, len(missing) // 10000)] - centre, full_matrices=False)
        direction = Vt[0]
        if direction[np.argmax(abs(direction))] < 0:
            direction = -direction
    direction = np.asarray(direction, float)
    projection = (missing - centre) @ direction
    edges = np.linspace(np.quantile(projection, 0.01), np.quantile(projection, 0.99), 42)
    indices = np.digitize(projection, edges) - 1
    points = []
    rows = []
    for b in range(41):
        group = missing[indices == b]
        if len(group) < 20:
            continue
        points.append(np.median(group, axis=0))
        rows.append(dict(bin=b, n=len(group)))
    if len(points) < 3:
        raise ValueError('UNKNOWN_TOO_FEW_CENTROID_STATIONS')
    points = np.array(points)
    return (points, dict(semantics='3D median centres of mirrored missing-surface PCA slabs; not a true medial skeleton', projection_direction=direction.tolist(), origin_proxy_mm=centre.tolist(), bins=rows, source_added_points=len(missing), arc_length_mm=float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum()), possible_failures='multimodal PCA slices, missing bilateral geometry, false mirrored additions; no anatomy labels'))

def freeze_r4():
    r3 = json.loads((ROOT / 'RESULTS_R3.json').read_text())
    freeze('R4', dict(prediction=['Pre only', 'same skull-centre sagittal mirror proxy'], evaluation=['Post virtual reference', 'Original defective geometry']), 'mirror added surface -> PCA-slab 3D median-centre rail -> thickness-constrained partitions -> corrected straight cylinder transport -> STL', dict(all_cases_available=True, median_primary_p95_mm_max=3.0, paired_gain_vs_raw_mirror_mm_min=1.0, solver_max_absolute_disagreement_mm=1e-08, all_miters_non_crossing=True, every_gate_fault_injection=True))
    path = ROOT / 'PREREG_R4.json'
    reg = json.loads(path.read_text())
    if 'centroid_contract' not in reg:
        reg['centroid_contract'] = dict(stations=41, projection_quantiles=[0.01, 0.99], minimum_points_per_station=20, centre='coordinate-wise median of whole 3D added-surface slab, no target Post', radius_mm=6.0, strongest_control='independent cross-product chord costs and exhaustive partitions under same finite miter constraints', baseline='all118 frozen raw Pre mirrors from FROZEN_MIRROR_BASELINE.json; also actual published OsteoOpt simplifier on same rail', external_arc_metric='same slab direction applied to expert Post-added surface; different from inferior-rail metric in earlier rounds', novelty='not claimed; PCA medians and finite partition search are conventional; new output information/representation is tested against source practice proxy')
        write(path, reg)
        (ROOT / 'PREREG_R4.sha256').write_text(sha(path) + '  ' + path.name + '\n')
        d = json.loads((ROOT / 'DECOMPOSITION_R4.json').read_text())
        d['leaves'].append(dict(name='3D centroid rail', status='CONSTITUTIVE_CLOSURE', basis='first PCA direction of Pre-mirror addition; median3D centres in 41 fixed slabs', stopping_argument='not a unique medial axis or anatomical border; slice multimodality and bilateral absence can refute it'))
        write(ROOT / 'DECOMPOSITION_R4.json', d)
    state('R4_preregistered', 'freeze valid Pre-only 3D centroid-rail plans before Post evaluation', r3['outcome'])
if __name__ == '__main__':
    freeze_r4()
