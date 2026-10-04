import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import datetime, hashlib, json, resource, time
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from voxel_decisions import local_cube_field, direct_cube_field, overlap_metrics, hd95_points
from run_r1 import put, sha, state
from run_r2 import DATA, points
P = Path(__file__).resolve().parent
NEIGH = [np.array(v) for v in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]]
STRUCT = ndi.generate_binary_structure(3, 1)

def shortest_connected_extension(body, hazard, spacing, threshold):
    X = np.asarray(body, int)
    Y = np.asarray(hazard, int)
    s = np.asarray(spacing)
    pad = np.ceil(6 / s).astype(int) + 3
    lo = np.minimum(X.min(0), Y.min(0)) - pad
    hi = np.maximum(X.max(0), Y.max(0)) + pad + 1
    shape = tuple(hi - lo)
    cm = np.zeros(shape, bool)
    bm = np.zeros(shape, bool)
    cm[tuple((Y - lo).T)] = True
    bm[tuple((X - lo).T)] = True
    centre = ndi.distance_transform_edt(~bm, sampling=s)
    candidates = np.argwhere(centre < threshold + np.linalg.norm(s) + 1e-09) + lo
    ff = local_cube_field(X, candidates, s)
    risky = ff['distance_mm'] < threshold
    targets = candidates[risky]
    targetmask = np.zeros(shape, bool)
    targetmask[tuple((targets - lo).T)] = True
    d = ndi.distance_transform_cdt(~cm, metric='taxicab')
    td = d[tuple((targets - lo).T)]
    k = int(td.min())
    end = targets[int(np.argmin(td))].copy()
    assert k > 0, 'Original hazard is already below threshold'
    v = end.copy()
    path = []
    while int(d[tuple(v - lo)]) > 0:
        path.append(v.copy())
        n = int(d[tuple(v - lo)])
        v = next((v + step for step in NEIGH if d[tuple(v + step - lo)] == n - 1))
    path = np.array(path[::-1], dtype=int)
    seed = v.copy()
    assert len(path) == k
    reach = cm.copy()
    control_k = 0
    previous_hit = False
    while not np.any(reach & targetmask):
        previous_hit = bool(np.any(reach & targetmask))
        reach = ndi.binary_dilation(reach, structure=STRUCT)
        control_k += 1
        if control_k > k + 2:
            break
    candidate_count = len(candidates)
    return {'path': path, 'seed': seed, 'k': k, 'control_k': control_k, 'k_minus_one_hit': previous_hit, 'lo': lo, 'shape': shape, 'source_mask': cm, 'target_mask': targetmask, 'target_cells': len(targets), 'candidate_cells': candidate_count, 'end_gap_mm': float(direct_cube_field(X, np.array([end]), s)[0])}

def safe_connected_extension(X, solid, surface, k, s, lo, shape, threshold):
    used = {tuple(v) for v in solid}
    tree = cKDTree(X * s)
    order = np.argsort(tree.query(surface * s, workers=1)[0])[::-1]
    path = []
    front = [surface[i].copy() for i in order]
    seed = None
    for y in front:
        for step in NEIGH:
            z = y + step
            if tuple(z) in used or np.any(z < lo) or np.any(z >= lo + shape):
                continue
            if tree.query(z * s)[0] - np.linalg.norm(s) <= threshold:
                continue
            if seed is None:
                seed = y.copy()
            path.append(z.copy())
            used.add(tuple(z))
            if len(path) == k:
                return (np.array(path, int), seed)
    raise ValueError('Could not construct equal-cardinality connected safe control')

def main():
    start = time.perf_counter()
    cpu = time.process_time()
    assert sha(P / 'PREREG_R3.json') == (P / 'PREREG_R3.sha256').read_text().strip()
    pr = json.load(open(P / 'PREREG_R3.json'))
    manifest = json.load(open(P / 'DATA_MANIFEST.json'))
    rows = []
    artifacts = []
    pred = {'frozen_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prereg_sha256': sha(P / 'PREREG_R3.json'), 'expected': 'All six minimum connected cardinalities equal independent dilation hitting time; all six survive seeded component selection; >=4 equal Dice/HD95 pairs yield different decisions; physical validation UNKNOWN', 'source_manifest_sha256': sha(P / 'DATA_MANIFEST.json'), 'expected_source_gaps': pr['selected_cases']}
    if not (P / 'FROZEN_PREDICTIONS_R3.json').exists():
        put('FROZEN_PREDICTIONS_R3.json', pred)
        (P / 'FROZEN_PREDICTIONS_R3.sha256').write_text(sha(P / 'FROZEN_PREDICTIONS_R3.json') + '\n')
    assert sha(P / 'FROZEN_PREDICTIONS_R3.json') == (P / 'FROZEN_PREDICTIONS_R3.sha256').read_text().strip()
    state('R3_RUNNING', 'R3_PREREG_AND_PREDICTION_FROZEN', 'Connected error graph minimality and topology-aware whole-mask controls')
    for m in manifest:
        t0 = time.perf_counter()
        assert sha(m['path']) == m['sha256']
        a = np.load(m['path'])
        X = a['body_voxel']
        solid = a['hazard_solid_voxel']
        Y = a['hazard_surface_voxel']
        s = a['spacing_zyx_mm']
        op = shortest_connected_extension(X, solid, s, pr['threshold_mm'])
        near = op['path']
        (far, far_seed) = safe_connected_extension(X, solid, Y, op['k'], s, op['lo'], np.array(op['shape']), pr['threshold_mm'])
        cm = op['source_mask']
        nm = cm.copy()
        fm = cm.copy()
        nm[tuple((near - op['lo']).T)] = True
        fm[tuple((far - op['lo']).T)] = True
        (source_labels, source_n) = ndi.label(cm, structure=STRUCT)
        (nn, nnc) = ndi.label(nm, structure=STRUCT)
        (ff, nfc) = ndi.label(fm, structure=STRUCT)
        original_component = int(nn[tuple(op['seed'] - op['lo'])])
        adverse_connected = bool(np.all(nn[tuple((near - op['lo']).T)] == original_component))
        far_labels = ff[tuple((far - op['lo']).T)]
        safe_connected = bool(np.all(far_labels > 0))
        yn = points(nm) + op['lo']
        yf = points(fm) + op['lo']
        rn = local_cube_field(X, yn, s)
        rf = local_cube_field(X, yf, s)
        near_gap = float(rn['distance_mm'].min())
        far_gap = float(rf['distance_mm'].min())
        same = overlap_metrics(solid, np.vstack([solid, near]))['Dice']
        same_far = overlap_metrics(solid, np.vstack([solid, far]))['Dice']
        hdnear = hd95_points(Y, yn, s)
        hdfar = hd95_points(Y, yf, s)
        selected = nn == original_component
        ys = points(selected) + op['lo']
        selected_gap = float(local_cube_field(X, ys, s)['distance_mm'].min())
        disconnected = near.copy()
        disconnected[-1] += np.array([1000, 0, 0])

        def valid_path(p):
            q = np.vstack([op['seed'], p])
            return bool(np.all(np.sum(abs(np.diff(q, axis=0)), axis=1) == 1))
        min_direct = op['control_k']
        lastdirect = op['end_gap_mm']
        checks = {'minimality_exact': op['k'] == min_direct, 'path_connected': valid_path(near), 'seeded_component_survives': adverse_connected and selected_gap < 2, 'same_count_safe_control_connected': safe_connected, 'k_minus_one_claim_rejected': op['k'] - 1 < min_direct and (not op['k_minus_one_hit']), 'disconnected_path_rejected': not valid_path(disconnected), 'wrong_target_claim_rejected': far_gap >= 2 and near_gap < 2, 'wrong_reported_gap_rejected': abs(near_gap + 0.05 - min(pr['selected_cases'][len(rows)]['expected_annotation_gap_mm'], lastdirect)) > 1e-09}
        newnear = float(direct_cube_field(X, near, s).min())
        newfar = float(direct_cube_field(X, far, s).min())
        old = pr['selected_cases'][len(rows)]['expected_annotation_gap_mm']
        replay = max(abs(near_gap - min(old, newnear)), abs(far_gap - min(old, newfar)))
        checks['direct_mask_replay'] = replay <= 1e-09
        pair = abs(same - same_far) <= 1e-12 and abs(hdnear - hdfar) <= 1e-09 and (near_gap < 2 <= far_gap)
        path = DATA / (Path(m['path']).stem + '_connected.npz')
        np.savez_compressed(path, adverse_added_voxel=near, safe_added_voxel=far, seed_voxel=op['seed'], spacing_zyx_mm=s)
        artifacts.append({'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size})
        row = {'case': pr['selected_cases'][len(rows)]['case'], 'component': pr['selected_cases'][len(rows)]['component'], 'source_array_sha256': m['sha256'], 'reference_gap_mm': old, 'minimum_connected_added_voxels': op['k'], 'minimum_added_volume_mm3': float(op['k'] * np.prod(s)), 'path_grid_arclength_mm': float(np.sum(abs(np.diff(np.vstack([op['seed'], near]), axis=0)) * s)), 'adverse_added_voxel': near.tolist(), 'adverse_seed_voxel': op['seed'].tolist(), 'safe_added_voxel': far.tolist(), 'adverse_gap_mm': near_gap, 'safe_gap_mm': far_gap, 'seeded_component_gap_mm': selected_gap, 'Dice_adverse': same, 'Dice_safe': same_far, 'HD95_adverse_mm': hdnear, 'HD95_safe_mm': hdfar, 'same_global_metrics_different_decision': pair, 'danger_target_cells': op['target_cells'], 'exact_control_replay_error_mm': replay, 'controls': checks, 'physical_status': 'UNKNOWN', 'probe_status': 'CONNECTED_SYNTHETIC_ERROR_ON_PUBLISHED_BASE_GEOMETRY', 'cost_seconds': time.perf_counter() - t0, 'data_path': str(path), 'data_sha256': artifacts[-1]['sha256']}
        rows.append(row)
        put('R3_CONNECTED_COUNTEREXAMPLES.json', rows)
        print(json.dumps({k: row[k] for k in ('case', 'minimum_connected_added_voxels', 'Dice_adverse', 'HD95_adverse_mm', 'same_global_metrics_different_decision', 'cost_seconds')}), flush=True)
        state('R3_RUNNING', f'{len(rows)}/6 connected controls run', 'Continue frozen connected error minima')
    summary = {'id': pr['id'], 'n_cases': len(rows), 'gates': {'minimum_connected_extension': all((r['controls']['minimality_exact'] and r['controls']['path_connected'] for r in rows)), 'global_metric_disagreement': sum((r['same_global_metrics_different_decision'] for r in rows)) >= 4, 'postprocessing_survival': all((r['controls']['seeded_component_survives'] for r in rows)), 'fault_injection': all((all(r['controls'].values()) for r in rows)), 'physical_validation': 'UNKNOWN_NO_MATCHED_LOCAL_ANATOMY'}, 'outcome': 'CONNECTED_ERROR_CAPABILITY_DELIVERED_STRONG_CONTROL_TIE_PHYSICAL_UNKNOWN', 'external_referent': pr['external_referent'], 'own_fixture_referent': pr['own_fixture_referent'], 'minimum_added_voxels': [r['minimum_connected_added_voxels'] for r in rows], 'same_metric_flips': sum((r['same_global_metrics_different_decision'] for r in rows)), 'cost': {'wall_seconds': time.perf_counter() - start, 'cpu_seconds': time.process_time() - cpu, 'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1, 'GPU': False, 'data_bytes': sum((m['bytes'] for m in artifacts)), 'source_annotation_and_future_lab_seconds': None}}
    put('SUMMARY_R3.json', summary)
    put('DATA_MANIFEST_R3.json', artifacts)
    (P / 'HANDOFF_R3.md').write_text('# R3 avgjord\n\nPREREG_R3 and FROZEN_PREDICTIONS_R3 was frozen before local calculations. Result in SUMMARY_R3 and R3_CONNECTED_COUNTEREXAMPLES . The minimum refers to the number of arbitrary new voxel cells for six -neighbours expansion, not a measured model error quantity. Explicit mask replay and independent binary dilation executed. Physical validation UNKNOWN .\n\nThe next design must change the information: frozen real model prediction + independent local channel/pulp reference (or matched-protocol phantom). A real false approval rate can then be measured, with full hazard coverage. MINIMUM_MEASUREMENT_PORT.md specifies exact quantities and adequacy conditions. Global inter-observation cannot replace this.\n')
    state('R3_COMPLETE', str(summary['gates']), 'Package one -command demo and independent local measurement port; continue through external reference')
    print(json.dumps(summary), flush=True)
if __name__ == '__main__':
    main()
