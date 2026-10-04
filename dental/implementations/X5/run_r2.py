from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import datetime, hashlib, json, resource, sys, time
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from voxel_decisions import local_cube_field, direct_cube_field, decision_report, overlap_metrics, hd95_points, verify_patch_selection, witness_subset_report, verify_witness_subset_claim
from run_r1 import put, sha, state
P = Path(__file__).resolve().parent
D = Path(__import__('os').environ.get('DENTAL_PROJECT_ROOT', str(P.parent.parent))).resolve()
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X5_decision_seg_error'))
DATA.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(D / 'cells/geometry'))
import tf2_io

def points(mask):
    return np.argwhere(mask & ~ndi.binary_erosion(mask))

def adversarial_additions(X, solid, surface, spacing, threshold, shape):
    """One-cell addition either near fixed design or far from it.

    Equal FP counts ensure exact equal Dice. The spurious island is an explicit
    topology failure probe, not a claim about realistic model-error prevalence.
    """
    occupied = {tuple(v) for v in solid}
    body = {tuple(v) for v in X}
    near = None
    for x in X:
        if tuple(x) not in occupied:
            near = x.copy()
            break
    if near is None:
        raise ValueError('No body point outside hazard to inject')
    tree = cKDTree(X * np.array(spacing))
    far = None
    for y in surface[::-1]:
        for delta in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]:
            z = y + np.array(delta)
            if np.any(z < 0) or np.any(z >= shape) or tuple(z) in occupied:
                continue
            dc = float(tree.query(z * np.array(spacing))[0])
            if dc > threshold + np.linalg.norm(spacing) + 2:
                far = z
                break
        if far is not None:
            break
    if far is None:
        raise ValueError('No safe remote addition found')
    return (near, far)

def main():
    st = time.perf_counter()
    cpu = time.process_time()
    assert sha(P / 'PREREG_R2.json') == (P / 'PREREG_R2.sha256').read_text().strip()
    assert sha(P / 'FROZEN_PREDICTIONS.json') == (P / 'FROZEN_PREDICTIONS.sha256').read_text().strip()
    pr = json.load(open(P / 'PREREG_R2.json'))
    rows = []
    manifest = []
    globalpairs = []
    state('R2_RUNNING', 'PREREG_AND_PREDICTIONS_HASH_OK', 'Full local fields, direct all-pairs controls and equal-Dice topology error probes')
    for target in pr['selected_cases']:
        t0 = time.perf_counter()
        case = target['case']
        (lab, sp, hdr) = tf2_io.load(case)
        assert hdr['TransformMatrix'] == '1 0 0 0 1 0 0 0 1' and hdr['Offset'] == '0 0 0'
        labelsha = hashlib.sha256(lab.tobytes()).hexdigest()
        assert labelsha == target['source_label_sha256']
        (cc, n) = ndi.label(lab == 10)
        X = np.argwhere(cc == target['component'])
        cm = np.isin(lab, [3, 4])
        solid = np.argwhere(cm)
        Y = points(cm)
        tfield = time.perf_counter()
        f = local_cube_field(X, Y, sp)
        fieldsec = time.perf_counter() - tfield
        tcontrol = time.perf_counter()
        direct = direct_cube_field(X, Y, sp)
        directsec = time.perf_counter() - tcontrol
        err = float(abs(direct - f['distance_mm']).max())
        r = decision_report(f, Y, pr['threshold_mm'], pr['scenario_patch_budget_mm'], sp, pr['patch_size_mm'])
        replay = abs(r['gap_mm'] - target['expected_annotation_gap_mm'])
        dvec = (f['body_witness_voxel'] - Y) * np.array(sp)
        separate = np.sign(dvec) * np.maximum(abs(dvec) - sp, 0)
        nz = f['distance_mm'] > 0
        direct_dir = separate[nz] / np.linalg.norm(separate[nz], axis=1)[:, None]
        direrr = float(abs(direct_dir - f['approach_direction_zyx'][nz]).max()) if nz.any() else 0.0
        wrongdir_rejected = bool(np.max(abs(-f['approach_direction_zyx'][nz] - direct_dir)) > 1e-06) if nz.any() else True
        dangerous_direct = np.array([d - pr['threshold_mm'] <= pr['scenario_patch_budget_mm'] for d in direct])
        dangerous = f['distance_mm'] - pr['threshold_mm'] <= pr['scenario_patch_budget_mm']
        patchcheck = bool(np.array_equal(dangerous, dangerous_direct))
        removed = dangerous.copy()
        if dangerous.any():
            removed[np.flatnonzero(dangerous)[0]] = False
        omit_rejected = not verify_patch_selection(direct, pr['threshold_mm'], pr['scenario_patch_budget_mm'], removed) if dangerous.any() else None
        i = int(np.argmin(f['distance_mm']))
        delta = f['distance_mm'][i] - pr['threshold_mm'] + 0.05
        moved = Y[i] * np.array(sp) + f['approach_direction_zyx'][i] * delta
        witness = f['body_witness_voxel'][i] * np.array(sp)
        newdist = float(np.linalg.norm(np.maximum(abs(moved - witness) - sp, 0)))
        direction_crossed = bool(newdist < pr['threshold_mm'])
        (near, far) = adversarial_additions(X, solid, Y, sp, pr['threshold_mm'], lab.shape)
        near_mask = cm.copy()
        near_mask[tuple(near)] = True
        far_mask = cm.copy()
        far_mask[tuple(far)] = True
        yn = points(near_mask)
        yf = points(far_mask)
        fn = local_cube_field(X, yn, sp)
        ff = local_cube_field(X, yf, sp)
        rn = decision_report(fn, yn, pr['threshold_mm'], spacing=sp)
        rf = decision_report(ff, yf, pr['threshold_mm'], spacing=sp)
        near_check = float(direct_cube_field(X, np.array([near]), sp)[0])
        far_check = float(direct_cube_field(X, np.array([far]), sp)[0])
        n_expected = min(r['gap_mm'], near_check)
        f_expected = min(r['gap_mm'], far_check)
        metrics_n = overlap_metrics(solid, np.vstack([solid, near]), sp)
        metrics_f = overlap_metrics(solid, np.vstack([solid, far]), sp)
        hn = hd95_points(Y, yn, sp)
        hf = hd95_points(Y, yf, sp)
        pair = {'case': case, 'component': target['component'], 'Dice_near': metrics_n['Dice'], 'Dice_far': metrics_f['Dice'], 'HD95_near_mm': hn, 'HD95_far_mm': hf, 'reference_gap_mm': r['gap_mm'], 'near_gap_mm': rn['gap_mm'], 'far_gap_mm': rf['gap_mm'], 'reference_pass': r['nominal_pass'], 'near_pass': rn['nominal_pass'], 'far_pass': rf['nominal_pass'], 'near_added_voxel': near.tolist(), 'far_added_voxel': far.tolist(), 'fault_kind': 'ONE_FALSE_POSITIVE_VOXEL_TOPOLOGY_PROBE_ON_REAL_DATA', 'external_truth': 'INJECTED_MASKS_ARE_OUR_OWN_FIXTURE; base geometry is published annotation', 'same_global_metrics_different_decision': abs(metrics_n['Dice'] - metrics_f['Dice']) <= 1e-12 and abs(hn - hf) <= 1e-09 and (rn['nominal_pass'] != rf['nominal_pass']), 'direct_prediction_control_error_mm': max(abs(rn['gap_mm'] - n_expected), abs(rf['gap_mm'] - f_expected))}
        globalpairs.append(pair)
        fpath = DATA / f"{case}_cc{target['component']}_local.npz"
        np.savez_compressed(fpath, body_voxel=X.astype(np.int32), hazard_surface_voxel=Y.astype(np.int32), hazard_solid_voxel=solid.astype(np.int32), spacing_zyx_mm=sp, distance_mm=f['distance_mm'], body_witness_voxel=f['body_witness_voxel'].astype(np.int32), approach_direction_zyx=f['approach_direction_zyx'], dangerous=dangerous, near_added_voxel=near, far_added_voxel=far)
        manifest.append({'path': str(fpath), 'bytes': fpath.stat().st_size, 'sha256': sha(fpath), 'source_label_array_sha256': labelsha, 'zip_member': f'{tf2_io.ROOT}/labelsTr/{case}.mha', 'zip_path': tf2_io.ZIP})
        cells = Y.shape[0]
        budget_sweep = []
        for B in [0.045, 0.08, 0.22, 0.3, 0.52, 0.77, 1.99]:
            q = decision_report(f, Y, pr['threshold_mm'], B, sp, pr['patch_size_mm'])
            budget_sweep.append({k: q[k] for k in ('scenario_budget_mm', 'dangerous_cells', 'measurement_cells', 'measurement_cell_fraction', 'measurement_patches')})
        out = {'case': case, 'component': target['component'], 'target': target, 'field': r, 'budget_sweep': budget_sweep, 'controls': {'source_replay_error_mm': replay, 'direct_allpairs_max_error_mm': err, 'direction_error': direrr, 'direction_crossed_threshold': direction_crossed, 'wrong_direction_rejected': wrongdir_rejected, 'dangerous_set_direct_replay': patchcheck, 'omitted_dangerous_cell_rejected': omit_rejected if dangerous.any() else 'NOT_APPLICABLE_EMPTY_SET', 'corrupted_gap_plus0.05_rejected': abs(r['gap_mm'] + 0.05 - target['expected_annotation_gap_mm']) > 1e-09, 'corrupted_distance_field_plus0.05_rejected': float(abs(direct - f['distance_mm'] - 0.05).max()) > 1e-09}, 'cost': {'field_seconds': fieldsec, 'direct_allpairs_seconds': directsec, 'total_case_seconds': time.perf_counter() - t0}, 'data_path': str(fpath), 'data_sha256': manifest[-1]['sha256']}
        rows.append(out)
        put('RAW_R2_FIELDS.json', rows)
        put('R2_GLOBAL_METRIC_COUNTEREXAMPLES.json', globalpairs)
        state('R2_RUNNING', f"{len(rows)}/6 local sources processed, replay max error {max((x['controls']['source_replay_error_mm'] for x in rows)):.3g} mm", 'Continue complete fields and equal-global-metric adversaries')
        print(json.dumps({'case': case, 'cells': cells, 'gap_mm': r['gap_mm'], 'direct_error_mm': err, 'measurement_fraction': r['measurement_cell_fraction'], 'same_metric_flip': pair['same_global_metrics_different_decision'], 'seconds': out['cost']['total_case_seconds']}), flush=True)
        del lab, cc, cm, near_mask, far_mask, X, Y, solid, f, direct, fn, ff
    nrows = list(map(json.loads, open(D / 'results/LANE_NEXT_N_CANAL_PROXIMITY/RAW_N2.jsonl')))
    nw = []
    for x in nrows:
        for method in ('SEG', 'HM'):
            ws = [w for w in x['witnesses'] if w['canal_method'] == method]
            if not ws:
                continue
            ds = [w['gap_mm'] for w in ws]
            limit = max(ds)
            correct = witness_subset_report(ds, 2.0)
            forged = {**correct, 'can_certify_full_mask_above_threshold': True}
            actual_rejected = not verify_witness_subset_claim(ds, 2.0, forged)
            nw.append({'case': x['case'], 'component': x['component'], 'method': method, 'witness_variants': len(ws), 'witness_min_mm': min(ds), 'witness_max_mm': limit, 'upper_witness_below_threshold': limit < 2, 'can_certify_full_mask_above_threshold': False, 'physical_status': 'UNKNOWN', 'full_safe_injection_rejected': actual_rejected, 'coverage': 'SUPPORTED_POINT_SUBSET; cannot bound unknown minimum from below'})
    put('R2_N_WITNESS_PORT.json', nw)
    put('DATA_MANIFEST.json', manifest)
    gates = {'source_replay': all((r['controls']['source_replay_error_mm'] <= 1e-09 for r in rows)), 'direct_allpairs': all((r['controls']['direct_allpairs_max_error_mm'] <= 1e-09 for r in rows)), 'equal_metrics_decision_disagreement': sum((r['same_global_metrics_different_decision'] for r in globalpairs)) >= 4, 'patch_compression': sum((r['field']['measurement_cell_fraction'] <= 0.1 for r in rows)) >= 4, 'direct_patch_replay': all((r['controls']['dangerous_set_direct_replay'] for r in rows)), 'direction_crossing': all((r['controls']['direction_crossed_threshold'] for r in rows)), 'wrong_direction_rejected': all((r['controls']['wrong_direction_rejected'] for r in rows)), 'corrupted_gap_rejected': all((r['controls']['corrupted_gap_plus0.05_rejected'] for r in rows)), 'corrupted_field_rejected': all((r['controls']['corrupted_distance_field_plus0.05_rejected'] for r in rows)), 'omitted_dangerous_cell_rejected': all((r['controls']['omitted_dangerous_cell_rejected'] is True for r in rows if r['field']['dangerous_cells'] > 0)), 'forged_full_surface_witness_rejected': all((r['full_safe_injection_rejected'] for r in nw)), 'prediction_control': all((r['direct_prediction_control_error_mm'] <= 1e-09 for r in globalpairs)), 'physical_facit': 'UNKNOWN_NO_LOCAL_INDEPENDENT_ANATOMY'}
    output = {'id': pr['id'], 'n_scans': len(rows), 'gates': gates, 'adversarial_same_metrics_flips': sum((r['same_global_metrics_different_decision'] for r in globalpairs)), 'patch_fractions': [r['field']['measurement_cell_fraction'] for r in rows], 'outcome': 'LOCAL_CAPABILITY_DELIVERED_DIRECT_CONTROL_TIE_PHYSICAL_VALIDATION_UNKNOWN', 'external_referent': pr['external_referent'], 'external_physical_referent': pr['external_physical_referent'], 'cost': {'wall_seconds': time.perf_counter() - st, 'cpu_seconds': time.process_time() - cpu, 'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1, 'GPU': False, 'output_bytes': sum((r['bytes'] for r in manifest)), 'source_preparation_and_independent_measurement_seconds': None}}
    put('SUMMARY_R2.json', output)
    (P / 'HANDOFF_R2.md').write_text('# R2 avgjord\n\nAll six frozen TF2 case executed. SUMMARY_R2.json preserves each gate , RAW_R2_FIELDS.json local values , DATA_MANIFEST.json arrayhashar. Complete submitted mask means full coverage of the mask, not full anatomy. Global-measurement probes are injected errors, not external measurements.\n\nR2 changes the representation to all targeted local channelytehot and provides a model-independent prediction mask-consumer as well as a conditional minimum measurement patch set. Strong direct control receives the same data and gives TIE . Independent physical edge calibration continues to be lacking.\n\nNext R3 : replace independent per-cell error with measured image edge profiles from N and a local review/measurement port contract. It must be able to separate strong alternative channel/metal contours and cannot be readjusted to the outcome.\n')
    state('R2_COMPLETE', str(gates), 'Freeze a genuinely new information construction: local independent measurement specification or N raw contour decision sensitivity')
    print(json.dumps(output), flush=True)
if __name__ == '__main__':
    main()
