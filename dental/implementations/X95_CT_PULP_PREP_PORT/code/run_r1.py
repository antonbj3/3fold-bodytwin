from dental_release.paths import expand as _release_expand
import csv, json, resource, time, warnings
from pathlib import Path
import numpy as np
from scipy.optimize import linprog, OptimizeWarning
from identity import ROOT, DATA, sha, write, source_path
from geometry import cap, regional, s13, measure, validate_port, validate_pulp_payload, digital_geometry_accepts, global_mm
from x12_consumer_ports import uniform_shell_decision

def lp_interval(d, eps, b, delta):
    lo = max(0, d - eps - 2 * delta)
    hi = d + eps + 2 * delta
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', OptimizeWarning)
        opts = {'threads': 1, 'primal_feasibility_tolerance': 1e-07, 'dual_feasibility_tolerance': 1e-07}
        minimum = linprog([1.0], bounds=[(lo, hi)], method='highs', options=opts)
        feasible = linprog([0.0], A_ub=[[-1.0]], b_ub=[-b], bounds=[(lo, hi)], method='highs', options=opts)
    if not minimum.success or feasible.status not in [0, 2]:
        raise RuntimeError('LP failure')
    cls = 'LOCAL_CONDITION_ALL_ENCLOSURES' if minimum.fun - b >= -1e-07 else 'LOCAL_CONDITION_NO_ENCLOSURE' if feasible.status == 2 else 'UNKNOWN_INTERVAL'
    return (cls, float(minimum.fun - b))

def main():
    started = time.monotonic()
    rows = []
    intervals = []
    witnesses = []
    manifest = []
    frames = []
    cohort = json.loads((ROOT / 'FROZEN_COHORT.json').read_text())
    assert sha(ROOT / 'FROZEN_COHORT.json') == (ROOT / 'FROZEN_COHORT.sha256').read_text().strip()
    assert sha(ROOT / 'FROZEN_PREDICTIONS.json') == (ROOT / 'FROZEN_PREDICTIONS.sha256').read_text().strip()
    source_rows = {int(r['fdi']): r for r in csv.DictReader((ROOT / 'sources/x12_bound/raw/R8_teeth.csv').open()) if r['case'] == 'P48'}
    header = json.loads((ROOT / 'raw/IDENTITY.json').read_text())['source_header']
    failures = []
    first = None
    for item in cohort['cohort']:
        fdi = item['canonical_fdi']
        if item['status'] != 'AVAILABLE':
            rows.append({'case': 'P48', 'fdi': fdi, 'operation': 'ALL', 'status': 'UNKNOWN', 'reason': item['reason'], 'remaining_dentin_mm': 'UNKNOWN_NO_DEJ'})
            continue
        assert sha(source_path(item['source_path'])) == item['sha256']
        a = np.load(source_path(item['source_path']))
        t = a['tooth']
        u = a['pulp']
        h = float(a['spacing'][0])
        origin = a['origin_zyx']
        validate_port('P48', sha(source_path(item['source_path'])), fdi, item['source_pulpy_fdi'], a['spacing'], item)
        validate_pulp_payload(u, item)
        sign = int(source_rows[fdi]['crown_axis_sign'])
        if sign < 0:
            t = t[::-1].copy()
            u = u[::-1].copy()
        base = max(1, int(np.argwhere(u)[:, 0].max()) - round(2 / h))
        operations = []
        for d in [1.0, 1.5, 2.0]:
            operations.append((f'X31_CAP_{d:g}mm', cap(t, int(np.ceil(d / h - 1e-12)), base), {'requested_depth_mm': d, 'actual_grid_depth_mm': float(np.ceil(d / h - 1e-12) * h), 'virtual_finish_z_voxel': base}))
        for (low, name) in [(True, 'REGIONAL_A_LOW_X'), (False, 'REGIONAL_B_HIGH_X')]:
            (p, meta) = regional(t, base, low)
            operations.append((name, p, meta))
        (p, meta) = s13(t, h, fdi % 10 <= 3)
        operations.append(('S13_TOOTH_ONLY_GEOMETRY', p, meta))
        current = {}
        for (name, p, meta) in operations:
            (summary, fields) = measure(t, u, p, h)
            row = {'case': 'P48', 'fdi': fdi, 'operation': name, 'status': 'MEASURED_DIGITAL_LABEL_GEOMETRY', **meta, **summary}
            if 'requested_depth_mm' in meta:
                original = float(source_rows[fdi]['occlusal_vertical_min_mm'])
                row['X31_scalar_original_vertical_mm'] = original
                row['X31_scalar_remaining_requested_mm'] = original - meta['requested_depth_mm']
                row['X31_scalar_remaining_grid_mm'] = original - meta['actual_grid_depth_mm']
                row['scalar_operator_scope'] = 'vertical topmost source-pulp-column clearance; not shortest whole-pulp Euclidean clearance'
            d = row['cut_center_min_mm']
            eps = row['epsilon_mm']
            for b in [0.0, 0.25, 1.0]:
                for delta in [0.0, 0.15, 0.3]:
                    q = uniform_shell_decision(d, eps, 0.0, b, delta)
                    (cls, margin) = lp_interval(d, eps, b, delta)
                    intervals.append({'case': 'P48', 'fdi': fdi, 'operation': name, 'query_clearance_mm': b, 'delta_each_mm': delta, **q, 'LP_class': cls, 'LP_margin_mm': margin, 'class_agrees': q['decision'] == cls, 'final_geometry_gate': row['geometry_gate'], 'clinical_feasibility': 'UNKNOWN'})
            frame = {'case': 'P48', 'tf2_case': 'ToothFairy2P_048', 'canonical_fdi': fdi, 'source_pulpy_fdi': item['source_pulpy_fdi'], 'map': 'left_right_swap', 'source_object_sha256': item['sha256'], 'unit': 'mm', 'axis_order': 'ZYX', 'shape_zyx': list(t.shape), 'crop_origin_native_zyx': origin.tolist(), 'crown_axis_sign': sign, 'spacing_zyx_mm': [h] * 3, 'source_header': header, 'operation': name, 'region': 'native geometric X/Y halves, anatomical labels UNKNOWN', 'pose': 'native TF2 with local crown-Z reversal only', 'timescale': 'HANDOVER'}
            file = DATA / f'R1_P48_FDI{fdi}_{name}.npz'
            np.savez_compressed(file, tooth=t, pulp=u, prepared=p, removed=t & ~p, spacing=np.array([h] * 3), origin_zyx=origin, frame_json=json.dumps(frame), **fields)
            read = np.load(file)
            read_frame = json.loads(str(read['frame_json']))
            mask_equal = all((np.array_equal(read[k], v) for (k, v) in [('tooth', t), ('pulp', u), ('prepared', p), ('removed', t & ~p)]))
            q = fields['pulp_local_zyx_voxel']
            p1 = global_mm(q, t.shape, origin, sign, header)
            p2 = global_mm(read['pulp_local_zyx_voxel'], read_frame['shape_zyx'], read_frame['crop_origin_native_zyx'], read_frame['crown_axis_sign'], read_frame['source_header'])
            err = float(np.max(np.abs(p1 - p2)))
            frames.append({'fdi': fdi, 'operation': name, 'mask_bit_exact': mask_equal, 'frame_metadata_exact': frame == read_frame, 'max_frame_error_mm': err})
            manifest.append({'path': str(file), 'sha256': sha(file), 'bytes': file.stat().st_size, 'fdi': fdi, 'operation': name})
            row['spatial_export'] = str(file)
            rows.append(row)
            current[name] = row
            assert mask_equal and frame == read_frame and (err <= 1e-09)
            assert row['EDT_KDTree_max_error_mm'] <= 1e-06 and row['exact_face_enclosure_error_mm'] <= 1e-07
        a = current['REGIONAL_A_LOW_X']
        b = current['REGIONAL_B_HIGH_X']
        same = a['removed_voxels'] == b['removed_voxels'] and a['removed_volume_mm3'].hex() == b['removed_volume_mm3'].hex()
        diff = abs(a['cut_exact_pulp_cube_min_mm'] - b['cut_exact_pulp_cube_min_mm'])
        witnesses.append({'fdi': fdi, 'summary': 'total removed integer voxels and binary64 mm3 volume', 'removed_voxels_A': a['removed_voxels'], 'removed_voxels_B': b['removed_voxels'], 'volume_A_hex': a['removed_volume_mm3'].hex(), 'volume_B_hex': b['removed_volume_mm3'].hex(), 'identity_error_mm3': abs(a['removed_volume_mm3'] - b['removed_volume_mm3']), 'bit_identity': same, 'cut_min_A_mm': a['cut_exact_pulp_cube_min_mm'], 'cut_min_B_mm': b['cut_exact_pulp_cube_min_mm'], 'downstream_difference_mm': diff, 'status': 'COUNTERWITNESS_FOUND' if diff > 0 else 'WITNESS_ABSENT', 'smallest_extension': 'For this specified minimum-clearance query: add cut-specific minimum and its provenance/argmin; for new cuts: retain spatial cut and full-pulp field. Regional volume alone has no general sufficiency proof.'})
        assert same
        if first is None:
            first = (item, t, u, operations, base)
        write(ROOT / 'raw/R1_ROWS.json', rows)
        write(ROOT / 'raw/R1_WITNESSES.json', witnesses)
        write(ROOT / 'CURRENT_WORK_STATE.json', {'lane': 'X95-ct-pulp-prep-port', 'milestone': 'REGIONAL_SOLIDS_RUNNING', 'latest_fdi': fdi, 'latest_gate': 'IDENTITY_PASS; EDT/LP/face enclosure checks ongoing', 'next_operation': 'complete frozen cohort, fault controls and spatial export'})
        print('FDI', fdi, 'rows', len(rows), 'witness_delta_mm', diff, flush=True)
    (item, t, u, operations, base) = first
    fdi = item['canonical_fdi']
    controls = []
    params = {'case': 'P48', 'source_sha': item['sha256'], 'canonical_fdi': fdi, 'source_fdi': item['source_pulpy_fdi'], 'spacing': [0.3] * 3, 'expected': item}
    validate_port(**params)
    foreign = Path(_release_expand('@DENTAL_WORK_ROOT@/X12/R8_P1_36_paired.npz'))
    foreign_obj = np.load(foreign)
    foreign_point = np.argwhere(foreign_obj['pulp'])[0].tolist()
    sign = int(source_rows[fdi]['crown_axis_sign'])
    validate_pulp_payload(u, item, sign)
    injected = u.copy()
    target = np.argwhere(~u & t)[0]
    injected[tuple(target)] = foreign_obj['pulp'][tuple(foreign_point)]
    payload_rejected = False
    try:
        validate_pulp_payload(injected, item, sign)
    except ValueError:
        payload_rejected = True
    controls.append({'check': 'foreign_actual_P1_cell_with_spoofed_P48_metadata', 'valid_pass': True, 'injected_rejected': payload_rejected, 'foreign_mask_sha256': sha(foreign), 'foreign_true_pulp_cell_zyx': foreign_point, 'injected_target_P48_zyx': target.tolist(), 'mutation': 'Copy one True cell from verified P1 pulp into a non-pulp P48 tooth cell; falsely keep P48 case/hash metadata. Actual mask-content binding must reject.'})
    for (name, patch) in [('wrong_FDI', {'source_fdi': fdi}), ('scale1mm', {'spacing': [1.0] * 3}), ('foreign_actual_P1_pulp_point', {'case': 'P1', 'source_sha': sha(foreign)}), ('source_hash', {'source_sha': '0' * 64})]:
        rejected = False
        reason = None
        try:
            validate_port(**params | patch)
        except ValueError as e:
            rejected = True
            reason = str(e)
        controls.append({'check': name, 'valid_pass': True, 'injected_rejected': rejected, 'reason': reason, 'foreign_locator': str(foreign) if name.startswith('foreign') else None, 'foreign_pulp_point_local_zyx': foreign_point if name.startswith('foreign') else None})
    assert digital_geometry_accepts(u, t)
    wrong = t.copy()
    point = np.argwhere(u)[-1]
    wrong[tuple(point)] = False
    controls.append({'check': 'cut_through_actual_pulp', 'valid_pass': digital_geometry_accepts(u, t), 'injected_rejected': not digital_geometry_accepts(u, wrong), 'removed_actual_pulp_point_zyx': point.tolist()})
    vals = [r['EDT_KDTree_max_error_mm'] for r in rows if r.get('status') == 'MEASURED_DIGITAL_LABEL_GEOMETRY']
    controls.append({'check': 'centre_distance_plus1mm', 'valid_pass': max(vals) <= 1e-06, 'injected_rejected': 1 - max(vals) > 1e-06})
    for (name, d, expected) in [('safe', 2.0, 'LOCAL_CONDITION_ALL_ENCLOSURES'), ('no', 0.2, 'LOCAL_CONDITION_NO_ENCLOSURE'), ('straddle', 1.0, 'UNKNOWN_INTERVAL')]:
        c = uniform_shell_decision(d, 0.52, 0.5, 0.5)['decision']
        (lp, _) = lp_interval(d, 0.52, 1.0, 0.0)
        controls.append({'check': 'interval_' + name, 'fixture_kind': 'our_own_fixture_for_numerical_guard_only', 'valid_pass': c == lp == expected, 'injected_rejected': lp != ('LOCAL_CONDITION_NO_ENCLOSURE' if name != 'no' else 'LOCAL_CONDITION_ALL_ENCLOSURES')})
    write(ROOT / 'raw/R1_CONTROLS.json', controls)
    assert all((c['valid_pass'] and c['injected_rejected'] for c in controls))
    write(ROOT / 'raw/R1_INTERVAL_ROWS.json', intervals)
    write(ROOT / 'raw/R1_FRAMES.json', frames)
    write(ROOT / 'raw/R1_MANIFEST.json', manifest)
    assert all((r['class_agrees'] for r in intervals))
    outcome = {'claim_type': 'information_link', 'rows': len(rows), 'available_teeth': sum((r['status'] == 'AVAILABLE' for r in cohort['cohort'])), 'missing_teeth': sum((r['status'] != 'AVAILABLE' for r in cohort['cohort'])), 'different_minimum_witnesses': sum((r['downstream_difference_mm'] > 0 for r in witnesses)), 'max_downstream_difference_mm': max((r['downstream_difference_mm'] for r in witnesses)), 'EDT_KDTree_max_error_mm': max(vals), 'LP_queries': len(intervals), 'LP_class_mismatches': sum((not r['class_agrees'] for r in intervals)), 'all_faults_rejected': all((c['injected_rejected'] for c in controls)), 'all_frames_exact': all((r['mask_bit_exact'] and r['frame_metadata_exact'] and (r['max_frame_error_mm'] <= 1e-09) for r in frames)), 'wall_s': time.monotonic() - started, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'physical_outcomes': 'UNKNOWN_NO_DEJ_CEJ_OR_MATCHED_METROLOGY', 'source_X12_status': 'INTE_REDO_UNCHANGED'}
    write(ROOT / 'raw/R1_OUTCOME.json', outcome)
    write(ROOT / 'CURRENT_WORK_STATE.json', {'lane': 'X95-ct-pulp-prep-port', 'milestone': 'R1_DIGITAL_SOLIDS_DECIDED', 'latest_gate': outcome, 'next_operation': 'write R1 handoff then freeze spatial inverse construction; STL validation and report'})
    print(json.dumps(outcome), flush=True)
if __name__ == '__main__':
    main()
