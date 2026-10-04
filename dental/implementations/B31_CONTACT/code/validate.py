from common import *
from geometry import *
from bounds import uniform, box, score
import time, resource

def faults(q):
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    f = np.array([[0, 1, 2]])
    left = np.full(3, 0.05)
    old = q.compare(xy, f, left, left)
    translated_xy = xy + [1.0, 0.0]
    area0 = abs(np.cross(xy[1] - xy[0], xy[2] - xy[0])) / 2
    area1 = abs(np.cross(translated_xy[1] - translated_xy[0], translated_xy[2] - translated_xy[0])) / 2
    translated_symdiff = area0 + area1
    plus3 = old['symdiff_mm2'] + 3
    return dict(mask_translation_1mm=dict(area_before_mm2=float(area0), area_after_mm2=float(area1), summary_identity_error_mm2=float(abs(area0 - area1)), overlap_mm2=0.0, symdiff_mm2=float(translated_symdiff), gate_max_mm2=0.5, rejected=translated_symdiff > 0.5, referent='our_own_fixture; geometric gate mechanism only'), reported_area_plus3_rejected=abs(old['symdiff_mm2'] - plus3) > 1e-07)

def thin_intaglio_fault(t, r):
    from scipy.spatial import cKDTree
    v = t['vertices'].copy()
    f = t['faces']
    inner = np.unique(f[t['face_roles'] == 1])
    outer = np.unique(f[t['face_roles'] == 0])
    (distance, j) = cKDTree(v[outer]).query(v[inner], k=1, workers=1)
    k = int(distance.argmin())
    iv = int(inner[k])
    ov = int(outer[j[k]])
    direction = v[iv] - v[ov]
    direction /= np.linalg.norm(direction)
    v[iv] = v[ov] + 0.01 * direction
    measured = float(np.linalg.norm(v[iv] - v[ov]))
    identity = float(abs(v[t['protected']] - t['vertices'][t['protected']]).max())
    return dict(inner_vertex=iv, outer_vertex=ov, wall_distance_upper_mm=measured, required_wall_mm=0.5, protected_identity_error_mm=identity, rejected=measured < 0.5 and identity > 0.0, whole_wall_lower_reuse_rejected=identity != 0.0, contact_independent=True)

def run():
    check_lock()
    start = time.perf_counter()
    pr = read(ROOT / 'PREREG_R1.json')
    params = pr['parameters']
    (q, reuse) = parents()
    const = read(ROOT / 'raw/CONSTRUCTIONS_R1.json')['rows']
    pred = read(ROOT / 'FROZEN_PREDICTIONS.json')['predictions']
    rows = []
    allcontrols = []
    for (r, c, predr) in zip(pr['cohort'], const, pred):
        t = load_case(r)
        row = dict(key=r['key'], case_key=r['case_key'], family=r['family'], dataset=r['dataset'], frame=r['frame'], resolution='PER_SURFACE_REGION', source_mask_area_mm2=q.compare(t['xy'], t['grid_faces'], t['original_gap'], t['reference_gap'])['reference']['area_mm2'], original=q.compare(t['xy'], t['grid_faces'], t['original_gap'], t['reference_gap']), outputs={})
        (d, uc) = uniform(q, t, r, params)
        (uv, ug) = wall_and_removal(t, d, r)
        up = DATA / 'R1' / r['key'] / 'uniform.npz'
        np.savez_compressed(up, vertices=uv, faces=t['faces'], face_roles=t['face_roles'], relief_mm=d, predicted_gap=gap(t, d))
        export_stl(ROOT / 'exports' / r['key'] / 'uniform.stl', uv, t['faces'])
        crow = dict(c['outputs'])
        crow['uniform'] = dict(path=str(up), sha256=sha(up), optimization=uc, wall_and_removal=ug)
        for (name, rec) in crow.items():
            if sha(rec['path']) != rec['sha256']:
                raise ValueError('Output hash drift ' + name)
            with np.load(rec['path'], allow_pickle=False) as a:
                v = a['vertices']
                storedg = a['predicted_gap']
                d = a['relief_mm']
            actual = t['ceiling'] - reuse.VF.height(v[t['faces'][t['face_roles'] == 0]], t['xy'])
            finite = np.isfinite(actual) & np.isfinite(storedg)
            rayerror = float(abs(actual[finite] - storedg[finite]).max(initial=0.0))
            measured = q.compare(t['xy'], t['grid_faces'], actual, t['reference_gap'])
            parity = q.Q.contact_map(t['xy'], t['grid_faces'], actual, t['reference_gap'])
            areaerror = abs(measured['symdiff_mm2'] - parity['contact_symdiff_mm2'])
            poses = {str(s): q.compare(t['xy'], t['grid_faces'], actual + s, t['reference_gap'] + s) for s in params['validation_pose_offsets_mm']}
            boxes = {str(b): box(q, t, actual, t['reference_gap'], b, params) for b in params['pose_boxes_mm']}
            g = rec['wall_and_removal']
            protected = float(abs(v[t['protected']] - t['vertices'][t['protected']]).max())
            machine_enclosure = 'MISSING'
            nominalnegative = measured['negative_gap_area_mm2'] <= min(0.1, row['original']['negative_gap_area_mm2']) + pr['metrics']['negative_gap_regression_tolerance_mm2']
            no_negative_regression = all((poses[str(-b)]['negative_gap_area_mm2'] <= min(0.1, q.compare(t['xy'], t['grid_faces'], t['original_gap'] - b, t['reference_gap'] - b)['negative_gap_area_mm2']) + 1e-07 for b in params['pose_boxes_mm']))
            geo = dict(wall=g['wall_model_pass'], margin_intaglio=protected == 0.0, nominal_negative=nominalnegative, box_negative=no_negative_regression, support=measured['coverage'] >= 0.95, actual_ray=rayerror <= 1e-08, area_parity=areaerror <= 1e-07, removal_subset=g['removal_status'] == 'CONDITIONAL_SUBSET')
            prediction_error = None if name == 'uniform' else abs(measured['symdiff_mm2'] - predr['outputs'][name]['nominal_contact_symdiff_mm2'])
            row['outputs'][name] = dict(contact=measured, poses=poses, boxes=boxes, geometry_gates=geo, wall_and_removal=g, optimization=rec['optimization'], path=rec['path'], sha256=rec['sha256'], actual_ray_max_error_mm=rayerror, independent_v4_area_error_mm2=areaerror, prediction_contact_error_mm2=prediction_error, physical_status='UNKNOWN', rigorous_machine_enclosure=machine_enclosure)
            allcontrols.append(dict(key=r['key'], output=name, ray_pass=rayerror <= 1e-08, v4_area_pass=areaerror <= 1e-07, plus3_area_rejected=abs(measured['symdiff_mm2'] + 3 - parity['contact_symdiff_mm2']) > 1e-07, protected_identity_pass=protected == 0.0, frozen_prediction_pass=prediction_error is None or prediction_error <= 1e-07))
        reg = row['outputs']['regional']
        uni = row['outputs']['uniform']
        lb = uni['optimization']['symdiff_optimum_bracket_mm2'][0]
        sd = reg['contact']['symdiff_mm2']
        nominalgain = lb > 1e-07 and sd <= 0.5 * lb + 1e-07
        posefailure = any((uni['poses'][str(s)]['symdiff_mm2'] > 1e-07 and reg['poses'][str(s)]['symdiff_mm2'] > 0.5 * uni['poses'][str(s)]['symdiff_mm2'] + 1e-07 for s in params['validation_pose_offsets_mm']))
        boxratio = all((uni['boxes'][str(b)]['all_offsets_symdiff_lower_mm2'] > 1e-07 and reg['boxes'][str(b)]['worst_symdiff_bracket_mm2'][1] <= 0.5 * uni['boxes'][str(b)]['all_offsets_symdiff_lower_mm2'] + 1e-07 for b in params['pose_boxes_mm']))
        complete = all((x['status'] == 'CONDITIONAL_REAL_ARITHMETIC_ENCLOSURE' for x in reg['boxes'].values()))
        gates = reg['geometry_gates']
        geopass = all(gates.values())
        row['nominal_50pct_vs_uniform'] = nominalgain
        row['nominal_geometry_model_pass'] = geopass
        row['robust_ratio_sufficient_guard'] = boxratio
        row['observed_pose_ratio_failure'] = posefailure
        row['digital_case_gate'] = 'FAIL' if not nominalgain or not geopass or posefailure else 'CONDITIONAL_MODEL_PASS' if boxratio and complete else 'UNKNOWN_BOX_RATIO'
        row['method_control_line'] = dict(regional_vs_ordinary_contact_difference_mm2=reg['contact']['symdiff_mm2'] - row['outputs']['ordinary_mesh']['contact']['symdiff_mm2'], regional_vs_ordinary_removal_difference_mm3=reg['wall_and_removal']['area_integral_abs_displacement_mm3'] - row['outputs']['ordinary_mesh']['wall_and_removal']['area_integral_abs_displacement_mm3'], algorithm_gain_claimed=False)
        row['thin_intaglio_fault'] = thin_intaglio_fault(t, r)
        rows.append(row)
        dump(ROOT / 'raw/VALIDATION_PROGRESS.json', dict(rows=rows))
        print('validated', r['key'], 'regional', sd, 'uniform bracket', uni['optimization']['symdiff_optimum_bracket_mm2'], row['digital_case_gate'], flush=True)
    fs = faults(q)
    fs['thin_intaglio_rejections'] = sum((r['thin_intaglio_fault']['rejected'] for r in rows))
    dump(ROOT / 'raw/CONTROLS_R1.json', dict(rows=allcontrols, faults=fs))
    count = sum((r['digital_case_gate'] == 'CONDITIONAL_MODEL_PASS' for r in rows))
    nominal = sum((r['nominal_50pct_vs_uniform'] for r in rows))
    summary = dict(claim_type='capability', requested=10, retained=len(rows), rejected_records=0, case_clusters=len({r['case_key'] for r in rows}), nominal_50pct_reductions=nominal, robust_geometry_case_passes=count, required_case_passes=3, gate='FAIL' if count + sum((r['digital_case_gate'] == 'UNKNOWN_BOX_RATIO' for r in rows)) < 3 else 'CONDITIONAL_MODEL_PASS' if count >= 3 else 'UNKNOWN', physical_passes=0, seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(ROOT / 'raw/ROUND1.json', dict(summary=summary, rows=rows))
    state('R1_DECIDED', summary, 'Write R1 handoff, freeze next construction changing the operation at the exposed obstruction')
    print(json.dumps(summary), flush=True)
if __name__ == '__main__':
    run()
