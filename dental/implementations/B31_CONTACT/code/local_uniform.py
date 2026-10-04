from common import *
from geometry import load_case, gap, wall_and_removal
from bounds import uniform, box, support
from scipy.sparse import csr_matrix
import time

def run():
    check_lock()
    start = time.perf_counter()
    (q, reuse) = parents()
    pr = read(ROOT / 'PREREG_R4.json')
    r3 = read(ROOT / 'raw/ROUND3.json')['rows']
    rows = []
    for (r, rr) in zip(pr['cohort'], r3):
        t = load_case(r)
        with np.load(rr['constraints_path']) as a:
            A = csr_matrix((a['A_data'], a['A_indices'], a['A_indptr']), shape=tuple(a['A_shape']))
            b = a['b']
        den = np.asarray(A @ t['taper']).ravel()
        positive = den > 0
        ratios = (b[positive] + 1e-08) / den[positive]
        cap = min(r['relief_cap_mm'], float(ratios.min()) if len(ratios) else r['relief_cap_mm'])
        local_r = dict(r, relief_cap_mm=max(0.0, cap))
        (d, uc) = uniform(q, t, local_r, pr['parameters'])
        (v, wall) = wall_and_removal(t, d, r)
        path = DATA / 'R4' / (r['key'] + '_local_uniform.npz')
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, vertices=v, faces=t['faces'], face_roles=t['face_roles'], relief_mm=d, predicted_gap=gap(t, d))
        fp = ROOT / ('FROZEN_PREDICTIONS_R4_' + r['key'] + '.json')
        prediction = dict(key=r['key'], artifact_sha256=sha(path), prereg_sha256=sha(ROOT / 'PREREG_R4.json'), local_cap_mm=cap, nominal_symdiff_mm2=q.compare(t['xy'], t['grid_faces'], gap(t, d), t['reference_gap'])['symdiff_mm2'], physical_measurement='NOT_RUN', source_used_in_design=True, before='Actual-mesh direct-height remeasurement and pose-box validation for this case')
        if not fp.exists():
            freeze(fp, prediction)
        elif read(fp)['artifact_sha256'] != prediction['artifact_sha256']:
            raise ValueError('R4 frozen artifact drift')
        actual = t['ceiling'] - reuse.VF.height(v[t['faces'][t['face_roles'] == 0]], t['xy'])
        c = q.compare(t['xy'], t['grid_faces'], actual, t['reference_gap'])
        v4 = q.Q.contact_map(t['xy'], t['grid_faces'], actual, t['reference_gap'])
        local_violation = float(np.max(A @ d - b, initial=0.0))
        ff = support(t)
        (R, _, _, _) = q.polygons(t['xy'], ff, t['reference_gap'])
        missing = 0.0
        for (j, p) in R.items():
            ids = ff[j]
            pp = q.Q.clip(p, t['original_gap'][ids] - 0.1)
            (ar, _) = q.measure(pp, t['xy'][ids])
            if np.all(t['original_gap'][ids] == 0.1):
                ar = 0.0
            missing += ar
        bounds = {str(x): box(q, t, actual, t['reference_gap'], x, pr['parameters']) for x in [0.01, 0.05]}
        poses = {str(s): q.compare(t['xy'], t['grid_faces'], actual + s, t['reference_gap'] + s) for s in [-0.05, -0.01, 0.0, 0.01, 0.05]}
        nominal = uc['symdiff_optimum_bracket_mm2'][0] > 1e-07 and rr['outputs']['regional']['contact']['symdiff_mm2'] <= 0.5 * uc['symdiff_optimum_bracket_mm2'][0] + 1e-07
        reg = rr['outputs']['regional']
        boxguard = all((bounds[str(x)]['all_offsets_symdiff_lower_mm2'] > 1e-07 and reg['boxes'][str(x)]['worst_symdiff_bracket_mm2'][1] <= 0.5 * bounds[str(x)]['all_offsets_symdiff_lower_mm2'] + 1e-07 for x in [0.01, 0.05]))
        pose_failure = False
        with np.load(reg['path']) as a:
            rgap = a['predicted_gap']
        for s in [-0.05, -0.01, 0.0, 0.01, 0.05]:
            rc = q.compare(t['xy'], t['grid_faces'], rgap + s, t['reference_gap'] + s)
            bc = poses[str(s)]['symdiff_mm2']
            if bc > 1e-07 and rc['symdiff_mm2'] > 0.5 * bc + 1e-07:
                pose_failure = True
        controls = dict(local_uniform_pass=local_violation <= 1e-08 + 1e-14, area_parity_pass=abs(c['symdiff_mm2'] - v4['contact_symdiff_mm2']) <= 1e-07, missing_floor_below_all_outputs=all((o['contact']['symdiff_mm2'] + 1e-07 >= missing for o in rr['outputs'].values())) and c['symdiff_mm2'] + 1e-07 >= missing, cap_fault_rejected=None)
        if len(ratios) and cap < r['relief_cap_mm']:
            badcap = cap + 0.001
            controls['cap_fault_rejected'] = bool(np.max(den * badcap - b, initial=0.0) > 1e-08)
        else:
            controls['cap_fault_rejected'] = bool(np.max(den * (r['relief_cap_mm'] + 0.001) - b, initial=0.0) > 1e-08) or r['relief_cap_mm'] + 0.001 > r['relief_cap_mm']
        assert all(controls.values())
        row = dict(key=r['key'], case_key=r['case_key'], uniform_local_cap_mm=cap, original_global_cap_mm=r['relief_cap_mm'], uniform_local_constraint_violation_mm=local_violation, uniform_profile=uc, uniform_contact=c, uniform_boxes=bounds, uniform_poses=poses, wall=wall, cap_independent_missing_source_area_mm2=missing, source_band_equality_answer='NEW_MANUFACTURE_REQUIRED_FOR_ANY_FIXED_XY_DOWNWARD_OPERATION' if missing > 1e-07 else 'NO_POSITIVE_HEIGHT_WITNESS', regional_R3_contact_mm2=reg['contact']['symdiff_mm2'], nominal50pct_vs_local_proxy=nominal, observed_pose_ratio_failure=pose_failure, robust_sufficient_ratio_guard=boxguard, case_gate='FAIL' if not nominal or pose_failure else 'UNKNOWN_MACHINE_SOLID_AND_RATIO' if not boxguard else 'CONDITIONAL_DIGITAL_RATIO', controls=controls, artifact_path=str(path), artifact_sha256=sha(path), resolution=dict(gap='PER_POINT', mask='PER_SURFACE_REGION'), physical_status='UNKNOWN')
        rows.append(row)
        dump(ROOT / 'raw/R4_PROGRESS.json', dict(rows=rows))
        print('local uniform', r['key'], 'cap', cap, 'nominal50%', nominal, 'posefail', pose_failure, flush=True)
    nominal_count = sum((r['nominal50pct_vs_local_proxy'] for r in rows))
    possible_count = sum((r['case_gate'] != 'FAIL' for r in rows))
    summary = dict(claim_type='capability', retained=10, case_clusters=4, nominal_50pct_vs_new_local_proxy=nominal_count, cases_not_felled_by_observed_poses=possible_count, cap_independent_positive_height_witnesses=sum((r['cap_independent_missing_source_area_mm2'] > 1e-07 for r in rows)), inherited50pct3of10='FAIL' if possible_count < 3 else 'UNKNOWN_MACHINE_SOLID_AND_ROBUST_RATIO', old_R1_R2_gate='FAIL_PRESERVED', seconds=time.perf_counter() - start, rigorous_machine_enclosure='MISSING', subset_original_validity='UNKNOWN')
    dump(ROOT / 'raw/ROUND4.json', dict(summary=summary, rows=rows))
    state('R4_DECIDED', summary, 'Seal complete demo with copied replay; next physical/source-overlay construction named in HANDOFF')
    print(json.dumps(summary), flush=True)
if __name__ == '__main__':
    run()
