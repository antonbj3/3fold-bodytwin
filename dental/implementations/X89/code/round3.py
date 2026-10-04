from dental_release.paths import expand as _release_expand
from codesign import *
from milling_gate import evaluate, parent_check
from sliver_scene import SliverScene
from calibration import displacement_bound
import copy

def build_contracts():
    pre = lock_prereg('R1')
    card = PARENT / 'NOMINAL_TIP_REFERENCE_CARD.json'
    rows = intake()
    out = []
    for (ix, row) in enumerate(rows):
        c = dict(units='mm', crown_sha256=row['source_sha256'], geometry_file=row['file'], geometry_sha256=row['sha256'], source_link_status='Inherited X53 normalized mesh or R4 export npz linked by parent source hash', region_faces=row['regions'], region_status=row['region_status'], tool_card=card, tool_card_sha256=sha(card), scenario=pre['metrics'], seed=pre['metrics']['seed'] + ix, samples_per_region=pre['metrics']['n_per_region'], measurement_locators={}, physical_scope='NOT_CALIBRATED; no commercial CAM fixture or observed process bounds')
        p = ROOT / 'contracts' / f"{row['id']}.json"
        dump(p, c)
        out.append(dict(id=row['id'], crown=row['source'], contract=p))
    dump(ROOT / 'contracts/INDEX.json', out)
    return out

def port_controls(m):
    import trimesh
    box = trimesh.creation.box(extents=[4, 4, 1])
    scene = SliverScene(box.vertices, box.faces)
    t = read(PARENT / 'NOMINAL_TIP_REFERENCE_CARD.json')['tools'][1]
    p = np.array([0.0, 0.0, 0.5])
    n = np.array([0.0, 0.0, 1.0])
    nominal = robust_search(scene, p, n, t, m)
    pose_id = t['id'] + ':12'
    port = dict(pose_id=pose_id, force_levels_N=[0.5, 1.0], tip_displacement_lower_mm=[0.0015, 0.0035], tip_displacement_upper_mm=[0.0025, 0.0045], curvature_bound_mm_N2=0.0008, curvature_source='our_own_fixture qualified by construction', data_kind='our_own_fixture')
    positive = robust_search(scene, p, n, t, m, {pose_id: port})
    missing = robust_search(scene, p, n, t, m, {pose_id: dict(port, curvature_bound_mm_N2=None)})
    wrong = robust_search(scene, p, n, t, m, {pose_id: dict(port, pose_id='wrong')})
    high = robust_search(scene, p, n, t, m, {pose_id: dict(port, tip_displacement_lower_mm=[0.05, 0.05], tip_displacement_upper_mm=[0.05, 0.05])})
    out = dict(nominal=nominal, positive=positive, missing_curvature=missing, wrong_pose=wrong, injected50um=high, nominal_rejected=nominal['status'] == 'MECHANICS_REJECTED', positive_pass=positive['status'] == 'CONDITIONAL_LOCAL_WITNESS', upper_stock_error_mm=positive['stock_interval_final_mm'][1] if positive['stock_interval_final_mm'] else None, missing_curvature_rejected=missing['status'] != 'CONDITIONAL_LOCAL_WITNESS', wrong_pose_rejected=wrong['status'] != 'CONDITIONAL_LOCAL_WITNESS', injected50um_rejected=high['status'] != 'CONDITIONAL_LOCAL_WITNESS', data_kind='our_own_fixture', physical_validation='NOT_ACQUIRED')
    assert all((out[k] for k in ('nominal_rejected', 'positive_pass', 'missing_curvature_rejected', 'wrong_pose_rejected', 'injected50um_rejected')))
    assert abs(out['upper_stock_error_mm'] - 0.017201) < 1e-12
    left = trimesh.creation.box(extents=[1.9, 4, 4])
    left.apply_translation([-1.05, 0, 0])
    right = trimesh.creation.box(extents=[1.9, 4, 4])
    right.apply_translation([1.05, 0, 0])
    slit = trimesh.util.concatenate([left, right])
    ss = SliverScene(slit.vertices, slit.faces)
    sr = [robust_search(ss, np.array([-0.1, 0, 0]), np.array([1.0, 0, 0]), tt, m) for tt in read(PARENT / 'NOMINAL_TIP_REFERENCE_CARD.json')['tools']]
    out['injected_0p2mm_slit'] = sr
    out['injected_0p2mm_slit_rejected'] = all((x['status'] != 'CONDITIONAL_LOCAL_WITNESS' for x in sr))
    assert out['injected_0p2mm_slit_rejected']
    return out

def run():
    lock_prereg('R3')
    start = time.perf_counter()
    rows = build_contracts()
    r2 = read(ROOT / 'RESULTS_R2.json')
    old = {x['id']: x for x in r2['designs']}
    reports = []
    for row in rows:
        state('R3_RUNNING', f'{len(reports)}/8 X49 rule queries', 'Verify data contracts and local rule; probe calibrated response coupling')
        report = evaluate(row['crown'], row['contract'])
        p = ROOT / 'raw/gate' / f"{row['id']}.json"
        dump(p, report)
        assert report['found'] == old[row['id']]['found'] and report['points'] == old[row['id']]['points']
        assert report['status'] == 'UNKNOWN'
        reports.append(dict(id=row['id'], status=report['status'], scenario_verdict=report['scenario_verdict'], points=report['points'], found=report['found'], report_file=p))
    m = r2['metrics']
    port = port_controls(m)
    faults = []
    row = rows[0]
    c = read(row['contract'])
    for (name, changes) in [('force100N', dict(force_N=100)), ('holder1mmN', dict(holder_compliance_mm_N=1)), ('reserve20um', dict(other_geometry_error_final_mm=0.02))]:
        cc = copy.deepcopy(c)
        cc['scenario'].update(changes)
        p = ROOT / 'raw/injections' / f'{name}.json'
        dump(p, cc)
        rr = evaluate(row['crown'], p)
        dump(ROOT / 'raw/injections' / f'{name}_result.json', rr)
        faults.append(dict(id=name, rejected=rr['status'] == 'FAIL', actual_status=rr['status'], found=rr['found'], points=rr['points']))
    cc = copy.deepcopy(c)
    cc['crown_sha256'] = '0' * 64
    p = ROOT / 'raw/injections/source_hash.json'
    dump(p, cc)
    try:
        evaluate(row['crown'], p)
        rejected = False
    except ValueError as ex:
        rejected = 'mismatch' in str(ex)
    faults.append(dict(id='source_hash', rejected=rejected))
    original = np.load(c['geometry_file'], allow_pickle=False)
    mutant = DATA / 'injections/shifted50um.npz'
    mutant.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(mutant, vertices=original['vertices'] + np.array([0.05, 0, 0]), faces=original['faces'])
    cc = copy.deepcopy(c)
    cc.update(geometry_file=str(mutant), geometry_sha256=sha(mutant))
    p = ROOT / 'raw/injections/coherent_geometry_mismatch.json'
    dump(p, cc)
    try:
        evaluate(row['crown'], p)
        rejected = False
    except ValueError as ex:
        rejected = 'coordinates mismatch' in str(ex)
    faults.append(dict(id='coherent_geometry_shift50um', rejected=rejected, mutant_sha256=sha(mutant)))
    src = read(ROOT / 'raw/EXTERNAL_FACIT.json')['rows']
    fake = [dict(x, mean_RMS_um=x['mean_RMS_um'] / 10) for x in src]
    check = lambda x: all((a['mean_RMS_um'] == b['mean_RMS_um'] for (a, b) in zip(x, src)))
    faults.append(dict(id='external_table_underreported_10x', baseline_pass=check(src), rejected=not check(fake)))
    parent = []
    ant = Path(_release_expand('@DENTAL_WORK_ROOT@/X49_design_gate/X1b_missing_antagonist_surrogate.stl'))
    for case in ['D1', 'M1']:
        d = DENTAL / 'results/LANE_X1B_CROWN_LOOP/exports' / case
        try:
            rr = parent_check(dict(prep=str(d / 'preparation.stl'), antagonist=str(ant), crown=str(d / 'crown.stl')), '3Y', None, 'mm', 1.25)
            parent.append(dict(id=case, parent_verdict=rr['verdict'], combined_verdict='FAIL' if rr['verdict'] == 'FAIL' else 'UNKNOWN', rules=rr['rules']))
        except Exception as ex:
            parent.append(dict(id=case, status='PARENT_EXECUTION_FAILED', error=str(ex)))
    dump(ROOT / 'raw/PARENT_INTEGRATION.json', parent)
    assert all((x['rejected'] for x in faults))
    assert any((x.get('parent_verdict') == 'FAIL' and x['combined_verdict'] == 'FAIL' for x in parent)), 'Parent FAIL preservation untested'
    r3 = dict(claim_type='capability', round='R3', status='INSTALLED_RESPONSE_PORT_CONSUMED_BY_SWEEP; NO_PHYSICAL_ADMISSION', reports=reports, calibrated_port_controls=port, injected_controls=faults, parent_integration=[{k: v for (k, v) in x.items() if k != 'rules'} for x in parent], prereg_sha256=sha(ROOT / 'PREREG_R3.json'), external_referent=r2['external_referent'], cost=dict(total_wall_s=time.perf_counter() - start, fit_s=0, questions_s=0, fallback_s=None, physical_measurements=0, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, gpu=0), resolution='PER_POINT', time_scale='HANDOVER')
    dump(ROOT / 'RESULTS_R3.json', r3)
    final = copy.deepcopy(r2)
    final['rounds'] = dict(R1='RESULTS_R1.json', R2='RESULTS_R2.json', R3='RESULTS_R3.json')
    final['integrated_rule'] = r3
    final['status'] = 'RUNNABLE_CROWN_TO_LOCAL_TOOL_AND_STIFFNESS_REQUIREMENTS; PHYSICAL_UNKNOWN'
    dump(ROOT / 'results.json', final)
    state('R3_COMPLETE', '8 rule replays match; installed response affects sweep; impossible/faults reject and parent FAIL preserved', 'Finalize one-command demo, figure, pending graph feedback and exact next lab measurement')
    return final
if __name__ == '__main__':
    run()
