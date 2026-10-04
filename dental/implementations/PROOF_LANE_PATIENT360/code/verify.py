from dental_release.paths import expand as _release_expand
import sys, copy
import numpy as np
from shared import *
run = Path(sys.argv[1])
checks = []

def record(name, good, bad, scope):
    checks.append(dict(check=name, valid_pass=bool(good), fault_rejected=bool(bad), scope=scope))
    if not (good and bad):
        raise AssertionError(name)
ios = run / _release_expand('@DENTAL_CASE_C@')
frame = load(ios / 'FRAME.json')
rr = np.array(frame['local_to_native_rotation'])
c = np.array(frame['native_origin_xyz_mm'])
p = np.array([[2.0, 3.0, 7.0], [-10.0, 7.0, -1.0]])
error = lambda tr: float(np.max(np.abs((p - c) @ rr @ tr.T + c - p)))
record('frame_inverse', error(rr) < 1e-08, error(rr + np.eye(3) * 0.01) > 1e-08, 'PER_POINT; algebra, not scanner accuracy')
new = load(ios / 'X21_REPLAY.json')['pairs']
old = load(RESULTS / _release_expand('LANE_X21_CONTACT_MAP/raw/cases/@DENTAL_SURFACE_ID_C@.json'))['pairs']
lookup = {(r['upper_fdi'], r['lower_fdi']): r for r in old}
err = max((abs(r['minimum_projected_gap_mm'] - lookup[r['upper_fdi'], r['lower_fdi']]['minimum_projected_gap_mm']) for r in new))
record('X21_original_cached_geometry_replay', err <= 0.005, err + 0.1 > 0.005, 'PER_SURFACE_REGION; source-mesh replay, not external physical contact')
cr = load(ios / 'CROWN_RESULTS.json')
ff = cr['FE']
accept = lambda r: r['force_balance_N'] <= 1e-06 and r['equilibrium_residual_N'] <= 1e-06
record('FE_equilibrium', all((accept(r) for r in ff.values())), not accept({**ff['100'], 'force_balance_N': 1.0}), 'PHENOMENOLOGICAL load, PER_POINT numerical response')
record('FE_linear_response', cr['FE_linearity_error'] <= 1e-08, abs((ff['100']['maximum_tensile_MPa'] + 1) / (2 * ff['50']['maximum_tensile_MPa']) - 1) > 1e-08, 'numerical consistency only')
for case in [_release_expand('@DENTAL_CASE_A@'), _release_expand('@DENTAL_CASE_B@')]:
    x = load(run / case / 'CBCT_RESULTS.json')
    co = x['controls']
    record(case + '_CT_identity', x['identity']['exact_CT_pair'], co['CT_corruption_rejected'], 'Full voxel equality; external paired published CT')
    record(case + '_pulp_distance', co['EDT_KDTree_max_error_mm'] < 1e-08, co['EDT_injected_1mm_rejected'], 'PER_POINT external label geometry')
    pred = load(run / case / 'FROZEN_PREDICTIONS.json')
    g = pred['virtual_implant']['geometry']
    op = module('full_geometry', RESULTS / 'LANE_X8_GUIDE_NERVE_RISK/full_geometry.py')
    vp = pred['virtual_implant']
    w = np.array(g['witness']['witness_zyx_mm'])
    q = op.project_cylinder(w, np.array(vp['entry_zyx_mm']), np.array(vp['axis_zyx']), vp['length_mm'], vp['radius_mm'])
    wd = float(np.linalg.norm(w - q))
    measure_error = abs(wd - g['upper_mm'])
    record(case + '_cylinder_witness', measure_error < 1e-08, abs(wd - (g['upper_mm'] + 1)) > 1e-08, 'PER_SURFACE_REGION actual nearest witness to finite cylinder')
    record(case + '_guide_numeric_control', co['guide_scalar_max_error'] < 1e-08, co['guide_scalar_max_error'] + 0.1 > 1e-08, 'POPULATION closure; not patient-risk validation')
r2 = load(run / 'R2/RESULTS_R2.json')
e = r2['export']
record('X38_units', e['export']['status'] == 'PASS', e['controls']['unit_fault_rejected'], 'geometry transport')
record('X38_vertex', e['facet_coordinate_change_mm'] < 1e-05, e['controls']['changed_vertex_rejected'], 'geometry transport')
for t in r2['thermal']:
    if t['status'] == 'SIMULATED_CONDITIONAL':
        ok = max((s['relative_energy_balance'] for s in t['scenarios'])) <= 1e-07
        record(t['patient_id'] + '_thermal_conservation', ok, not 0.01 <= 1e-07, 'Energy ledger, not physical thermal validation')
        v = t['numerical_controls']
        record(t['patient_id'] + '_thermal_slab', v['slab_pass'], v['flux_x1p1_rejected'], 'External closed-form slab, not anatomy')
r3 = load(run / 'R3/RESULTS_R3.json')
record('shared_pose_direct_control', r3['maximum_direct_control_error_mm'] < 1e-08, r3['maximum_direct_control_error_mm'] + 0.1 > 1e-08, 'Directly translated original source triangles')
for k in ['patient_id', 'frame_sha256', 'unit', 'source_sha256']:
    record('port_' + k, r3['ports']['valid_port_pass'], r3['ports'][k + '_fault_rejected'], 'Actual guarded edge; no clinical accuracy assertion')
exp = load(DATA / run.name / 'R3/export/FROZEN_PREDICTIONS.json')
expected = load(run / 'R3/FROZEN_PREDICTIONS.json')
accept_export = lambda x: x['patient_id'] == expected['patient_id'] and x['theta_rows'] == expected['theta_rows'] and (x['frame'] == expected['frame'])
bad = copy.deepcopy(exp)
bad['theta_rows'][0]['shared_theta']['IOS_pose_delta_mm'] += 0.1
record('export_shared_theta_identity', accept_export(exp), not accept_export(bad), 'Actual exported frozen theta and patient frame')
foreign = load(RESULTS / 'LANE_X1B_CROWN_LOOP/FROZEN_PREDICTIONS.json')
accept_patient_prediction = lambda x: x.get('patient_id') == _release_expand('@DENTAL_CASE_C@') and x.get('frame') == expected['frame']
record('foreign_X1b_blocked', accept_patient_prediction(exp), not accept_patient_prediction(foreign), 'Actual frozen X1b lacks this patient/frame; excluded')
for b in load(run / 'R4/RESULTS_R4.json')['patients']:
    record(b['patient_id'] + '_bone_exact_intervals', b['exact_gate'] == 'PASS', b['fault_injections_rejected'], 'Independent cube slab vs full grid-face label intervals; original sampled gate FAIL retained')
dump(run / 'VERIFICATION.json', dict(status='PASS', checks=checks, controls=len(checks), X21_replay_error_mm=err, scope='Computational consistency and fault rejection, not independent scientific review'))
print('Verified', len(checks), 'controls', flush=True)
