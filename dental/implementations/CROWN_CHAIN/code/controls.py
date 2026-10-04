from dental_release.paths import expand as _release_expand
from common_chain import *
import copy, time
import numpy as np
from fractions import Fraction as Q
from contracts import join, field_query, compare_observation

def run():
    rows = []

    def check(name, good, reject, scope):
        rows.append(dict(name=name, valid_case_accepted=bool(good), injected_fault_rejected=bool(reject), scope=scope))
    r1 = read(ROOT / 'raw/R1_ALL.json')['rows']
    pr = read(ROOT / 'PREREG_R1.json')
    sys.path.insert(0, _release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_FULL_CROWN_R5/deps'))
    import igl, trimesh
    from scipy.optimize import linprog
    oracle = []
    for r in r1 + read(ROOT / 'raw/R3_ALL.json')['rows']:
        with np.load(r['mesh_path'], allow_pickle=False) as a:
            (v, f) = (a['vertices'], a['faces'])
            roles = a['roles'] if 'roles' in a else a['face_roles']
        ext_ids = np.flatnonzero(roles == 0)
        inner = f[roles == 1]
        selected = ext_ids[np.linspace(0, len(ext_ids) - 1, min(64, len(ext_ids))).astype(int)]
        points = v[f[selected]].mean(1)
        (squared, _, _) = igl.point_mesh_squared_distance(points, v, inner)
        mesh = trimesh.Trimesh(v, inner, process=False)
        (_, td, _) = trimesh.proximity.closest_point(mesh, points)
        error = float(np.max(np.abs(np.sqrt(np.maximum(0, squared)) - td)))
        check('independent_triangle_distance_' + r['key'], error <= 1e-07, np.max(np.abs(np.sqrt(np.maximum(0, squared)) - (td + 1))) > 1e-07, 'libigl versus trimesh actual64 source points,1mm injection')
        entry = dict(key=r['key'], distance_oracle_max_difference_mm=error, points=len(points))
        if r in r1:
            normals = -trimesh.Trimesh(v, f, process=False).face_normals[roles == 1]
            normals = normals[np.linalg.norm(normals, axis=1) > 0]
            N = np.unique(np.rint(normals * 10 ** 6).astype(np.int64), axis=0) / 10 ** 6
            infeasible = []
            for axis in range(3):
                for sign in (-1, 1):
                    bounds = [(-1, 1)] * 3
                    bounds[axis] = (sign, sign)
                    lp = linprog(np.zeros(3), A_ub=-N, b_ub=np.full(len(N), 2e-06), bounds=bounds, method='highs')
                    infeasible.append(lp.status == 2)
            check('six_chart_LP_' + r['key'], all(infeasible), not all(infeasible[:-1] + [False]), 'six independent numerical primal charts versus rational Farkas infeasibility; removed rejection chart cannot prove no direction')
            entry['infeasible_six_charts'] = infeasible
        oracle.append(entry)
    dump(ROOT / 'raw/INDEPENDENT_CONTROLS.json', dict(records=oracle, libigl_module=str(igl.__file__), tolerance_mm=1e-07))
    for r in r1:
        error = r['replay_delta_mm']['shape']
        check('shape_' + r['key'], error <= 1e-07, abs(r['shape']['p95_mm'] + 1 - (r['shape']['p95_mm'] - error)) > 1e-07, 'actual1mm corrupted replay scalar')
        good = r['margin_max_mm']
        check('margin_' + r['key'], good <= 0.025, good + 0.1 > 0.025, 'actual0.1mm corrupted margin')
    cone = module('controls_cone', BASE / 'LANE_NEXT_D_INSERTION_PROOF/code/cone.py')
    for r in r1:
        bad = copy.deepcopy(r['cone'])
        bad['certificates'][0]['y'] = ['0'] * len(bad['certificates'][0]['y'])
        check('Farkas_' + r['key'], cone.recheck(r['cone']), not cone.recheck(bad), 'exact positive support certificate changed to zero')
    (_, collision, _) = scoped_modules(BASE / 'LANE_X53_MILLING_TOOLS/code', ['common', 'collision', 'run_r1'])
    source = next((r for r in r1 if any((x['has_witness'] for x in r['tool_points']))))
    with np.load(source['mesh_path'], allow_pickle=False) as a:
        (v, f) = (a['vertices'], a['faces'])
    valid = np.linalg.norm(np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]]), axis=1) / 2 >= 1e-14
    scene = collision.Scene(v, f[valid])
    (co, collision, xrun) = scoped_modules(BASE / 'LANE_X53_MILLING_TOOLS/code', ['common', 'collision', 'run_r1'])
    tools = xrun.lib('vhf', 5)
    point = next((p for p in source['tool_points'] if p['has_witness']))
    w = next((t for t in point['tools'] if t['status'] == 'FOUND'))
    tool = next((t for t in tools if t['id'] == w['tool_id']))
    approach = float(np.linalg.norm(scene.bounds[1] - scene.bounds[0]) + 2 * tool['holder_diameter_mm'] + 2)
    good = collision.pose_clearance(scene, np.array(w['centre']), np.array(point['normal']), np.array(w['direction']), tool, approach)[0]
    bad = collision.pose_clearance(scene, np.array(w['centre']), np.array(point['normal']), -np.array(w['direction']), tool, approach)[0]
    check('actual_tool_approach', good, not bad, 'same mesh/point/envelope, reversed noncutting approach')
    port = dict(case_key='A', geometry_sha256='a' * 64, frame='A_tooth31_mm', unit='N', quantity='axial_tooth_force', resolution='PER_TOOTH', observation_state='loaded_after_seating', time_scale='SIMULTANEOUS', epistemic='EXTERNALLY_CALIBRATED')
    for (k, v) in [('case_key', 'B'), ('geometry_sha256', 'b' * 64), ('frame', 'other_jaw'), ('unit', 'MPa'), ('quantity', 'strain_signal'), ('resolution', 'POPULATION'), ('observation_state', 'dry_unseated'), ('time_scale', 'HANDOVER')]:
        bad = dict(port, **{k: v})
        check('join_' + k, join(port, port)['status'] == 'PASS', join(bad, port)['status'] == 'FAIL', 'typed join fixture; no fixture counts as external measurement')
    model = module('chain_field_model', PKG / 'batch18/demos/FALT_TANDLAST/code/model.py')
    crown = dict(case_key='fixture', geometry_sha256='a' * 64, frame='fixture', fdi=11)
    for (cls, edges) in [('MUST', [(11, 41, Q(0), Q(0)), (12, 42, Q(1), Q(1))]), ('CAN', [(11, 41, Q(2), Q(2)), (12, 42, Q(1), Q(1))]), ('NEVER', [(12, 42, Q(0), Q(0))])]:
        graph = dict(crown, complete_declared_pair_graph=True, edges=edges)
        ans = field_query(crown, graph, model)
        bad = dict(graph, geometry_sha256='b' * 64)
        try:
            field_query(crown, bad, model)
            rejected = False
        except ValueError:
            rejected = True
        check('Field_' + cls, ans['classification'] == cls, rejected, 'reviewed exact classifier on three rational test graphs; absent patient graph remains UNKNOWN')
    check('Field_missing_is_unknown', field_query(crown)['classification'] == 'UNKNOWN', field_query(crown)['classification'] != 'NEVER', 'no missing-data-to-zero-load promotion')
    freeze = dict(key='fixture', design_mesh_sha256='a' * 64, quantity='axial_force', unit='N', frame='fixture', measurement_state='after_seating', prediction_interval=[9, 11])
    obs = dict(freeze, value=10, absolute_error_bound=0.1, raw_measurement_sha256='c' * 64, calibration_locator='own fixture')
    check('lab_interval_alarm', compare_observation(freeze, obs)['status'] == 'CONSISTENT_AT_DECLARED_ERRORS', compare_observation(freeze, dict(obs, value=20))['status'] == 'ALARM_DISJOINT_INTERVALS', 'synthetic comparator self-check only')
    check('lab_wrong_unit', compare_observation(freeze, obs)['status'] == 'CONSISTENT_AT_DECLARED_ERRORS', compare_observation(freeze, dict(obs, unit='kN'))['status'] == 'REJECTED_BINDING', 'units reject')
    check('lab_unpredicted_force', compare_observation(dict(freeze, prediction_interval=None), obs)['status'] == 'UNKNOWN_NO_PHYSICAL_PREDICTION', compare_observation(dict(freeze, prediction_interval=None), obs)['status'] != 'CONSISTENT_AT_DECLARED_ERRORS', 'null prediction cannot be passed')
    assert all((x['valid_case_accepted'] and x['injected_fault_rejected'] for x in rows))
    out = dict(controls=rows, all_pass=True, external_facit=False, scope='Contract and numerical fault controls; physical validation is separate')
    dump(ROOT / 'raw/CONTROLS.json', out)
    return out
if __name__ == '__main__':
    run()
