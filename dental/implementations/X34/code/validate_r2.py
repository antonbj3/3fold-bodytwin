"""Replay frozen cases and independently falsify every new witness checker."""
import copy
import json
import time
import sys
from pathlib import Path
import numpy as np
import trimesh
sys.path.insert(0, str(Path(__file__).resolve().parent))
from designgate.gate import check
from designgate.geometry import load, faces, continuous
from designgate.distance import nearest
from designgate.vendor.geometry import export_stl
from initialize import H, dump, sha
from prepare_demo import D

def paths(name):
    if name in ('rotation', 'micrometers', 'wrong_scale'):
        stem = 'rot_' if name == 'rotation' else 'um_'
        return {k: D / (stem + k + '.stl') for k in ('prep', 'crown', 'antagonist')}
    keys = ['prep', 'crown', 'antagonist'] + (['neighbor'] if name.startswith('proximal') else [])
    return {k: D / name / (k + '.stl') for k in keys}

def contract(name):
    return D / (name + '_R2.json') if name in ('rotation', 'micrometers', 'wrong_scale') else D / name / 'contract_R2.json'

def naive_distance(mesh, p):
    return float(trimesh.proximity.closest_point_naive(mesh, np.asarray(p)[None, :])[1][0])

def witness_valid(w, mesh, signed=False):
    p = np.asarray(w['source_point_mm'])
    q = np.asarray(w['target_point_mm'])
    d = w['distance_mm']
    if not np.isfinite(p).all() or not np.isfinite(q).all() or (not np.isfinite(d)):
        return False
    got = naive_distance(mesh, p)
    return abs(np.linalg.norm(p - q) - abs(d)) <= 1e-07 and abs(got - abs(d)) <= 1e-07

def solid_angle(point, tri):
    v = tri - np.asarray(point)
    (a, b, c) = (v[:, 0], v[:, 1], v[:, 2])
    (la, lb, lc) = (np.linalg.norm(x, axis=1) for x in (a, b, c))
    num = np.sum(a * np.cross(b, c), axis=1)
    den = la * lb * lc + np.sum(a * b, axis=1) * lc + np.sum(b * c, axis=1) * la + np.sum(c * a, axis=1) * lb
    return float(np.sum(2 * np.arctan2(num, den)))

def main():
    freeze = json.loads((H / 'FROZEN_PREDICTIONS_R2.json').read_text())
    assert sha(H / 'PREREG_R2.json') == freeze['prereg_sha256']
    for row in freeze['files']:
        assert sha(row['path']) == row['sha256']
    for (p, h) in freeze['code'].items():
        assert sha(H / p) == h
    start = time.perf_counter()
    tests = []
    reports = {}

    def record(n, b, detail):
        tests.append(dict(name=n, passed=bool(b), detail=detail))
    for (name, expected) in freeze['expected'].items():
        report = check(paths(name), contract=contract(name), round_version='R2')
        reports[name] = report
        for (k, v) in expected.items():
            record(name + '_' + k, report['rules'][k]['status'] == v, dict(expected=v, observed=report['rules'][k]))
    from designgate.distance import enclosure
    tri = np.array([[[-2, -2, 0.8], [2, -2, 0.8], [-2, 2, 0.8]], [[2, -2, 0.8], [2, 2, 0.8], [-2, 2, 0.8]]], float)
    target = trimesh.Trimesh([[-2, -2, 0], [2, -2, 0], [2, 2, 0], [-2, 2, 0], [0, 0, 0]], [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]], process=False)
    d = enclosure(tri, target, minimum_required=0.5)
    record('unmatched_tessellation', d['minimum_status'] == 'PASS', d)
    bad = copy.deepcopy(d['minimum_witness'])
    bad['distance_mm'] += 0.2
    record('distance_verifier_positive', witness_valid(d['minimum_witness'], target), d['minimum_witness'])
    record('distance_verifier_injected_bad_value', not witness_valid(bad, target), bad)
    bad['distance_mm'] = float('nan')
    record('distance_verifier_injected_NaN', not witness_valid(bad, target), 'Nonfinite distance rejected')
    points = np.array([[0, 0, 0.8], [1.1, 0.3, 0.8], [-2, -2, 0.8], [0, 0, -0.2]])
    (_, fast, _) = nearest(target, points)
    control = trimesh.proximity.closest_point_naive(target, points)[1]
    record('published_code_full_face_parity', max(abs(fast - control)) <= 1e-07, dict(fast_mm=fast.tolist(), all_face_control_mm=control.tolist()))
    ball = reports['stepped_cavity']['rules']['ball_milling']
    w = ball.get('witness', {})
    mesh = load(paths('stepped_cavity')['crown'], 1)['mesh']

    def ball_valid(w):
        try:
            center = np.array(w['ball_center_mm'])
            obs = np.array(w['obstacle_point_mm'])
            r = ball['radius_final_design_mm']
            actual = naive_distance(mesh, center)
            return np.isfinite(center).all() and abs(np.linalg.norm(center - obs) - w['obstacle_distance_to_center_mm']) <= 1e-07 and (abs(actual - w['obstacle_distance_to_center_mm']) <= 1e-07) and (actual < r - 1e-05)
        except (KeyError, ValueError):
            return False
    record('ball_witness_independent_naive', ball_valid(w), w)
    bad = copy.deepcopy(w)
    bad['obstacle_distance_to_center_mm'] += 0.5
    record('ball_verifier_injected_bad_value', not ball_valid(bad), bad)
    ins = reports['stepped_cavity']['rules']['insertion']
    w = ins.get('witness', {})

    def insertion_valid(w):
        try:
            p = np.array(w['point_in_crown_frame_mm'])
            q = np.array(w['preparation_point_mm'])
            t = np.array(w['crown_translation_mm'])
            ang = solid_angle(p, mesh.triangles)
            return np.linalg.norm(p - (q - t)) <= 1e-07 and abs(abs(ang) - 4 * np.pi) <= 1e-06 and (abs(naive_distance(mesh, p) - w['interior_depth_mm']) <= 1e-07)
        except (KeyError, ValueError):
            return False
    record('insertion_witness_independent_solid_angle', insertion_valid(w), dict(witness=w, solid_angle_sr=solid_angle(w['point_in_crown_frame_mm'], mesh.triangles) if w else None, locator='doi:10.1109/TBME.1983.325207'))
    bad = copy.deepcopy(w)
    bad['interior_depth_mm'] += 0.2
    record('insertion_verifier_injected_bad_value', not insertion_valid(bad), bad)
    meta = json.loads(contract('valid').read_text())
    meta.pop('preparation_outward_axis')
    meta['geometry_scope'] = 'general_3D'
    tmp = D / 'attack_missing_sign.json'
    dump(tmp, meta)
    r = check(paths('valid'), contract=tmp, round_version='R2')
    record('open_prep_unknown_sign', r['rules']['film_min']['status'] == 'UNKNOWN' and r['rules']['film_max']['status'] == 'UNKNOWN', r['rules']['film_min'])
    meta = json.loads(contract('valid').read_text())
    meta['geometry_scope'] = 'general_3D'
    meta['crown_solid_integrity_asserted'] = True
    meta['extraction_axis'] = [0, 0, 1]
    meta['extraction_travel_mm'] = 1
    meta['ball_radius_final_mm'] = 0.3
    base = load(paths('valid')['crown'], 1)
    dest = D / 'open_crown.stl'
    export_stl(dest, base['vertices'], base['faces'][:-1])
    pp = paths('valid')
    pp['crown'] = dest
    meta['input_sha256']['crown'] = sha(dest)
    tmp = D / 'attack_open_crown.json'
    dump(tmp, meta)
    r = check(pp, contract=tmp, round_version='R2')
    record('open_crown_no_milling_insertion_PASS', all((r['rules'][k]['status'] == 'UNKNOWN' for k in ['ball_milling', 'insertion'])), {k: r['rules'][k] for k in ['ball_milling', 'insertion']})
    meta['input_sha256']['crown'] = 'f' * 64
    tmp = D / 'attack_wrong_hash_R2.json'
    dump(tmp, meta)
    try:
        check(pp, contract=tmp, round_version='R2')
        rejected = False
    except ValueError:
        rejected = True
    record('R2_wrong_hash', rejected, 'Changed STL cannot reuse region mapping')
    nv = base['vertices'].copy()
    nv[0, 0] = np.nan
    dest = D / 'nonfinite.stl'
    export_stl(dest, nv, base['faces'])
    pp['crown'] = dest
    try:
        check(pp, units='mm', common_frame=True, round_version='R2')
        rejected = False
    except ValueError:
        rejected = True
    record('nonfinite_STL_rejected', rejected, 'NaN coordinates do not enter numerical checks')
    real = check(paths('real_X18_110_34'), contract=contract('real_X18_110_34'), round_version='R2')
    reports['real'] = real
    rr = real['rules']['material_wall']
    wi = rr.get('witness')
    cc = load(paths('real_X18_110_34')['crown'], 1)
    ids = json.loads(contract('real_X18_110_34').read_text())['regions']['exterior']
    (ot, oi) = faces(cc, ids)
    tar = trimesh.Trimesh(ot.reshape(-1, 3), np.arange(3 * len(ot)).reshape(-1, 3), process=False)
    record('real_wall_independent_full_face', wi is not None and witness_valid(wi, tar), dict(status=rr['status'], witness=wi))
    cg = real['rules']['occlusal_contact']['measurement']
    record('real_contact_frozen_target', abs(cg['minimum_gap_mm'] - 0.025) <= 1e-05, cg)
    (u, ui) = faces(load(paths('real_X18_110_34')['antagonist'], 1))
    (l, li) = faces(cc, ids)
    full = continuous(u, l, [0, 0, 1], False)
    record('real_continuous_unpruned_parity', abs(cg['minimum_gap_mm'] - full['minimum_gap_mm']) <= 1e-07, full)
    m = json.loads(contract('valid').read_text())
    m.pop('ifu_profile')
    tmp = D / 'attack_class_only_R2.json'
    dump(tmp, m)
    r = check(paths('valid'), contract=tmp, material='3Y', round_version='R2')
    record('R2_generic_3Y_UNKNOWN', r['rules']['material_wall']['status'] == 'UNKNOWN', r['rules']['material_wall'])
    result = dict(round='R2', claim_type='capability', tests=tests, reports=reports, total=len(tests), passed=sum((x['passed'] for x in tests)), all_gates_pass=all((x['passed'] for x in tests)), cost_wall_s=time.perf_counter() - start, physical_measurement_performed=False, new_metrology_calibration=False)
    dump(H / 'rounds/R2.json', result)
    print(json.dumps(dict(passed=result['passed'], total=result['total'], failed=[x['name'] for x in tests if not x['passed']], real_rules={k: v['status'] for (k, v) in real['rules'].items()})), flush=True)
    return 0 if result['all_gates_pass'] else 1
if __name__ == '__main__':
    raise SystemExit(main())
