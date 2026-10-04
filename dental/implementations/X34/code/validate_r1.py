from dental_release.paths import expand as _release_expand
import json
import time
import sys
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
sys.path.insert(0, str(Path(__file__).resolve().parent))
from designgate.gate import check
from designgate.geometry import continuous, load, faces
from designgate.vendor.geometry import export_stl
from initialize import H, dump, sha
D = Path(_release_expand('@DENTAL_WORK_ROOT@/X34_stl_design_gate'))

def paths(name):
    return {k: D / name / (k + '.stl') for k in ('prep', 'crown', 'antagonist')}

def main():
    frozen = json.loads((H / 'FROZEN_PREDICTIONS.json').read_text())
    assert sha(H / 'PREREG_R1.json') == frozen['prereg_sha256']
    for x in frozen['files']:
        assert sha(x['path']) == x['sha256']
    for (p, h) in frozen['code'].items():
        assert sha(H / p) == h
    start = time.perf_counter()
    rows = {}
    tests = []

    def record(name, condition, detail):
        tests.append(dict(name=name, passed=bool(condition), detail=detail))
    for (case, expect) in frozen['expected'].items():
        r = check(paths(case), contract=D / case / 'contract.json')
        rows[case] = r
        for (k, status) in expect.items():
            record(case + '_' + k, r['rules'][k]['status'] == status, dict(expected=status, observed=r['rules'][k]['status']))
    v = rows['valid']
    record('closed_form_wall', max((abs(x - 0.8) for x in v['rules']['material_wall']['minimum_interval_mm'])) <= 1e-05, v['rules']['material_wall']['minimum_interval_mm'])
    record('closed_form_film', abs(v['rules']['film_min']['value_mm'] - 0.05) <= 1e-05, v['rules']['film_min']['value_mm'])
    lower = np.array([[[-1, 0, 0], [1, 0, 0], [0, 2, 0]]], float)
    upper = np.array([[[-1, 1, 0.1], [1, 1, 0.1], [0, -1, -0.3]]], float)
    a = continuous(upper, lower, [0, 0, 1], True)
    b = continuous(upper, lower, [0, 0, 1], False)
    record('continuous_crossing_analytic', abs(a['minimum_gap_mm'] - -0.1) <= 1e-05, a)
    record('continuous_full_pair_control', abs(a['minimum_gap_mm'] - b['minimum_gap_mm']) <= 1e-07, b)
    for (reason, kw) in [('missing_units', dict()), ('missing_frame', dict(units='mm'))]:
        r = check(paths('valid'), **kw)
        record(reason, all((x['status'] == 'UNKNOWN' for x in r['rules'].values())), r['verdict'])
    meta = json.loads((D / 'valid/contract.json').read_text())
    meta.pop('ifu_profile')
    meta['indication'] = 'posterior'
    p = D / 'attack_3Y.json'
    dump(p, meta)
    r = check(paths('valid'), contract=p, material='3Y')
    record('class_3Y_no_IFU', r['rules']['material_wall']['status'] == 'UNKNOWN', r['rules']['material_wall'])
    meta = json.loads((D / 'valid/contract.json').read_text())
    meta['input_sha256']['crown'] = '0' * 64
    p = D / 'attack_hash.json'
    dump(p, meta)
    try:
        check(paths('valid'), contract=p)
        rejected = False
    except ValueError:
        rejected = True
    record('changed_input_hash', rejected, 'bad hash must reject')
    rr = Rotation.from_rotvec([0.35, -0.6, 0.2]).as_matrix()
    tt = np.array([10, 20, -7])
    rotpaths = {}
    upaths = {}
    for (k, p) in paths('valid').items():
        m = load(p, 1)
        dest = D / ('rot_' + k + '.stl')
        export_stl(dest, m['vertices'] @ rr.T + tt, m['faces'])
        rotpaths[k] = dest
        dest = D / ('um_' + k + '.stl')
        export_stl(dest, m['vertices'] * 1000, m['faces'])
        upaths[k] = dest
    meta = json.loads((D / 'valid/contract.json').read_text())
    meta['contact_axis'] = rr[:, 2].tolist()
    meta['input_sha256'] = {k: sha(p) for (k, p) in rotpaths.items()}
    p = D / 'rotation.json'
    dump(p, meta)
    rot = check(rotpaths, contract=p)
    rows['rotation'] = rot
    record('rotation_contact_invariant', abs(rot['rules']['occlusal_contact']['measurement']['minimum_gap_mm'] - 0.025) <= 0.0001, rot['rules']['occlusal_contact'])
    base_wit = np.array(v['rules']['occlusal_contact']['measurement']['witness']['lower_point_mm'])
    rp = np.array(rot['rules']['occlusal_contact']['measurement']['witness']['lower_point_mm'])
    record('rotation_witness_plane', abs(((rp - tt) @ rr)[2] - 0.85) <= 0.0001, rp.tolist())
    meta = json.loads((D / 'valid/contract.json').read_text())
    meta['units'] = 'um'
    meta['input_sha256'] = {k: sha(p) for (k, p) in upaths.items()}
    p = D / 'micrometers.json'
    dump(p, meta)
    um = check(upaths, contract=p)
    record('explicit_um_conversion', um['rules']['material_wall']['status'] == 'PASS' and abs(um['rules']['film_min']['value_mm'] - 0.05) <= 1e-05, um['rules']['material_wall'])
    meta['units'] = 'mm'
    p = D / 'wrong_scale.json'
    dump(p, meta)
    wrong = check(upaths, contract=p)
    record('injected_wrong_scale', wrong['rules']['film_max']['status'] == 'FAIL', wrong['rules']['film_max'])
    real = check(paths('real_X18_110_34'), contract=D / 'real_X18_110_34/contract.json')
    rows['real'] = real
    g = real['rules']['occlusal_contact']['measurement']
    record('real_X18_frozen_gap', g['minimum_gap_mm'] is not None and abs(g['minimum_gap_mm'] - 0.025) <= 1e-05, g)
    import trimesh
    ma = load(paths('real_X18_110_34')['antagonist'], 1)['mesh']
    w = np.array(g['witness']['upper_point_mm'])
    origin = w + [0, 0, 100]
    (loc, ray, tri) = ma.ray.intersects_location([origin], [[0, 0, -1]], multiple_hits=True)
    ray_z = float(loc[:, 2].min()) if len(loc) else None
    record('independent_trimesh_ray', ray_z is not None and abs(ray_z - w[2]) <= 1e-05, dict(published_code='https://trimesh.org/trimesh.ray.ray_triangle.html', witness_z_mm=float(w[2]), independent_lowest_z_mm=ray_z))
    raw = dict(round='R1', claim_type='capability', tests=tests, reports=rows, passed=sum((x['passed'] for x in tests)), total=len(tests), cost_wall_s=time.perf_counter() - start, all_gates_pass=all((x['passed'] for x in tests)), full_restoration_validated=False)
    dump(H / 'rounds/R1.json', raw)
    print(json.dumps(dict(passed=raw['passed'], total=raw['total'], failed=[x['name'] for x in tests if not x['passed']], real_rules={k: v['status'] for (k, v) in real['rules'].items()})))
    return 0 if raw['all_gates_pass'] else 1
if __name__ == '__main__':
    raise SystemExit(main())
