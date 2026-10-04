import time, copy
from common import *
from collision import *
from run_r1 import lib

def run():
    pre = read(ROOT / 'PREREG_R3.json')
    assert sha(ROOT / 'PREREG_R3.json') == (ROOT / 'PREREG_R3.sha256').read_text().strip()
    assert sha(ROOT / 'FROZEN_PREDICTIONS.json') == pre['frozen_R1_sha256']
    start = time.perf_counter()
    rows = []
    manifest = {r['id']: r for r in read(ROOT / 'INPUT_MANIFEST.json')['rows']}
    scale = 0.8
    for (file, expected) in read(ROOT / 'FROZEN_PREDICTIONS.json')['payload']['point_files'].items():
        assert sha(file) == expected
        source = read(file)
        v = next((v for v in source['points_and_tools'] if v['library'] == 'vhf' and v['axes'] == 5))
        selected = [(i, r) for (i, r) in enumerate(v['records']) if r['region'] == 'intaglio' and r['tools'][0]['status'] != 'FOUND' and any((t['status'] == 'FOUND' for t in r['tools'][1:]))]
        if not selected:
            continue
        src = manifest[source['id']]
        a = np.load(src['file'], allow_pickle=False)
        scene = Scene(a['vertices'], a['faces'])
        base = lib('vhf', 5)[0]
        for (i, r) in selected:
            p = np.array(r['point'])
            n = np.array(r['normal'])
            lo = 3.0
            hi = 22.0
            t = dict(base, neck_reach_mm=hi * scale)
            upper = search(scene, p, n, t, 5)
            if upper['status'] != 'FOUND':
                rows.append(dict(id=source['id'], point_index=i, status='UNKNOWN_NO22MM_POSE', lower_green_mm=lo, upper_green_mm=None))
                continue
            while hi - lo > pre['metrics']['bisection_bracket_green_mm_max']:
                mid = (hi + lo) / 2
                t = dict(base, neck_reach_mm=mid * scale)
                w = search(scene, p, n, t, 5)
                if w['status'] == 'FOUND':
                    hi = mid
                    upper = w
                else:
                    lo = mid
            t = dict(base, neck_reach_mm=hi * scale)
            approach = float(np.linalg.norm(scene.bounds[1] - scene.bounds[0]) + 2 * t['holder_diameter_mm'] + 2.0)
            replay = pose_clearance(scene, np.array(upper['centre']), n, np.array(upper['direction']), t, approach)[0]
            wrong = dict(base, neck_reach_mm=3 * scale)
            corrupt = pose_clearance(scene, np.array(upper['centre']), n, np.array(upper['direction']), wrong, approach)[0]
            assert replay and (not corrupt)
            rows.append(dict(id=source['id'], point_index=i, point=p, normal=n, status='BRACKETED_GRID_NECK_REQUIREMENT', lower_green_mm=lo, upper_green_mm=hi, bracket_green_mm=hi - lo, nominal_neck_green_mm=3.0, upper_pose=upper, replayed_upper_pass=replay, injected3mm_neck_rejected=not corrupt, tool_diameter_green_mm=0.6, zero_stock_offset=True, physical_status='UNKNOWN_NECK_PROFILE_HOLDER_AND_COMPATIBILITY'))
        dump(ROOT / 'raw/R3_PARTIAL.json', rows)
        state('R3_INVERSE_NECK_RUNNING', str(len(rows)) + ' point requirements saved', 'Continue inverse tool specification without changing intaglio')
    controls = []
    valid = [r for r in rows if r['status'] == 'BRACKETED_GRID_NECK_REQUIREMENT']
    for k in sorted(set([0, len(valid) // 2, len(valid) - 1])):
        r = valid[k]
        a = np.load(manifest[r['id']]['file'], allow_pickle=False)
        scene = Scene(a['vertices'], a['faces'])
        p = np.array(r['point'])
        n = np.array(r['normal'])
        d = np.array(r['upper_pose']['direction'])
        c = np.array(r['upper_pose']['centre'])
        t = lib('vhf', 5)[0]
        vals = []
        for ell in [3.0, r['lower_green_mm'], r['upper_green_mm'], 22.0]:
            t1 = dict(t, neck_reach_mm=ell * 0.8)
            (ok, _, _) = pose_clearance(scene, c, n, d, t1, 40.0)
            shank_a = c + ell * 0.8 * d
            shank_b = c + 40 * d
            rad = t1['shank_mm'] / 2
            lo = np.minimum(shank_a, shank_b) - rad - 1
            hi = np.maximum(shank_a, shank_b) + rad + 1
            keep = np.all(scene.hi >= lo, 1) & np.all(scene.lo <= hi, 1)
            if not keep.any():
                enclosure = [0.0, np.inf]
                model = scene.clearance(shank_a, shank_b, rad) + rad
            else:
                tri = scene.tri[keep]
                mesh = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=False)
                pts = shank_a + (shank_b - shank_a) * np.linspace(0, 1, 801)[:, None]
                (_, dist, _) = trimesh.proximity.closest_point_naive(mesh, pts)
                up = dist.min()
                low = max(0.0, up - np.linalg.norm(shank_b - shank_a) / 1600)
                model = segment_triangle(shank_a, shank_b, tri).min()
                enclosure = [low, up]
                assert low - 1e-06 <= model <= up + 1e-06
            vals.append(dict(neck_green_mm=ell, fixed_pose_found=ok, independent_distance_enclosure_mm=enclosure, segment_distance_mm=model))
        seq = [x['fixed_pose_found'] for x in vals]
        assert seq == sorted(seq)
        controls.append(dict(id=r['id'], point_index=r['point_index'], lengths=vals, monotone=True, scope='fixed-direction conventional control on same scene; no method win'))
    upper = [r['upper_green_mm'] for r in valid]
    result = dict(round='R3', claim_type='capability', rows=rows, selected=len(rows), bracketed=len(valid), unknown=len(rows) - len(valid), required_neck_green_mm_range=[min(upper), max(upper)] if upper else None, required_neck_green_mm_median=float(np.median(upper)) if upper else None, max_bracket_mm=max((r['bracket_green_mm'] for r in valid)), all_upper_witnesses_replayed=all((r['replayed_upper_pass'] for r in valid)), all_short_neck_injections_rejected=all((r['injected3mm_neck_rejected'] for r in valid)), controls=controls, seconds=time.perf_counter() - start, construction_gate='PASS' if len(rows) == 111 and len(valid) == 111 else 'PARTIAL', physical_gate='UNKNOWN', external_referent=pre['external_referent'], resolution='PER_POINT')
    dump(ROOT / 'RESULTS_R3.json', result)
    freeze(ROOT / 'FROZEN_PREDICTIONS_R3.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R3.json'), results_sha256=sha(ROOT / 'RESULTS_R3.json'), code_sha256=sha(Path(__file__)), physical_measurement_performed=False))
    return result
if __name__ == '__main__':
    run()
