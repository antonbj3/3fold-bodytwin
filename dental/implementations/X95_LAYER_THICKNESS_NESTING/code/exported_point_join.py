from pose_join import *
from measure_geometry import closest
import trimesh

def run():
    st = time.perf_counter()
    pr = json.loads((ROOT / 'PREREG_R5.json').read_text())
    rule = json.loads((ROOT / 'RULE_CONTRACT.json').read_text())
    assert source_gate(rule)
    rows = []
    for e in json.loads((ROOT / 'raw/R4_RESULTS.json').read_text())['rows']:
        a = dict(np.load(ROOT / e['array_path'], allow_pickle=False))
        wall = dict(np.load(ROOT / 'raw' / f"{e['family']}_wall.npz", allow_pickle=False))
        m = trimesh.Trimesh(a['vertices'], a['faces'], process=False)
        p = wall['points']
        (near, t, ids) = closest(m, p)
        pose = e['source_best_pose']
        n = axis(pose['tilt_deg'], pose['azimuth_deg'])
        z = pose['top_depth_mm'] + np.max(wall['vertices'] @ n) - p @ n
        (high, low) = (0.8, 0.4) if e['family'] == 'anterior' else (1.0, 0.5)
        req = np.where(z <= pose['height_mm'] / 2, high, low)
        slack = t - req
        k = int(slack.argmin())
        assert slack.min() >= 0
        js = np.unique(np.r_[k, np.linspace(0, len(p) - 1, 32, dtype=int)])
        ind = np.array([point_triangle_distance(p[j], m.triangles) for j in js])
        error = float(abs(ind - t[js]).max())
        assert error <= 1e-10
        fault = m.copy()
        fault.vertices = fault.vertices + (p[k] - fault.vertices[0])
        faultdist = closest(fault, p[k:k + 1])[1][0]
        faultslack = float(faultdist - req[k])
        assert faultslack < 0
        path = ROOT / 'exports' / f"{e['family']}_final_point_slack.npz"
        np.savez_compressed(path, points_mm=p, thickness_mm=t, required_mm=req, slack_mm=slack, disc_depth_mm=z, nearest_intaglio_mm=near, area_mm2=wall['area'], original_intaglio_thickness_mm=wall['distance'])
        with (ROOT / 'exports' / f"{e['family']}_final_point_slack.vtk").open('w') as f:
            f.write('# vtk DataFile Version 3.0\nX95 virtual geometry PER_POINT mm NOT FABRICATION QUALIFIED\nASCII\nDATASET POLYDATA\nPOINTS %d double\n' % len(p))
            np.savetxt(f, p, fmt='%.17g')
            f.write('VERTICES %d %d\n' % (len(p), 2 * len(p)))
            np.savetxt(f, np.c_[np.ones(len(p), int), np.arange(len(p))], fmt='%d')
            f.write('POINT_DATA %d\n' % len(p))
            for (name, data) in [('thickness_mm', t), ('required_mm', req), ('slack_mm', slack), ('disc_depth_mm', z)]:
                f.write('SCALARS %s double 1\nLOOKUP_TABLE default\n' % name)
                np.savetxt(f, data, fmt='%.17g')
        row = {'family': e['family'], 'resolution': 'PER_POINT', 'points': len(p), 'point_min_thickness_mm': float(t.min()), 'point_min_slack_mm': float(slack.min()), 'negative_point_slacks': int(np.sum(slack < 0)), 'source_pose': pose, 'independent_point_triangle_max_error_mm': error, 'control_points': len(js), 'fault_core_translation_min_slack_mm': faultslack, 'fault_rejected': True, 'array_path': str(path.relative_to(ROOT)), 'array_sha256': sha(path), 'rigorous_float_source_enclosure': 'MISSING', 'qualification': 'CONDITIONAL_VIRTUAL_ONLY_NOT_FABRICATION_QUALIFIED'}
        rows.append(row)
        print(e['family'], 'localminimumslack', row['point_min_slack_mm'], 'fault', faultslack, flush=True)
    out = {'round': 'R5', 'claim_type': 'capability', 'primary_gate': 'PASS_SCOPED_EXPORTED_GEOMETRY_POINT_QUERY', 'rows': rows, 'seconds': time.perf_counter() - st, 'resolution': 'PER_POINT', 'timescale': 'HANDOVER', 'physical_measurement': 'NOT_PERFORMED'}
    dump(ROOT / 'raw/R5_RESULTS.json', out)
    return out
if __name__ == '__main__':
    run()
