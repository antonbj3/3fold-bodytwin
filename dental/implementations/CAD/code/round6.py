import thread_guard
from cadlib import *
from readonly import parents
from pipeline import stitch
import trimesh, three_mf, resource
from shapely.geometry import Polygon
P = parents()

def run():
    st = time.perf_counter()
    pr = read(ROOT / 'PREREG_R6.json')
    old = {r['uid']: r for r in read(ROOT / 'raw/R5_PREDICTIONS.json')}
    inputs = {r['key']: r for r in read(ROOT / 'INPUT_LOCK_R3.json')['rows']}
    rows = []
    from obstacle import constraints
    for uid in pr['cohort']:
        r = old[uid]
        a = np.load(r['mesh_path'])
        inp = inputs[r['key']]
        dat = np.load(inp['public_geometry'])
        t = P['tasks'].load_task(inp['task_path'])
        rec = dict(uid=uid, key=r['key'], family=r['family'], status='FAILED')
        try:
            tri = dat['donor'].copy()
            tri[:, :, 0] *= -1
            tri = tri[:, ::-1]
            center = np.mean(tri.reshape(-1, 3), axis=0)
            tri[:, :, :2] -= center[:2]
            dm = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(tri.size // 3).reshape(-1, 3), process=True)
            fit = r['donor_fit']
            dm.vertices[:, :2] *= fit['scale_xy']
            dm.vertices = dm.vertices @ np.array(fit['R']).T + fit['translation_mm']
            (ext, drop, nc) = P['r4'].component(dm.triangles)
            z0 = float(a['iv'][:, 2].min())
            ext = ext.slice_plane([0, 0, z0], [0, 0, 1], cap=False)
            ext.merge_vertices(digits_vertex=10)
            loops = P['r4'].loops(ext)
            if drop > 0.05 or len(loops) != 1:
                raise ValueError('UNKNOWN_DONOR_TOPOLOGY_OR_COMPONENT_LOSS')
            ids = loops[0]
            boundary = ext.vertices[ids]
            phi = np.mod(np.arctan2(boundary[:, 1], boundary[:, 0]), 2 * np.pi)
            walk = np.mod(np.roll(phi, -1) - phi, 2 * np.pi).sum()
            if abs(walk - 2 * np.pi) > 1e-06 and abs(walk - 2 * np.pi * (len(ids) - 1)) > 1e-06:
                raise ValueError('UNKNOWN_NONSTAR_BOUNDARY')
            if not Polygon(boundary[:, :2]).is_valid:
                raise ValueError('UNKNOWN_CROSSING_BOUNDARY')
            v = ext.vertices.copy()
            f = ext.faces.copy()
            N = len(v)
            bottom = boundary.copy()
            bottom[:, 2] = z0
            if np.min(boundary[:, 2] - z0) < 1e-05:
                raise ValueError('UNKNOWN_MIXED_CUT_AND_NATIVE_BOUNDARY')
            vv = np.r_[v, bottom]
            ff = f.tolist()
            for (j, i) in enumerate(ids):
                k = (j + 1) % len(ids)
                b = ids[k]
                ff.extend([[int(i), int(b), N + k], [int(i), N + k, N + j]])
            ff = np.array(ff)
            basis = np.clip((vv[:, 2] - z0) / 2, 0, 1)
            (A, b) = constraints(t['antagonist'], vv[:, :2], ff, 0.0)
            slack = b - 0.02 - A @ vv[:, 2]
            w = A @ basis
            zero = w < 1e-12
            if np.any(slack[zero] < -1e-08):
                raise ValueError('INFEASIBLE_FIXED_CERVICAL_OCCLUSION')
            alpha = min(0.0, float(np.min(slack[~zero] / w[~zero])) - 1e-08) if np.any(~zero) else 0.0
            vv[:, 2] += alpha * basis
            (m, roles) = stitch(vv, ff, a['iv'], a['inf'], np.arange(24), np.arange(N, len(vv)))
            dest = DATA / 'R6' / uid
            dest.mkdir(parents=True, exist_ok=True)
            path = dest / 'crown.npz'
            np.savez_compressed(path, vertices=m.vertices, faces=m.faces, roles=roles, ov=vv, of=ff, iv=a['iv'], inf=a['inf'])
            three_mf.write(dest / 'crown.3mf', m.vertices, m.faces, roles, dict(material=r['material'], physical_release='UNKNOWN'))
            rec.update(status='GENERATED', mesh_path=path, mesh_sha256=sha(path), height_coefficient_mm=alpha, source_component_drop=drop, donor_source_components=nc, occlusion_linear_max_violation_mm=float(np.max(A @ vv[:, 2] - (b - 0.02))), outer_faces=len(ff), topology=dict(watertight=bool(m.is_watertight), winding=bool(m.is_winding_consistent), volume_mm3=float(m.volume)))
        except Exception as ex:
            rec['reason'] = repr(ex)
        rows.append(rec)
        print(uid, rec['status'], rec.get('reason', ''), flush=True)
    dump(ROOT / 'raw/R6_PREDICTIONS.json', rows)
    freeze(ROOT / 'FROZEN_PREDICTIONS_R6.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R6.json'), rows_sha256=sha(ROOT / 'raw/R6_PREDICTIONS.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, reference_scored=False))
    for rec in rows:
        if rec['status'] != 'GENERATED':
            continue
        a = np.load(rec['mesh_path'])
        ext = a['ov'][a['of']]
        inp = inputs[rec['key']]
        t = P['tasks'].load_task(inp['task_path'])
        ref = np.load(inp['reference_path'])
        target = ref['target_tri']
        target = target[target.mean(1)[:, 2] >= a['iv'][:, 2].min()]
        rec['shape'] = P['r4'].metrics(ext, target, n=2048)
        rec['shape_gate'] = 'PASS' if rec['shape']['p95_mm'] <= pr['metrics']['shape_p95_mm_by_family'][rec['family']] else 'FAIL'
        rec['wall'] = P['cap'].wall_check(a['ov'], a['of'], a['iv'], a['inf'], 0.5)
        native = P['v2g'].height(target, t['xy'])
        pred = P['v2g'].height(ext, t['xy'])
        cm = P['v6'].compare(t['xy'], t['faces'], t['ceiling'] - pred, t['ceiling'] - native)
        rec['V6_contact'] = cm
        rec['V6_gates'] = P['v6'].gates(cm, pr['metrics']['V6'])
        vv = a['ov']
        ff = a['of']
        em = trimesh.Trimesh(vv, ff, process=False)
        bd = P['r4'].loops(em)[0]
        cc = vv[bd].mean(0)
        cv = np.r_[vv, [cc]]
        cf = np.r_[ff, [[int(i), len(vv), int(j)] for (i, j) in zip(bd, np.roll(bd, -1))]]
        solid = trimesh.Trimesh(cv, cf, process=False)
        trimesh.repair.fix_normals(solid)
        pts = a['iv'].copy()
        pts[pts[:, 2] == pts[:, 2].min(), 2] += 1e-05
        rec['intaglio_vertex_containment'] = dict(inside_count=int(solid.contains(pts).sum()), points=len(pts), scope='Point sample, not continuous containment proof')
        rec['complete_status'] = 'UNKNOWN_PHYSICAL_AND_CONTINUOUS_CONTAINMENT'
        print('R6_SCORE', rec['uid'], rec['shape_gate'], rec['shape']['p95_mm'], rec['wall']['status'], flush=True)
    out = dict(claim_type='capability', requested=3, generated=sum((r['status'] == 'GENERATED' for r in rows)), shape_pass=sum((r.get('shape_gate') == 'PASS' for r in rows)), rows=rows, external_referent=dict(pr['external_referent'], refutes_us=any((r.get('shape_gate') == 'FAIL' for r in rows))), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(ROOT / 'RESULTS_R6.json', out)
    (ROOT / 'HANDOFF_R6.md').write_text(json.dumps(clean(out), indent=2) + '\nNext: source collar/pose and all-surface containment; no threshold change or clinical claim.\n')
    print('R6 complete', out['shape_pass'])
if __name__ == '__main__':
    run()
