from bounded_outer import repair_outer
from score6 import *

def generate():
    rows = []
    st = time.perf_counter()
    for (a, b) in inputs():
        t = time.perf_counter()
        r = dict(key=a['key'], family=a['family'], status='REJECTED')
        old = npz(D / 'B/distance_local_thickening' / (a['key'] + '.npz'))
        try:
            kr = next((r for r in read(R6 / 'FROZEN_PREDICTIONS_K_ALL.json')['rows'] if r['key'] == a['key'] and r['status'] == 'GENERATED'), None)
            if kr:
                frozen = npz(kr['mesh_path'])
                (v, f) = compact(frozen['vertices'], frozen['faces'][frozen['roles'] == 0])
                info = dict(reused_frozen_K_exterior=True, source=kr['mesh_path'], source_sha256=kr['mesh_sha256'], note='Same D operation already executed, no second fit or native-target choice')
            else:
                (v, f, info) = repair_outer(old, a['key'] + '_D_BATCH')
            dest = D6 / (a['key'] + '_D_BATCH.npz')
            np.savez_compressed(dest, vertices=v, faces=f)
            r.update(status='GENERATED', mesh_path=dest, mesh_sha256=sha(dest), repair=info)
        except Exception as e:
            r['reason'] = repr(e)
        r['seconds'] = time.perf_counter() - t
        rows.append(r)
        save(R6 / 'raw/D_BATCH_GENERATION.json', rows)
        print('D_BUILD', r['key'], r['status'], flush=True)
    freeze(R6 / 'FROZEN_PREDICTIONS_D_BATCH.json', dict(rows=rows, prereg_sha256=sha(R6 / 'PREREG_D_BATCH.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))

def evaluate():
    st = time.perf_counter()
    idx = {b['key']: b for (a, b) in inputs()}
    rows = []
    geometry_r4.closest = fast_nearest
    scorer.closest = fast_nearest
    for rec in read(R6 / 'FROZEN_PREDICTIONS_D_BATCH.json')['rows']:
        r = dict(rec)
        t = time.perf_counter()
        if rec['status'] == 'GENERATED':
            try:
                m = npz(rec['mesh_path'])
                assert sha(rec['mesh_path']) == rec['mesh_sha256']
                old = npz(D / 'B/distance_local_thickening' / (rec['key'] + '.npz'))
                p = npz(idx[rec['key']]['public_path'])
                target = npz(idx[rec['key']]['private_path'])['target']
                (v, f) = (m['vertices'], m['faces'])
                mm = trimesh.Trimesh(v, f, process=False)
                ll = loops(mm)
                shape = scorer.metrics(v[f], target)
                wall = scorer.wall(v[f], old['vertices'][old['faces'][old['roles'] == 1]])
                dd = cKDTree(v).query(p['margin_curve'], workers=1)[0]
                si = intersect(v, f, rec['key'] + '_D_BATCH_check')
                pairs = np.array(si.pop('pairs', []), int)
                dest = D6 / (rec['key'] + '_D_BATCH_si.npz')
                np.savez_compressed(dest, pairs=pairs)
                si.update(pairs_path=dest, pairs_sha256=sha(dest))
                gate = dict(shape=shape['p95_mm'] <= read(R6 / 'PREREG_A.json')['metrics']['shape_p95_mm'][rec['family']], margin_identity=bool((dd == 0).all()), wall=wall['continuous_lower_mm'] >= 0.5, simple_disk=si['status'] == 'PASS' and len(ll) == 1 and bool(mm.is_winding_consistent) and (mm.euler_number == 1) and (len(trimesh.graph.connected_components(mm.face_adjacency, nodes=np.arange(len(f)))) == 1))
                r.update(shape=shape, wall=wall, margin_identity_error_mm=float(dd.max()), self_intersection=si, boundary_loop_lengths=[len(x) for x in ll], gates=gate, scope='Exterior only against unchanged R5 B intaglio. No complete-crown claim.')
            except Exception as e:
                r['score_error'] = repr(e)
        r['score_seconds'] = time.perf_counter() - t
        rows.append(r)
        save(R6 / 'raw/D_BATCH_SCORE.json', rows)
        print('D_SCORE', r['key'], r.get('gates'), flush=True)
    save(R6 / 'RESULTS_D_BATCH.json', dict(claim_type='capability', rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    generate()
    evaluate()
