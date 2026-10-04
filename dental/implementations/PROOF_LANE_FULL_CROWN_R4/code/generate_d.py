from generate_c import exact_shell
from geometry import *

def adapt(ext, p):
    base = ext.vertices.copy()
    out = ext.copy()
    ids = loops(ext)[0]
    free = np.ones(len(base), bool)
    free[ids] = False
    info = {}
    for side in ['mesial', 'distal']:
        tri = p[side]
        if not len(tri):
            info[side] = dict(status='ABSTAIN_NO_NEIGHBOUR')
            continue
        (near, dist, fi) = closest(tri, out.vertices)
        poss = free & (dist < 0.75)
        if not poss.any():
            info[side] = dict(status='ABSTAIN_NO_PATCH_WITHIN_0.75MM')
            continue
        cand = np.flatnonzero(poss)
        idx = cand[np.argmin(dist[cand])]
        center = out.vertices[idx].copy()
        rad = np.linalg.norm(out.vertices - center, axis=1)
        w = np.where(rad < 0.8, 0.5 * (1 + np.cos(np.pi * np.minimum(rad / 0.8, 1))), 0)
        w[~free] = 0
        initial = float(dist[idx])
        for _ in range(4):
            (q, d, _) = closest(tri, out.vertices)
            direction = (out.vertices - q) / np.maximum(d[:, None], 1e-12)
            target = q + 0.05 * direction
            want = out.vertices + w[:, None] * (target - out.vertices)
            u = want - base
            nn = np.linalg.norm(u, axis=1)
            u *= np.minimum(1, 0.5 / np.maximum(nn, 1e-20))[:, None]
            u[ids] = 0
            out.vertices = base + u
        info[side] = dict(status='UNSIGNED_NEAR_PATCH_ADAPTED', initial_nearest_mm=initial, patch_center_mm=center, moved_vertices=int(np.sum(w > 0)))
    info['antagonist'] = dict(status='ABSTAIN_UNVERIFIED_BITE_POSE', transformed=False)
    info['additional_displacement_bound_mm'] = float(np.linalg.norm(out.vertices - base, axis=1).max())
    info['boundary_coordinate_error_mm'] = float(np.linalg.norm(out.vertices[ids] - base[ids], axis=1).max())
    return (out, info)

def run(tag='D2'):
    st = time.perf_counter()
    inputs = {r['key']: r for r in read(ROOT / 'FROZEN_INPUTS_B.json')['records']}
    prior = read(ROOT / 'FROZEN_PREDICTIONS_C.json')
    rows = []
    previous = {r['key']: r for r in read(ROOT / 'FROZEN_PREDICTIONS_D.json')['rows']} if tag == 'D2' else {}
    for rec in prior['rows']:
        if rec['method'] != 'fixed_prep_projection':
            continue
        if rec['key'] in previous and previous[rec['key']]['status'] == 'GENERATED':
            rows.append(previous[rec['key']])
            continue
        r = {k: rec[k] for k in ['key', 'case_key', 'family']}
        r.update(method='measured_neighbours', status='FAILED')
        t0 = time.perf_counter()
        try:
            p = npz(inputs[r['key']]['public_path'])
            raw = npz(rec['mesh_path'])
            ext = trimesh.Trimesh(raw['vertices'], raw['faces'][raw['roles'] == 0], process=False)
            ext.remove_unreferenced_vertices()
            (e, info) = adapt(ext, p)
            (m, roles) = exact_shell(e, p)
            path = DATA / tag / r['method'] / (r['key'] + '.npz')
            path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(path, vertices=m.vertices, faces=m.faces, roles=roles)
            r.update(status='GENERATED', mesh_path=path, mesh_sha256=sha(path), context_adaptation=info, build=rec['build'])
        except Exception as e:
            r['reason'] = repr(e)
        r['seconds'] = time.perf_counter() - t0
        rows.append(r)
        dump(ROOT / 'raw' / f'{tag}_GENERATION.json', rows)
        print(tag, r['key'], r['status'], r.get('reason', ''), flush=True)
    freeze(ROOT / f'FROZEN_PREDICTIONS_{tag}.json', dict(claim_type='capability', prereg_sha256=sha(ROOT / 'PREREG_D.json'), parent_sha256=sha(ROOT / 'FROZEN_PREDICTIONS_C.json'), code_sha256=sha(Path(__file__)), rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, reused_previously_generated=16 if tag == 'D2' else 0))
if __name__ == '__main__':
    run()
