from construct_f import *

def inner(m, key):
    cache = D6 / (key + '_I_inner_raw.npz')
    if not cache.exists():
        raise ValueError('No I offset with verified source support')
    attempts = read(R6 / 'raw' / ('I_CLOSURE_' + key + '.json'))
    best = next((x for x in attempts if x.get('watertight') and x.get('intersections') == 0 and (x.get('original_source_faces_retained') == x.get('original_source_faces'))))
    name = key + '_I_closure_' + str(best['factor']) + '_after.mesh'
    with (D6 / name).open() as h:
        (n, nf) = map(int, h.readline().split())
        v = np.loadtxt(h, max_rows=n)
        f = np.loadtxt(h, dtype=int)
    dat = npz(cache)
    (iv, inf, pv, pf) = (dat['vertices'], dat['faces'], dat['prep_vertices'], dat['prep_faces'])
    rim = pv[loops(trimesh.Trimesh(pv, pf, process=False))[0]]
    cm = trimesh.Trimesh(iv, inf, process=False)
    cc = trimesh.graph.connected_components(cm.face_adjacency, nodes=np.arange(len(inf)))
    area = [cm.area_faces[c].sum() for c in cc]
    ix = int(np.argmax(area))
    (iv, inf) = compact(iv, inf[cc[ix]])
    cm = trimesh.Trimesh(iv, inf, process=False)
    info = dict(initial_components=len(cc), discarded_offset_area_mm2=float(sum(area) - area[ix]), closure_path=D6 / name, closure_sha256=sha(D6 / name))
    ll = loops(cm)
    if len(ll) != 1:
        raise ValueError('Principal I component has ' + str(len(ll)) + ' boundaries')
    history = []
    for level in range(4):
        (sd, _, q, n) = igl.signed_distance(iv, v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
        iv += (boundary_gap(q, rim) - sd)[:, None] * n
        (q, g, _) = fast_nearest(pv[pf], iv[inf].mean(1))
        err = abs(g - boundary_gap(q, rim))
        bad = err > 0.005
        history.append(dict(level=level, residual_max_mm=float(err.max()), bad_facets=int(bad.sum())))
        if not bad.any() or level == 3:
            break
        (nv, nf) = subdiv(iv, inf, bad)
        if len(nv) > 300000:
            raise ValueError('K vertex budget')
        (iv, inf) = (nv, nf)
    cap = trimesh.Trimesh(iv, inf, process=False)
    trimesh.repair.fix_normals(cap, multibody=True)
    (iv, inf) = (cap.vertices, cap.faces)
    (_, _, _, nn) = igl.signed_distance(iv[inf].mean(1), v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
    if np.einsum('ij,ij->i', np.cross(iv[inf][:, 1] - iv[inf][:, 0], iv[inf][:, 2] - iv[inf][:, 0]), nn).sum() > 0:
        inf = inf[:, ::-1]
    return (iv, inf, dict(info, history=history))

def run():
    st = time.perf_counter()
    rows = []
    for (a, b) in inputs():
        t = time.perf_counter()
        r = dict(key=a['key'], family=a['family'], status='REJECTED')
        m = npz(a['mesh_path'])
        try:
            (iv, inf, info) = inner(m, a['key'])
            old = npz(D / 'B/distance_local_thickening' / (a['key'] + '.npz'))
            (ov, of, ri) = repair_outer(old, a['key'] + '_K')
            (v, f, roles) = join(ov, of, iv, inf)
            dest = D6 / (a['key'] + '_K.npz')
            np.savez_compressed(dest, vertices=v, faces=f, roles=roles, prep_triangles=m['vertices'][m['faces'][m['roles'] == 1]])
            r.update(status='GENERATED', mesh_path=dest, mesh_sha256=sha(dest), offset=info, outer_repair=ri)
        except Exception as e:
            r['reason'] = repr(e)
        r['seconds'] = time.perf_counter() - t
        rows.append(r)
        save(R6 / 'raw/K_ALL_GENERATION.json', rows)
        print(r['key'], r['status'], r.get('reason'), flush=True)
    freeze(R6 / 'FROZEN_PREDICTIONS_K_ALL.json', dict(claim_type='capability', rows=rows, prereg_sha256=sha(R6 / 'PREREG_K.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run()
