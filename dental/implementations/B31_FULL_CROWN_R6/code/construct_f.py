from local import *
from check_intersections import intersect
from clip import continuous_cut
from skimage.measure import marching_cubes
import subprocess

def continuation(m):
    (pv, pf) = compact(m['vertices'], m['faces'][m['roles'] == 1][:, ::-1])
    ids = loops(trimesh.Trimesh(pv, pf, process=False))[0]
    rim = pv[ids]
    apex = pv[np.setdiff1d(np.arange(len(pv)), ids)[0]]
    N = len(pv)
    v = np.r_[pv, apex + 1.5 * (rim - apex)]
    f = pf.tolist()
    edges = {tuple(e) for face in pf for e in zip(face, np.roll(face, -1))}
    for (j, a) in enumerate(ids):
        k = (j + 1) % len(ids)
        b = ids[k]
        fs = [[int(a), N + j, N + k], [int(a), N + k, int(b)]]
        if (int(a), int(b)) not in edges:
            fs = [q[::-1] for q in fs]
        f.extend(fs)
    return (v, np.array(f, int), pv, pf, rim)

def raw_offset(m, key):
    (v, f, pv, pf, rim) = continuation(m)
    cache = D6 / (key + '_C_inner_raw.npz')
    if cache.exists():
        x = npz(cache)
        return (x['vertices'], x['faces'], (v, f, pv, pf, rim))
    si = intersect(v, f, key + '_continued_F')
    save(R6 / 'raw' / ('F_CONTINUED_' + key + '.json'), si)
    if si['status'] != 'PASS':
        raise ValueError('continued source self intersections ' + str(si.get('count')))
    step = 0.05
    lo = pv.min(0) - 0.12
    hi = pv.max(0) + 0.12
    axes = [np.arange(lo[i], hi[i] + step, step) for i in range(3)]
    shape = tuple(map(len, axes))
    size = int(np.prod(shape))
    if size > 12000000:
        raise ValueError('grid budget')
    field = np.empty(size, np.float32)
    for k in range(0, size, 65536):
        j = np.arange(k, min(k + 65536, size))
        ijk = np.array(np.unravel_index(j, shape)).T
        points = np.column_stack([axes[l][ijk[:, l]] for l in range(3)])
        (sd, _, q, n) = igl.signed_distance(points, v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
        field[k:k + len(j)] = sd - boundary_gap(q, rim)
    (iv, inf, _, _) = marching_cubes(field.reshape(shape), level=0.0, spacing=(step,) * 3, allow_degenerate=False)
    iv = iv.astype(float) + lo
    (_, dc, _) = fast_nearest(pv[pf], iv)
    (_, de, _) = fast_nearest(v[f[len(pf):]], iv)
    (iv, inf) = continuous_cut(iv, inf, de - dc - 1e-07)
    np.savez_compressed(cache, vertices=iv, faces=inf, prep_vertices=pv, prep_faces=pf)
    return (iv, inf, (v, f, pv, pf, rim))

def reconstruct(m, key):
    (iv, inf, ct) = raw_offset(m, key)
    (v, f, pv, pf, rim) = ct
    cap = trimesh.Trimesh(iv, inf, process=False)
    cc = trimesh.graph.connected_components(cap.face_adjacency, nodes=np.arange(len(inf)))
    areas = np.array([cap.area_faces[c].sum() for c in cc])
    best = int(areas.argmax())
    (iv, inf) = compact(iv, inf[cc[best]])
    info = dict(initial_components=len(cc), discarded_offset_area_mm2=float(areas.sum() - areas[best]), retained_offset_area_mm2=float(areas[best]), source_triangles_retained=len(pf))
    ll = loops(trimesh.Trimesh(iv, inf, process=False))
    perimeters = [np.linalg.norm(iv[l] - np.roll(iv[l], 1, axis=0), axis=1).sum() for l in ll]
    main = int(np.argmax(perimeters))
    info['initial_loops'] = len(ll)
    info['filled_loops'] = len(ll) - 1
    verts = iv.tolist()
    faces = inf.tolist()
    edges = {tuple(e) for face in inf for e in zip(face, np.roll(face, -1))}
    for (j, l) in enumerate(ll):
        if j == main:
            continue
        x = iv[l].mean(0)[None]
        for _ in range(5):
            (sd, _, q, n) = igl.signed_distance(x, v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
            x += (boundary_gap(q, rim) - sd)[:, None] * n
        c = len(verts)
        verts.append(x[0].tolist())
        for (a, b) in zip(l, np.roll(l, -1)):
            faces.append([int(b), int(a), c] if (a, b) in edges else [int(a), int(b), c])
    (iv, inf) = (np.array(verts), np.array(faces, int))
    history = []
    for level in range(4):
        (sd, _, q, n) = igl.signed_distance(iv, v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
        iv += (boundary_gap(q, rim) - sd)[:, None] * n
        (q, g, _) = fast_nearest(pv[pf], iv[inf].mean(1))
        err = abs(g - boundary_gap(q, rim))
        bad = err > 0.005
        history.append(dict(level=level, max_centroid_residual_mm=float(err.max()), bad_facets=int(bad.sum())))
        if not bad.any() or level == 3:
            break
        (iv, inf) = subdiv(iv, inf, bad)
    cap = trimesh.Trimesh(iv, inf, process=False)
    trimesh.repair.fix_normals(cap, multibody=True)
    (iv, inf) = (cap.vertices, cap.faces)
    (_, _, _, nn) = igl.signed_distance(iv[inf].mean(1), v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
    if np.einsum('ij,ij->i', np.cross(iv[inf][:, 1] - iv[inf][:, 0], iv[inf][:, 2] - iv[inf][:, 0]), nn).sum() > 0:
        inf = inf[:, ::-1]
    info.update(history=history, final_loops=[len(l) for l in loops(cap)], final_components=len(trimesh.graph.connected_components(cap.face_adjacency, nodes=np.arange(len(inf)))))
    return (iv, inf, info)

def repair_outer(old, key):
    (v, f) = compact(old['vertices'], old['faces'][old['roles'] == 0])
    intersect(v, f, key + '_outer_before')
    p = subprocess.run([str(D6 / 'repair'), str(D6 / (key + '_outer_before.mesh')), str(D6 / (key + '_outer_after.mesh'))], capture_output=True, text=True)
    if p.returncode:
        raise ValueError('repair aborted ' + p.stderr[:250])
    with (D6 / (key + '_outer_after.mesh')).open() as h:
        (n, nf) = map(int, h.readline().split())
        v = np.loadtxt(h, max_rows=n)
        f = np.loadtxt(h, dtype=int)
    return (v, f, dict(returncode=p.returncode, backend_output=p.stdout, stderr=p.stderr))

def run(allcases=False):
    cases = inputs()
    cases.sort(key=lambda q: (q[0]['key'] != '079905ebf9504544_molar', q[0]['key']))
    rows = []
    st = time.perf_counter()
    tag = 'F_ALL' if allcases else 'F_FIRST'
    if not allcases:
        cases = cases[:1]
    for (a, b) in cases:
        t = time.perf_counter()
        r = dict(key=a['key'], family=a['family'], status='REJECTED')
        m = npz(a['mesh_path'])
        try:
            (iv, inf, info) = reconstruct(m, a['key'])
            old = npz(D / 'B/distance_local_thickening' / (a['key'] + '.npz'))
            (ov, of, ri) = repair_outer(old, a['key'])
            (v, f, roles) = join(ov, of, iv, inf)
            dest = D6 / (a['key'] + '_F.npz')
            np.savez_compressed(dest, vertices=v, faces=f, roles=roles, prep_triangles=m['vertices'][m['faces'][m['roles'] == 1]])
            r.update(status='GENERATED', mesh_path=dest, mesh_sha256=sha(dest), offset=info, outer_repair=ri)
        except Exception as e:
            r['reason'] = repr(e)
        r['seconds'] = time.perf_counter() - t
        rows.append(r)
        save(R6 / f'raw/{tag}_GENERATION.json', rows)
        state6(tag + '_GENERATING', r, 'Freeze all generated candidates and rejections then independent geometry scoring')
        print(clean(r), flush=True)
    freeze(R6 / f'FROZEN_PREDICTIONS_{tag}.json', dict(claim_type='capability', rows=rows, prereg_sha256=sha(R6 / 'PREREG_F.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run('--all' in sys.argv)
