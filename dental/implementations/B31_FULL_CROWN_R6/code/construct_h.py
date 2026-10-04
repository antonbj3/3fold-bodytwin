from construct_f import *

def retained_faces(v, f, pv, pf):

    def key(t):
        return tuple(sorted((tuple((float(x) for x in p)) for p in t)))
    s = {key(t) for t in v[f]}
    return np.array([key(t) in s for t in pv[pf]])

def closure(m, key):
    (pv, pf) = compact(m['vertices'], m['faces'][m['roles'] == 1][:, ::-1])
    ids = loops(trimesh.Trimesh(pv, pf, process=False))[0]
    rim = pv[ids]
    apex = pv[np.setdiff1d(np.arange(len(pv)), ids)[0]]
    attempts = []
    for factor in [1.5, 2.0, 3.0]:
        N = len(pv)
        n = len(ids)
        outer = apex + factor * (rim - apex)
        low = outer.copy()
        low[:, 2] = min(outer[:, 2].min(), pv[:, 2].min()) - 2.0
        v = np.r_[pv, outer, low, [[0.0, 0.0, low[0, 2]]]]
        f = pf.tolist()
        edges = {tuple(e) for face in pf for e in zip(face, np.roll(face, -1))}
        for (j, a) in enumerate(ids):
            k = (j + 1) % n
            b = ids[k]
            fs = [[int(a), N + j, N + k], [int(a), N + k, int(b)], [N + j, N + n + j, N + n + k], [N + j, N + n + k, N + k], [N + n + j, N + 2 * n, N + n + k]]
            if (int(a), int(b)) not in edges:
                fs = [q[::-1] for q in fs]
            f.extend(fs)
        cm = trimesh.Trimesh(v, np.array(f, int), process=False)
        trimesh.repair.fix_normals(cm, multibody=True)
        if cm.volume < 0:
            cm.invert()
        name = key + '_H_closure_' + str(factor)
        pre = intersect(cm.vertices, cm.faces, name + '_before')
        dest = D6 / (name + '_after.mesh')
        process = subprocess.run([str(D6 / 'repair_local'), str(D6 / (name + '_before.mesh')), str(dest)], capture_output=True, text=True)
        row = dict(factor=factor, before_intersections=pre.get('count'), returncode=process.returncode, stderr=process.stderr[:1000])
        if process.returncode:
            attempts.append(row)
            continue
        with dest.open() as h:
            (nv, nf) = map(int, h.readline().split())
            v = np.loadtxt(h, max_rows=nv)
            f = np.loadtxt(h, dtype=int)
        cm = trimesh.Trimesh(v, f, process=False)
        si = intersect(v, f, name + '_check')
        keep = retained_faces(v, f, pv, pf)
        row.update(watertight=cm.is_watertight, winding=cm.is_winding_consistent, intersections=si.get('count'), original_source_faces_retained=int(keep.sum()), original_source_faces=len(pf))
        attempts.append(row)
        save(R6 / 'raw' / ('H_CLOSURE_' + key + '.json'), attempts)
        if cm.is_watertight and cm.is_winding_consistent and (si['status'] == 'PASS') and keep.all():

            def fk(t):
                return tuple(sorted((tuple((float(x) for x in p)) for p in t)))
            keys = {fk(t) for t in pv[pf]}
            isorig = np.array([fk(t) in keys for t in v[f]])
            return (v, f, pv, pf, rim, isorig, attempts)
    raise ValueError('No admissible remote closure retaining all original preparation triangles: ' + str(attempts))

def construct(m, key):
    (v, f, pv, pf, rim, isorig, attempts) = closure(m, key)
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
    (_, de, _) = fast_nearest(v[f[~isorig]], iv)
    (iv, inf) = continuous_cut(iv, inf, de - dc - 1e-07)
    cache = D6 / (key + '_H_inner_raw.npz')
    np.savez_compressed(cache, vertices=iv, faces=inf, prep_vertices=pv, prep_faces=pf)
    cm = trimesh.Trimesh(iv, inf, process=False)
    ll = loops(cm)
    info = dict(closure_attempts=attempts, loops=[len(l) for l in ll], components=len(trimesh.graph.connected_components(cm.face_adjacency, nodes=np.arange(len(inf)))))
    if len(ll) != 1 or info['components'] != 1:
        raise ValueError('H raw cap non-disk ' + str(info))
    for level in range(4):
        (sd, _, q, n) = igl.signed_distance(iv, v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
        iv += (boundary_gap(q, rim) - sd)[:, None] * n
        (q, g, _) = fast_nearest(pv[pf], iv[inf].mean(1))
        err = abs(g - boundary_gap(q, rim))
        bad = err > 0.005
        if not bad.any() or level == 3:
            break
        (iv, inf) = subdiv(iv, inf, bad)
    cap = trimesh.Trimesh(iv, inf, process=False)
    trimesh.repair.fix_normals(cap, multibody=True)
    (iv, inf) = (cap.vertices, cap.faces)
    (_, _, _, nn) = igl.signed_distance(iv[inf].mean(1), v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
    if np.einsum('ij,ij->i', np.cross(iv[inf][:, 1] - iv[inf][:, 0], iv[inf][:, 2] - iv[inf][:, 0]), nn).sum() > 0:
        inf = inf[:, ::-1]
    info['centroid_residual_max_mm'] = float(err.max())
    return (iv, inf, info)

def run(allcases=False):
    cases = inputs()
    cases.sort(key=lambda q: (q[0]['key'] != '079905ebf9504544_molar', q[0]['key']))
    tag = 'H_ALL' if allcases else 'H_FIRST'
    rows = []
    st = time.perf_counter()
    if not allcases:
        cases = cases[:1]
    for (a, b) in cases:
        t = time.perf_counter()
        r = dict(key=a['key'], family=a['family'], status='REJECTED')
        m = npz(a['mesh_path'])
        try:
            (iv, inf, info) = construct(m, a['key'])
            old = npz(D / 'B/distance_local_thickening' / (a['key'] + '.npz'))
            (ov, of, ri) = repair_outer(old, a['key'])
            (v, f, roles) = join(ov, of, iv, inf)
            dest = D6 / (a['key'] + '_H.npz')
            np.savez_compressed(dest, vertices=v, faces=f, roles=roles, prep_triangles=m['vertices'][m['faces'][m['roles'] == 1]])
            r.update(status='GENERATED', mesh_path=dest, mesh_sha256=sha(dest), offset=info, outer_repair=ri)
        except Exception as e:
            r['reason'] = repr(e)
        r['seconds'] = time.perf_counter() - t
        rows.append(r)
        save(R6 / f'raw/{tag}_GENERATION.json', rows)
        state6(tag + '_GENERATING', r, 'Freeze all generated candidates and rejections before scoring')
        print(clean(r), flush=True)
    freeze(R6 / f'FROZEN_PREDICTIONS_{tag}.json', dict(claim_type='capability', rows=rows, prereg_sha256=sha(R6 / 'PREREG_H.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run('--all' in sys.argv)
