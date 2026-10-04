from local import *
from clip import continuous_cut
from skimage.measure import marching_cubes
from check_intersections import intersect

def make_inner_c(m, key):
    (pv, pf) = compact(m['vertices'], m['faces'][m['roles'] == 1][:, ::-1])
    top = trimesh.Trimesh(pv, pf, process=False)
    ids = loops(top)[0]
    rim = pv[ids]
    apex = pv[np.setdiff1d(np.arange(len(pv)), ids)[0]]
    ex = apex + 1.5 * (rim - apex)
    N = len(pv)
    v = np.r_[pv, ex]
    f = pf.tolist()
    edges = {tuple(e) for face in pf for e in zip(face, np.roll(face, -1))}
    for (j, a) in enumerate(ids):
        k = (j + 1) % len(ids)
        b = ids[k]
        fs = [[int(a), N + j, N + k], [int(a), N + k, int(b)]]
        if (int(a), int(b)) not in edges:
            fs = [q[::-1] for q in fs]
        f.extend(fs)
    f = np.array(f, int)
    ext = v[f[len(pf):]]
    orig = pv[pf]
    si = intersect(v, f, key + '_continued')
    save(R6 / 'raw' / ('CONTINUED_' + key + '.json'), si)
    if si['status'] != 'PASS':
        raise ValueError('continued source intersects: ' + str(si.get('count')))
    step = 0.025
    lo = pv.min(0) - 0.12
    hi = pv.max(0) + 0.12
    axes = [np.arange(lo[i], hi[i] + step, step) for i in range(3)]
    shape = tuple(map(len, axes))
    size = int(np.prod(shape))
    if size > 12000000:
        raise ValueError('grid point budget ' + str(size))
    field = np.empty(size, np.float32)
    for k in range(0, size, 65536):
        j = np.arange(k, min(k + 65536, size))
        ijk = np.array(np.unravel_index(j, shape)).T
        points = np.column_stack([axes[l][ijk[:, l]] for l in range(3)])
        (sd, _, q, n) = igl.signed_distance(points, v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
        field[k:k + len(j)] = sd - boundary_gap(q, rim)
    (iv, inf, _, _) = marching_cubes(field.reshape(shape), level=0.0, spacing=(step,) * 3, allow_degenerate=False)
    iv = iv.astype(float) + lo
    (_, dc, _) = fast_nearest(orig, iv)
    (_, de, _) = fast_nearest(ext, iv)
    (iv, inf) = continuous_cut(iv, inf, de - dc - 1e-07)
    cap = trimesh.Trimesh(iv, inf, process=False)
    ll = loops(cap)
    ncomp = len(cap.split(only_watertight=False))
    info = dict(boundary_loops=[len(l) for l in ll], components=ncomp, grid_points=size, grid_step_mm=step)
    cache = D6 / (key + '_C_inner_raw.npz')
    np.savez_compressed(cache, vertices=iv, faces=inf, prep_vertices=pv, prep_faces=pf)
    save(R6 / 'raw' / ('C_TOPOLOGY_' + key + '.json'), info)
    if len(ll) != 1 or ncomp != 1:
        raise ValueError('non-disk cap ' + str(info))
    for level in range(4):
        for _ in range(5):
            (sd, _, q, n) = igl.signed_distance(iv, v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
            iv += (boundary_gap(q, rim) - sd)[:, None] * n
        tri = iv[inf]
        (q, g, _) = fast_nearest(orig, tri.mean(1))
        err = abs(g - boundary_gap(q, rim))
        bad = err > 0.005
        if not bad.any() or level == 3:
            break
        (iv, inf) = subdiv(iv, inf, bad)
    (sd, _, q, n) = igl.signed_distance(iv[inf].mean(1), v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
    if np.einsum('ij,ij->i', np.cross(iv[inf][:, 1] - iv[inf][:, 0], iv[inf][:, 2] - iv[inf][:, 0]), n).sum() > 0:
        inf = inf[:, ::-1]
    info['gap_centroid_residual_mm'] = float(err.max())
    return (iv, inf, info)

def run(allcases=False):
    cases = inputs()
    cases.sort(key=lambda q: (q[0]['key'] != '079905ebf9504544_molar', q[0]['key']))
    rows = []
    st = time.perf_counter()
    if not allcases:
        cases = cases[:1]
    for (a, b) in cases:
        t = time.perf_counter()
        r = dict(key=a['key'], family=a['family'], status='REJECTED')
        m = npz(a['mesh_path'])
        try:
            (iv, inf, info) = make_inner_c(m, a['key'])
            old = npz(D / 'B/distance_local_thickening' / (a['key'] + '.npz'))
            (ov, of) = compact(old['vertices'], old['faces'][old['roles'] == 0])
            (v, f, roles) = join(ov, of, iv, inf)
            dest = D6 / (a['key'] + '_C.npz')
            np.savez_compressed(dest, vertices=v, faces=f, roles=roles, prep_triangles=m['vertices'][m['faces'][m['roles'] == 1]])
            r.update(status='GENERATED', mesh_path=dest, mesh_sha256=sha(dest), offset=info)
        except Exception as e:
            r['reason'] = repr(e)
        r['seconds'] = time.perf_counter() - t
        rows.append(r)
        save(R6 / 'raw/C_GENERATION.json', rows)
        state6('C_GENERATING', r, 'Freeze candidate before native scoring; no rejected tooth replacement')
        print(clean(r), flush=True)
    tag = 'C_ALL' if allcases else 'C_FIRST'
    freeze(R6 / f'FROZEN_PREDICTIONS_{tag}.json', dict(claim_type='capability', rows=rows, prereg_sha256=sha(R6 / 'PREREG_C.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run('--all' in sys.argv)
