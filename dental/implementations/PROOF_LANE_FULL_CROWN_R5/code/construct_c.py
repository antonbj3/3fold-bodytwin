from construct_b import *
from skimage.measure import marching_cubes
from trimesh.ray.ray_pyembree import RayMeshIntersector

def implicit_inner(m, key):
    (pv, pf) = compact(m['vertices'], m['faces'][m['roles'] == 1][:, ::-1])
    top = trimesh.Trimesh(pv, pf, process=False)
    ids = loops(top)[0]
    rim = pv[ids]
    N = len(pv)
    low = rim.copy()
    low[:, 2] = rim[:, 2].min() - 1.0
    bc = np.r_[rim[:, :2].mean(0), low[0, 2]]
    vv = np.r_[pv, low, bc[None]]
    faces = pf.tolist()
    nt = len(pf)
    for (i, a) in enumerate(ids):
        k = (i + 1) % len(ids)
        b = ids[k]
        faces.extend([[int(a), N + i, N + k], [int(a), N + k, int(b)], [N + i, N + len(ids), N + k]])
    closed = trimesh.Trimesh(vv, np.array(faces), process=False)
    trimesh.repair.fix_normals(closed, multibody=True)
    if closed.volume < 0:
        closed.invert()
    (cv, cf) = (closed.vertices, closed.faces)
    step = 0.035
    lo = pv.min(0) - 0.12
    hi = pv.max(0) + 0.12
    axes = [np.arange(lo[i], hi[i] + step, step) for i in range(3)]
    shape = tuple(map(len, axes))
    size = int(np.prod(shape))
    if size > 16000000:
        raise ValueError('grid point cap')
    proxy = []
    for t in cv[cf]:
        n = max(1, int(np.ceil(np.linalg.norm(t - np.roll(t, 1, axis=0), axis=1).max() / 0.2)))
        weights = np.array([(i / n, j / n) for i in range(n + 1) for j in range(n + 1 - i)])
        proxy.append(t[0] + weights[:, 0, None] * (t[1] - t[0]) + weights[:, 1, None] * (t[2] - t[0]))
    proxy = np.concatenate(proxy)
    tree = cKDTree(proxy)
    ray = RayMeshIntersector(closed)
    field = np.empty(size, np.float32)
    near_count = 0
    sign_mismatches = 0
    sign_probes = 0
    (xx, yy) = np.meshgrid(axes[0], axes[1], indexing='ij')
    origins = np.c_[xx.ravel(), yy.ravel(), np.full(xx.size, cv[:, 2].min() - 0.1)]
    directions = np.tile([0.0, 0.0, 1.0], (len(origins), 1))
    (locations, rayids, _) = ray.intersects_location(origins, directions, multiple_hits=True)
    hits = {}
    for (ri, z) in zip(rayids, locations[:, 2]):
        hits.setdefault(int(ri), []).append(float(z))
    column_sign = np.ones(shape, np.int8)
    for (ri, zs) in hits.items():
        (ix, iy) = divmod(ri, len(axes[1]))
        column_sign[ix, iy] = np.where(np.searchsorted(np.unique(zs), axes[2], side='right') % 2, -1, 1)
    column_sign = column_sign.ravel()
    print('C_GRID', key, shape, 'column rays', len(origins), flush=True)
    for i in range(0, size, 65536):
        j = np.arange(i, min(i + 65536, size))
        ijk = np.array(np.unravel_index(j, shape)).T
        points = np.column_stack([axes[k][ijk[:, k]] for k in range(3)])
        dd = tree.query(points, distance_upper_bound=0.26, workers=1)[0]
        near = dd <= 0.26
        inside = column_sign[j] < 0
        values = column_sign[j].astype(float)
        if near.any():
            (sd, fi, q, n) = igl.signed_distance(points[near], cv, cf, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
            values[near] = sd - boundary_gap(q, rim)
            near_count += int(near.sum())
        checks = np.flatnonzero(~near)[::4096]
        if len(checks):
            (ss, _, _, _) = igl.signed_distance(points[checks], cv, cf, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
            sign_mismatches += int(np.sum((ss < 0) != inside[checks]))
            sign_probes += len(checks)
        field[i:i + len(j)] = values
    (iv, inf, _, _) = marching_cubes(field.reshape(shape), level=0.0, spacing=(step,) * 3, allow_degenerate=False)
    iv = iv.astype(float) + lo
    cache = D / 'C_IMPLICIT_CACHE' / (key + '.npz')
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, vertices=iv, faces=inf, prep_vertices=cv, prep_faces=cf, original_cap_faces=np.array(nt), rim=rim)
    print('C_ISO', key, len(iv), len(inf), flush=True)
    (sd, j, q, n) = igl.signed_distance(iv[inf].mean(1), cv, cf, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
    sel = j < nt
    (iv, inf) = compact(iv, inf[sel])
    cap = trimesh.Trimesh(iv, inf, process=False)
    if len(loops(cap)) != 1:
        raise ValueError('implicit cap boundary count ' + str(len(loops(cap))))
    if len(cap.split(only_watertight=False)) != 1:
        raise ValueError('implicit cap multiple components')
    for _ in range(8):
        (sd, idx, q, n) = igl.signed_distance(iv, cv, cf, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
        res = boundary_gap(q, rim) - sd
        iv += res[:, None] * n
    tri = iv[inf]
    (sd, idx, q, n) = igl.signed_distance(tri.mean(1), cv, cf, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
    gg = boundary_gap(q, rim)
    res = abs(sd - gg)
    norm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    agreement = np.einsum('ij,ij->i', norm, n).sum()
    if agreement > 0:
        inf = inf[:, ::-1]
    return (iv, inf, dict(max_plane_offset_residual_mm=float(res.max()), metric='centroid implicit offset residual against computational closed prep; independent original-cap distance follows', implicit_grid_mm=step, grid_points=size, near_band_queries=near_count, far_sign_probe_mismatches=sign_mismatches, far_sign_probes=sign_probes, selected_fraction=float(sel.mean()), closure_watertight=closed.is_watertight, closure_self_intersection='UNKNOWN', vertices=len(iv), faces=len(inf)))

def run():
    st = time.perf_counter()
    rows = read(R / 'raw/C_GENERATION.json') if (R / 'raw/C_GENERATION.json').exists() else []
    for (rec, rr) in inputs():
        if sum((r['key'] == rec['key'] for r in rows)) == 2:
            continue
        rows = [r for r in rows if r['key'] != rec['key']]
        m = npz(rec['mesh_path'])
        t = time.perf_counter()
        try:
            (iv, inf, offset) = implicit_inner(m, rec['key'])
            error = None
        except Exception as exc:
            error = repr(exc)
        for method in ['implicit_offset_only', 'implicit_local_thickening']:
            row = dict(key=rec['key'], family=rec['family'], method=method, resolution='PER_POINT', status='REJECTED')
            if error:
                row.update(reason=error)
            else:
                try:
                    (ov, of) = compact(m['vertices'], m['faces'][m['roles'] == 0])
                    info = {}
                    if method == 'implicit_local_thickening':
                        (ov, of, info) = outer_repair(m, iv[inf])
                    (v, f, roles) = join(ov, of, iv, inf)
                    dest = D / 'C' / method / (rec['key'] + '.npz')
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    np.savez_compressed(dest, vertices=v, faces=f, roles=roles, prep_triangles=m['vertices'][m['faces'][m['roles'] == 1]])
                    row.update(mesh_path=dest, mesh_sha256=sha(dest), offset=offset, repair=info, status='GENERATED')
                except Exception as exc:
                    row.update(reason=repr(exc))
            row['construction_seconds'] = time.perf_counter() - t
            rows.append(row)
            save(R / 'raw/C_GENERATION.json', rows)
            print(rec['key'], method, row['status'], row.get('reason', offset if not error else ''), flush=True)
        checkpoint('C_GENERATION', len(rows), 'Freeze all implicit-offset candidates and rejected cases before scoring')
    freeze(R / 'FROZEN_PREDICTIONS_C.json', dict(claim_type='capability', prereg_sha256=sha(R / 'PREREG_C.json'), code_sha256=sha(Path(__file__)), rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run()
