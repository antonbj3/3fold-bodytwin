import os
os.environ['IGL_NUM_THREADS'] = '1'
from construct_a import *
from threadpoolctl import threadpool_limits
threadpool_limits(limits=1)
sys.path.insert(0, str(D / 'deps'))
import igl
NEAREST_CACHE = {}

def fast_nearest(tri, points):
    key = id(tri)
    if key not in NEAREST_CACHE:
        if len(NEAREST_CACHE) >= 3:
            NEAREST_CACHE.clear()
        vv = np.ascontiguousarray(tri.reshape(-1, 3))
        ff = np.arange(len(vv)).reshape(-1, 3)
        tree = igl.AABB()
        tree.init(vv, ff)
        NEAREST_CACHE[key] = (tri, vv, ff, tree)
    (_, vv, ff, tree) = NEAREST_CACHE[key]
    (d, j, q) = tree.squared_distance(vv, ff, np.ascontiguousarray(points))
    return (q, np.sqrt(np.maximum(0, d)), j)

def subdiv(v, f, marked):
    ff = f[marked]
    length = np.linalg.norm(v[ff] - np.roll(v[ff], -1, axis=1), axis=2)
    k = length.argmax(1)
    chosen = set((tuple(sorted((int(face[j]), int(face[(j + 1) % 3])))) for (face, j) in zip(ff, k)))
    if not chosen:
        return (v, f)
    verts = v.tolist()
    mid = {}
    for (a, b) in sorted(chosen):
        mid[a, b] = len(verts)
        verts.append(((v[a] + v[b]) / 2).tolist())
    out = []
    for face in f:
        mids = [mid.get(tuple(sorted((int(face[j]), int(face[(j + 1) % 3]))))) for j in range(3)]
        n = sum((x is not None for x in mids))
        if n == 0:
            out.append(face.tolist())
        elif n == 1:
            j = next((j for j in range(3) if mids[j] is not None))
            (a, b, c) = np.roll(face, -j)
            m = mids[j]
            out.extend([[a, m, c], [m, b, c]])
        elif n == 2:
            j = next((j for j in range(3) if mids[j] is not None and mids[(j - 1) % 3] is not None))
            a = face[j]
            b = face[(j + 1) % 3]
            c = face[(j - 1) % 3]
            u = mids[j]
            w = mids[(j - 1) % 3]
            out.extend([[a, u, w], [w, u, b], [w, b, c]])
        else:
            (a, b, c) = face
            (u, w, z) = mids
            out.extend([[a, u, z], [u, b, w], [z, w, c], [u, w, z]])
    return (np.array(verts), np.array(out, int))

def compact(v, f):
    ids = np.unique(f)
    lookup = np.full(len(v), -1, int)
    lookup[ids] = np.arange(len(ids))
    return (v[ids].copy(), lookup[f])

def boundary_gap(p, boundary):
    ed = np.c_[np.arange(len(boundary)), np.roll(np.arange(len(boundary)), -1)]
    (dd, _, _) = igl.point_mesh_squared_distance(np.asarray(p), boundary, ed)
    return 0.025 + 0.025 * np.clip(np.sqrt(np.maximum(0, dd)), 0, 1)

def make_inner(m):
    (v, f) = (m['vertices'], m['faces'])
    roles = m['roles']
    (pv, pf) = compact(v, f[roles == 1][:, ::-1])
    orig = trimesh.Trimesh(pv, pf, process=False)
    rim = pv[loops(orig)[0]]
    (iv, inf) = (pv.copy(), pf.copy())
    history = []
    for _ in range(5):
        tri = iv[inf]
        long = np.linalg.norm(tri - np.roll(tri, 1, axis=1), axis=2).max(1) > 0.3
        if not long.any():
            break
        (nv, nf) = subdiv(iv, inf, long)
        if len(nv) > 150000:
            raise ValueError('inner initialization vertex cap')
        (iv, inf) = (nv, nf)
    normal = igl.per_vertex_normals(iv, inf)
    iv += boundary_gap(iv, rim)[:, None] * normal
    for level in range(6):
        for it in range(8):
            (sd, idx, q, n) = igl.signed_distance(iv, pv, pf, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
            gg = boundary_gap(q, rim)
            delta = gg - sd
            iv += np.clip(delta, -0.1, 0.1)[:, None] * n
            if np.max(abs(delta)) < 1e-07:
                break
        tri = iv[inf]
        cen = tri.mean(1)
        (d, _, q, _) = igl.signed_distance(cen, pv, pf, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
        err = abs(d - boundary_gap(q, rim))
        edge = np.linalg.norm(tri - np.roll(tri, 1, axis=1), axis=2).max(1)
        bad = (err > 0.005) & (edge > 0.05)
        history.append(dict(level=level, vertices=len(iv), faces=len(inf), max_centroid_residual_mm=float(err.max()), refine=int(bad.sum())))
        if not bad.any() or level == 5:
            break
        (nv, nf) = subdiv(iv, inf, bad)
        if len(nv) > 150000:
            history[-1]['vertex_budget_exhausted'] = True
            break
        (iv, inf) = (nv, nf)
    return (iv, inf[:, ::-1], dict(history=history, max_plane_offset_residual_mm=float(err.max()), metric='centroid oriented-distance field residual, vertex/probe checks follow in scorer', maximum_vertex_offset_mm=None), pv, pf)

def outer_repair(m, inner):
    (v, f) = compact(m['vertices'], m['faces'][m['roles'] == 0])
    base_mesh = trimesh.Trimesh(v, f, process=False)
    base = v.copy()
    base_faces = f.copy()
    history = []
    for level in range(6):
        obj = trimesh.Trimesh(v, f, process=False)
        bound = loops(obj)[0]
        free = np.ones(len(v), bool)
        free[bound] = False
        (_, idx, _) = igl.point_mesh_squared_distance(v, base, base_faces)
        normal = base_mesh.face_normals[idx]
        for iteration in range(3):
            (q, d, j) = fast_nearest(inner, v)
            active = free & (d < 0.6 - 1e-06)
            if not active.any():
                break
            x = v[active]
            n = normal[active]
            low = np.zeros(len(x))
            high = np.full(len(x), 0.6)
            for k in range(5):
                (_, dd, _) = fast_nearest(inner, x + high[:, None] * n)
                short = dd < 0.6
                if not short.any():
                    break
                high[short] = np.minimum(5.0, 2 * high[short])
            for k in range(20):
                mid = (low + high) / 2
                (_, dd, _) = fast_nearest(inner, x + mid[:, None] * n)
                short = dd < 0.6
                low[short] = mid[short]
                high[~short] = mid[~short]
            v[active] = x + high[:, None] * n
        tri = v[f]
        cen = tri.mean(1)
        (_, dc, _) = fast_nearest(inner, cen)
        rad = np.linalg.norm(tri - cen[:, None], axis=2).max(1)
        lb = dc - rad
        bad = lb < 0.5
        history.append(dict(level=level, vertices=len(v), faces=len(f), sampled_min_mm=float(dc.min()), cover_min_mm=float(lb.min()), refine=int(bad.sum())))
        if not bad.any() or level == 5:
            break
        (nv, nf) = subdiv(v, f, bad)
        if len(nv) > 150000:
            history[-1]['vertex_budget_exhausted'] = True
            break
        (v, f) = (nv, nf)
    (q, disp, _) = fast_nearest(base[base_faces], v)
    return (v, f, dict(history=history, max_displacement_mm=float(disp.max()), p95_displacement_mm=float(np.quantile(disp, 0.95)), boundary_identity_error_mm=0.0, displacement_metric='distance to original exterior after refinement; exact initial boundary vertices retained'))

def join(ov, of, iv, inf):
    from generate import curve_parameter
    out = trimesh.Trimesh(ov, of, process=False)
    inn = trimesh.Trimesh(iv, inf, process=False)
    a = loops(out)[0]
    b = loops(inn)[0]
    from generate import correspond
    (_, info) = correspond(ov[a], iv[b])
    b = b[::-1] if info['flipped'] else b
    start = int(np.argmin(np.linalg.norm(iv[b] - ov[a[0]], axis=1)))
    b = np.roll(b, -start)
    (sa, la) = curve_parameter(ov[a])
    (sb, lb) = curve_parameter(iv[b])
    sa /= la
    sb /= lb
    N = len(ov)
    faces = of.tolist() + (inf + N).tolist()
    roles = [0] * len(of) + [1] * len(inf)
    i = j = 0
    while i < len(a) or j < len(b):
        na = sa[i + 1] if i < len(a) else np.inf
        nb = sb[j + 1] if j < len(b) else np.inf
        if na <= nb:
            faces.append([int(a[i % len(a)]), int(a[(i + 1) % len(a)]), N + int(b[j % len(b)])])
            i += 1
        else:
            faces.append([int(a[i % len(a)]), N + int(b[(j + 1) % len(b)]), N + int(b[j % len(b)])])
            j += 1
        roles.append(2)
    mesh = trimesh.Trimesh(np.r_[ov, iv], np.array(faces), process=False)
    trimesh.repair.fix_normals(mesh, multibody=True)
    if mesh.volume < 0:
        mesh.invert()
    return (mesh.vertices, mesh.faces, np.array(roles))

def run():
    st = time.perf_counter()
    rows = read(R / 'raw/B_GENERATION.json') if (R / 'raw/B_GENERATION.json').exists() else []
    for (rec, rr) in inputs():
        if sum((r['key'] == rec['key'] for r in rows)) == 2:
            continue
        rows = [r for r in rows if r['key'] != rec['key']]
        m = npz(rec['mesh_path'])
        (iv, inf, offset, pv, pf) = make_inner(m)
        inner = iv[inf]
        for method in ['distance_offset_only', 'distance_local_thickening']:
            t0 = time.perf_counter()
            (ov, of) = compact(m['vertices'], m['faces'][m['roles'] == 0])
            info = {}
            if method == 'distance_local_thickening':
                (ov, of, info) = outer_repair(m, inner)
            (v, f, roles) = join(ov, of, iv, inf)
            dest = D / 'B' / method / (rec['key'] + '.npz')
            dest.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(dest, vertices=v, faces=f, roles=roles, prep_triangles=m['vertices'][m['faces'][m['roles'] == 1]])
            row = dict(key=rec['key'], family=rec['family'], method=method, mesh_path=dest, mesh_sha256=sha(dest), offset=offset, repair=info, construction_seconds=time.perf_counter() - t0, resolution='PER_POINT')
            rows.append(row)
            save(R / 'raw/B_GENERATION.json', rows)
            print(rec['key'], method, 'faces', len(f), 'offset', offset['max_plane_offset_residual_mm'], 'repair', info.get('history', [{}])[-1], flush=True)
        checkpoint('B_GENERATING', len(rows), 'Complete fixed18 teeth; freeze predictions before independent scoring')
    freeze(R / 'FROZEN_PREDICTIONS_B.json', dict(claim_type='capability', prereg_sha256=sha(R / 'PREREG_B.json'), code_sha256=sha(Path(__file__)), rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run()
