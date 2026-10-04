"""Research-only straight graft planner. Source geometry is read-only.

Geometry supports, donor radius, coordinate orientation and occlusal landmarks
have no independent anatomical calibration. Post is a virtual expert target.
"""
from dental_release.paths import expand as _release_expand
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_k] = '2'
import argparse, itertools, json, resource, time
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from freeze import ROOT, OLD, DATA, sha, write, state, now
CACHE = Path(_release_expand('@DENTAL_WORK_ROOT@/G_mandible_postop/R2/components'))
PARAM = json.loads((ROOT / 'PREREG_R1.json').read_text())['parameters']
MAN = json.loads((OLD / 'CASE_MANIFEST.json').read_text())

def dist(p, q):
    return cKDTree(q).query(p, workers=1)[0]

def transform(q, T):
    return q @ T[:3, :3].T + T[:3, 3]

def stats(d):
    if not len(d):
        return dict(n=0, mean_mm=None, p95_mm=None, max_mm=None)
    return dict(n=len(d), mean_mm=float(np.mean(d)), p95_mm=float(np.quantile(d, 0.95)), max_mm=float(np.max(d)))

def load(cid, role):
    pp = CACHE / f'{cid}_{role}.npz'
    if not pp.exists():
        pp = DATA / 'components' / f'{cid}_{role}.npz'
        if not pp.exists():
            prepare_component(cid, role, pp)
    meta = json.loads(pp.with_suffix('.json').read_text())
    if sha(pp) != meta['cache_sha256']:
        raise ValueError('cache SHA mismatch')
    if meta['source_md5'] != MAN['roles'][cid][role]['source_md5']:
        raise ValueError('source association mismatch')
    with np.load(pp) as z:
        q = z['jaw_points'].astype(float)
        a = z['anchor_points'].astype(float)
    return (q, a, dict(path=str(pp), sha256=meta['cache_sha256'], source_md5=meta['source_md5'], support='GEOMETRY_PROXY_PENDING_ANATOMICAL_REVIEW'))

def prepare_component(cid, role, out):
    """Reuse predecessor's component operation, write only this lane's cache."""
    import hashlib
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    src = Path(MAN['roles'][cid][role]['path'])
    tick = time.perf_counter()
    h = hashlib.md5()
    with src.open('rb') as ff:
        for b in iter(lambda : ff.read(1048576), b''):
            h.update(b)
    if h.hexdigest() != MAN['roles'][cid][role]['source_md5']:
        raise ValueError('source MD5 mismatch')
    dt = np.dtype([('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')])
    soup = np.memmap(src, dtype=dt, offset=84, mode='r')
    (v, ix) = np.unique(soup['v'].reshape(-1, 3), axis=0, return_inverse=True)
    f = ix.reshape(-1, 3)
    edges = np.r_[f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]
    g = coo_matrix((np.ones(len(edges), 'i1'), (edges[:, 0], edges[:, 1])), shape=(len(v), len(v))).tocsr()
    (_, labels) = connected_components(g, directed=False)
    counts = np.bincount(labels)
    ai = int(np.argmax(counts))
    anchor = v[labels == ai]
    ac = anchor.mean(0)
    az = np.median(anchor[:, 2])
    jaw = []
    for k in np.argsort(counts)[::-1]:
        if counts[k] < 1000:
            continue
        points = v[labels == k]
        if k != ai and points.mean(0)[2] <= ac[2] - 10 and (points[:, 2].min() <= az - 15):
            jaw.append(k)
    if not jaw:
        raise ValueError('UNKNOWN_COMPONENT_PROXY')
    fl = labels[f[:, 0]]
    jf = f[np.isin(fl, jaw)]
    af = f[fl == ai]
    points = surface_samples([(v, jf)], 60000, 20261002 + int(cid))
    ap = surface_samples([(v, af)], 60000, 20261001 + int(cid))
    out.parent.mkdir(exist_ok=True)
    np.savez_compressed(out, jaw_points=points.astype('f4'), anchor_points=ap.astype('f4'))
    write(out.with_suffix('.json'), dict(case=cid, role=role, cache_sha256=sha(out), source_md5=h.hexdigest(), source=str(src), wall_seconds=time.perf_counter() - tick, component_operation='same connected-component selection as predecessor experiment_r2.py', jaw_component_ids=[int(k) for k in jaw], anchor_component_id=ai, anatomical_support='PENDING_INDEPENDENT_REVIEW'))

def rigid_fit(p, q):
    """Same rigid SVD/trimmed ICP mechanism as predecessor, no scale fit."""
    tree = cKDTree(q)
    T = np.eye(4)
    prev = np.inf
    for _ in range(PARAM['register_iterations']):
        x = transform(p, T)
        (d, ix) = tree.query(x, workers=1)
        keep = d <= np.quantile(d, 0.9)
        aa = x[keep]
        bb = q[ix[keep]]
        ac = aa.mean(0)
        bc = bb.mean(0)
        (U, _, Vt) = np.linalg.svd((aa - ac).T @ (bb - bc))
        rot = Vt.T @ U.T
        if np.linalg.det(rot) < 0:
            Vt[-1] *= -1
            rot = Vt.T @ U.T
        D = np.eye(4)
        D[:3, :3] = rot
        D[:3, 3] = bc - rot @ ac
        T = D @ T
        loss = np.mean(d[keep] ** 2)
        if abs(prev - loss) < 1e-08:
            break
        prev = loss
    d = dist(transform(p, T), q)
    return (T, stats(d[d <= np.quantile(d, 0.9)]))

def center(a):
    return (np.quantile(a, 0.005, axis=0) + np.quantile(a, 0.995, axis=0)) / 2

def reflection(q, a):
    c = center(a)
    out = q.copy()
    out[:, 0] = 2 * c[0] - out[:, 0]
    return (np.r_[q, out], dict(normal=[1.0, 0.0, 0.0], offset_mm=float(c[0]), method='skull bounding-centre sagittal proxy'))

def extract_rail(q, pre, a, radius=6.0):
    """Fixed angular inferior contour operator; no manual per-case edits."""
    c = center(a)
    theta = np.arctan2(q[:, 0] - c[0], -(q[:, 1] - c[1]))
    edges = np.deg2rad(np.linspace(-100, 100, 74))
    bi = np.digitize(theta, edges) - 1
    dd = dist(q, pre)
    vals = []
    rows = []
    for b in range(73):
        ix = np.flatnonzero(bi == b)
        if len(ix) < 20:
            continue
        low = ix[q[ix, 2] <= np.quantile(q[ix, 2], 0.22)]
        if len(low) < 4:
            continue
        point = np.r_[np.median(q[low, :2], axis=0), np.quantile(q[ix, 2], 0.08) + radius]
        frac = float(np.mean(dd[low] > 3))
        rows.append(dict(bin=b, n=len(ix), inferior_n=len(low), added_fraction=frac, missing=frac > 0.5))
        vals.append(point)
    if not vals:
        raise ValueError('UNKNOWN_NO_ANGULAR_SUPPORT')
    vals = np.array(vals)
    missing = np.flatnonzero([r['missing'] for r in rows])
    if not len(missing):
        raise ValueError('UNKNOWN_NO_SUPPORTED_MISSING_RAIL')
    lo = max(0, int(missing[0]) - 1)
    hi = min(len(vals) - 1, int(missing[-1]) + 1)
    curve = vals[lo:hi + 1]
    used = rows[lo:hi + 1]
    if len(curve) < 3:
        raise ValueError('UNKNOWN_TOO_SHORT_RAIL')
    return (curve, dict(origin_proxy_mm=c.tolist(), bins=used, all_bins=rows, angular_support_gaps=int(sum((used[i + 1]['bin'] - used[i]['bin'] > 1 for i in range(len(used) - 1)))), selection='bounding interval of missing inferior bins; may cross retained anatomy', radius_mm=radius, arc_length_mm=float(np.linalg.norm(np.diff(curve, axis=0), axis=1).sum())))

def chord_errors(curve):
    n = len(curve)
    E = np.full((n, n), np.inf)
    L = np.zeros((n, n))
    for i in range(n - 1):
        for j in range(i + 1, n):
            v = curve[j] - curve[i]
            l = np.linalg.norm(v)
            L[i, j] = l
            if l < 10:
                continue
            t = np.clip((curve[i:j + 1] - curve[i]) @ v / (l * l), 0, 1)
            d = np.linalg.norm(curve[i:j + 1] - (curve[i] + t[:, None] * v), axis=1)
            E[i, j] = d.max()
    return (E, L)

def dp_plan(curve, force_k=None):
    (E, L) = chord_errors(curve)
    n = len(curve)
    D = np.full((4, n), np.inf)
    D[0, 0] = 0.0
    parents = np.full((4, n), -1, int)
    for k in range(1, 4):
        for j in range(1, n):
            v = np.maximum(D[k - 1, :j], E[:j, j])
            i = int(np.argmin(v))
            D[k, j] = v[i]
            parents[k, j] = i
    feasible = [k for k in range(1, 4) if np.isfinite(D[k, -1])]
    if not feasible:
        raise ValueError('UNKNOWN_NO_10MM_SEGMENT_PARTITION')
    ks = [k for k in feasible if D[k, -1] <= 2.0]
    k = force_k if force_k is not None else min(ks) if ks else min(feasible, key=lambda k: D[k, -1])
    if k not in feasible:
        raise ValueError('requested count infeasible')
    j = n - 1
    ids = [j]
    for kk in range(k, 0, -1):
        j = int(parents[kk, j])
        ids.append(j)
    ids = ids[::-1]
    return (curve[ids], dict(count=k, indices=ids, max_chord_error_mm=float(D[k, -1]), errors_by_count={str(kk): float(D[kk, -1]) if kk in feasible else None for kk in range(1, 4)}, minimum_count_reaches_2mm=bool(ks), length_mm=float(sum((L[ids[i], ids[i + 1]] for i in range(k))))))

def exhaustive(curve, k, radius=None):
    """Independent chord distances using cross products and endpoint Voronoi regions."""
    n = len(curve)
    E = {}
    for i in range(n - 1):
        for j in range(i + 1, n):
            v = curve[j] - curve[i]
            l = float(np.linalg.norm(v))
            if l < 10:
                E[i, j] = np.inf
                continue
            x = curve[i:j + 1] - curve[i]
            proj = x @ v / l
            d = np.linalg.norm(np.cross(x, v), axis=1) / l
            d[proj < 0] = np.linalg.norm(x[proj < 0], axis=1)
            d[proj > l] = np.linalg.norm(curve[i:j + 1][proj > l] - curve[j], axis=1)
            E[i, j] = float(d.max())
    best = np.inf
    ids = None
    count = 0
    for inside in itertools.combinations(range(1, n - 1), k - 1):
        idx = (0,) + inside + (n - 1,)
        cost = max((E[idx[i], idx[i + 1]] for i in range(k)))
        count += 1
        if radius is not None:
            from constrained import cut_feasible
            if not cut_feasible(curve[list(idx)], radius)[0]:
                continue
        if cost < best:
            best = cost
            ids = idx
    return dict(max_chord_error_mm=float(best), indices=list(ids) if ids else None, partitions_examined=count)

def uniform(curve, k):
    arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(curve, axis=0), axis=1))]
    ids = [0] + [int(np.argmin(abs(arc - arc[-1] * i / k))) for i in range(1, k)] + [len(curve) - 1]
    if len(set(ids)) != len(ids):
        raise ValueError('UNKNOWN_DUPLICATE_UNIFORM_NODE')
    (E, L) = chord_errors(curve)
    if not all((np.isfinite(E[ids[i], ids[i + 1]]) for i in range(k))):
        raise ValueError('UNKNOWN_UNIFORM_SHORT_SEGMENT')
    return (curve[ids], dict(count=k, indices=ids, max_chord_error_mm=float(max((E[ids[i], ids[i + 1]] for i in range(k)))), length_mm=float(sum((L[ids[i], ids[i + 1]] for i in range(k))))))

def unit(v):
    d = np.linalg.norm(v)
    if d < 1e-10:
        raise ValueError('degenerate direction')
    return v / d

def planes(nodes):
    u = np.diff(nodes, axis=0)
    u = u / np.linalg.norm(u, axis=1)[:, None]
    normals = [u[0]] + [unit(u[i - 1] + u[i]) for i in range(1, len(u))] + [u[-1]]
    angles = [float(np.rad2deg(np.arccos(np.clip(u[i - 1] @ u[i], -1, 1)))) for i in range(1, len(u))]
    return (u, np.array(normals), angles)

def graft_mesh_legacy(nodes, radius=6.0, sides=48):
    """Exact circular-cylinder cut rings; adjacent rings share one ellipse."""
    (u, nn, angles) = planes(nodes)
    rings = []
    residual = []
    theta = np.arange(sides) * 2 * np.pi / sides
    for (j, (p, n)) in enumerate(zip(nodes, nn)):
        axis = u[max(0, j - 1)]
        if 0 < j < len(nodes) - 1 and angles[j - 1] > 1e-06:
            e = unit(u[j - 1] - u[j])
            other = unit(np.cross(n, e))
            major = radius / abs(axis @ n)
        else:
            ref = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
            e = unit(np.cross(n, ref))
            other = unit(np.cross(n, e))
            major = radius
        ring = p + major * np.cos(theta)[:, None] * e + radius * np.sin(theta)[:, None] * other
        rings.append(ring)
        residual.append(float(np.max(abs((ring - p) @ n))))
    meshes = []
    for i in range(len(u)):
        v = np.r_[rings[i], rings[i + 1], nodes[i:i + 2]]
        f = []
        for j in range(sides):
            k = (j + 1) % sides
            f.extend([[j, k, sides + k], [j, sides + k, sides + j], [2 * sides, k, j], [2 * sides + 1, sides + j, sides + k]])
        meshes.append((v, np.array(f, int)))
    detail = dict(radius_mm=radius, count=len(u), nodes_mm=nodes.tolist(), lengths_mm=np.linalg.norm(np.diff(nodes, axis=0), axis=1).tolist(), bend_angles_deg=angles, miter_obliquities_deg=[a / 2 for a in angles], cut_plane_normals=nn.tolist(), cut_plane_points_mm=nodes.tolist(), cut_plane_residual_max_mm=max(residual), segment_axes=u.tolist(), donor_anatomy='circular surrogate; not a patient fibula', host_contact_angles='UNKNOWN; host cut planes and cortex not labelled')
    return (meshes, detail)

def graft_mesh(nodes, radius=6.0, sides=48):
    """Transport one circular cross-section along each segment's actual axis.

The ring shared by neighboring axes lies on their equal-radius intersection
ellipse. Project it onto the next axis-normal plane; carry identical radial
vectors to the next cut. Every corresponding side edge is parallel to the axis.
Legacy independently phased rings are deliberately retained in R1/R2 snapshots.
"""
    (u, nn, angles) = planes(nodes)
    theta = np.arange(sides) * 2 * np.pi / sides
    ref = np.array([0.0, 0.0, 1.0]) if abs(u[0, 2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e = unit(np.cross(u[0], ref))
    other = unit(np.cross(u[0], e))
    rings = [nodes[0] + radius * (np.cos(theta)[:, None] * e + np.sin(theta)[:, None] * other)]
    for i in range(len(u)):
        d = rings[-1] - nodes[i]
        radial = d - (d @ u[i])[:, None] * u[i]
        denominator = u[i] @ nn[i + 1]
        if denominator <= 1e-10:
            raise ValueError('UNKNOWN_DEGENERATE_MITER')
        ring = nodes[i + 1] + radial - u[i] * (radial @ nn[i + 1] / denominator)[:, None]
        rings.append(ring)
    meshes = []
    edge_errors = []
    radius_errors = []
    plane_errors = []
    for i in range(len(u)):
        v = np.r_[rings[i], rings[i + 1], nodes[i:i + 2]]
        f = []
        for j in range(sides):
            k = (j + 1) % sides
            f.extend([[j, k, sides + k], [j, sides + k, sides + j], [2 * sides, k, j], [2 * sides + 1, sides + j, sides + k]])
        edges = rings[i + 1] - rings[i]
        perp = edges - (edges @ u[i])[:, None] * u[i]
        edge_errors.append(float(np.linalg.norm(perp, axis=1).max()))
        radial = v[:2 * sides] - nodes[i]
        radial -= (radial @ u[i])[:, None] * u[i]
        radius_errors.append(float(np.max(abs(np.linalg.norm(radial, axis=1) - radius))))
        plane_errors.extend([float(np.max(abs((rings[i] - nodes[i]) @ nn[i]))), float(np.max(abs((rings[i + 1] - nodes[i + 1]) @ nn[i + 1])))])
        meshes.append((v, np.array(f, int)))
    detail = dict(radius_mm=radius, count=len(u), nodes_mm=nodes.tolist(), lengths_mm=np.linalg.norm(np.diff(nodes, axis=0), axis=1).tolist(), bend_angles_deg=angles, miter_obliquities_deg=[a / 2 for a in angles], cut_plane_normals=nn.tolist(), cut_plane_points_mm=nodes.tolist(), cut_plane_residual_max_mm=max(plane_errors), longitudinal_edge_transverse_motion_max_mm=max(edge_errors), cylinder_ring_radius_error_max_mm=max(radius_errors), segment_axes=u.tolist(), donor_anatomy='circular surrogate; not patient fibula', host_contact_angles='UNKNOWN; host cut planes and cortex not labelled', mesh_operation='axis transport of common miter rings, straight faceted cylinder')
    return (meshes, detail)

def surface_samples(meshes, n=3000, seed=42):
    rng = np.random.default_rng(seed)
    parts = []
    for (v, f) in meshes:
        tri = np.asarray(v[f], dtype=np.float64)
        ar = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
        cumulative = np.cumsum(ar, dtype=np.float64)
        ix = np.searchsorted(cumulative, rng.random(n) * cumulative[-1])
        uv = rng.random((n, 2))
        uv[uv.sum(1) > 1] = 1 - uv[uv.sum(1) > 1]
        parts.append(tri[ix, 0] + uv[:, :1] * (tri[ix, 1] - tri[ix, 0]) + uv[:, 1:] * (tri[ix, 2] - tri[ix, 0]))
    return np.concatenate(parts)

def stl(path, mesh):
    (v, f) = mesh
    tri = v[f]
    normal = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normal /= np.maximum(np.linalg.norm(normal, axis=1)[:, None], 1e-12)
    dt = np.dtype([('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')])
    out = np.zeros(len(f), dt)
    out['n'] = normal
    out['v'] = tri
    with Path(path).open('wb') as ff:
        ff.write(b'X4 research graft surrogate'.ljust(80, b' '))
        ff.write(np.array([len(f)], '<u4').tobytes())
        ff.write(out.tobytes())

def graft_levelset(points, nodes, radius):
    """Exact sign of cylinder/halfspace CSG. Not Euclidean SDF at cut edges."""
    (u, nn, _) = planes(nodes)
    phi = np.full(len(points), np.inf)
    for i in range(len(u)):
        d = points - nodes[i]
        radial = np.linalg.norm(d - (d @ u[i])[:, None] * u[i], axis=1) - radius
        part = np.maximum.reduce([radial, -d @ nn[i], (points - nodes[i + 1]) @ nn[i + 1]])
        phi = np.minimum(phi, part)
    return phi

def predict(round_id):
    reg = json.loads((ROOT / f'PREREG_{round_id}.json').read_text())
    frozen = ROOT / f'FROZEN_PREDICTIONS_{round_id}.json'
    if frozen.exists():
        for r in json.loads(frozen.read_text())['predictions']:
            if sha(r['path']) != r['sha256']:
                raise ValueError('frozen prediction modified')
        print(round_id, 'FROZEN_REUSE', flush=True)
        return
    tick = time.perf_counter()
    records = []
    code = sha(__file__)
    (ROOT / 'sources' / f'planner_{round_id}.py').write_bytes(Path(__file__).read_bytes())
    for cid in reg['population']:
        start = time.perf_counter()
        out = DATA / round_id / cid
        out.mkdir(parents=True, exist_ok=True)
        row = dict(case=cid, round=round_id, created_utc=now(), prereg_sha256=sha(ROOT / f'PREREG_{round_id}.json'), code_sha256=code)
        arrays = {}
        inputs = []
        try:
            (pre, anchor, pmeta) = load(cid, 'Pre')
            inputs.append(pmeta)
            if round_id in ('R1', 'R4'):
                (target, plane) = reflection(pre, anchor)
                row['reflection_plane'] = plane
            else:
                (post, pa, pm) = load(cid, 'Post')
                inputs.append(pm)
                (T, fit) = rigid_fit(anchor[:4000], pa[:30000])
                target = transform(post, np.linalg.inv(T))
                row['supplied_information'] = 'Full virtual target used only to acquire inferior rail and common cranial alignment; assisted CAD planning, not blind reconstruction'
                row['registration'] = dict(T_Pre_to_Post=T.tolist(), fit=fit)
            if round_id == 'R4':
                from centroid_rail import rail as centroid
                (curve, rail) = centroid(target, pre, anchor)
            else:
                (curve, rail) = extract_rail(target, pre, anchor)
            row['rail'] = rail
            if round_id == 'R1':
                (nodes, plan) = dp_plan(curve)
            else:
                from constrained import constrained_plan
                (E, _) = chord_errors(curve)
                (nodes, plan) = constrained_plan(curve, E)
            row['dp'] = plan
            tcontrol = time.perf_counter()
            control = exhaustive(curve, plan['count'], None if round_id == 'R1' else 6.0)
            row['exhaustive_control'] = control
            row['control_seconds'] = time.perf_counter() - tcontrol
            row['solver_disagreement_mm'] = abs(control['max_chord_error_mm'] - plan['max_chord_error_mm'])
            methods = {'dp': nodes}
            try:
                (un, ud) = uniform(curve, plan['count'])
                methods['uniform'] = un
                row['uniform'] = ud
            except ValueError as exc:
                row['uniform_unknown'] = str(exc)
            if round_id != 'R1':
                from competitor import simplify
                try:
                    (on, od) = simplify(curve, 3)
                    methods['osteoopt'] = on
                    row['osteoopt'] = od
                    from constrained import cut_feasible
                    row['osteoopt']['non_crossing_miters'] = cut_feasible(on)[0]
                except ValueError as exc:
                    row['osteoopt_unknown'] = str(exc)
            row['grafts'] = {}
            for (name, nd) in methods.items():
                (meshes, detail) = (graft_mesh if round_id == 'R4' else graft_mesh_legacy)(nd)
                row['grafts'][name] = detail
                points = surface_samples(meshes, 3000, 20261002 + int(cid))
                arrays[name + '_nodes'] = nd.astype('f4')
                arrays[name + '_points'] = points.astype('f4')
                arrays[name + '_points_dense'] = surface_samples(meshes, 12000, 20261002 + int(cid)).astype('f4')
                for (si, mesh) in enumerate(meshes):
                    stl(out / f'{name}_segment_{si + 1}.stl', mesh)
            arrays['curve'] = curve.astype('f8')
            if round_id == 'R1':
                arrays['mirror_points'] = target.astype('f4')
            row['status'] = 'PREDICTED'
        except (ValueError, FileNotFoundError) as exc:
            row['status'] = 'UNKNOWN'
            row['error'] = str(exc)
        row['inputs'] = inputs
        row['target_files_opened'] = ['Pre'] if round_id in ('R1', 'R4') else ['Pre', 'Post_assisted_target']
        row['prediction_seconds'] = time.perf_counter() - start
        if arrays:
            ap = out / 'prediction.npz'
            np.savez_compressed(ap, **arrays)
            row['arrays'] = dict(path=str(ap), sha256=sha(ap))
        path = out / 'plan.json'
        write(path, row)
        records.append(dict(case=cid, path=str(path), sha256=sha(path)))
        state(f'{round_id}_predicted_{cid}', 'finish and hash all predictions before evaluation')
        if len(records) % 10 == 0:
            print(round_id, 'PREDICTED', len(records), row['status'], flush=True)
    write(frozen, dict(round=round_id, frozen_utc=now(), prereg_sha256=sha(ROOT / f'PREREG_{round_id}.json'), predictor_code_sha256=code, predictions=records, wall_seconds=time.perf_counter() - tick, cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime + resource.getrusage(resource.RUSAGE_SELF).ru_stime, peak_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    (ROOT / f'FROZEN_PREDICTIONS_{round_id}.sha256').write_text(sha(frozen) + '\n')
    if round_id == 'R1':
        write(ROOT / 'FROZEN_PREDICTIONS.json', json.loads(frozen.read_text()))
    state(round_id + '_predictions_frozen', 'evaluate independent virtual target with frozen prediction hashes')
    print(round_id, 'FROZEN', len(records), flush=True)

def evaluate(round_id, shard=None):
    frozen = json.loads((ROOT / f'FROZEN_PREDICTIONS_{round_id}.json').read_text())
    rows = []
    tick = time.perf_counter()
    for (record_index, rr) in enumerate(frozen['predictions']):
        if shard is not None and record_index % shard[1] != shard[0]:
            continue
        if sha(rr['path']) != rr['sha256']:
            raise ValueError('frozen plan hash mismatch')
        pred = json.loads(Path(rr['path']).read_text())
        cid = rr['case']
        start = time.perf_counter()
        saved = ROOT / round_id / 'cases' / f'{cid}.json'
        if saved.exists():
            existing = json.loads(saved.read_text())
            if existing.get('prediction_manifest') == rr and ('raw' not in existing or sha(existing['raw']['path']) == existing['raw']['sha256']):
                rows.append(existing)
                continue
        row = dict(case=cid, HCL=MAN['records'][cid]['HCL classification'], quality=MAN['records'][cid]['Data Quality'], prediction_status=pred['status'], prediction_manifest=rr, methods={}, solver_disagreement_mm=pred.get('solver_disagreement_mm'))
        if pred['status'] == 'PREDICTED':
            assert sha(pred['arrays']['path']) == pred['arrays']['sha256']
            (pre, a, _) = load(cid, 'Pre')
            (post, pa, _) = load(cid, 'Post')
            (original, oa, _) = load(cid, 'Original')
            (T, fit) = rigid_fit(a[:4000], pa[:30000])
            q = transform(post, np.linalg.inv(T))
            (O, ofit) = rigid_fit(oa[:4000], a[:30000])
            orig = transform(original, O)
            qd = dist(q, pre)
            roi = q[qd > 3.0]
            base_distance = qd[qd > 3.0]
            row['roi_n'] = len(roi)
            row['registration'] = dict(T_Pre_to_Post=T.tolist(), anchor_trimmed=fit)
            row['original_preservation'] = stats(dist(orig, pre))
            row['original_semantics'] = 'defective; measures preparation/repositioning, not intact completion'
            row['condyle_landmark_error_mm'] = None
            row['condyle_status'] = 'UNKNOWN_NO_LABELS; retained Pre is unchanged by construction'
            row['occlusal_landmark_error_mm'] = None
            raw = {'expert_roi': roi.astype('f4'), 'post_aligned': q.astype('f4')}
            with np.load(pred['arrays']['path']) as z:
                curve = z['curve'].astype(float)
                try:
                    if round_id == 'R4':
                        from centroid_rail import rail as centroid
                        (erail, erm) = centroid(q, pre, a, pred['rail']['projection_direction'])
                    else:
                        (erail, erm) = extract_rail(q, pre, a)
                    row['expert_arc_length_mm'] = erm['arc_length_mm']
                except ValueError:
                    row['expert_arc_length_mm'] = None
                for method in ('dp', 'uniform', 'mirror', 'osteoopt', 'shape_only'):
                    key = method + '_points'
                    if key not in z:
                        continue
                    graft = z[key].astype(float)
                    candidate = graft[dist(graft, pre) > 3.0]
                    fd = np.minimum(base_distance, dist(roi, graft))
                    rd = dist(candidate, q)
                    value = max(np.quantile(fd, 0.95), np.quantile(rd, 0.95)) if len(fd) and len(rd) else None
                    dense = z[method + '_points_dense'].astype(float) if method + '_points_dense' in z else graft
                    if method + '_points_dense' in z:
                        fdd = np.minimum(base_distance, dist(roi, dense))
                        dense_c = dense[dist(dense, pre) > 3.0]
                        rdd = dist(dense_c, q)
                    else:
                        fdd = fd
                        rdd = rd
                    vd = max(np.quantile(fdd, 0.95), np.quantile(rdd, 0.95)) if len(fdd) and len(rdd) else None
                    detail = pred.get('grafts', {}).get(method, {})
                    length = sum(detail.get('lengths_mm', [])) or None
                    row['methods'][method] = dict(primary_p95_mm=float(value) if value is not None else None, forward=stats(fd), reverse=stats(rd), density_probe_change_mm=abs(float(value) - float(vd)) if value is not None and vd is not None else None, length_mm=length, arc_length_error_mm=abs(length - row['expert_arc_length_mm']) if length and row['expert_arc_length_mm'] else None, count=detail.get('count'), bend_angles_deg=detail.get('bend_angles_deg'), cut_plane_residual_max_mm=detail.get('cut_plane_residual_max_mm'), endpoint_envelope_gap_mm=np.maximum(dist(np.array(detail['nodes_mm'])[[0, -1]], pre) - detail['radius_mm'], 0).tolist() if detail else None)
                    if detail:
                        from constrained import cut_feasible
                        row['methods'][method]['non_crossing_miters'] = cut_feasible(np.array(detail['nodes_mm']), detail['radius_mm'])[0]
                    raw[method + '_forward'] = fd.astype('f4')
                    raw[method + '_reverse'] = rd.astype('f4')
            rawpath = DATA / round_id / cid / 'evaluation_raw.npz'
            np.savez_compressed(rawpath, **raw)
            row['raw'] = dict(path=str(rawpath), sha256=sha(rawpath))
            row['numerically_available'] = len(roi) >= 200 and row['methods'].get('dp', {}).get('primary_p95_mm') is not None
        else:
            row['numerically_available'] = False
            row['unknown_reason'] = pred['error']
        row['evaluation_seconds'] = time.perf_counter() - start
        rows.append(row)
        write(ROOT / round_id / 'cases' / f'{cid}.json', row)
        state(round_id + '_evaluated_' + cid, 'aggregate frozen gates; preserve unknown cases')
        if len(rows) % 20 == 0:
            print(round_id, 'EVALUATED', len(rows), flush=True)
    if shard is not None:
        print(round_id, 'SHARD_COMPLETE', shard, 'cases', len(rows), flush=True)
        return
    available = [r for r in rows if r['numerically_available']]
    methods = {}
    for name in ('dp', 'uniform', 'mirror', 'osteoopt', 'shape_only'):
        vals = [r['methods'][name]['primary_p95_mm'] for r in rows if name in r['methods'] and r['methods'][name]['primary_p95_mm'] is not None]
        methods[name] = dict(n=len(vals), diagnostic_median_p95_mm=float(np.median(vals)) if vals else None, complete_cohort_median_p95_mm=float(np.median(vals)) if len(vals) == len(rows) else None)
    gains = [r['methods']['uniform']['primary_p95_mm'] - r['methods']['dp']['primary_p95_mm'] for r in available if 'uniform' in r['methods'] and r['methods']['uniform']['primary_p95_mm'] is not None]
    reg = json.loads((ROOT / f'PREREG_{round_id}.json').read_text())
    gates = dict(all_cases_available=len(available) == len(rows), median_primary_p95_mm_max=methods['dp']['complete_cohort_median_p95_mm'] is not None and methods['dp']['complete_cohort_median_p95_mm'] <= 3.0, paired_gain_vs_uniform_mirror_mm_min=len(gains) == len(rows) and np.median(gains) >= 0.5, solver_matches_exhaustive=all((r.get('solver_disagreement_mm') is not None and r['solver_disagreement_mm'] <= 1e-08 for r in rows)))
    extra = {}
    if round_id == 'R2':
        del gates['paired_gain_vs_uniform_mirror_mm_min']
        r1 = {r['case']: r for r in json.loads((ROOT / 'RESULTS_R1.json').read_text())['cases']}
        paired = [r1[r['case']]['methods']['dp']['primary_p95_mm'] - r['methods']['dp']['primary_p95_mm'] for r in available if r1[r['case']]['numerically_available']]
        expected_pairs = len(reg['additional_contract']['paired_R1_case_ids'])
        gates['paired_gain_vs_R1_mm_min'] = len(paired) == expected_pairs and np.median(paired) >= 1.0
        gates['all_miters_non_crossing'] = all((r['methods'].get('dp', {}).get('non_crossing_miters', False) for r in rows))
        extra = dict(paired_R1_diagnostic_n=len(paired), paired_gain_vs_R1_diagnostic_mm=float(np.median(paired)) if paired else None)
    if round_id == 'R3':
        penalties = [r['methods']['dp']['primary_p95_mm'] - r['methods']['shape_only']['primary_p95_mm'] for r in available if r['methods'].get('shape_only', {}).get('primary_p95_mm') is not None]
        plans = [json.loads(Path(rr['path']).read_text()) for rr in frozen['predictions']]
        ncontain = sum((p.get('all_envelopes_contained', False) for p in plans))
        gates = dict(all_cases_available=len(available) == len(rows), median_primary_p95_mm_max=methods['dp']['complete_cohort_median_p95_mm'] is not None and methods['dp']['complete_cohort_median_p95_mm'] <= 3.0, fraction_all_segment_envelopes_contained_min=ncontain / len(rows) >= 0.8, median_surface_penalty_vs_shape_only_mm_max=len(penalties) == len(rows) and np.median(penalties) <= 2.0, all_miters_non_crossing=all((r['methods'].get('dp', {}).get('non_crossing_miters', False) for r in rows)))
        extra = dict(n_all_envelopes_contained=ncontain, fraction_all_envelopes_contained=ncontain / len(rows), surface_penalty_vs_shape_only_diagnostic_mm=float(np.median(penalties)) if penalties else None, equal_information_controls_agree=sum((p.get('equal_information_control_agrees', False) for p in plans)), clinical_implant_feasibility='UNKNOWN: donor/cortex/occlusion and other anatomy absent; containment is a scenario calculation')
    if round_id == 'R4':
        mirrors = {r['case']: r for r in json.loads((ROOT / 'RESULTS_MIRROR_ALL.json').read_text())['cases']}
        paired = [mirrors[r['case']]['primary_p95_mm'] - r['methods']['dp']['primary_p95_mm'] for r in available if mirrors[r['case']]['primary_p95_mm'] is not None]
        gates = dict(all_cases_available=len(available) == len(rows), median_primary_p95_mm_max=methods['dp']['complete_cohort_median_p95_mm'] is not None and methods['dp']['complete_cohort_median_p95_mm'] <= 3.0, paired_gain_vs_raw_mirror_mm_min=len(paired) == len(rows) and np.median(paired) >= 1.0, solver_matches_exhaustive=all((r.get('solver_disagreement_mm') is not None and r['solver_disagreement_mm'] <= 1e-08 for r in rows)), all_miters_non_crossing=all((r['methods'].get('dp', {}).get('non_crossing_miters', False) for r in rows)))
        extra = dict(paired_gain_vs_raw_mirror_diagnostic_mm=float(np.median(paired)) if paired else None, blind_information_contract='Pre only; no Post until all plans frozen')
    out = dict(round=round_id, prereg_sha256=sha(ROOT / f'PREREG_{round_id}.json'), prediction_manifest_sha256=sha(ROOT / f'FROZEN_PREDICTIONS_{round_id}.json'), n_designated=147, n_complete_triplets=len(rows), n_available=len(available), methods=methods, gates={k: bool(v) for (k, v) in gates.items()}, paired_gain_uniform_diagnostic_mm=float(np.median(gains)) if gains else None, outcome='PASS_PENDING_REVIEW' if all(gates.values()) else 'FAIL_FROZEN_GATES', external_referent=dict(reg['external_referent'], refutes_us=not all(gates.values())), scientific_admission=False, review_state='PENDING_INDEPENDENT_REVIEW', evaluation_wall_seconds=time.perf_counter() - tick, peak_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, cases=rows, **extra)
    write(ROOT / f'RESULTS_{round_id}.json', out)
    state(round_id + '_decided', 'write handoff and change load-bearing information/representation', out['outcome'])
    print(json.dumps({k: v for (k, v) in out.items() if k != 'cases'}, indent=2), flush=True)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['predict', 'evaluate'])
    ap.add_argument('--round', default='R1')
    ap.add_argument('--shard')
    args = ap.parse_args()
    if args.stage == 'predict':
        predict(args.round)
    else:
        evaluate(args.round, tuple(map(int, args.shard.split('/'))) if args.shard else None)
