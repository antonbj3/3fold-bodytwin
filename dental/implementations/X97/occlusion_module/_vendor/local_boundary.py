import numpy as np
import time
from scipy.sparse import csr_matrix
from shapely import polygons, STRtree

def overlap_bary(q, upper, lower):
    xy = lower[:, :2]
    det = (xy[1, 0] - xy[0, 0]) * (xy[2, 1] - xy[0, 1]) - (xy[1, 1] - xy[0, 1]) * (xy[2, 0] - xy[0, 0])
    if det < 0:
        lower = lower[[0, 2, 1]]
        xy = lower[:, :2]
    p = np.eye(3)
    for (a, b) in zip(xy, np.roll(xy, -1, axis=0)):
        u = upper[:, :2] - a
        edge = b - a
        values = edge[0] * u[:, 1] - edge[1] * u[:, 0]
        p = q.Q.clip(p, values)
        if len(p) < 3:
            return (None, None)
    points = p @ upper[:, :2]
    (ar, _) = q.measure(p, upper[:, :2])
    if ar <= 1e-14:
        return (None, None)
    mat = np.c_[lower[:, :2], np.ones(3)]
    coefficients = np.linalg.solve(mat, lower[:, 2])
    zl = np.c_[points, np.ones(len(points))] @ coefficients
    p = np.maximum(p, 0)
    p /= p.sum(1)[:, None]
    return (p, p @ upper[:, 2] - zl)

def build(q, t, r):
    start = time.perf_counter()
    v = t['vertices']
    f = t['faces']
    tri = v[f]
    nz = np.cross(tri[:, 1, :2] - tri[:, 0, :2], tri[:, 2, :2] - tri[:, 0, :2])
    upids = np.flatnonzero((t['face_roles'] == 0) & (nz > 1e-12) & np.any(t['taper'][f] > 0, axis=1))
    downids = np.flatnonzero(nz < -1e-12)
    tree = STRtree(polygons(tri[downids, :, :2]))
    rows = []
    cols = []
    vals = []
    rhs = []
    pairs = 0
    vertices = 0
    crossings = []
    pruned = 0
    pruned_lb = np.inf
    fault = None
    for chunk in np.array_split(upids, max(1, int(np.ceil(len(upids) / 512)))):
        candidates = tree.query(polygons(tri[chunk, :, :2]))
        for (ui, di) in candidates.T:
            uf = int(chunk[ui])
            df = int(downids[di])
            (p, separation) = overlap_bary(q, tri[uf], tri[df])
            if p is None:
                continue
            pairs += 1
            vertices += len(p)
            if separation.min() < -1e-09 and separation.max() > 1e-09:
                crossings.append(dict(upper_face=uf, lower_face=df, separation_min_mm=float(separation.min()), separation_max_mm=float(separation.max())))
                continue
            if separation.max() <= 1e-09 and separation.min() < -1e-09:
                continue
            for (a, s) in zip(p, separation):
                s = float(s)
                maxmotion = r['relief_cap_mm'] * float(a @ t['taper'][f[uf]])
                if fault is None and a.max() > 0.1 and (s > 1e-06) and np.any(t['taper'][f[uf]] > 0):
                    fault = dict(upper_face=uf, lower_face=df, upper_vertex_ids=f[uf], barycentric=a, separation_mm=s)
                if s >= maxmotion:
                    pruned += 1
                    pruned_lb = min(pruned_lb, s - maxmotion)
                    continue
                rid = len(rhs)
                rows.extend([rid] * 3)
                cols.extend(f[uf])
                vals.extend(a)
                rhs.append(max(0.0, s))
    A = csr_matrix((vals, (rows, cols)), shape=(len(rhs), len(v)))
    b = np.array(rhs)
    stats = dict(upward_potentially_moved_faces=len(upids), fixed_downward_faces=len(downids), positive_area_overlap_pairs=pairs, overlap_polygon_vertices_tested=vertices, cap_implied_constraints=pruned, cap_implied_minimum_slack_lower_mm=None if not np.isfinite(pruned_lb) else float(pruned_lb), explicit_local_constraints=len(rhs), original_mixed_order_overlap_pairs=len(crossings), original_crossing_witnesses=crossings[:5], seconds=time.perf_counter() - start, full_original_self_intersection_proof='UNKNOWN', rigorous_machine_enclosure='MISSING')
    return (A, b, stats, fault)
