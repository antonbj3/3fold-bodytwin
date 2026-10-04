"""Threshold-directed 3D surface covering; published nearest-facet query reused."""
import time
import numpy as np
import trimesh
EPS = 1e-06

def subdivide(t):
    (a, b, c) = (t[:, 0], t[:, 1], t[:, 2])
    ab = (a + b) / 2
    bc = (b + c) / 2
    ca = (c + a) / 2
    return np.stack([np.stack([a, ab, ca], 1), np.stack([ab, b, bc], 1), np.stack([ca, bc, c], 1), np.stack([ab, bc, ca], 1)], 1).reshape(-1, 3, 3)

def nearest(mesh, points):
    chunks = []
    for i in range(0, len(points), 1024):
        chunks.append(trimesh.proximity.closest_point(mesh, points[i:i + 1024]))
    return tuple((np.concatenate([x[j] for x in chunks]) for j in range(3)))

def enclosure(source, target, minimum_required=None, maximum_allowed=None, signed=False, source_face_ids=None, budget=65536):
    """Numeric interval over the entire labelled source, conditional on nearest-query accuracy.

    Early failure returns a conservative interval and an actual violating point.
    Budget exhaustion returns UNKNOWN; every unrefined source cell is retained.
    """
    start = time.perf_counter()
    cells = source.copy()
    owners = np.arange(len(source)) if source_face_ids is None else np.asarray(source_face_ids)
    queries = 0
    level = 0

    def query(cells):
        centers = cells.mean(1)
        radius = np.linalg.norm(cells - centers[:, None, :], axis=2).max(1)
        (closest, dist, ids) = nearest(target, centers)
        if signed:
            dist = np.where(target.contains(centers), -dist, dist)
        return (centers, radius, closest, dist, ids)
    (cen, rad, cl, d, ids) = query(cells)
    queries += len(cells)
    while True:
        lows = d - rad - EPS
        highs = d + rad + EPS
        j = int(np.argmin(d))
        k = int(np.argmax(d))
        min_interval = [float(lows.min()), float(d[j] + EPS)]
        max_interval = [float(d[k] - EPS), float(highs.max())]

        def witness(i):
            return dict(source_point_mm=cen[i].tolist(), target_point_mm=cl[i].tolist(), distance_mm=float(d[i]), source_STL_face=int(owners[i]), target_region_face=int(ids[i]), source_subtriangle_mm=cells[i].tolist())
        lo_status = 'UNKNOWN' if minimum_required is None else 'PASS' if min_interval[0] >= minimum_required - EPS else 'FAIL' if min_interval[1] < minimum_required - EPS else 'UNKNOWN'
        hi_status = 'UNKNOWN' if maximum_allowed is None else 'PASS' if max_interval[1] <= maximum_allowed + EPS else 'FAIL' if max_interval[0] > maximum_allowed + EPS else 'UNKNOWN'
        need_lo = minimum_required is not None and lo_status == 'UNKNOWN'
        need_hi = maximum_allowed is not None and hi_status == 'UNKNOWN'
        ambiguous = (lows < minimum_required if need_lo else np.zeros(len(cells), bool)) | (highs > maximum_allowed if need_hi else np.zeros(len(cells), bool))
        if not ambiguous.any() or len(cells) + 3 * int(ambiguous.sum()) > budget or level >= 16:
            return dict(minimum_interval_mm=min_interval, maximum_interval_mm=max_interval, minimum_status=lo_status, maximum_status=hi_status, minimum_witness=witness(j), maximum_witness=witness(k), leaf_cells=len(cells), nearest_queries=queries, refinement_levels=level, budget_exhausted=bool(ambiguous.any()), query_wall_s=time.perf_counter() - start, numeric_error_allowance_mm=EPS, formal_outward_rounding=False, scope='Euclidean separation to labelled target facets over complete labelled source region; conditional on facet/region semantics and nearest-query accuracy')
        old = ~ambiguous
        new = subdivide(cells[ambiguous])
        newowners = np.repeat(owners[ambiguous], 4)
        (nc, nr, ncl, nd, ni) = query(new)
        queries += len(new)
        cells = np.concatenate([cells[old], new])
        owners = np.r_[owners[old], newowners]
        cen = np.r_[cen[old], nc]
        rad = np.r_[rad[old], nr]
        cl = np.r_[cl[old], ncl]
        d = np.r_[d[old], nd]
        ids = np.r_[ids[old], ni]
        level += 1

def planar_exterior(source, target, normal):
    n = np.asarray(normal, float)
    if n.shape != (3,) or not np.isfinite(n).all() or np.linalg.norm(n) < 1e-12:
        return False
    n /= np.linalg.norm(n)
    v = np.asarray(target.vertices)
    h = v @ n
    return bool(np.ptp(h) <= EPS and np.min(source.reshape(-1, 3) @ n - h.mean()) >= -EPS)

def tangent_ball_obstruction(crown, ids, r):
    if type(r) not in (int, float) or not np.isfinite(r) or r <= 0:
        raise ValueError('Ball radius must be finite, positive, in final-design mm')
    m = crown['mesh']
    ids = np.asarray(ids, int)
    sel = ids[::max(1, int(np.ceil(len(ids) / 4096)))]
    points = m.triangles_center[sel]
    normals = m.face_normals[sel]
    centers = points + r * normals
    (closest, d, fi) = nearest(m, centers)
    j = int(np.argmin(d - r))
    if d[j] < r - 1e-05:
        return dict(status='FAIL', reason='Tangent ideal ball intersects another crown boundary, so this face-interior point cannot be machined by that ball', radius_final_design_mm=r, collision_deficit_mm=float(r - d[j]), sampled_faces=len(sel), source_faces=len(ids), witness=dict(surface_point_mm=points[j].tolist(), surface_STL_face=int(sel[j]), normal_out_of_material=normals[j].tolist(), ball_center_mm=centers[j].tolist(), obstacle_point_mm=closest[j].tolist(), obstacle_STL_face=int(fi[j]), obstacle_distance_to_center_mm=float(d[j])), scope='Necessary local ideal-ball fit obstruction; no shaft, fixture, green-body conversion or global tool path certified')
    return dict(status='UNKNOWN', reason='No local tangent-ball obstruction sampled; global path/shaft/fixture remain unknown', radius_final_design_mm=r, sampled_faces=len(sel), source_faces=len(ids))

def insertion_obstruction(prep, crown, axis, travel):
    n = np.asarray(axis, float)
    if n.shape != (3,) or not np.isfinite(n).all() or np.linalg.norm(n) < 1e-12:
        raise ValueError('Invalid extraction axis')
    if type(travel) not in (float, int) or not np.isfinite(travel) or travel <= 0:
        raise ValueError('Extraction travel required in mm')
    n /= np.linalg.norm(n)
    points = prep['vertices']
    steps = np.linspace(0, travel, 41)
    for s in steps:
        relative = points - s * n
        inside = crown['mesh'].contains(relative)
        if not inside.any():
            continue
        ii = np.flatnonzero(inside)
        (closest, d, ids) = nearest(crown['mesh'], relative[ii])
        j = int(np.argmax(d))
        if d[j] > 1e-05:
            i = int(ii[j])
            return dict(status='FAIL', reason='Preparation point enters moving crown material along declared rigid straight extraction', travel_mm=float(s), witness=dict(preparation_vertex=i, preparation_point_mm=points[i].tolist(), point_in_crown_frame_mm=relative[i].tolist(), crown_translation_mm=(s * n).tolist(), nearest_crown_point_mm=closest[j].tolist(), interior_depth_mm=float(d[j]), nearest_crown_STL_face=int(ids[j])), scope='One actual point collision refutes this declared straight path; not a search over possible paths', sampled_steps=41)
    return dict(status='UNKNOWN', reason='No sampled preparation point entered material; continuous triangle/edge swept volume is not certified', sampled_steps=41, travel_mm=travel)
