"""STL/region adapter. Coordinates stay in the original common frame (mm)."""
import hashlib
from pathlib import Path
import numpy as np
import trimesh
from .vendor.geometry import parse_stl
from .vendor import continuous as C
TOL = 1e-06

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def load(path, scale):
    p = Path(path)
    if p.stat().st_size > 150000000:
        raise ValueError('STL exceeds 150 MB input budget')
    (v, f) = parse_stl(p.read_bytes())
    if not len(f) or len(f) > 1000000 or (not np.isfinite(v).all()):
        raise ValueError('empty, nonfinite or oversized STL')
    v = v * scale
    m = trimesh.Trimesh(v, f, process=False)
    deg = np.linalg.norm(np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]]), axis=1) <= 1e-12
    sf = np.sort(f, axis=1)
    duplicates = len(f) - len(np.unique(sf, axis=0))
    return dict(path=str(p.resolve()), sha256=sha(p), vertices=v, faces=f, mesh=m, diagnostics=dict(faces=len(f), canonical_vertices=len(v), exact_duplicate_vertex_occurrences_merged=int(3 * len(f) - len(v)), degenerate_faces=int(deg.sum()), duplicate_faces=int(duplicates), watertight=bool(m.is_watertight), winding_consistent=bool(m.is_winding_consistent), extent_mm=np.ptp(v, axis=0).tolist(), repair='exact vertex dedup only; no filling, rounding, independent alignment or normal flip'))

def basis(axis):
    a = np.asarray(axis, float)
    if a.shape != (3,) or not np.isfinite(a).all() or np.linalg.norm(a) < 1e-12:
        raise ValueError('invalid axis')
    z = a / np.linalg.norm(a)
    ref = np.eye(3)[int(np.argmin(abs(z)))]
    x = np.cross(ref, z)
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    return np.array([x, y, z])

def faces(data, ids=None):
    if ids is None:
        ids = np.arange(len(data['faces']))
    a = np.asarray(ids)
    if a.ndim != 1 or a.dtype.kind not in 'iu' or len(a) == 0 or (len(np.unique(a)) != len(a)) or (a.min() < 0) or (a.max() >= len(data['faces'])):
        raise ValueError('invalid/empty/duplicate region face IDs')
    return (data['vertices'][data['faces'][a]], a)

def continuous(upper, lower, axis, bounded=True):
    R = basis(axis)
    U = upper @ R.T
    L = lower @ R.T
    um = np.abs(C.cross(U[:, 1, :2] - U[:, 0, :2], U[:, 2, :2] - U[:, 0, :2])) > 1e-12
    lm = np.abs(C.cross(L[:, 1, :2] - L[:, 0, :2], L[:, 2, :2] - L[:, 0, :2])) > 1e-12
    lost = int((~um).sum() + (~lm).sum())
    U = U[um]
    L = L[lm]
    result = dict(excluded_projected_faces=lost, arithmetic_tolerance_mm=TOL, formal_outward_rounding=False, scope='continuous signed projected facet gap along declared axis; not contact pressure or a general solid collision certificate')
    if not len(U) or not len(L):
        return dict(result, minimum_gap_mm=None, reason='No nondegenerate projection')
    try:
        pairs = C.broad_phase(U, L)
    except ValueError:
        return dict(result, minimum_gap_mm=None, reason='No overlapping projected bounding boxes')
    if not len(pairs):
        return dict(result, minimum_gap_mm=None, reason='No projected pairs')
    up = C.plane(U)
    lp = C.plane(L)
    if not any((np.isfinite(C.exact_batch(U, L, up, lp, pairs[k:k + 4096])[0]) for k in range(0, len(pairs), 4096))):
        return dict(result, minimum_gap_mm=None, reason='No actual projected intersections')
    out = C.extremum(U, L, pairs, bounded)
    (ui, li) = out['witness_triangle_indices']
    xy = np.asarray(out['witness_xy_mm'])
    uz = float(up[ui, :2] @ xy + up[ui, 2])
    lz = float(lp[li, :2] @ xy + lp[li, 2])
    out['witness'] = dict(upper_point_mm=(np.r_[xy, uz] @ R).tolist(), lower_point_mm=(np.r_[xy, lz] @ R).tolist(), upper_region_face=int(np.flatnonzero(um)[ui]), lower_region_face=int(np.flatnonzero(lm)[li]))
    return dict(result, **out, broad_phase_pairs=len(pairs))

def graph_region(tri):
    """Return XY-sorted unique graph nodes and topology, or refuse height folds."""
    (v, inv) = np.unique(tri.reshape(-1, 3), axis=0, return_inverse=True)
    (xy, ix) = np.unique(v[:, :2], axis=0, return_inverse=True)
    if len(xy) != len(v):
        raise ValueError('Region has multiple heights at same projected point')
    f = ix[inv.reshape(-1, 3)]
    z = np.empty(len(v))
    z[ix] = v[:, 2]
    if np.any(np.abs(C.cross(xy[f[:, 1]] - xy[f[:, 0]], xy[f[:, 2]] - xy[f[:, 0]])) < 1e-12):
        raise ValueError('Region has vertical/degenerate projections')
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    polygons = [Polygon(t) for t in xy[f]]
    area = sum((p.area for p in polygons))
    union = unary_union(polygons)
    if abs(area - union.area) > 1e-08 * max(1, area):
        raise ValueError('Overlapping projected region triangles')
    canon = np.sort(f, axis=1)
    order = np.lexsort(canon.T[::-1])
    canon = canon[order]
    return (xy, z, canon)

def graph_task(prep, inner, outer, axis, film_min, film_max, wall):
    R = basis(axis)
    (xi, zi, fi) = graph_region(inner @ R.T)
    (xo, zo, fo) = graph_region(outer @ R.T)
    if xi.shape != xo.shape or not np.allclose(xi, xo, rtol=0, atol=1e-09) or (not np.array_equal(fi, fo)):
        raise ValueError('Graph domains/topology are unmatched')
    P = prep @ R.T
    if np.ptp(P[:, :, 2]) > 1e-06:
        raise ValueError('Inherited preparation must be horizontal')
    from .vendor.geometry import height
    pz = height(P, xi)
    if not np.isfinite(pz).all():
        raise ValueError('Preparation does not cover intaglio graph')
    task = dict(task_id='STL', frame='declared-common-frame', xy=xi, faces=fi, preparation_z=pz, A=np.empty((0, len(xi))), obstacle_b=np.empty(0), family='single_crown', requirements=dict(film_min_mm=film_min, film_max_mm=film_max, wall_mm=wall))
    design = dict(task_id='STL', status='DESIGN', units='mm', frame=task['frame'], outer_vertices=np.c_[xi, zo], inner_vertices=np.c_[xi, zi], faces=fi)
    j = int(np.argmin(zo - zi))
    return (task, design, dict(inner_point_mm=(np.r_[xi[j], zi[j]] @ R).tolist(), outer_point_mm=(np.r_[xo[j], zo[j]] @ R).tolist()))
