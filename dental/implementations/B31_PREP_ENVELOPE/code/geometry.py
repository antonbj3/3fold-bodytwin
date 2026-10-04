from common import *
import trimesh
from shapely.geometry import Polygon
from shapely.ops import triangulate

def source_mesh(t):
    (v, ii) = np.unique(t.reshape(-1, 3), axis=0, return_inverse=True)
    return trimesh.Trimesh(v, ii.reshape(-1, 3), process=False)

def boundary(m):
    e = np.vstack([m.faces[:, [0, 1]], m.faces[:, [1, 2]], m.faces[:, [2, 0]]])
    (u, c) = np.unique(np.sort(e, axis=1), axis=0, return_counts=True)
    if np.any(c > 2):
        raise ValueError('source nonmanifold edges')
    adj = {}
    for (a, b) in u[c == 1]:
        adj.setdefault(int(a), []).append(int(b))
        adj.setdefault(int(b), []).append(int(a))
    if not adj or any((len(q) != 2 for q in adj.values())):
        raise ValueError('source branching/absent boundary')
    start = min(adj)
    loop = [start]
    prev = -1
    cur = start
    while True:
        nxt = next((a for a in adj[cur] if a != prev))
        if nxt == start:
            break
        if nxt in loop:
            raise ValueError('source repeats boundary')
        loop.append(nxt)
        (prev, cur) = (cur, nxt)
    if len(loop) != len(adj):
        raise ValueError('source multiple loops')
    return np.array(loop)

def planar_disk(poly, z):
    tt = []
    for tr in triangulate(poly):
        if poly.covers(tr):
            tt.append(np.c_[np.array(tr.exterior.coords)[:3], np.full(3, z)])
    covered = sum((Polygon(t[:, :2]).area for t in tt))
    if abs(covered - poly.area) > 1e-09:
        raise ValueError('unresolved projected polygon triangulation')
    return np.array(tt)

def closed_source(t, base):
    m = source_mesh(t)
    ids = boundary(m)
    rim = m.vertices[ids].copy()
    p = Polygon(rim[:, :2])
    if not p.is_valid or p.area < 1e-06:
        raise ValueError('NON_SIMPLE_NATIVE_MARGIN_XY_PROJECTION')
    floor = rim.copy()
    floor[:, 2] = base
    skirts = []
    for i in range(len(rim)):
        j = (i + 1) % len(rim)
        skirts.extend([[rim[i], rim[j], floor[j]], [rim[i], floor[j], floor[i]]])
    disk = planar_disk(p, base)
    out = source_mesh(np.r_[t, np.array(skirts), disk])
    trimesh.repair.fix_normals(out)
    if out.volume < 0:
        out.invert()
    if not out.is_watertight or not out.is_winding_consistent:
        raise ValueError('VIRTUAL_SOURCE_CLOSURE_INVALID')
    return (out, rim, p)

def distances(m, pts):
    parts = []
    for p in np.array_split(pts, max(1, int(np.ceil(len(pts) / 128)))):
        if len(p):
            parts.append(trimesh.proximity.closest_point(m, p)[1])
    return np.concatenate(parts) if parts else np.array([])

def wall_bound(outer, inner, threshold=0.5, depth=6):
    mesh = source_mesh(inner)
    t = outer.copy()
    lower = np.inf
    sample = np.inf
    facets = 0
    unresolved = 0
    witness = None
    for level in range(depth + 1):
        c = t.mean(1)
        rad = np.linalg.norm(t - c[:, None], axis=2).max(1)
        d = distances(mesh, c)
        lb = d - rad
        jj = int(np.argmin(d))
        if d[jj] < sample:
            sample = float(d[jj])
            witness = {'point_mm': c[jj], 'facet_mm': t[jj], 'distance_mm': sample}
        refine = (lb < threshold) & (d >= threshold) & (level < depth)
        if np.any(~refine):
            lower = min(lower, float(lb[~refine].min()))
        if level == depth:
            unresolved = int(np.sum((lb < threshold) & (d >= threshold)))
        facets += len(t)
        if not refine.any():
            break
        (a, b, c) = t[refine].transpose(1, 0, 2)
        ab = (a + b) / 2
        bc = (b + c) / 2
        ca = (c + a) / 2
        t = np.concatenate([np.stack([a, ab, ca], 1), np.stack([ab, b, bc], 1), np.stack([ca, bc, c], 1), np.stack([ab, bc, ca], 1)])
    return {'sampled_min_mm': sample, 'geometric_lower_mm': lower, 'unresolved_facets': unresolved, 'facets_evaluated': facets, 'witness': witness, 'whole_original_exterior_covered': True, 'rigorous_roundoff_enclosure': 'MISSING', 'source_tissue_error': 'UNKNOWN', 'resolution': 'PER_POINT'}
