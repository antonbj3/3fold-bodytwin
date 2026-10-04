"""Research IPR geometry. All lengths in mm. No calibrated clinical allowance."""
import hashlib
import json
from pathlib import Path
import numpy as np
import trimesh

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_case(path):
    with np.load(path, allow_pickle=False) as d:
        meta = json.loads(str(d['meta']))
        center = d['molar_centroid']
        m = d['mesial_centroid'] - center
        m[2] = 0
        m /= np.linalg.norm(m)
        z = np.array([0.0, 0.0, 1.0 if meta['jaw'] == 'lower' else -1.0])
        b = np.cross(z, m)
        R = np.stack([m, b, z])
        meshes = {}
        for name in ['Dentin', 'Enamel']:
            meshes['Tooth' if name == 'Dentin' else name] = trimesh.Trimesh((d[name + '_V'] - center) @ R.T, d[name + '_F'], process=False)
        return (meshes, meta, {'center_mm': center.tolist(), 'rotation': R.tolist()})

def ray_hits(mesh, yz, side):
    """Return unique x-coordinates sorted outside to inside for each yz ray."""
    yz = np.asarray(yz)
    start = float(np.max(np.abs(mesh.vertices[:, 0])) + 5)
    origins = np.column_stack([np.full(len(yz), side * start), yz])
    directions = np.tile([-side, 0.0, 0.0], (len(yz), 1))
    (loc, idx, _) = mesh.ray.intersects_location(origins, directions, multiple_hits=True)
    out = [[] for _ in yz]
    if len(loc) == 0:
        return out
    for (x, i) in zip(loc[:, 0], idx):
        out[i].append(float(x))
    return [np.unique(np.round(xs, 9))[np.argsort(-side * np.unique(np.round(xs, 9)))].tolist() for xs in out]

def capacity(thickness, boundary=0.3, residual=0.5, fraction=0.5, overcut=0.05):
    lower = np.maximum(0.0, np.asarray(thickness) - boundary)
    return np.maximum(0.0, np.minimum(fraction * lower, lower - residual) - overcut)

def control_bisection(thickness, boundary=0.3, residual=0.5, fraction=0.5, overcut=0.05):
    """Independent feasibility query, returning zero for an infeasible zero cut."""
    t = max(0.0, float(thickness) - boundary)
    (lo, hi) = (0.0, max(0.0, t))
    for _ in range(60):
        r = (lo + hi) / 2
        if t - r - overcut >= residual and r + overcut <= fraction * t:
            lo = r
        else:
            hi = r
    return lo

def measure_patch(meshes, side, height_fraction):
    (enamel, tooth) = (meshes['Enamel'], meshes['Tooth'])
    (z0, z1) = enamel.bounds[:, 2]
    zc = float(z0 + height_fraction * (z1 - z0))
    vs = tooth.vertices
    band = vs[np.abs(vs[:, 2] - zc) <= 0.4]
    sx = side * band[:, 0]
    yc = float(np.median(band[sx >= np.quantile(sx, 0.95), 1]))
    (yy, zz) = np.meshgrid(np.linspace(yc - 0.5, yc + 0.5, 9), np.linspace(zc - 0.4, zc + 0.4, 9))
    yz = np.column_stack([yy.ravel(), zz.ravel()])
    (eh, th) = (ray_hits(enamel, yz, side), ray_hits(tooth, yz, side))
    (t, gap, outer, inner) = ([], [], [], [])
    for (es, ts) in zip(eh, th):
        valid = len(es) >= 2 and len(ts) >= 2
        t.append(side * (es[0] - es[1]) if valid else np.nan)
        gap.append(side * (ts[0] - es[0]) if valid else np.nan)
        outer.append(es[0] if valid else np.nan)
        inner.append(es[1] if valid else np.nan)
    return {'side': 'mesial' if side == 1 else 'distal', 'height_fraction': height_fraction, 'height_center_mm': zc, 'buccolingual_center_mm': yc, 'yz_mm': yz, 'thickness_mm': np.array(t), 'surface_gap_mm': np.array(gap), 'enamel_outer_x_mm': np.array(outer), 'dej_x_mm': np.array(inner)}

def jsonable(obj):
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, dict):
        return {k: jsonable(v) for (k, v) in obj.items()}
    if isinstance(obj, (tuple, list)):
        return [jsonable(v) for v in obj]
    return obj

def write_json(path, obj):
    Path(path).write_text(json.dumps(jsonable(obj), indent=2, allow_nan=False) + '\n')
