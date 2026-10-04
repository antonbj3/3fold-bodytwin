"""Exact nearest-triangle search using centroid radius enclosures.
The bound follows triangle convexity; no k-nearest approximation survives into
an answer. Every triangle that could improve the initial distance is projected.
"""
import numpy as np
from scipy.spatial import cKDTree
from trimesh.triangles import closest_point
_CACHE = {}

def build(v, f):
    key = (id(v), id(f))
    value = _CACHE.get(key)
    if value is not None and value['v'] is v and (value['f'] is f):
        return value
    t = v[f]
    c = t.mean(1)
    r = np.max(np.linalg.norm(t - c[:, None, :], axis=2), axis=1)
    bins = np.ceil(np.log2(np.maximum(r, 2 ** (-20)))).astype(int)
    groups = []
    for exp in np.unique(bins):
        ids = np.flatnonzero(bins == exp)
        groups.append((float(2.0 ** exp), ids, cKDTree(c[ids])))
    value = {'v': v, 'f': f, 'triangles': t, 'centers': c, 'radii': r, 'tree': cKDTree(c), 'groups': groups}
    if len(_CACHE) >= 12:
        _CACHE.pop(next(iter(_CACHE)))
    _CACHE[key] = value
    return value

def exact_distances(points, v, f):
    if len(f) == 0:
        raise ValueError('EMPTY_SURFACE')
    m = build(v, f)
    tri = m['triangles']
    tree = m['tree']
    centers = m['centers']
    radii = m['radii']
    out = []
    k = min(8, len(f))
    for start in range(0, len(points), 128):
        pts = points[start:start + 128]
        (_, ii) = tree.query(pts, k=k, workers=1)
        if k == 1:
            ii = ii[:, None]
        query = np.repeat(pts, k, axis=0)
        cc = closest_point(tri[ii.ravel()], query)
        best = np.linalg.norm(cc - query, axis=1).reshape(len(pts), k).min(1)
        pools = [[] for _ in pts]
        for (radius, ids, kt) in m['groups']:
            choices = kt.query_ball_point(pts, best + radius + 1e-12, workers=1)
            for (i, items) in enumerate(choices):
                if items:
                    pools[i].extend(ids[items].tolist())
        sizes = np.array([len(a) for a in pools])
        pointidx = np.repeat(np.arange(len(pts)), sizes)
        idx = np.concatenate([np.asarray(a, dtype=np.int64) for a in pools])
        bound = np.linalg.norm(centers[idx] - pts[pointidx], axis=1) <= best[pointidx] + radii[idx] + 1e-12
        idx = idx[bound]
        pointidx = pointidx[bound]
        for a in range(0, len(idx), 20000):
            inds = idx[a:a + 20000]
            pi = pointidx[a:a + 20000]
            qp = pts[pi]
            cp = closest_point(tri[inds], qp)
            d = np.linalg.norm(cp - qp, axis=1)
            np.minimum.at(best, pi, d)
        out.append(best)
    return np.concatenate(out)
