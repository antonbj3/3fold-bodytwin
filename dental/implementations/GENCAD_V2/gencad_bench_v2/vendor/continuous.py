"""Reused projected continuous-gap primitive; provenance in SOURCE_REUSE.json."""
import time
import numpy as np
from scipy.spatial import cKDTree

def cross(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]

def plane(tri):
    a = tri[:, 1] - tri[:, 0]
    b = tri[:, 2] - tri[:, 0]
    d = cross(a[:, :2], b[:, :2])
    x = (a[:, 2] * b[:, 1] - b[:, 2] * a[:, 1]) / d
    y = (a[:, 0] * b[:, 2] - b[:, 0] * a[:, 2]) / d
    return np.c_[x, y, tri[:, 0, 2] - x * tri[:, 0, 0] - y * tri[:, 0, 1]]

def inside(points, tri):
    edge = np.roll(tri, -1, axis=1) - tri
    side = cross(edge[:, None, :, :], points[:, :, None, :] - tri[:, None, :, :])
    orientation = np.sign(cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]))
    return np.all(side * orientation[:, None, None] >= -1e-11, axis=2)

def exact_batch(U, L, up, lp, ids):
    """All vertices of projected triangle intersections, then affine gaps."""
    u = U[ids[:, 0], :, :2]
    l = L[ids[:, 1], :, :2]
    eu = np.roll(u, -1, axis=1) - u
    el = np.roll(l, -1, axis=1) - l
    delta = l[:, None, :, :] - u[:, :, None, :]
    den = cross(eu[:, :, None, :], el[:, None, :, :])
    nonparallel = np.abs(den) > 1e-14
    t = np.divide(cross(delta, el[:, None, :, :]), den, out=np.zeros_like(den), where=nonparallel)
    s = np.divide(cross(delta, eu[:, :, None, :]), den, out=np.zeros_like(den), where=nonparallel)
    valid_edges = nonparallel & (t >= 0) & (t <= 1) & (s >= 0) & (s <= 1)
    intersections = (u[:, :, None, :] + t[:, :, :, None] * eu[:, :, None, :]).reshape(-1, 9, 2)
    points = np.concatenate([u, l, intersections], axis=1)
    valid = np.c_[inside(u, l), inside(l, u), valid_edges.reshape(-1, 9)]
    coeff = up[ids[:, 0]] - lp[ids[:, 1]]
    gap = points[:, :, 0] * coeff[:, None, 0] + points[:, :, 1] * coeff[:, None, 1] + coeff[:, None, 2]
    gap = np.where(valid, gap, np.inf)
    flat = int(np.argmin(gap))
    (row, column) = np.unravel_index(flat, gap.shape)
    return (float(gap[row, column]), ids[row].copy(), points[row, column].copy(), int(valid.any(axis=1).sum()))

def broad_phase(U, L):
    uc = U[:, :, :2].mean(axis=1)
    lc = L[:, :, :2].mean(axis=1)
    ur = np.linalg.norm(U[:, :, :2] - uc[:, None, :], axis=2).max(axis=1)
    lr = np.linalg.norm(L[:, :, :2] - lc[:, None, :], axis=2).max(axis=1)
    ub0 = U[:, :, :2].min(axis=1)
    ub1 = U[:, :, :2].max(axis=1)
    lb0 = L[:, :, :2].min(axis=1)
    lb1 = L[:, :, :2].max(axis=1)
    tree = cKDTree(lc)
    parts = []
    for first in range(0, len(U), 1024):
        neighbors = tree.query_ball_point(uc[first:first + 1024], ur[first:first + 1024] + lr.max(), workers=1)
        for (offset, js) in enumerate(neighbors):
            if not js:
                continue
            i = first + offset
            j = np.asarray(js, int)
            good = np.all((lb1[j] >= ub0[i]) & (lb0[j] <= ub1[i]), axis=1)
            if good.any():
                parts.append(np.c_[np.full(good.sum(), i, int), j[good]].astype('int32'))
    return np.concatenate(parts)

def extremum(U, L, pairs, bounded):
    start = time.perf_counter()
    up = plane(U)
    lp = plane(L)
    bounds = U[pairs[:, 0], :, 2].min(axis=1) - L[pairs[:, 1], :, 2].max(axis=1)
    if bounded:
        order = np.argsort(bounds, kind='stable')
        pairs = pairs[order]
        bounds = bounds[order]
    best = np.inf
    witness = None
    evaluated = 0
    overlaps = 0
    for first in range(0, len(pairs), 4096):
        ids = pairs[first:first + 4096]
        if bounded:
            ids = ids[bounds[first:first + 4096] < best]
            if len(ids) == 0:
                break
        (value, pair, xy, n) = exact_batch(U, L, up, lp, ids)
        evaluated += len(ids)
        overlaps += n
        if value < best:
            best = value
            witness = (pair, xy)
    return dict(minimum_gap_mm=best, maximum_penetration_mm=max(0, -best), query_s=time.perf_counter() - start, evaluated_pairs=evaluated, overlapping_pairs_evaluated=overlaps, witness_triangle_indices=witness[0].tolist(), witness_xy_mm=witness[1].tolist())
