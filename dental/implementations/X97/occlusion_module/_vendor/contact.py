import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from . import quality as Q

def measure(poly, xy):
    if len(poly) < 3:
        return (0.0, np.zeros(2))
    q = poly @ xy
    tri = np.stack([np.repeat(q[:1], len(q) - 2, axis=0), q[1:-1], q[2:]], axis=1)
    a = abs(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])) / 2
    return (float(a.sum()), np.sum(a[:, None] * tri.mean(1), axis=0))

def polygons(xy, faces, g, lo=0.0, hi=0.1):
    gg = g[faces]
    finite = np.isfinite(gg).all(1)
    active = finite & (gg.max(1) >= lo) & (gg.min(1) <= hi)
    poly = {}
    areas = []
    moments = []
    idx = []
    for fi in np.flatnonzero(active):
        ids = faces[fi]
        v = g[ids]
        p = np.eye(3)
        if np.isfinite(lo) and v.min() < lo:
            p = Q.clip(p, v - lo)
        if np.isfinite(hi) and v.max() > hi:
            p = Q.clip(p, hi - v)
        (a, m) = measure(p, xy[ids])
        if a > 1e-12:
            poly[int(fi)] = p
            areas.append(a)
            moments.append(m)
            idx.append(fi)
    return (poly, np.asarray(areas), np.asarray(moments).reshape(-1, 2), idx)

def spatial(xy, faces, g, lo=0.0, hi=0.1):
    (poly, a, m, idx) = polygons(xy, faces, g, lo, hi)
    links = []
    edge = {}
    for (n, fi) in enumerate(idx):
        ids = faces[fi]
        for (i, j) in zip(ids, np.roll(ids, -1)):
            (ga, gb) = (g[i], g[j])
            if ga == gb:
                positive = lo <= ga <= hi
            else:
                (x, y) = sorted([(lo - ga) / (gb - ga), (hi - ga) / (gb - ga)])
                positive = min(1.0, y) > max(0.0, x)
            if positive:
                k = tuple(sorted((int(i), int(j))))
                if k in edge:
                    links.append((n, edge[k]))
                else:
                    edge[k] = n
    if len(idx):
        ij = np.array(links, int).reshape(-1, 2)
        (n, labels) = connected_components(coo_matrix((np.ones(len(ij)), (ij[:, 0], ij[:, 1])), shape=(len(idx), len(idx))).tocsr(), directed=False)
        aa = np.bincount(labels, weights=a)
        mm = np.array([np.bincount(labels, weights=m[:, k]) for k in range(2)]).T
        regs = [dict(area_mm2=float(ar), centroid_xy_mm=mm[k] / ar) for (k, ar) in enumerate(aa)]
    else:
        n = 0
        regs = []
    return (dict(count=int(n), area_mm2=float(a.sum()), centroid_xy_mm=m.sum(0) / a.sum() if len(a) else None, regions=regs), poly)

def intersection_area(xy, faces, poly, g, lo=0.0, hi=0.1):
    total = 0.0
    for (fi, p) in poly.items():
        ids = faces[fi]
        v = g[ids]
        if np.isfinite(v).all():
            total += Q.area(Q.clip(Q.clip(p, v - lo), hi - v), xy[ids])
    return total

def compare(xy, faces, pg, rg):
    finite = np.isfinite(pg[faces]).all(1) & np.isfinite(rg[faces]).all(1)
    ff = faces[finite]
    support = float((abs(np.cross(xy[ff[:, 1]] - xy[ff[:, 0]], xy[ff[:, 2]] - xy[ff[:, 0]])) / 2).sum())
    reference_finite = np.isfinite(rg[faces]).all(1)
    rf = faces[reference_finite]
    reference_support = float((abs(np.cross(xy[rf[:, 1]] - xy[rf[:, 0]], xy[rf[:, 2]] - xy[rf[:, 0]])) / 2).sum())
    if support <= 0:
        return dict(status='UNKNOWN_NO_COMMON_SUPPORT', support_mm2=0.0, reference_support_mm2=reference_support, coverage=0.0)
    (a, poly) = spatial(xy, ff, pg)
    (b, _) = spatial(xy, ff, rg)
    intersection = intersection_area(xy, ff, poly, rg)
    (_, na, _, _) = polygons(xy, ff, pg, -np.inf, 0.0)
    flat_zero = np.all(pg[ff] == 0, axis=1)
    zero_area = float((abs(np.cross(xy[ff[flat_zero, 1]] - xy[ff[flat_zero, 0]], xy[ff[flat_zero, 2]] - xy[ff[flat_zero, 0]])) / 2).sum())
    negative = max(0.0, float(na.sum()) - zero_area)
    dist = float(np.linalg.norm(np.asarray(a['centroid_xy_mm']) - b['centroid_xy_mm'])) if a['count'] and b['count'] else 0.0 if a['count'] == b['count'] else None
    return dict(status='SCORED', predicted=a, reference=b, intersection_mm2=intersection, symdiff_mm2=max(0.0, a['area_mm2'] + b['area_mm2'] - 2 * intersection), count_error=abs(a['count'] - b['count']), centroid_distance_mm=dist, negative_gap_area_mm2=negative, support_mm2=support, reference_support_mm2=reference_support, coverage=support / reference_support if reference_support else 0.0)

def gates(m, pr):
    if m['status'] != 'SCORED':
        return {'support': False}
    tol = pr['floating_comparison_tolerance']
    return dict(support=m['coverage'] + tol >= pr['common_reference_support_fraction_min'], pattern=m['symdiff_mm2'] <= pr['symdiff_max_mm2'] + tol, count=m['count_error'] <= pr['component_count_error_max'], location=m['centroid_distance_mm'] is not None and m['centroid_distance_mm'] <= pr['centroid_distance_max_mm'] + tol, interference=m['negative_gap_area_mm2'] <= pr['negative_gap_area_max_mm2'] + tol)

def enclosure(xy, faces, pg, rg, delta=0.05):
    """Conservative set bounds for ALL offsets/errors in +/-delta on the PL field.
 Even permits independent errors; shared-offset worlds are a subset.
 Arithmetic is floating: set-theoretic enclosure, no rigorous rounding bound.
 """
    ff = faces[np.isfinite(pg[faces]).all(1) & np.isfinite(rg[faces]).all(1)]
    if not len(ff):
        return dict(status='UNKNOWN_NO_SUPPORT')
    (pd, pa, _, _) = polygons(xy, ff, pg, delta, 0.1 - delta)
    (rd, ra, _, _) = polygons(xy, ff, rg, delta, 0.1 - delta)
    (pp, pb, _, _) = polygons(xy, ff, pg, -delta, 0.1 + delta)
    (rp, rb, _, _) = polygons(xy, ff, rg, -delta, 0.1 + delta)
    lower = pa.sum() + ra.sum() - intersection_area(xy, ff, pd, rg, -delta, 0.1 + delta) - intersection_area(xy, ff, rd, pg, -delta, 0.1 + delta)
    support = float((abs(np.cross(xy[ff[:, 1]] - xy[ff[:, 0]], xy[ff[:, 2]] - xy[ff[:, 0]])) / 2).sum())
    upper = min(support, pb.sum() + rb.sum() - 2 * intersection_area(xy, ff, pd, rg, delta, 0.1 - delta))
    return dict(status='CONDITIONAL_SET_ENCLOSURE', symdiff_lower_mm2=max(0.0, float(lower)), symdiff_upper_mm2=max(0.0, float(upper)), predicted_area_bounds_mm2=[float(pa.sum()), float(pb.sum())], reference_area_bounds_mm2=[float(ra.sum()), float(rb.sum())], delta_mm=delta, roundoff_enclosure='MISSING', source_geometry_enclosure='MISSING')
