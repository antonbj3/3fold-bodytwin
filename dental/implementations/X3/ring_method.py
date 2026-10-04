"""R2 adds an observed gingival hole rim to the positional information.
Only the surviving mask and a visible homolog are read. No removed-tooth
border/center/normal enters the fit. This depends on the synthetic deletion
leaving a cervical footprint; healed extraction sites need a different model.
"""
from crownbench import *

def edge_components(edges):
    adj = {}
    for (a, b) in edges:
        adj.setdefault(int(a), set()).add(int(b))
        adj.setdefault(int(b), set()).add(int(a))
    seen = set()
    out = []
    for root in adj:
        if root in seen:
            continue
        comp = set()
        stack = [root]
        while stack:
            a = stack.pop()
            if a in comp:
                continue
            comp.add(a)
            stack.extend(adj[a] - comp)
        seen |= comp
        out.append(np.array(sorted(comp), dtype=np.int32))
    return out

def hole_ring(v, f, labels, mirror_center, F):
    (edges, _) = boundary(v, f)
    candidates = []
    for ids in edge_components(edges):
        if len(ids) < 20:
            continue
        points = v[ids]
        local = (points - mirror_center) @ F
        span = np.ptp(local[:, :2], axis=0)
        gingiva = np.mean(labels[ids] == 0)
        center_dist = np.linalg.norm(local[:, :2].mean(0))
        radius = np.quantile(np.linalg.norm(local[:, :2], axis=1), 0.95)
        if gingiva < 0.5 or span.max() > 25 or center_dist > 6 or (radius > 12):
            continue
        candidates.append((float(center_dist), -len(ids), ids, {'boundary_vertices': len(ids), 'gingiva_fraction': float(gingiva), 'planar_centroid_distance_mm': float(center_dist), 'planar_span_mm': span.tolist(), 'radius_q95_mm': float(radius)}))
    if not candidates:
        raise ValueError('NO_VISIBLE_LOCAL_GINGIVAL_RING')
    (_, _, ids, detail) = min(candidates, key=lambda a: (a[0], a[1]))
    points = v[ids]
    if len(points) > 256:
        points = points[np.linspace(0, len(points) - 1, 256).astype(int)]
    return (points, detail)

def ring_fit(mv, mf, v, f, labels, k, F, teeth):
    c = mv.mean(0)
    (rim, detail) = hole_ring(v, f, labels, c, F)
    (edges, _) = boundary(mv, mf)
    parts = edge_components(edges)
    if not parts:
        raise ValueError('NO_SOURCE_CERVICAL_BOUNDARY')
    ids = max(parts, key=len)
    border = mv[ids]
    if len(border) > 256:
        border = border[np.linspace(0, len(border) - 1, 256).astype(int)]
    local = (mv - c) @ F
    bl = (border - c) @ F
    target = (rim - c) @ F
    tt = cKDTree(target)

    def warp(p, z):
        out = z.copy()
        out[:, 0] *= p[4]
        out[:, 1] *= p[5]
        co = np.cos(p[3])
        si = np.sin(p[3])
        xy = out[:, :2].copy()
        out[:, 0] = co * xy[:, 0] - si * xy[:, 1]
        out[:, 1] = si * xy[:, 0] + co * xy[:, 1]
        return out + p[:3]

    def residual(p):
        b = warp(p, bl)
        near = target[tt.query(b, workers=1)[1]]
        bt = cKDTree(b)
        back = b[bt.query(target, workers=1)[1]]
        return np.r_[(b - near).ravel() / 0.3, (target - back).ravel() / 0.3, 0.35 * p[:3], 0.35 * p[3] / 0.2, 0.35 * (p[4:] - 1) / 0.2]
    shift = np.clip(target.mean(0) - bl.mean(0), -2.5, 2.5)
    p0 = np.r_[shift, 0, 1, 1]
    fit = least_squares(residual, p0, bounds=([-3, -3, -3, -0.4, 0.7, 0.7], [3, 3, 3, 0.4, 1.3, 1.3]), loss='soft_l1', f_scale=1, max_nfev=80, diff_step=0.001, ftol=1e-06, xtol=1e-06, gtol=1e-06)
    new = warp(fit.x, local) @ F.T + c
    detail.update(params=fit.x.tolist(), success=bool(fit.success), nfev=int(fit.nfev), cost=float(fit.cost), source_border_vertices=len(ids), sampled_rim_vertices=len(rim), observed_rim_native=rim.tolist(), information='remaining-mask hole rim, not removed tooth boundary', synthetic_mask_footprint=True)
    return (new, mf, detail)
