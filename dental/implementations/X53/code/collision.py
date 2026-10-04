"""Continuous swept capsule vs retained triangle surface; millimetres.

Floating computations, not an interval-arithmetic certificate. Segments are
continuous, not discrete motion samples. AABB pruning uses the same inflation
as the clearance test. Closed solid and an exterior start are required.
"""
import numpy as np
import trimesh
EPS = 1e-07

def point_triangle(p, tri):
    if not len(tri):
        return np.empty(0)
    cp = trimesh.triangles.closest_point(tri, np.broadcast_to(p, (len(tri), 3)))
    return np.linalg.norm(cp - p, axis=1)

def segment_segment(a, b, c, d):
    """One segment against an array; convex constrained quadratic solution."""
    u = b - a
    v = d - c
    w = a - c
    aa = np.dot(u, u)
    bb = v @ u
    cc = np.sum(v * v, axis=1)
    dd = w @ u
    ee = np.sum(v * w, axis=1)
    den = aa * cc - bb * bb
    s = np.divide(bb * ee - cc * dd, den, out=np.zeros_like(cc), where=den > 1e-25)
    t = np.divide(aa * ee - bb * dd, den, out=np.zeros_like(cc), where=den > 1e-25)
    good = (den > 1e-25) & (s >= 0) & (s <= 1) & (t >= 0) & (t <= 1)
    dist = np.full(len(c), np.inf)
    dist[good] = np.linalg.norm(w[good] + s[good, None] * u - t[good, None] * v[good], axis=1)
    for ss in [0.0, 1.0]:
        tt = np.clip(np.divide(ee + ss * bb, cc, out=np.zeros_like(cc), where=cc > 1e-25), 0, 1)
        dist = np.minimum(dist, np.linalg.norm(w + ss * u - tt[:, None] * v, axis=1))
    for tt in [0.0, 1.0]:
        ss = np.clip((tt * bb - dd) / max(aa, 1e-25), 0, 1)
        dist = np.minimum(dist, np.linalg.norm(w + ss[:, None] * u - tt * v, axis=1))
    return dist

def segment_triangle(a, b, tri):
    if not len(tri):
        return np.empty(0)
    dist = np.minimum(point_triangle(a, tri), point_triangle(b, tri))
    for k in range(3):
        dist = np.minimum(dist, segment_segment(a, b, tri[:, k], tri[:, (k + 1) % 3]))
    e1 = tri[:, 1] - tri[:, 0]
    e2 = tri[:, 2] - tri[:, 0]
    u = b - a
    h = np.cross(np.broadcast_to(u, e2.shape), e2)
    det = np.sum(e1 * h, axis=1)
    inv = np.divide(1.0, det, out=np.zeros_like(det), where=abs(det) > 1e-20)
    s = a - tri[:, 0]
    uu = inv * np.sum(s * h, axis=1)
    q = np.cross(s, e1)
    vv = inv * (q @ u)
    tt = inv * np.sum(e2 * q, axis=1)
    hit = (abs(det) > 1e-20) & (uu >= -EPS) & (vv >= -EPS) & (uu + vv <= 1 + EPS) & (tt >= -EPS) & (tt <= 1 + EPS)
    dist[hit] = 0.0
    return dist

class Scene:

    def __init__(self, vertices, faces):
        self.v = np.asarray(vertices, float)
        self.f = np.asarray(faces, int)
        if self.v.ndim != 2 or self.v.shape[1] != 3 or (not np.isfinite(self.v).all()):
            raise ValueError('bad vertices')
        if self.f.ndim != 2 or self.f.shape[1] != 3 or (not len(self.f)) or (self.f.min() < 0) or (self.f.max() >= len(self.v)):
            raise ValueError('bad faces')
        self.tri = self.v[self.f]
        self.lo = self.tri.min(1)
        self.hi = self.tri.max(1)
        area = np.linalg.norm(np.cross(self.tri[:, 1] - self.tri[:, 0], self.tri[:, 2] - self.tri[:, 0]), axis=1) / 2
        if np.any(area < 1e-14):
            raise ValueError('degenerate triangles')
        self.area = area
        self.bounds = np.array([self.v.min(0), self.v.max(0)])

    def clearance(self, a, b, r):
        lo = np.minimum(a, b) - r - EPS
        hi = np.maximum(a, b) + r + EPS
        keep = np.all(self.hi >= lo, axis=1) & np.all(self.lo <= hi, axis=1)
        if not keep.any():
            return np.inf
        return float(segment_triangle(np.asarray(a), np.asarray(b), self.tri[keep]).min() - r)

    def ball_clearance(self, c, r):
        return self.clearance(c, c, r)

def directions(axes):
    A = np.radians(np.arange(0, 360, 15.0))
    Bs = [0.0] if axes == 4 else np.radians([0, -20, 20, -35, 35])
    rows = []
    for b in Bs:
        for a in A:
            rows.append([-np.sin(b), np.sin(a) * np.cos(b), -np.cos(a) * np.cos(b)])
    return np.asarray(rows)

def tool_capsules(c, d, tool, approach, holder=True):
    r = tool['diameter_mm'] / 2
    neck = tool['neck_reach_mm']
    gauge = tool['gauge_mm']
    parts = [('cutting_envelope', c, c + (neck + approach) * d, r)]
    parts.append(('shank', c + neck * d, c + (gauge + approach) * d, tool['shank_mm'] / 2))
    if holder:
        parts.append(('holder', c + gauge * d, c + (gauge + tool['holder_length_mm'] + approach) * d, tool['holder_diameter_mm'] / 2))
    return parts

def pose_clearance(scene, c, n, d, tool, approach, holder=True):
    if np.dot(n, d) < 0.001:
        return (False, 'noncutting_hemisphere', None)
    for (name, a, b, r) in tool_capsules(c, d, tool, approach, holder):
        val = scene.clearance(a, b, r)
        if val < -EPS:
            return (False, name, val)
    return (True, 'witness', None)

def search(scene, p, n, tool, axes, offsets=(0.0,), holder=True):
    ds = directions(axes)
    scores = ds @ n
    order = np.argsort(-scores, kind='stable')
    approach = float(np.linalg.norm(scene.bounds[1] - scene.bounds[0]) + 2 * tool['holder_diameter_mm'] + 2.0)
    reasons = {}
    r = tool['diameter_mm'] / 2
    for offset in offsets:
        c = p + (r + offset + 5 * EPS) * n
        if scene.ball_clearance(c, r) < -EPS:
            reasons['ball'] = reasons.get('ball', 0) + 1
            continue
        for k in order:
            if scores[k] < 0.001:
                continue
            (ok, why, _) = pose_clearance(scene, c, n, ds[k], tool, approach, holder)
            if ok:
                return dict(status='FOUND', offset_mm=offset, centre=c, direction=ds[k], pose_index=int(k), tested_rejections=reasons)
            reasons[why] = reasons.get(why, 0) + 1
    return dict(status='NOT_FOUND_ON_POSE_GRID', offset_mm=None, centre=None, direction=None, tested_rejections=reasons)

def first_union_ball_entry(points, normals, centres, radii, window=(-0.4, 0.8)):
    """Entry of the union on each signed normal ray, exact for spheres.

    Surface material is t<0, cavity air t>0. Entry>0 means stock; entry<0
    means overcut. This is a ball-witness set, not a complete CAM stock model.
    """
    if not len(centres):
        return np.full(len(points), np.nan)
    v = np.asarray(centres)[None, :, :] - np.asarray(points)[:, None, :]
    a = np.sum(v * np.asarray(normals)[:, None, :], axis=2)
    rho = np.asarray(radii)[None, :] ** 2 - (np.sum(v * v, axis=2) - a * a)
    valid = rho >= 0
    root = np.sqrt(np.maximum(rho, 0))
    lo = a - root
    hi = a + root
    valid &= (hi >= window[0]) & (lo <= window[1])
    lo = np.where(valid, np.maximum(lo, window[0]), np.inf)
    out = lo.min(1)
    out[~np.isfinite(out)] = np.nan
    return out
