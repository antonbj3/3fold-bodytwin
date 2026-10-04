"Generic cylinder fitting and surface-coordinate helpers. Wrap evaluation requires an independently supplied evaluator."
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Source: results/N7b/h_wrap.py WRAP (from LegTLEM/TLEM2.2/ModelParameters.any:244-266).
FIT_BAND_MM = 10.0
TRIM_MM = 3.0


def _perp_basis(a):
    """results/X1b/geomgr_ll2/joints.py _perp_basis (copied)."""
    a = a / np.linalg.norm(a)
    t = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(a, t)
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(a, e1)




def cyl_coords(X, cyl):
    rel = X - cyl['c']
    ax = rel @ cyl['a']
    r = np.linalg.norm(rel - np.outer(ax, cyl['a']), axis=1)
    return r, (X - cyl['s0']) @ cyl['a']


def circle_fit_on_axis(P, axis, c0, trim=2.0, iters=5):
    """results/X1b/geomgr_ll2/wrap.py circle_fit_on_axis (copied): circle in the plane normal to a fixed axis,
    Kasa start, geometric Gauss-Newton, trimming |res| > trim."""
    a = np.asarray(axis, float) / np.linalg.norm(axis)
    e1, e2 = _perp_basis(a)
    P = np.asarray(P, float)
    keep = np.ones(len(P), bool)
    for _ in range(iters + 1):
        X = P[keep] - c0
        u, v = X @ e1, X @ e2
        A = np.c_[2 * u, 2 * v, np.ones_like(u)]
        sol, *_ = np.linalg.lstsq(A, u ** 2 + v ** 2, rcond=None)
        cu, cv = sol[0], sol[1]
        R = np.sqrt(sol[2] + cu ** 2 + cv ** 2)
        for _ in range(20):
            du, dv = u - cu, v - cv
            r = np.hypot(du, dv)
            res = r - R
            J = np.c_[-du / r, -dv / r, -np.ones_like(r)]
            step, *_ = np.linalg.lstsq(J, -res, rcond=None)
            cu, cv, R = cu + step[0], cv + step[1], R + step[2]
            if np.linalg.norm(step) < 1e-12:
                break
        Xall = P - c0
        res_all = np.hypot(Xall @ e1 - cu, Xall @ e2 - cv) - R
        new_keep = np.abs(res_all) < trim
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    c = c0 + cu * e1 + cv * e2
    res = res_all[keep]
    return dict(c=c, a=a, R=float(R), rms=float(np.sqrt(np.mean(res ** 2))), n=int(keep.sum()))




def tps_weights(src_cp, targets):
    """Linear map W (len(targets), len(src_cp)) of scipy RBFInterpolator(kernel='linear', degree=1): the value of the
    interpolant of data d at `targets` is W @ d. Affine reproduction: W @ src_cp == targets, rows sum to 1."""
    from scipy.interpolate import RBFInterpolator
    n = len(src_cp)
    f = RBFInterpolator(np.asarray(src_cp, float), np.eye(n), kernel='linear', degree=1)
    return f(np.atleast_2d(np.asarray(targets, float)))


def fps(X, n, seed=0):
    """results/P1/p1_common.py fps (copied): farthest point sampling."""
    idx = [int(np.random.default_rng(seed).integers(len(X)))]
    d = np.linalg.norm(X - X[idx[0]], axis=1)
    for _ in range(n - 1):
        i = int(np.argmax(d))
        idx.append(i)
        d = np.minimum(d, np.linalg.norm(X - X[i], axis=1))
    return np.array(idx)


def construct_entities(part, src_V, points_mm, names, source_note, n_cp=400, seed=0):
    """Entities (kind 'construct') for points given on a source geometry in correspondence (src_V, same topology)."""
    from .identity import Entity
    cp = fps(np.asarray(src_V, float), n_cp, seed=seed)
    W = tps_weights(np.asarray(src_V, float)[cp], points_mm)
    out = []
    for nm, w in zip(names, W):
        w = w / w.sum()
        out.append(Entity(part, 'construct', nm, indices=[int(i) for i in cp], weights=[float(x) for x in w],
                          rule='tps_linear_deg1_400fps', source='atlas_internal', license='internal', note=source_note))
    return out, cp, W




# ---------------------------------------------------------------- geodesic wrap (X1b wrap.py, copied)
@dataclass
class Cylinder:
    name: str
    c: np.ndarray
    a: np.ndarray
    R: float

    def frame(self):
        a = np.asarray(self.a, float) / np.linalg.norm(self.a)
        e1, e2 = _perp_basis(a)
        return a, e1, e2

    def transformed(self, Rot, t):
        return Cylinder(self.name, Rot @ self.c + t, Rot @ self.a, self.R)


def wrap_length(P, Q, cyl, side=None):
    """results/X1b/geomgr_ll2/wrap.py wrap_length (copied): shortest path P->Q over an infinite cylinder if the
    straight line cuts it; `side` (3D) forces the arc whose midpoint points to `side` (BT-LV2 anatomical side)."""
    a, e1, e2 = cyl.frame()
    R = cyl.R
    p = np.asarray(P, float) - cyl.c
    q = np.asarray(Q, float) - cyl.c
    p2 = np.array([p @ e1, p @ e2])
    q2 = np.array([q @ e1, q @ e2])
    zp, zq = p @ a, q @ a
    rp, rq = np.linalg.norm(p2), np.linalg.norm(q2)
    straight = float(np.linalg.norm(q - p))
    if rp <= R or rq <= R:
        return straight, 'inside'
    d = q2 - p2
    dd = d @ d
    s = np.clip(-(p2 @ d) / dd, 0.0, 1.0) if dd > 0 else 0.0
    if np.linalg.norm(p2 + s * d) >= R:
        return straight, 'straight'
    tp, tq = np.sqrt(rp ** 2 - R ** 2), np.sqrt(rq ** 2 - R ** 2)
    ap, aq = np.arctan2(p2[1], p2[0]), np.arctan2(q2[1], q2[0])
    bp, bq = np.arccos(R / rp), np.arccos(R / rq)
    best = np.inf
    for sgn in (1.0, -1.0):
        phi_p = ap + sgn * bp
        phi_q = aq - sgn * bq
        dphi = (sgn * (phi_q - phi_p)) % (2 * np.pi)
        if side is not None:
            mid = phi_p + sgn * dphi / 2
            if (np.cos(mid) * e1 + np.sin(mid) * e2) @ np.asarray(side, float) <= 0:
                continue
        best = min(best, tp + tq + R * dphi)
    if not np.isfinite(best):
        return straight, 'straight'
    return float(np.sqrt(best ** 2 + (zq - zp) ** 2)), 'wrap'


def evaluate_wrap_entity(e, V, Q, scale=1.0):
    raise NotImplementedError("Wrapping requires an independently supplied surface evaluator.")
