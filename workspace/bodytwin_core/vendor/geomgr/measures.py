"""geomgr.measures: CCD / anteversion definitions on a femur in template correspondence (N7c, BT-B42 choice).

Two automatic definitions disagree (N7b: r = 0.45 for AV, 0.59 for CCD on 70 Imperial femurs):
  * `p1_smooth`  -- P1 measures() with N7b's continuous shaft axis (Gauss-weighted slice centroids, sigma 1.5 mm;
                    results/N7b/n7b_common.py:measures_smooth). Neck axis = head sphere centre (6 head landmarks)
                    -> centre of the 4 neck-isthmus landmarks (SNI/ANI/INI/PNI); AV = neck vs posterior condylar
                    line (PMC-PLC) in the plane normal to the mechanical axis, SIGNED (retroversion < 0).
  * `b42`        -- BT-B42 surface method (results/BT-B42/bt42_runner.py fit_shaft/fit_collum/angle_values, copied
                    unchanged except inputs): shaft = trimmed PCA of the 25-75 % diaphysis, collum = HJC -> shaft
                    centre-line point 1.25 head radii below the HJC (a proxy line, not the neck), AV from the
                    posterior condyle line in the plane normal to the shaft; B42 reports |AV| only. `b42_signed`
                    adds the sign (neck anterior of the condylar line = +).
  * `neck_pca`   -- independent reference for the neck axis only: first principal axis of the registry's
                    'neck' region vertices (|v - neck centre| <= 20 mm), oriented towards the head.

CHOSEN = 'p1_smooth' (decision rule and numbers in results/N7c/k5_definition.json; PREREG K5). The geometry
manager reports and edits CCD/AV in this definition, versioned as MEASURE_VERSION.

numpy only.
"""
from __future__ import annotations

import math

import numpy as np

from . import core as C

MEASURE_VERSION = 'ccd_av@p1_smooth-1'
CHOSEN = 'p1_smooth'
NECK4 = ('SNI', 'ANI', 'INI', 'PNI')


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else None


def p1_smooth(L, V, sigma=1.5):
    """N7b measures_smooth (copied, results/N7b/n7b_common.py): P1 measures with a continuous shaft axis."""
    m = C.measures(L, V)
    knee_c, mech, head_c = m['_knee_c'], m['_mech'], m['_head_c']
    proj = (V - knee_c) @ mech
    lo, hi = proj.min(), proj.max()
    cents = []
    for f in np.linspace(0.30, 0.70, 9):
        h = lo + f * (hi - lo)
        w = np.exp(-0.5 * ((proj - h) / sigma) ** 2)
        cents.append((w[:, None] * V).sum(0) / w.sum())
    Cc = np.array(cents)
    _, _, vt = np.linalg.svd(Cc - Cc.mean(0))
    shaft = C.unit(vt[0] * np.sign(vt[0] @ mech))
    neck_c = np.mean([L[k] for k in NECK4], 0)
    neck = C.unit(head_c - neck_c)
    m['CCD'] = float(np.degrees(np.arccos(np.clip(neck @ -shaft, -1, 1))))
    m['_shaft'] = shaft
    m['_neck'] = neck
    return m


# ---------------------------------------------------------------- BT-B42 (copied, results/BT-B42/bt42_runner.py)
def _pca_axis(points, hint):
    points = np.asarray(points, dtype=float)
    if len(points) < 3:
        return None
    centered = points - points.mean(axis=0)
    vals, vecs = np.linalg.eigh(centered.T @ centered / float(len(centered)))
    axis = vecs[:, int(np.argmax(vals))]
    if hint is not None and np.dot(axis, hint) < 0:
        axis = -axis
    return _unit(axis)


def _fit_shaft(points, hjc, condyle_mid):
    u0 = _unit(hjc - condyle_mid)
    rel = (points - condyle_mid) @ u0
    span = float(np.dot(hjc - condyle_mid, u0))
    base = (rel >= 0.25 * span) & (rel <= 0.75 * span)
    radial = np.linalg.norm((points - condyle_mid) - rel[:, None] * u0, axis=1)
    base &= radial <= 40.0
    axis = u0.copy()
    for _ in range(2):
        d = (points - condyle_mid) - ((points - condyle_mid) @ axis)[:, None] * axis
        r = np.linalg.norm(d, axis=1)
        selected = base & (r <= 40.0)
        med = float(np.median(r[selected]))
        mad = float(np.median(np.abs(r[selected] - med)))
        selected &= r <= med + 3.0 * 1.4826 * mad
        axis = _pca_axis(points[selected], u0)
    return axis


def _fit_collum(points, hjc, condyle_mid, shaft, head_radius):
    span = float(np.dot(hjc - condyle_mid, shaft))
    t = (points - condyle_mid) @ shaft
    offs = []
    for lo in np.arange(0.30 * span, 0.76 * span, 5.0):
        q = points[(t >= lo) & (t < lo + 5.0)]
        d = q - condyle_mid - ((q - condyle_mid) @ shaft)[:, None] * shaft
        sel = np.linalg.norm(d, axis=1) <= 40.0
        if int(sel.sum()) >= 20:
            offs.append(np.median(d[sel], axis=0))
    offset = np.median(np.asarray(offs), axis=0)
    endpoint = condyle_mid + (float(np.dot(hjc - condyle_mid, shaft)) - 1.25 * head_radius) * shaft + offset
    return _unit(endpoint - hjc)


def surface_points(V, F):
    """Vertices + face centroids (deterministic stand-in for B42's 100k area-stratified samples)."""
    V = np.asarray(V, float)
    return np.r_[V, V[np.asarray(F)].mean(1)]


def b42(L, V, F):
    """BT-B42 CCD and AV on the template mesh; returns |AV| (as B42) and the signed AV."""
    hjc, r_head = C.sphere_fit(np.array([L[k] for k in C.HEAD6]))
    # B42 used posterior-condyle sphere centres (B25) as the distal origin; on the template the nearest stand-in
    # inside the bone is the epicondyle midpoint (the posterior condyle SURFACE points lie > 40 mm off the shaft line).
    cmid = 0.5 * (L['MEC'] + L['LEC'])
    P = surface_points(V, F)
    shaft = _fit_shaft(P, hjc, cmid)
    collum = _fit_collum(P, hjc, cmid, shaft, r_head)
    cond = L['PLC'] - L['PMC']
    cond = _unit(cond - (cond @ shaft) * shaft)
    ccd = 180.0 - math.degrees(math.acos(float(np.clip(abs(collum @ shaft), 0, 1))))
    n_perp = _unit(-collum - (-collum @ shaft) * shaft)        # head-ward neck direction (collum points away from head)
    av_abs = math.degrees(math.acos(float(np.clip(abs(n_perp @ cond), 0, 1))))
    ant = cmid - 0.5 * (L['PMC'] + L['PLC'])
    ant = _unit(ant - (ant @ shaft) * shaft - (ant @ cond) * cond)
    av_signed = math.degrees(math.atan2(float(n_perp @ ant), abs(float(n_perp @ cond))))
    return dict(CCD=ccd, AV_abs=av_abs, AV=av_signed, _shaft=shaft, _neck=-collum, _hjc=hjc)


def neck_cylinder_axis(V, region_idx, seeds):
    """Least-squares cylinder axis of the neck region (X1b joints.cylinder_fit, copied logic), best of several
    seed axes (lowest RMS). Oriented towards the head by the caller. Independent of both CCD/AV definitions
    (the seeds only start the optimiser; the best-RMS solution is kept)."""
    from scipy.optimize import least_squares
    P = np.asarray(V, float)[np.asarray(region_idx)]
    p0 = P.mean(0)
    best = None
    for a0 in seeds:
        a0 = np.asarray(a0, float) / np.linalg.norm(a0)
        t = np.array([1.0, 0, 0]) if abs(a0[0]) < 0.9 else np.array([0, 1.0, 0])
        e1 = np.cross(a0, t)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(a0, e1)

        def unpack(x):
            c = p0 + x[0] * e1 + x[1] * e2
            a = a0 + x[2] * e1 + x[3] * e2
            return c, a / np.linalg.norm(a), x[4]

        def f(x):
            c, a, r = unpack(x)
            w = P - c
            return np.linalg.norm(w - np.outer(w @ a, a), axis=1) - r
        w = P - p0
        r0 = float(np.median(np.linalg.norm(w - np.outer(w @ a0, a0), axis=1)))
        sol = least_squares(f, np.r_[0, 0, 0, 0, r0], method='lm')
        c, a, r = unpack(sol.x)
        rms = float(np.sqrt(np.mean(sol.fun ** 2)))
        if best is None or rms < best['rms']:
            best = dict(point=c, axis=a, radius=float(r), rms=rms)
    return best


def neck_pca_axis(V, region_idx, head_c):
    P = np.asarray(V, float)[np.asarray(region_idx)]
    ax = _pca_axis(P, head_c - P.mean(0))
    return ax


def ccd_av(L, V, F=None, definition=CHOSEN):
    """CCD and AV (deg) in the named definition. L: landmark dict (registry landmark names)."""
    if definition == 'p1_smooth':
        m = p1_smooth(L, V)
        return dict(CCD=m['CCD'], AV=m['AV'], definition=MEASURE_VERSION)
    if definition == 'p1':
        m = C.measures(L, V)
        return dict(CCD=m['CCD'], AV=m['AV'], definition='ccd_av@p1-1')
    if definition == 'b42':
        m = b42(L, V, F)
        return dict(CCD=m['CCD'], AV=m['AV'], AV_abs=m['AV_abs'], definition='ccd_av@b42_signed-1')
    raise ValueError(definition)


def landmarks_of(values, part='femur_r'):
    """Registry values -> {short landmark name: xyz}."""
    pre = f'{part}/landmark/'
    return {k[len(pre):]: v for k, v in values.items() if k.startswith(pre)}
