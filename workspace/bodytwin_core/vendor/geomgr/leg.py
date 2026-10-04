"""geomgr.leg: moment arms through the individual's own femur and wrap cylinders (N7c).

Model (copied from results/X1b/geomgr_ll2/paths.py LegModel + leg2.py LegModel2, positions only; X1b docstrings
apply): pelvis = world (TLEM 2.0 pelvis, generic), ball hip, hinge knee (TLEM axis carried to the individual),
patella hinge with constant patellar-ligament length, ankle locked at neutral (foot points expressed in the tibia
frame). Muscle paths = TLEM 2.0 Table A3 elements; every FEMUR-fixed point comes from the individual instance:
  * origins/insertions on the femur -> registry attachments `femur_r/attachment/<Element>|<Origin|Insertion>` (N7a);
  * femoral via points, knee-axis and patella-axis points -> registry constructs `femur_r/construct/TLEM.*`
    (TLEM table points carried from the TLEM femur in the population by the same TPS correspondence as wraps);
  * hip centre = sphere fitted to the registry 'head' region; femur orientation from the registry ISB frame
    (R0_f maps the individual's ISB axes onto the TLEM femur's ISB axes);
  * wrapping: rectus femoris over `femur_r/wrap/RectusVastii` (anterior side), hamstrings over
    `femur_r/wrap/Hamstring` and gastrocnemius over `femur_r/wrap/Gastro` (posterior side) -- BT-LV2 forced
    anatomical side, X1b geodesic cylinder wrap on the femoral->distal crossing segment.
Pelvis, tibia, patella and foot are TLEM 2.0 (not individualised: UNKNOWN for the individual).
Moment arm r = -dL/dq (central differences, h = 1e-4 rad; positive = pulls towards +q; knee flexion positive,
hip adduction positive, hip flexion positive). Group value = PCSA-weighted mean of elements (TLEM Table A7).

API:  leg = LegFrozen.load(path)                       # TLEM tables frozen to JSON (tlem_leg.json)
      arms = moment_arms(values, leg)                  # values = Instance.values() in mm -> {quantity: mm}
      band = arm_band(gm, res, leg, n=50, kappa=1.0)   # posterior draws -> mean/p05/p95/sd per quantity
numpy only.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import wrap as W

COORDS = ('hip_flex', 'hip_add', 'hip_rot', 'knee_flex')
GROUPS = {
    'gmed': [f'Gluteus Medius {p} {k}' for p in ('Anterior', 'Posterior') for k in range(1, 7)],
    'addmag': [f'Adductor Magnus Distal {k}' for k in range(1, 4)] + [f'Adductor Magnus Mid {k}' for k in range(1, 7)]
              + [f'Adductor Magnus Proximal {k}' for k in range(1, 5)],
    'hamstrings': ['Semimembranosus 1', 'Semimembranosus 2', 'Semimembranosus 3', 'Semitendinosus 1',
                   'Biceps Femoris Caput Longum 1'],
    'rectus': ['Rectus Femoris 1', 'Rectus Femoris 2'],
    'gastro': ['Gastrocnemius Lateralis 1', 'Gastrocnemius Medialis 1'],
}
WRAP_OF = {'rectus': ('RectusVastii', +1.0), 'hamstrings': ('Hamstring', -1.0), 'gastro': ('Gastro', -1.0)}
KNEE = (0.0, 30.0, 60.0, 90.0)
QUANTITIES = ([('gmed', 'hip_add', 0.0), ('addmag', 'hip_add', 0.0), ('addmag', 'hip_flex', 0.0)]
              + [('hamstrings', 'knee_flex', k) for k in KNEE] + [('hamstrings', 'hip_flex', 0.0)]
              + [('rectus', 'knee_flex', k) for k in KNEE] + [('rectus', 'hip_flex', 0.0)]
              + [('gastro', 'knee_flex', k) for k in KNEE])


def qname(g, c, k):
    return f'{g}:{c}@knee{int(k)}'


def rot(axis, ang):
    a = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


@dataclass
class LegModel:
    """results/X1b/geomgr_ll2/paths.py LegModel (copied; positions only) + leg2.py wrapping (copied)."""
    H: np.ndarray
    c_p: np.ndarray
    c_f: np.ndarray
    R0_f: np.ndarray
    knee_c: np.ndarray
    knee_a: np.ndarray
    tib_c: np.ndarray
    R0_t: np.ndarray
    pat_c: np.ndarray
    pat_a: np.ndarray
    patc_p: np.ndarray
    R0_pa: np.ndarray
    lig_pat: np.ndarray
    lig_tib: np.ndarray
    ankle_t: np.ndarray
    muscles: dict = field(default_factory=dict)
    wraps: dict = field(default_factory=dict)      # element -> wrap.Cylinder (femur coordinates)
    wrap_side: dict = field(default_factory=dict)  # element -> side direction (femur coordinates)
    s_k: float = 1.0
    s_p: float = 1.0
    L0: float = float('nan')
    h: float = 1e-4

    def __post_init__(self):
        self.knee_a = self.knee_a / np.linalg.norm(self.knee_a)
        self.pat_a = self.pat_a / np.linalg.norm(self.pat_a)
        ex_f = self.R0_f.T @ self.H[:, 0]
        a0 = self._tib_to_fem(self.ankle_t, 0.0)
        a1 = self._tib_to_fem(self.ankle_t, 0.2)
        if (a1 - a0) @ ex_f > 0:
            self.s_k = -1.0
        self.L0 = float(np.linalg.norm(self._pat_to_fem(self.lig_pat, 0.0) - self._tib_to_fem(self.lig_tib, 0.0)))
        self.s_p = self.s_k * float(np.sign(self.knee_a @ self.pat_a))

    def _tib_to_fem(self, x, th):
        return (rot(self.knee_a, self.s_k * th) @ (self.R0_t @ (np.asarray(x) - self.tib_c).T)).T + self.knee_c

    def _pat_to_fem(self, x, ph):
        return (rot(self.pat_a, self.s_p * ph) @ (self.R0_pa @ (np.asarray(x) - self.patc_p).T)).T + self.pat_c

    def R_hip(self, q):
        H = self.H
        return H @ rot([0, 0, 1], q['hip_flex']) @ rot([1, 0, 0], q['hip_add']) @ rot([0, 1, 0], q['hip_rot']) @ H.T

    def patella_angle(self, th, tol=1e-12):
        T = self._tib_to_fem(self.lig_tib, th)
        ph = 0.7 * th
        for _ in range(60):
            P = self._pat_to_fem(self.lig_pat, ph)
            w = P - T
            g = w @ w - self.L0 ** 2
            dP = self.s_p * np.cross(self.pat_a, P - self.pat_c)
            step = g / (2 * w @ dP)
            ph -= step
            if abs(step) < tol:
                break
        return ph

    def positions(self, pts, q):
        q = {k: float(q.get(k, 0.0)) for k in COORDS}
        th = q['knee_flex']
        ph = self.patella_angle(th)
        Rw = self.R_hip(q) @ self.R0_f
        X = np.zeros((len(pts), 3))
        for i, (seg, x) in enumerate(pts):
            x = np.asarray(x, float)
            if seg == 'pelvis':
                X[i] = x
                continue
            xf = x if seg == 'femur' else self._tib_to_fem(x, th) if seg == 'tibia' else self._pat_to_fem(x, ph)
            X[i] = Rw @ (xf - self.c_f) + self.c_p
        return X

    def crossing_segment(self, element):
        pts = self.muscles[element]
        for k in range(len(pts) - 1):
            if pts[k][0] in ('pelvis', 'femur') and pts[k + 1][0] in ('tibia', 'patella'):
                return k
        return None

    def length(self, element, q, status=None):
        pts = self.muscles[element]
        X = self.positions(pts, q)
        seg = np.linalg.norm(np.diff(X, axis=0), axis=1)
        cyl = self.wraps.get(element)
        if cyl is not None:
            k = self.crossing_segment(element)
            qq = {c: float(q.get(c, 0.0)) for c in COORDS}
            Rw = self.R_hip(qq) @ self.R0_f
            cw = cyl.transformed(Rw, self.c_p - Rw @ self.c_f)
            side = self.wrap_side.get(element)
            L, st = W.wrap_length(X[k], X[k + 1], cw, None if side is None else Rw @ np.asarray(side, float))
            seg[k] = L
            if status is not None:
                status.append(st)
        return float(seg.sum())

    def arm(self, element, q, coord):
        qp, qm = dict(q), dict(q)
        qp[coord] = q.get(coord, 0.0) + self.h
        qm[coord] = q.get(coord, 0.0) - self.h
        return -(self.length(element, qp) - self.length(element, qm)) / (2 * self.h)


class LegFrozen:
    """TLEM 2.0 tables frozen by n7c_freeze.py (tlem_leg.json) + TLEM femur ISB axes and the TLEM->native map."""

    def __init__(self, d):
        self.d = d
        self.J = {k: dict(center=np.asarray(v['center'], float), axis=None if v['axis'] is None else np.asarray(v['axis'], float))
                  for k, v in d['joints'].items()}
        self.muscles = d['muscles']
        self.lig = d['patellar_ligament']
        self.A_tlem = np.asarray(d['tlem_isb_axes'], float)        # rows, TLEM femur coordinates

    @classmethod
    def load(cls, path):
        return cls(json.loads(Path(path).read_text()))

    def foot_to_tibia(self, x):
        J = self.J
        return (np.asarray(x, float) + J['Subtalar|Talus']['center'] - J['Subtalar|Foot']['center']
                + J['Talocrural|Tibia']['center'] - J['Talocrural|Talus']['center'])


def tlem_construct_points(leg):
    """TLEM femur-fixed points that are not registry attachments (TLEM femur frame, mm) -> construct names."""
    J = leg.J
    pts = {'TLEM.knee.c': J['Knee|Femur']['center'], 'TLEM.knee.a30': J['Knee|Femur']['center'] + 30 * J['Knee|Femur']['axis'],
           'TLEM.patella.c': J['Patella|Femur']['center'],
           'TLEM.patella.a30': J['Patella|Femur']['center'] + 30 * J['Patella|Femur']['axis']}
    for el in sorted(set(sum(GROUPS.values(), []))):
        for k, (kind, seg, xyz) in enumerate(leg.muscles[el]['points']):
            if seg == 'femur' and kind == 'Via':
                pts[f'TLEM.via.{el.replace(" ", "_")}.{k}'] = np.asarray(xyz, float)
    return pts


def _sphere(P):
    A = np.c_[2 * P, np.ones(len(P))]
    c, *_ = np.linalg.lstsq(A, (P ** 2).sum(1), rcond=None)
    return c[:3]


def build_leg(values, leg, part='femur_r', V=None, wraps=True, femur_points=None):
    """LegModel for one instance. values = Instance.values() (mm). V = instance vertices (for the head sphere).
    femur_points: optional override {name: xyz} (used by the landmark-scaling baseline baseline)."""
    J = leg.J
    fp = femur_points
    if fp is None:
        g = lambda k, n: values[f'{part}/{k}/{n}']  # noqa: E731
        hj = _sphere(np.asarray(V, float)[values[f'{part}/region/head']['indices']])
        A_i = values[f'{part}/frame/ISB']['axes']
        R0_f = leg.A_tlem.T @ A_i
        kc, ka = g('construct', 'TLEM.knee.c'), g('construct', 'TLEM.knee.a30') - g('construct', 'TLEM.knee.c')
        pc, pa = g('construct', 'TLEM.patella.c'), g('construct', 'TLEM.patella.a30') - g('construct', 'TLEM.patella.c')
    else:
        g = None
        hj, R0_f = fp['_hjc'], np.eye(3)
        kc, ka = fp['TLEM.knee.c'], fp['TLEM.knee.a30'] - fp['TLEM.knee.c']
        pc, pa = fp['TLEM.patella.c'], fp['TLEM.patella.a30'] - fp['TLEM.patella.c']
    muscles, wr, side = {}, {}, {}
    for grp, els in GROUPS.items():
        for el in els:
            pts = []
            for k, (kind, seg, xyz) in enumerate(leg.muscles[el]['points']):
                if seg == 'femur':
                    key = (f'{el.replace(" ", "_")}|{kind}' if kind in ('Origin', 'Insertion')
                           else f'TLEM.via.{el.replace(" ", "_")}.{k}')
                    if fp is not None:
                        x = fp[key]
                    else:
                        x = g('attachment', key) if kind in ('Origin', 'Insertion') else g('construct', key)
                    pts.append(('femur', np.asarray(x, float)))
                elif seg == 'foot':
                    pts.append(('tibia', leg.foot_to_tibia(xyz)))
                else:
                    pts.append((seg, np.asarray(xyz, float)))
            muscles[el] = pts
            if wraps and grp in WRAP_OF:
                cname, sgn = WRAP_OF[grp]
                cy = values[f'{part}/wrap/{cname}'] if fp is None else fp[f'wrap:{cname}']
                wr[el] = W.Cylinder(cname, np.asarray(cy['c'], float), np.asarray(cy['a'], float), float(cy['R']))
                side[el] = R0_f.T @ np.array([sgn, 0.0, 0.0])      # anterior (+x TLEM) / posterior in femur coords
    RT = R0_f.T
    return LegModel(H=np.eye(3), c_p=J['Hip|Pelvis']['center'], c_f=np.asarray(hj, float), R0_f=R0_f,
                    knee_c=np.asarray(kc, float), knee_a=np.asarray(ka, float), tib_c=J['Knee|Tibia']['center'],
                    R0_t=RT, pat_c=np.asarray(pc, float), pat_a=np.asarray(pa, float),
                    patc_p=J['Patella|Patella']['center'], R0_pa=RT,
                    lig_pat=np.asarray([p[2] for p in leg.lig if p[1] == 'patella'][0], float),
                    lig_tib=np.asarray([p[2] for p in leg.lig if p[1] == 'tibia'][0], float),
                    ankle_t=J['Talocrural|Tibia']['center'], muscles=muscles, wraps=wr, wrap_side=side)


def moment_arms(values, leg, V=None, wraps=True, femur_points=None, per_element=False):
    """Group moment arms (mm) for QUANTITIES. Returns {qname: value} (+ per-element dict if asked)."""
    m = build_leg(values, leg, V=V, wraps=wraps, femur_points=femur_points)
    out, per = {}, {}
    for g, c, k in QUANTITIES:
        q = {'knee_flex': np.radians(k)}
        els = GROUPS[g]
        r = np.array([m.arm(e, q, c) for e in els])
        w = np.array([leg.muscles[e]['pcsa_cm2'] for e in els], float)
        w = np.where(np.isfinite(w), w, np.nanmean(w))
        out[qname(g, c, k)] = float(w @ r / w.sum())
        if per_element:
            per[qname(g, c, k)] = dict(zip(els, r.tolist()))
    return (out, per) if per_element else out


def posterior_draws(gm, res, n, seed=0, kappa=1.0):
    """Shape-frame posterior draws of the vertices (b ~ N(b_post, kappa^2 Sbb)); pose is irrelevant for arms."""
    from . import core as C
    rng = np.random.default_rng(seed)
    L = C._psd_sqrt(res.post.Sbb)
    for _ in range(n):
        yield gm.model.shape(res.post.b + kappa * (L @ rng.standard_normal(L.shape[1])))


def arm_band(gm, res, leg, n=50, seed=0, kappa=1.0):
    """Posterior band of the moment arms: arms at the posterior-mean shape + mean/sd/p05/p95 over n draws."""
    from .identity import evaluate
    reg = res.instance.reg
    point = moment_arms(evaluate(reg, res.X_shape), leg, V=res.X_shape)
    D = []
    for Vd in posterior_draws(gm, res, n, seed=seed, kappa=kappa):
        D.append(moment_arms(evaluate(reg, Vd), leg, V=Vd))
    keys = list(point)
    A = np.array([[d[k] for k in keys] for d in D])
    return {k: dict(point=point[k], mean=float(A[:, j].mean()), sd=float(A[:, j].std(ddof=1)),
                    p05=float(np.percentile(A[:, j], 5)), p95=float(np.percentile(A[:, j], 95)))
            for j, k in enumerate(keys)}, A
