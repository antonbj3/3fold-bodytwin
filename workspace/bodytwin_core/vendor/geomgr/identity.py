"""geomgr.identity: versioned, stable identities for landmarks, derived points, attachments, regions and frames.

An entity is `<part>/<kind>/<name>` with an integer version. Its definition (anchor = linear functional of the
part's vertices, or a rule over other entity ids, or a frozen vertex-id set) is hashed; the registry manifest
hash covers all entities. Operations on an Instance never touch the definitions except where the topology
changes (retopologize), which bumps the version and records the displacement. Identity checks count lost ids.

numpy only (runs on the cloud box).
"""
from __future__ import annotations

import copy
import hashlib
import json
from dataclasses import dataclass, field

import numpy as np

KINDS = ('landmark', 'derived', 'attachment', 'region', 'frame', 'construct', 'wrap')
# carried by dense-correspondence TPS, which is linear in the vertices); 'wrap' = wrapping cylinder evaluated
SOURCES = ('measured', 'atlas', 'derived', 'synthetic', 'atlas_internal')


def _h(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


@dataclass
class Entity:
    part: str
    kind: str
    name: str
    version: int = 1
    indices: list = field(default_factory=list)     # anchor vertex indices (landmark/attachment) or region vertex set
    weights: list = field(default_factory=list)     # anchor weights (sum 1)
    rule: str = ''                                  # human/machine rule (derived, frame)
    deps: list = field(default_factory=list)        # entity base ids used by a rule
    source: str = 'derived'
    license: str = 'open'
    note: str = ''

    @property
    def base(self):
        return f'{self.part}/{self.kind}/{self.name}'

    @property
    def id(self):
        return f'{self.base}@{self.version}'

    def definition(self):
        return dict(part=self.part, kind=self.kind, name=self.name, version=self.version,
                    indices=[int(i) for i in self.indices], weights=[float(w) for w in self.weights],
                    rule=self.rule, deps=list(self.deps), source=self.source, license=self.license)

    def hash(self):
        return _h(self.definition())

    def validate(self, n_vertices):
        if self.kind not in KINDS:
            raise ValueError(f'{self.id}: kind {self.kind}')
        if self.source not in SOURCES:
            raise ValueError(f'{self.id}: source {self.source}')
        if self.kind in ('landmark', 'attachment'):
            if not self.indices or len(self.indices) != len(self.weights):
                raise ValueError(f'{self.id}: anchor needs indices+weights')
            if abs(sum(self.weights) - 1.0) > 1e-9 or min(self.weights) < -1e-12:
                raise ValueError(f'{self.id}: weights must be a partition of unity')
            if max(self.indices) >= n_vertices:
                raise ValueError(f'{self.id}: index outside part')
        if self.kind == 'region':
            if not self.indices or max(self.indices) >= n_vertices:
                raise ValueError(f'{self.id}: region vertex set invalid')
        if self.kind == 'construct':
            if not self.indices or len(self.indices) != len(self.weights) or max(self.indices) >= n_vertices:
                raise ValueError(f'{self.id}: construct needs indices+weights inside the part')
            if abs(sum(self.weights) - 1.0) > 1e-9:
                raise ValueError(f'{self.id}: construct weights must sum to 1 (affine combination)')
        if self.kind == 'wrap':
            if len(self.deps) != 5 or not self.indices or len(self.indices) != len(self.weights):
                raise ValueError(f'{self.id}: wrap needs 5 construct deps and fit region with offsets')


class Registry:
    """Entities of one part (bone), bound to a topology (n_vertices, faces hash)."""

    def __init__(self, part, n_vertices, faces, entities=(), version_note=''):
        self.part = part
        self.n_vertices = int(n_vertices)
        self.faces_sha256 = hashlib.sha256(np.ascontiguousarray(np.asarray(faces, np.int64)).tobytes()).hexdigest()
        self.entities = {}
        self.history = []
        for e in entities:
            self.add(e)
        self.version_note = version_note

    def add(self, e):
        if e.part != self.part:
            raise ValueError(f'{e.id}: part {e.part} != registry {self.part}')
        e.validate(self.n_vertices)
        if e.base in self.entities:
            raise ValueError(f'duplicate {e.base}')
        for d in e.deps:
            if d not in self.entities:
                raise ValueError(f'{e.id}: unknown dependency {d}')
        self.entities[e.base] = e

    def ids(self):
        return sorted(e.id for e in self.entities.values())

    def bases(self):
        return sorted(self.entities)

    def manifest(self):
        return {'part': self.part, 'n_vertices': self.n_vertices, 'faces_sha256': self.faces_sha256,
                'entities': {e.id: e.hash() for e in sorted(self.entities.values(), key=lambda x: x.base)}}

    def manifest_hash(self):
        return _h(self.manifest())

    def by_kind(self, kind):
        return [e for e in self.entities.values() if e.kind == kind]

    def to_json(self):
        return dict(part=self.part, n_vertices=self.n_vertices, faces_sha256=self.faces_sha256,
                    version_note=self.version_note, history=self.history,
                    entities=[e.definition() | {'note': e.note} for e in sorted(self.entities.values(), key=lambda x: x.base)])

    @classmethod
    def from_json(cls, d, faces=None):
        r = cls.__new__(cls)
        r.part, r.n_vertices, r.faces_sha256 = d['part'], d['n_vertices'], d['faces_sha256']
        r.version_note, r.history, r.entities = d.get('version_note', ''), list(d.get('history', [])), {}
        for x in d['entities']:
            e = Entity(part=x['part'], kind=x['kind'], name=x['name'], version=x['version'], indices=x['indices'],
                       weights=x['weights'], rule=x['rule'], deps=x['deps'], source=x['source'],
                       license=x['license'], note=x.get('note', ''))
            r.entities[e.base] = e
        if faces is not None:
            fs = hashlib.sha256(np.ascontiguousarray(np.asarray(faces, np.int64)).tobytes()).hexdigest()
            if fs != r.faces_sha256:
                raise ValueError('registry bound to a different topology')
        return r


# ------------------------------------------------------------------ evaluation of entities on vertices
def _sphere_center(P):
    A = np.c_[2 * P, np.ones(len(P))]
    c, *_ = np.linalg.lstsq(A, (P ** 2).sum(1), rcond=None)
    return c[:3]


def _unit(v):
    return v / np.linalg.norm(v)


def evaluate(reg, V, scale=1.0):
    """All entities on vertices V: points (landmark, attachment, construct, derived), regions (index arrays +
    centroid), frames (origin, rows=axes), wraps (cylinder dict c, a, R, s0, L). `scale` = V units per mm
    (1 for mm, 1e-3 for m); only wraps need it (their offsets are stored in mm). Returns dict base -> value."""
    V = np.asarray(V, float)
    out = {}
    kinds = ('landmark', 'attachment', 'construct', 'region', 'derived', 'frame', 'wrap')
    order = sorted(reg.entities.values(), key=lambda e: (kinds.index(e.kind), e.base))
    for e in order:
        if e.kind in ('landmark', 'attachment', 'construct'):
            out[e.base] = np.asarray(e.weights) @ V[np.asarray(e.indices)]
        elif e.kind == 'region':
            idx = np.asarray(e.indices)
            out[e.base] = dict(indices=idx, centroid=V[idx].mean(0))
        elif e.kind == 'derived':
            P = [out[d] if not isinstance(out[d], dict) else out[d]['centroid'] for d in e.deps]
            if e.rule == 'sphere_center':
                out[e.base] = _sphere_center(np.array(P))
            elif e.rule == 'midpoint':
                out[e.base] = 0.5 * (P[0] + P[1])
            else:
                raise ValueError(e.rule)
    for e in order:
        if e.kind == 'frame':
            P = [out[d] if not isinstance(out[d], dict) else out[d]['centroid'] for d in e.deps]
            out[e.base] = frame_from_rule(e.rule, P)
    for e in order:
        if e.kind == 'wrap':
            from .wrap import evaluate_wrap_entity
            out[e.base] = evaluate_wrap_entity(e, V, [out[d] for d in e.deps], scale)
    return out


def frame_from_rule(rule, P):
    if rule in ('isb_femur_right', 'isb_femur_left'):
        # Wu 2002 femur: origin hip centre, y = hip centre - knee centre, z along MEC->LEC (right) ortho, x = y x z
        hc, kc, mec, lec = P
        y = _unit(hc - kc)
        zr = lec - mec if rule.endswith('right') else mec - lec
        z = _unit(zr - (zr @ y) * y)
        x = np.cross(y, z)
        return dict(origin=hc, axes=np.stack([x, y, z]))
    if rule == 'knee_p1':
        mec, lec, pmc, plc = P
        kc = 0.5 * (lec + mec)
        x = _unit(lec - mec)
        yy = kc - 0.5 * (plc + pmc)
        y = _unit(yy - (yy @ x) * x)
        return dict(origin=kc, axes=np.stack([x, y, np.cross(x, y)]))
    if rule in ('tibia_right', 'tibia_left'):
        mtp, ltp, mm = P
        o = 0.5 * (mtp + ltp)
        y = _unit(o - mm)
        zr = ltp - mtp if rule.endswith('right') else mtp - ltp
        z = _unit(zr - (zr @ y) * y)
        return dict(origin=o, axes=np.stack([np.cross(y, z), y, z]))
    raise ValueError(rule)


# ------------------------------------------------------------------ instance + operations
class Instance:
    """Geometry of one individual on a registry's topology plus an operation log."""

    def __init__(self, reg, V, F, *, units='mm', side='right', frame='model', subject='?', meta=None):
        self.reg = reg
        self.V = np.array(V, float)
        self.F = np.array(F, np.int64)
        self.units = units
        self.side = side
        self.frame = frame
        self.subject = subject
        self.meta = dict(meta or {})
        self.log = [dict(op='create', registry=reg.manifest_hash(), n_ids=len(reg.entities))]
        self.id_map = {b: b for b in reg.bases()}   # original base id -> current base id
        self.fields = {}   # N7c: per-vertex fields name -> dict(values (n,), kind 'length'|'scalar', units, source)

    def copy(self):
        c = copy.copy(self)
        c.V, c.F, c.log, c.id_map, c.meta = self.V.copy(), self.F.copy(), list(self.log), dict(self.id_map), dict(self.meta)
        c.fields = {k: dict(v, values=np.array(v['values'], float)) for k, v in getattr(self, 'fields', {}).items()}
        return c

    def unit_scale(self):
        return {'mm': 1.0, 'm': 1e3}[self.units]

    def values(self):
        return evaluate(self.reg, self.V, scale=1.0 / self.unit_scale())


def _log(inst, op, **kw):
    inst.log.append(dict(op=op, registry=inst.reg.manifest_hash(), n_ids=len(inst.reg.entities), **kw))


def op_rigid(inst, R, t, frame):
    out = inst.copy()
    out.V = inst.V @ np.asarray(R).T + np.asarray(t)
    out.frame = frame
    _log(out, 'rigid', frame=frame)
    return out


def op_units(inst, units):
    f = {('mm', 'm'): 1e-3, ('m', 'mm'): 1e3}.get((inst.units, units), 1.0 if inst.units == units else None)
    if f is None:
        raise ValueError(f'units {inst.units}->{units}')
    out = inst.copy()
    out.V = inst.V * f
    out.units = units
    for fd in out.fields.values():
        if fd.get('kind') == 'length':
            fd['values'] = fd['values'] * f
            fd['units'] = units
    _log(out, 'units', to=units, factor=f)
    return out


def rename_part(reg, new_part, frame_rule_map):
    r = Registry.__new__(Registry)
    r.part, r.n_vertices, r.faces_sha256, r.version_note = new_part, reg.n_vertices, reg.faces_sha256, reg.version_note
    r.history = list(reg.history) + [dict(op='rename_part', from_part=reg.part, to=new_part)]
    r.entities = {}
    for e in reg.entities.values():
        e2 = copy.deepcopy(e)
        e2.part = new_part
        e2.deps = [d.replace(reg.part + '/', new_part + '/', 1) for d in e.deps]
        e2.rule = frame_rule_map.get(e.rule, e.rule)
        r.entities[e2.base] = e2
    return r


def op_mirror(inst, new_part):
    """x -> -x, faces reversed (orientation kept outward). Entities keep name/kind/anchors; part renamed."""
    out = inst.copy()
    out.V = inst.V * np.array([-1.0, 1.0, 1.0])
    out.F = inst.F[:, ::-1].copy()
    rule_map = {'isb_femur_right': 'isb_femur_left', 'isb_femur_left': 'isb_femur_right',
                'tibia_right': 'tibia_left', 'tibia_left': 'tibia_right'}
    old_part = inst.reg.part
    out.reg = rename_part(inst.reg, new_part, rule_map)
    out.reg.faces_sha256 = hashlib.sha256(np.ascontiguousarray(out.F.astype(np.int64)).tobytes()).hexdigest()
    out.side = 'left' if inst.side == 'right' else 'right'
    out.id_map = {k: v.replace(old_part + '/', new_part + '/', 1) for k, v in inst.id_map.items()}
    _log(out, 'mirror', new_part=new_part)
    return out


def subdivide(V, F):
    """Midpoint subdivision (each triangle -> 4). Original vertices keep their indices."""
    V = np.asarray(V, float)
    F = np.asarray(F, np.int64)
    e = np.sort(np.r_[F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]], axis=1)
    uniq, inv = np.unique(e, axis=0, return_inverse=True)
    inv = inv.ravel()
    mid = len(V) + np.arange(len(uniq))
    V2 = np.r_[V, 0.5 * (V[uniq[:, 0]] + V[uniq[:, 1]])]
    nF = len(F)
    m01, m12, m20 = mid[inv[:nF]], mid[inv[nF:2 * nF]], mid[inv[2 * nF:]]
    F2 = np.concatenate([np.c_[F[:, 0], m01, m20], np.c_[F[:, 1], m12, m01], np.c_[F[:, 2], m20, m12],
                         np.c_[m01, m12, m20]], 0)
    return V2, F2, uniq


def op_subdivide(inst):
    """Nested remesh: anchors stay exact (original vertex indices unchanged); regions gain midpoints whose both
    edge ends are in the region (version bump, rule recorded)."""
    V2, F2, edges = subdivide(inst.V, inst.F)
    out = inst.copy()
    out.V, out.F = V2, F2
    reg = copy.deepcopy(inst.reg)
    reg.n_vertices = len(V2)
    reg.faces_sha256 = hashlib.sha256(np.ascontiguousarray(F2).tobytes()).hexdigest()
    for e in reg.entities.values():
        if e.kind == 'region':
            s = np.zeros(len(inst.V), bool)
            s[np.asarray(e.indices)] = True
            add = len(inst.V) + np.flatnonzero(s[edges[:, 0]] & s[edges[:, 1]])
            e.indices = sorted(list(map(int, e.indices)) + add.tolist())
            e.version += 1
            e.note = (e.note + ' | ' if e.note else '') + 'subdivision closure'
    reg.history.append(dict(op='subdivide', n_vertices=len(V2)))
    out.reg = reg
    for fd in out.fields.values():         # fields: midpoint average on new vertices
        v = fd['values']
        fd['values'] = np.r_[v, 0.5 * (v[edges[:, 0]] + v[edges[:, 1]])]
    _log(out, 'subdivide', n_vertices=len(V2))
    return out


def closest_on_mesh(V, F, P, k=24):
    """Closest point on triangle mesh (V,F) for points P: candidate faces by centroid kd-tree, exact
    point-triangle distance (Ericson). Returns face index, barycentric (3,), distance."""
    from scipy.spatial import cKDTree
    V = np.asarray(V, float)
    T = V[F]
    tree = cKDTree(T.mean(1))
    k = min(k, len(F))
    _, cand = tree.query(P, k=k)
    fi, bcs, ds = [], [], []
    for p, cs in zip(P, np.atleast_2d(cand)):
        best = (np.inf, None, None)
        for c in cs:
            q, bc = _closest_tri(p, *T[c])
            d = np.linalg.norm(p - q)
            if d < best[0]:
                best = (d, c, bc)
        ds.append(best[0])
        fi.append(best[1])
        bcs.append(best[2])
    return np.array(fi), np.array(bcs), np.array(ds)


def _closest_tri(p, a, b, c):
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = ab @ ap, ac @ ap
    if d1 <= 0 and d2 <= 0:
        return a, np.array([1., 0, 0])
    bp = p - b
    d3, d4 = ab @ bp, ac @ bp
    if d3 >= 0 and d4 <= d3:
        return b, np.array([0, 1., 0])
    vc = d1 * d4 - d3 * d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        v = d1 / (d1 - d3)
        return a + v * ab, np.array([1 - v, v, 0])
    cp = p - c
    d5, d6 = ab @ cp, ac @ cp
    if d6 >= 0 and d5 <= d6:
        return c, np.array([0, 0, 1.])
    vb = d5 * d2 - d1 * d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        w = d2 / (d2 - d6)
        return a + w * ac, np.array([1 - w, 0, w])
    va = d3 * d6 - d5 * d4
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return b + w * (c - b), np.array([0, 1 - w, w])
    den = 1.0 / (va + vb + vc)
    v, w = vb * den, vc * den
    return a + ab * v + ac * w, np.array([1 - v - w, v, w])


def op_retopologize(inst, V2, F2):
    """Non-nested remesh: anchors re-projected to the closest point on the new surface (version +1,
    displacement recorded), regions transferred by nearest new vertex -> old region label (version +1)."""
    from scipy.spatial import cKDTree
    old = inst.values()
    reg = copy.deepcopy(inst.reg)
    reg.n_vertices = len(V2)
    reg.faces_sha256 = hashlib.sha256(np.ascontiguousarray(np.asarray(F2, np.int64)).tobytes()).hexdigest()
    anchored = [e for e in reg.entities.values() if e.kind in ('landmark', 'attachment')]
    P = np.array([old[e.base] for e in anchored])
    fi, bc, d = closest_on_mesh(V2, F2, P)
    disp = {}
    for e, f, w, dd in zip(anchored, fi, bc, d):
        e.indices = [int(i) for i in np.asarray(F2)[f]]
        w = np.clip(w, 0, None)
        w = w / w.sum()
        e.weights = [float(x) for x in w]
        e.version += 1
        e.note = (e.note + ' | ' if e.note else '') + f'reprojected {dd:.3g} mm'
        disp[e.base] = float(dd)
    nearest = cKDTree(inst.V).query(V2)[1]
    back = cKDTree(np.asarray(V2, float)).query(inst.V)[1]        # old vertex -> nearest new vertex
    for e in reg.entities.values():
        if e.kind == 'construct':            # re-express on the new mesh: each old vertex -> closest point (bary)
            idx = np.asarray(e.indices)
            fi2, bc2, d2 = closest_on_mesh(V2, F2, inst.V[idx])
            acc = {}
            for w, f_, b_ in zip(e.weights, fi2, bc2):
                b_ = np.clip(b_, 0, None)
                b_ = b_ / b_.sum()
                for j, bj in zip(np.asarray(F2)[f_], b_):
                    acc[int(j)] = acc.get(int(j), 0.0) + float(w) * float(bj)
            e.indices, e.weights = list(acc), list(acc.values())
            e.version += 1
            e.note = (e.note + ' | ' if e.note else '') + f'reprojected max {float(d2.max()):.3g} mm'
            disp[e.base] = float(d2.max())
        elif e.kind == 'wrap':
            idx = np.asarray(e.indices)
            new = back[idx]
            _, first = np.unique(new, return_index=True)
            e.indices = new[first].astype(int).tolist()
            e.weights = np.asarray(e.weights, float)[first].tolist()
            e.version += 1
            e.note = (e.note + ' | ' if e.note else '') + 'fit region nearest-vertex transfer'
    for e in reg.entities.values():
        if e.kind == 'region':
            s = np.zeros(len(inst.V), bool)
            s[np.asarray(e.indices)] = True
            e.indices = np.flatnonzero(s[nearest]).astype(int).tolist()
            e.version += 1
            e.note = (e.note + ' | ' if e.note else '') + 'nearest-vertex transfer'
    reg.history.append(dict(op='retopologize', n_vertices=len(V2), max_anchor_displacement_mm=float(max(disp.values()))))
    out = inst.copy()
    out.V, out.F, out.reg = np.asarray(V2, float), np.asarray(F2, np.int64), reg
    for fd in out.fields.values():
        fd['values'] = fd['values'][nearest]
    _log(out, 'retopologize', max_anchor_displacement=float(max(disp.values())))
    return out, disp


def laplacian_smooth(V, F, iters=2, lam=0.3):
    import scipy.sparse as sp
    n = len(V)
    r = np.r_[F[:, 0], F[:, 1], F[:, 2], F[:, 1], F[:, 2], F[:, 0]]
    c = np.r_[F[:, 1], F[:, 2], F[:, 0], F[:, 0], F[:, 1], F[:, 2]]
    A = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(n, n))
    A.data[:] = 1.0
    deg = np.asarray(A.sum(1)).ravel()
    X = np.asarray(V, float).copy()
    for _ in range(iters):
        X = X + lam * (np.asarray(A @ X) / deg[:, None] - X)
    return X


# ------------------------------------------------------------------ identity check
def check_identity(before, after, *, exact=True, tol=1e-9, transform=None):
    """Compare two instances: every original base id must still exist (through id_map), definitions unchanged
    unless the op declares a version bump, and for exact ops |value(after) - transform(value(before))| <= tol.
    Returns dict(n_ids, lost, changed_definition, max_point_dev)."""
    lost = [ob for ob in before.id_map if after.id_map.get(ob) not in after.reg.entities]
    va, vb = after.values(), before.values()
    maxdev, changed = 0.0, []
    for ob, cur_b in before.id_map.items():
        cur_a = after.id_map.get(ob)
        if cur_a not in after.reg.entities or cur_b not in before.reg.entities:
            continue
        ea, eb = after.reg.entities[cur_a], before.reg.entities[cur_b]
        if ea.kind != eb.kind or ea.name != eb.name:
            changed.append(ob)
            continue
        if exact and ea.kind in ('landmark', 'attachment', 'derived', 'construct'):
            x = vb[cur_b]
            if transform is not None:
                x = transform(x[None])[0]
            maxdev = max(maxdev, float(np.linalg.norm(va[cur_a] - x)))
        if exact and ea.kind == 'wrap':      # cylinder: axis point distance, radius, axis direction (mm-equivalent)
            ca, cb = va[cur_a], vb[cur_b]
            c0, a0 = cb['c'], cb['a']
            if transform is not None:
                p0, p1 = transform(np.array([c0, c0 + a0]))
                s_ = np.linalg.norm(p1 - p0)
                c0, a0, R0 = p0, (p1 - p0) / s_, cb['R'] * s_
            else:
                R0 = cb['R']
            a1 = ca['a']
            off = (ca['c'] - c0) - ((ca['c'] - c0) @ a1) * a1          # centre distance to the other axis
            ang = np.linalg.norm(np.cross(a1, a0))
            maxdev = max(maxdev, float(np.linalg.norm(off)), abs(ca['R'] - R0), float(ang) * max(R0, 1.0))
    return dict(n_ids=len(before.id_map), lost=lost, n_lost=len(lost), changed=changed,
                max_point_dev=maxdev, ok=(not lost and not changed and (not exact or maxdev <= tol)))
