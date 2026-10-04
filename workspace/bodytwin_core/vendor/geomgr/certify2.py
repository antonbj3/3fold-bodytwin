"""geomgr.certify2: validity certificate v2 per geometry-changing operation (BT-N5 + unit contract), N7c.

certificate_v2(V_src, F_src, V_def, F_def, *, units, coordinate_units, attachments, joints, local=None)
  gates (all must pass):
    unit_contract          declared units == coordinate units == 'mm'                       (BT-B22 contract)
    v1                     BodyTwin geometry_certificate v1 (closed, winding, positive volume, conforming,
                           no self-intersection) -- LOCAL ONLY (trimesh). Without it (cloud) the numeric
                           stand-in `closed_oriented_positive` is used and self-intersection is UNKNOWN.
    no_negative_jacobian   BT-N5 per-element material-frame Jacobian (rotation removed) > -1e-9
    no_fold                BT-N5 / X1b rotation-invariant fold count == 0
    no_collapse            min deformed/source face-area ratio > 1e-6 (BT-B22 threshold)
    attachments_on_surface declared attachment points within 1e-6 mm of the deformed surface
    joints_on_head         declared joint centre within 10 mm of the sphere fitted to its head points
  (functions: geomgr/vendor/geometry_validity_v2_n5.py, copied from results/BT-N5/geometry_validity_v2.py)

instance_certificate(inst_src, inst_def, op, local=None) builds the declarations from the registry:
  attachments = every landmark/attachment entity (barycentric on the surface by construction, so this gate
  tests the geometry the ids point to), joint = HJC (registry derived, 4 head landmarks) vs the 'head' region.
numpy (+ scipy) only unless local=True.
"""
from __future__ import annotations

import numpy as np

from .vendor import geometry_validity_v2_n5 as N5


def _edges_closed_oriented(F):
    F = np.asarray(F)
    e = np.r_[F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]
    und = np.sort(e, 1)
    _, cnt = np.unique(und, axis=0, return_counts=True)
    _, dcnt = np.unique(e, axis=0, return_counts=True)
    return bool(np.all(cnt == 2)), bool(np.all(dcnt == 1))


def _signed_volume(V, F):
    t = np.asarray(V, float)[np.asarray(F)]
    return float(np.einsum('ij,ij->i', t[:, 0], np.cross(t[:, 1], t[:, 2])).sum() / 6.0)


def _local_available():
    try:
        import trimesh  # noqa: F401
        from ._local import _path, BT_SRC
        _path(BT_SRC)
        import bodytwin.geometry.geometry_validity_v1  # noqa: F401
        return True
    except Exception:
        return False


def certificate_v2(V_src, F_src, V_def, F_def, *, units='mm', coordinate_units='mm', attachments=(), joints=(),
                   local=None, op='?'):
    V0, V1 = np.asarray(V_src, float), np.asarray(V_def, float)
    F0, F1 = np.asarray(F_src), np.asarray(F_def)
    local = _local_available() if local is None else local
    gates, checks = {}, {}
    gates['unit_contract'] = units == 'mm' and coordinate_units == 'mm'
    checks['units'] = dict(units=units, coordinate_units=coordinate_units)
    closed, oriented = _edges_closed_oriented(F1)
    vol = _signed_volume(V1, F1)
    checks['closed_oriented_positive'] = dict(closed=closed, oriented=oriented, signed_volume=vol)
    if local:
        from ._local import _path, BT_SRC
        _path(BT_SRC)
        v1 = N5.geometry_certificate_v1(V1, F1, tolerance_mm=1e-6, declared_units='mm', evaluate=True)
        checks['v1'] = N5._scalars(dict(v1))
        gates['v1'] = bool(v1['accepted'])
    else:
        checks['v1'] = 'UNKNOWN: not run (BodyTwin v1 needs trimesh; numeric stand-in used, self-intersection unchecked)'
        gates['closed_oriented_positive'] = closed and oriented and vol > 0
    same_topology = V0.shape == V1.shape and np.array_equal(F0, F1)
    if same_topology:
        jac = N5.per_face_jacobian(V0, V1, F0)
        fold = N5.rotation_invariant_fold(V0, V1, F0)
        checks['jacobian'] = N5._scalars({k: v for k, v in jac.items() if not isinstance(v, np.ndarray)})
        checks['fold'] = N5._scalars(fold)
        gates['no_negative_jacobian'] = jac['n_negative_jacobian_faces'] == 0
        gates['no_fold'] = fold['n_folded_faces'] == 0
        a0 = np.linalg.norm(np.cross(V0[F0[:, 1]] - V0[F0[:, 0]], V0[F0[:, 2]] - V0[F0[:, 0]]), axis=1)
        a1 = np.linalg.norm(np.cross(V1[F0[:, 1]] - V1[F0[:, 0]], V1[F0[:, 2]] - V1[F0[:, 0]]), axis=1)
        ok = a0 > 0
        mar = float(np.min(a1[ok] / a0[ok]))
        checks['min_area_ratio'] = mar
        gates['no_collapse'] = mar > 1e-6      # BT-B22 registered threshold (deformed/source face area)
    else:
        checks['jacobian'] = checks['fold'] = 'not applicable: topology changed (retopologize/subdivide)'
    att = [dict(id=a['id'], part='p', kind=a.get('kind', 'barycentric'), pos=np.asarray(a['pos'], float)) for a in attachments]
    if att:
        ar = N5.attachment_residuals({'p': (V1, F1)}, att)
        checks['attachments'] = dict(max_dist_mm=ar['attachment_max_dist_mm'], worst=ar['attachment_worst_id'], n=ar['n_gating'])
        gates['attachments_on_surface'] = ar['attachment_max_dist_mm'] <= N5.ATTACH_TOL_MM * (1e-3 if units == 'm' else 1.0)
    sph = [j for j in joints if 'head_points' in j]
    if sph:            # BT-N5 rule: declared centre within 10 mm of the sphere fitted to the (spherical) head
        jr = N5.joint_residuals([dict(id=j['id'], center=j['center'], head_points=j['head_points']) for j in sph])
        checks['joints'] = N5._scalars(jr)
        gates['joints_on_head'] = jr['joint_max_residual_mm'] <= N5.JOINT_TOL_MM * (1e-3 if units == 'm' else 1.0)
    decl = [j for j in joints if 'head_center' in j]
    if decl:           # BT-B22 contract for a declared (non-spherical) head: |c - head centre| <= 1.05 head radius
        rr = {j['id']: float(np.linalg.norm(np.asarray(j['center'], float) - np.asarray(j['head_center'], float))
                             / float(j['head_radius'])) for j in decl}
        checks['joints_declared_head'] = dict(max_ratio=max(rr.values()), per_joint=rr, limit=1.05)
        gates['joints_in_declared_head'] = max(rr.values()) <= 1.05
    gates = {k: bool(v) for k, v in gates.items()}
    acc = all(gates.values())
    return dict(op=op, version='certificate_v2@N7c-1 (BT-N5 v2 + unit contract)', accepted=acc, gates=gates,
                reasons=[k for k, v in gates.items() if not v], checks=checks, local_v1=bool(local))


def instance_declarations(inst, max_points=None):
    vals = inst.values()
    att = [dict(id=b, pos=v, kind='barycentric') for b, v in vals.items()
           if b in inst.reg.entities and inst.reg.entities[b].kind in ('landmark', 'attachment')]
    if max_points:
        att = att[:max_points]
    part = inst.reg.part
    joints = []
    hj = f'{part}/derived/HJC'
    head = f'{part}/region/head'
    if hj in vals and head in vals:
        joints.append(dict(id=hj, center=vals[hj], head_points=inst.V[vals[head]['indices']]))
    return att, joints


def instance_certificate(src, inst, op, local=None):
    """v2 certificate of an operation's output `inst` against its input `src` (same registry lineage)."""
    att, joints = instance_declarations(inst)
    return certificate_v2(src.V, src.F, inst.V, inst.F, units=inst.units, coordinate_units=inst.units,
                          attachments=att, joints=joints, local=local, op=op)


def b22_case_certificate(lib, i, local=None):
    """Certificate for BT-B22 library case i. lib = load_b22(dir)."""
    c = lib['cases'][i]
    base = lib['base']
    V = lib['arr']['vertices'][i]
    F = lib['arr']['faces'][i]
    att = [dict(id=f'a{k}', pos=p) for k, p in enumerate(lib['arr']['attachments'][i])]
    joints = [dict(id='joint', center=lib['arr']['joint_centers'][i], head_center=lib['arr']['head_centers'][i],
                   head_radius=lib['arr']['head_radii'][i])]
    return certificate_v2(base['vertices'], base['faces'], V, F, units=c['units'], coordinate_units=c['coordinate_units'],
                          attachments=att, joints=joints, local=local, op=c['id'])


def load_b22(d):
    import json
    from pathlib import Path
    d = Path(d)
    cases = json.loads((d / 'cases.json').read_text())['cases']
    b = np.load(d / 'base_mesh.npz')
    a = np.load(d / 'cases.npz')
    return dict(cases=cases, base={k: b[k] for k in b.files}, arr={k: a[k] for k in a.files})
