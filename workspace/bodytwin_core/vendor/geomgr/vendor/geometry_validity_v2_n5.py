"BT-N5: extended validity certificate v2 for a *organized individual*.\n\nCopy/extension of the published `bodytwin.geometry.geometry_validity_v1` (source:\nreferences/current_bodytwin/src, byte-identical to inputs/baseline/src/geometry_validity_v1.py,\nsha256 47bba924…). The v1 check `geometry_certificate` is imported unchanged and supplemented with\nfour checks missing from v1 (see PREREG.md):\n\n  1. Jacobian sign per element (per triangle, in the source's own tangent basis).\n  2. Rotation-invariant fold measure (X1b's Kabsch method).\n  3. Attachment–surface distance (anchored attachments against the deformed part surface).\n  4. Joint-centre–sphere residual (declared joint centre against a sphere fit of the joint head).\n\nCopied functions with source references (no `import *`):\n  * `_patch_rotations`, `_fold_check`  <- results/X1b/geomgr_ll2/model.py (lines 128-185).\n  * `sphere_fit_algebraic`, `sphere_fit` <- results/X1b/geomgr_ll2/joints.py (lines 14-38).\n  * `distance_to_surface`              <- results/X1b/geomgr_ll2/landmarks.py (lines 92-96).\n\nThresholds (fixed in PREREG.md): JAC_TOL = -1e-9, FOLD = 0, ATTACH_TOL_MM = 1e-6, JOINT_TOL_MM = 10.0.\n"
from __future__ import annotations

import numpy as np

def geometry_certificate_v1(*a, **k):  # N7c: lazy (local-only: BodyTwin + trimesh)
    from bodytwin.geometry.geometry_validity_v1 import geometry_certificate
    return geometry_certificate(*a, **k)

__all__ = [
    'JAC_TOL', 'ATTACH_TOL_MM', 'JOINT_TOL_MM',
    'per_face_jacobian', 'rotation_invariant_fold', 'distance_to_surface',
    'attachment_residuals', 'sphere_fit', 'joint_residuals',
    'geometry_certificate_v2',
]

JAC_TOL = -1e-9
ATTACH_TOL_MM = 1e-6
JOINT_TOL_MM = 10.0
DEGEN_REL = 1e-3   # D1:s degeneracitetskriterium: hogd/langsta kant < DEGEN_REL


# ---------------------------------------------------------------- 1. Jacobian per element
def _patch_rotation_and_sign(V0, V1, F, rings=2):
    """As `_patch_rotations` but also returns the orientation sign d = sign(det(V U^T)) per face.

    d = +1: the local source->deformed map is (best) orientation-preserving (a rotation).
    d = -1: no rotation fits; the element is orientation-reversed (a reflection / inversion),
    which is the per-element Jacobian sign this certificate gates on.
    """
    import scipy.sparse as sp
    nV, nF = len(V0), len(F)
    r = np.r_[F[:, 0], F[:, 1], F[:, 2], F[:, 1], F[:, 2], F[:, 0]]
    c = np.r_[F[:, 1], F[:, 2], F[:, 0], F[:, 0], F[:, 1], F[:, 2]]
    A = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(nV, nV)) + sp.identity(nV, format='csr')
    A.data[:] = 1.0
    M = sp.csr_matrix((np.ones(3 * nF), (np.repeat(np.arange(nF), 3), F.ravel())), shape=(nF, nV))
    P = M
    for _ in range(rings):
        P = P @ A
        P.data[:] = 1.0
    n = np.asarray(P.sum(1)).ravel()
    s0 = P @ V0
    s1 = P @ V1
    outer = (V0[:, :, None] * V1[:, None, :]).reshape(nV, 9)
    H = (P @ outer).reshape(nF, 3, 3) - np.einsum('fi,fj->fij', s0, s1) / n[:, None, None]
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum('fji,fkj->fik', Vt, U)))
    D = np.ones((nF, 3))
    D[:, 2] = np.where(d == 0, 1.0, d)
    R = np.einsum('fji,fj,fkj->fik', Vt, D, U)
    return R, d


def per_face_jacobian(V0, V1, F, rings=2):
    """Per-element Jacobian of the source->deformed surface deformation.

    Two complementary per-element signs are combined:
      * **material-frame area Jacobian** det(J2): the deformed triangle is first rotated back into
        the source triangle's frame with the neighbourhood's best rigid rotation R (Kabsch, same
        R as the fold check), then both triangles are expressed in the source tangent basis
        (e1 = unit first edge, n = area normal, e2 = n x e1). det(J2) = signed 2-D area ratio of
        the rotation-removed map. This removes the global rigid part (TLEM->CT, or a deliberate
        rigid transform), so a valid morph gives det > 0 and a locally inverted/folded element
        gives det < 0. It is rotation-invariant because the same R removes any added rigid
        rotation of source and deformed together.

    Degenerate source triangles (height/longest-edge < DEGEN_REL) are excluded and counted.
    Returns dict with (n_faces,) arrays and the scalar gates used by v2.
    """
    V0 = np.asarray(V0, float)
    V1 = np.asarray(V1, float)
    F = np.asarray(F)
    t0, t1 = V0[F], V1[F]
    a0, b0, c0 = t0[:, 0], t0[:, 1], t0[:, 2]
    a1, b1, c1 = t1[:, 0], t1[:, 1], t1[:, 2]

    e = b0 - a0
    e_len = np.linalg.norm(e, axis=1)
    e1 = e / e_len[:, None]
    n = np.cross(b0 - a0, c0 - a0)
    area2 = np.linalg.norm(n, axis=1)
    n_hat = n / area2[:, None]
    e2 = np.cross(n_hat, e1)

    def coords(P, Q, R):
        u = np.einsum('ij,ij->i', P, e1)
        v = np.einsum('ij,ij->i', P, e2)
        w = np.einsum('ij,ij->i', Q, e1)
        z = np.einsum('ij,ij->i', Q, e2)
        return u, v, w, z

    u0, v0, w0, z0 = coords(b0 - a0, c0 - a0, None)
    src2 = u0 * z0 - w0 * v0          # > 0 by construction (right-handed basis)

    # neighbourhood best-fit rotation (same as the fold check): remove the local rigid part
    R, orientation = _patch_rotation_and_sign(V0, V1, F, rings=rings)
    d1 = np.einsum('fji,fj->fi', R, b1 - a1)   # R^T (b1-a1)
    d2 = np.einsum('fji,fj->fi', R, c1 - a1)   # R^T (c1-a1)
    u1, v1, w1, z1 = coords(d1, d2, None)
    def2 = u1 * z1 - w1 * v1
    with np.errstate(divide='ignore', invalid='ignore'):
        det = def2 / src2
    det = np.where(np.abs(src2) > 0.0, det, np.nan)

    side = [np.linalg.norm(b0 - a0, axis=1), np.linalg.norm(c0 - b0, axis=1),
            np.linalg.norm(a0 - c0, axis=1)]
    longest = np.max(np.stack(side), axis=0)
    height_rel = np.where(longest > 0, (area2 / longest) / longest, 0.0)
    degenerate = height_rel < DEGEN_REL

    ok = ~degenerate & np.isfinite(det)
    n_negative = int(np.sum(ok & (det < JAC_TOL)))
    min_det = float(np.min(det[ok])) if np.any(ok) else float('nan')
    min_face = int(np.argmin(np.where(ok, det, np.inf))) if np.any(ok) else -1
    return dict(jacobian_det=det, source_area2=src2, degenerate=degenerate,
                orientation_sign=orientation,
                n_degenerate=int(degenerate.sum()),
                n_negative_jacobian_faces=n_negative,
                n_negative_jacobian_elements=n_negative,
                n_orientation_flips=int(np.sum(~degenerate & (orientation < 0))),
                jacobian_min_det=min_det, jacobian_min_det_face=min_face,
                jacobian_median_det=float(np.median(det[ok])) if np.any(ok) else float('nan'))


# ---------------------------------------------------------------- 2. rotation-invariant fold check
def _patch_rotations(V0, V1, F, rings=2):
    "COPY with source reference: results/X1b/geomgr_ll2/model.py lines 128-157 (identical logic;\n    shares implementation with `_patch_rotation_and_sign`)."
    return _patch_rotation_and_sign(V0, V1, F, rings=rings)[0]


def rotation_invariant_fold(V0, V1, F, rings=2):
    "KOPIA/adaptation with source citation: results/X1b/geomgr_ll2/model.py rad 160-185.\n\n    Source normal rotated with the neighbourhood's best rigid rotation (Kabsch), then compared with\n    the deformed normal. Folded = cos <= 0. Degenerate source triangles counted separately (D1 rule).\n    "
    def normals(V):
        t = V[F]
        nn = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
        return nn, np.linalg.norm(nn, axis=1)
    V0 = np.asarray(V0, float)
    V1 = np.asarray(V1, float)
    F = np.asarray(F)
    n0, a0 = normals(V0)
    n1, a1 = normals(V1)
    R = _patch_rotations(V0, V1, F, rings=rings)
    n0r = np.einsum('fij,fj->fi', R, n0)
    with np.errstate(divide='ignore', invalid='ignore'):
        cosang = np.einsum('ij,ij->i', n0r, n1) / (a0 * a1)
    t = V0[F]
    longest = np.max(np.stack([np.linalg.norm(t[:, 1] - t[:, 0], axis=1),
                               np.linalg.norm(t[:, 2] - t[:, 1], axis=1),
                               np.linalg.norm(t[:, 0] - t[:, 2], axis=1)]), axis=0)
    degenerate = (a0 / longest) / longest < DEGEN_REL
    flipped = cosang <= 0.0
    return dict(n_folded_faces=int(np.sum(flipped & ~degenerate)),
                n_degenerate_source_faces=int(degenerate.sum()),
                min_cos_normal_change=float(np.nanmin(cosang)) if len(cosang) else float('nan'),
                rotation_free=True, rings=rings)


# ---------------------------------------------------------------- 3. attachment-to-surface
def distance_to_surface(V, F, points):
    "KOPIA with source citation: results/X1b/geomgr_ll2/landmarks.py rad 92-96 (identisk).\n    N7c: numpy exact point-triangle distance when trimesh is absent (cloud)."
    try:
        import trimesh  # noqa: F401
    except ImportError:
        from ..identity import closest_on_mesh
        return closest_on_mesh(np.asarray(V, float), np.asarray(F), np.atleast_2d(points), k=min(64, len(F)))[2]
    import trimesh
    mesh = trimesh.Trimesh(np.asarray(V), np.asarray(F), process=False)
    _cp, dist, _tri = trimesh.proximity.closest_point(mesh, np.atleast_2d(points))
    return np.asarray(dist)


def attachment_residuals(parts_def, attachments):
    """Distance from each attachment position to its deformed part surface.

    parts_def: {part: (V_def, F_def)}. attachments: list of dict(id, part, kind, pos).
    Only 'vertex'/'barycentric' gate; 'vertex_set' (region centroid, not on surface) is reported
    separately. Returns dict with per-attachment (id -> mm), max over gating kinds, and the id.
    """
    dists = {}
    for a in attachments:
        V, F = parts_def[a['part']]
        d = distance_to_surface(V, F, np.asarray(a['pos'], float)[None])[0]
        dists[a['id']] = dict(mm=float(d), kind=a['kind'], part=a['part'])
    gating = {k: v for k, v in dists.items() if v['kind'] in ('vertex', 'barycentric')}
    if gating:
        worst = max(gating, key=lambda k: gating[k]['mm'])
        mx, wid = gating[worst]['mm'], worst
    else:
        mx, wid = 0.0, None
    return dict(per_attachment=dists, attachment_max_dist_mm=float(mx), attachment_worst_id=wid,
                n_gating=len(gating))


# ---------------------------------------------------------------- 4. joint-centre sphere residual
def sphere_fit_algebraic(P):
    "COPY with source reference: results/X1b/geomgr_ll2/joints.py lines 14-19 (identical)."
    P = np.asarray(P, float)
    A = np.c_[2 * P, np.ones(len(P))]
    x = np.linalg.lstsq(A, (P ** 2).sum(1), rcond=None)[0]
    c = x[:3]
    return c, float(np.sqrt(x[3] + c @ c))


def sphere_fit(P):
    "COPY with source reference: results/X1b/geomgr_ll2/joints.py lines 22-38 (identical)."
    from scipy.optimize import least_squares
    P = np.asarray(P, float)
    c0, r0 = sphere_fit_algebraic(P)
    f = lambda x: np.linalg.norm(P - x[:3], axis=1) - x[3]  # noqa: E731
    sol = least_squares(f, np.r_[c0, r0], method='lm')
    res = sol.fun
    return dict(center=sol.x[:3], radius=float(sol.x[3]), rms=float(np.sqrt(np.mean(res ** 2))),
                max_abs=float(np.max(np.abs(res))), n=len(P))


def joint_residuals(joints):
    """Residual of a declared joint centre from the sphere fitted to the joint-head points.

    joints: list of dict(id, center (3,), head_points (k,3)). Returns per-joint residual (mm),
    the head RMS, and the max residual.
    """
    out = {}
    for j in joints:
        c = np.asarray(j['center'], float)
        head = np.asarray(j['head_points'], float)
        if len(head) < 4:
            out[j['id']] = dict(residual_mm=None, head_rms_mm=None, n=int(len(head)),
                                note='too few head points')
            continue
        fit = sphere_fit(head)
        out[j['id']] = dict(residual_mm=float(np.linalg.norm(c - fit['center'])),
                            head_rms_mm=fit['rms'], head_radius_mm=fit['radius'],
                            fitted_center=fit['center'].tolist(), n=fit['n'])
    vals = [v['residual_mm'] for v in out.values() if v['residual_mm'] is not None]
    return dict(per_joint=out, joint_max_residual_mm=float(max(vals)) if vals else 0.0)


# ---------------------------------------------------------------- certificate
def _scalars(obj):
    if isinstance(obj, dict):
        return {k: _scalars(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_scalars(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return None  # drop big arrays; per-face arrays stay in the npz if needed
    if isinstance(obj, (np.floating, np.integer, np.bool_)):
        return obj.item()
    return obj


def geometry_certificate_v2(vertices_src, faces_src, vertices_def, faces_def, *,
                            tolerance_mm=1e-6, declared_units='mm', attachments=None, joints=None,
                            evaluate=True, fold_rings=2):
    """v1 certificate on the deformed surface + the four BT-N5 controls.

    vertices_src/faces_src: source geometry (undeformed). vertices_def/faces_def: deformed.
    attachments: list of dict(id, part, kind, pos) with pos on the deformed part.
    joints: list of dict(id, center, head_points) in deformed coordinates.
    Returns a certificate dict; accepted = all gates pass. Arrays are dropped from the summary
    (kept only as scalars/minima); a per-face npz is written by the caller if needed.
    """
    V0 = np.asarray(vertices_src, float)
    V1 = np.asarray(vertices_def, float)
    F0 = np.asarray(faces_src)
    F1 = np.asarray(faces_def)

    cert = dict(geometry_certificate_v1(V1, F1, tolerance_mm=tolerance_mm,
                                        declared_units=declared_units, evaluate=evaluate))

    jac = per_face_jacobian(V0, V1, F0)
    fold = rotation_invariant_fold(V0, V1, F0, rings=fold_rings)
    # Attachments reference a part; each certified case is a single part mesh, so every
    # attachment's surface is this same (deformed) mesh.
    att = attachment_residuals({a['part']: (V1, F1) for a in (attachments or [])}, attachments or [])
    jnt = joint_residuals(joints or [])

    cert['jacobian'] = _scalars(jac)
    cert['fold'] = _scalars(fold)
    cert['attachments'] = _scalars(att)
    cert['joints'] = _scalars(jnt)
    cert['thresholds'] = dict(jac_tol=JAC_TOL, fold_faces=0, attach_tol_mm=ATTACH_TOL_MM,
                              joint_tol_mm=JOINT_TOL_MM)

    gates = dict(
        v1_accepted=bool(cert['accepted']),
        no_negative_jacobian=jac['n_negative_jacobian_faces'] == 0,
        no_fold=fold['n_folded_faces'] == 0,
        attachments_on_surface=att['attachment_max_dist_mm'] <= ATTACH_TOL_MM,
        joints_on_head=jnt['joint_max_residual_mm'] <= JOINT_TOL_MM,
    )
    cert['gates'] = gates
    cert['accepted'] = bool(all(gates.values()))
    if cert['accepted']:
        cert['reason_v2'] = 'accepted_v2: v1 accepted and 0 negative-Jacobian elements, 0 folds, ' \
                            'attachments on surface, joint centres inside head sphere'
    else:
        bad = [k for k, v in gates.items() if not v]
        cert['reason_v2'] = 'rejected_v2: ' + ', '.join(bad)
    return cert
