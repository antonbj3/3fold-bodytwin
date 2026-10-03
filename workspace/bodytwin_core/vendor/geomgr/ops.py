"geomgr.ops: N7b interpretations as operations on an Instance (stable ids, certificate v2 per op). N7c.\n\n  op_edit_osteotomy(inst, target, delta)  (b) proximal rigid osteotomy above 78 % of L_mech, +-10 mm smoothstep\n                                          band; Newton on (twist, varus, lengthening) so that `target` in\n                                          (CCD, AV, L_mech) moves by `delta` and the other two are held.\n                                          Measures = geomgr.measures CHOSEN definition (p1_smooth). Distal part\n                                          bit-identical. Source: results/N7b/b_edit.py (frame, osteotomy, solve).\n  op_tissue(inst, thickness_mm, source)   (e) per-vertex fields cortical_thickness_mm (CT-measured, N7b e) and\n                                          cartilage_thickness_mm (2 mm on tissue.head_articular, ASSUMED); layer\n                                          shell (inner surface V - t n) with exact prism volume per region.\n                                          Regions/footprints are registry entities (geomgr.extend).\n  op_register_surface(gm, dense, normals) (c) surface -> template correspondence: rigid ICP (4 PCA starts) ->\n                                          staged TPS-ICP -> projection onto the tangent plane of the nearest\n                                          surface sample (P1 projects to the closest mesh point; only a point\n                                          cloud is given here; snapping to the nearest sample folded triangles,\n                                          PREREG_ADD2). numpy/scipy copy of results/N7b/cloud_job3/cloud_c.py.\n  op_outlier(gm, landmarks, names)        BT-B35 leave-one-landmark-out test (no geometry change).\nEvery op returns a new Instance (log entry appended) and, where geometry changes, a v2 certificate."
from __future__ import annotations

import time

import numpy as np

from . import certify2 as C2
from . import core as C
from . import identity as I
from . import measures as MS

CUT_FRAC = 0.78
BAND_MM = 10.0
TARGETS = ('CCD', 'AV', 'L_mech')
MATERIALS = {   # results/N7b/e_regions.py MATERIALS (copied)
    'cortical_bone': dict(E_MPa=17000.0, status='literature', src='Reilly & Burstein 1975, longitudinal'),
    'trabecular_bone': dict(E_MPa='E = 6850*rho_app^1.49 (Morgan 2003, femoral neck)',
                            status='literature; rho from HU needs calibration (UNKNOWN here)'),
    'articular_cartilage': dict(E_MPa=10.0, thickness_mm=2.0, status='assumed; thickness UNKNOWN (not visible in CT)'),
    'tendon_attachment': dict(status='region only; footprint radius 8 mm assumed, material UNKNOWN'),
}


# ------------------------------------------------------------------ (b) edit
def _smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)


def _rodrigues(X, axis, ang):
    a = np.asarray(axis, float) / np.linalg.norm(axis)
    c, s = np.cos(ang)[:, None], np.sin(ang)[:, None]
    return X * c + np.cross(a, X) * s + np.outer(X @ a, a) * (1 - c)


def _measure(inst_reg, V, F):
    L = {e.name: np.asarray(e.weights) @ V[np.asarray(e.indices)] for e in inst_reg.by_kind('landmark')}
    m = MS.p1_smooth(L, V)
    return m, L


def edit_frame(reg, V, F):
    m, L = _measure(reg, V, F)
    h = (V - m['_knee_c']) @ m['_mech']
    hcut = CUT_FRAC * m['L_mech']
    c = V[np.abs(h - hcut) < 2.0].mean(0)
    neck_c = np.mean([L[k] for k in MS.NECK4], 0)
    neck = C.unit(m['_head_c'] - neck_c)
    a_v = C.unit(np.cross(neck, m['_shaft']))
    w = _smoothstep((h - (hcut - BAND_MM)) / (2 * BAND_MM))
    return dict(y={k: float(m[k]) for k in TARGETS}, c=c, shaft=m['_shaft'], a_v=a_v, w=w, hcut=hcut)


def osteotomy(V, fr, u):
    th, ph, s = u
    w = fr['w']
    act = w > 0
    Vn = V.copy()
    if th == 0 and ph == 0 and s == 0:
        return Vn
    X = V[act] - fr['c']
    wa = w[act]
    Xr = _rodrigues(_rodrigues(X, fr['shaft'], wa * th), fr['a_v'], wa * ph)
    Vn[act] = V[act] + (wa * s)[:, None] * fr['shaft'] + (Xr - X)
    return Vn


def op_edit_osteotomy(inst, target='AV', delta=10.0, iters=12, tol=1e-4):
    """Clinical-parameter edit (N7b b). Returns (new Instance, report incl. certificate v2)."""
    t0 = time.process_time()
    if inst.units != 'mm':
        raise ValueError('edit expects mm')
    V, F, reg = inst.V, inst.F, inst.reg
    fr = edit_frame(reg, V, F)
    k = TARGETS.index(target)
    y0 = np.array([fr['y'][t] for t in TARGETS])
    tgt = y0.copy()
    tgt[k] += delta

    def ym(u):
        m, _ = _measure(reg, osteotomy(V, fr, u), F)
        return np.array([m[t] for t in TARGETS])
    u = np.zeros(3)
    steps = np.array([np.radians(0.2), np.radians(0.2), 0.2])
    hist = []
    for _ in range(iters):
        err = ym(u) - tgt
        hist.append(float(np.abs(err).max()))
        if np.abs(err).max() < tol:
            break
        J = np.zeros((3, 3))
        for j in range(3):
            du = np.zeros(3)
            du[j] = steps[j]
            J[:, j] = (ym(u + du) - ym(u - du)) / (2 * steps[j])
        u = u - np.linalg.solve(J, err)
    V1 = osteotomy(V, fr, u)
    y1 = ym(u)
    out = inst.copy()
    out.V = V1
    I._log(out, 'edit_osteotomy', target=target, delta=float(delta), measure=MS.MEASURE_VERSION,
           u=[float(np.degrees(u[0])), float(np.degrees(u[1])), float(u[2])])
    cert = C2.instance_certificate(inst, out, f'edit_osteotomy:{target}{delta:+g}', local=False)
    rep = dict(target=target, delta=float(delta), measure=MS.MEASURE_VERSION, y0=dict(zip(TARGETS, y0.tolist())),
               y1=dict(zip(TARGETS, y1.tolist())), target_err=float(y1[k] - tgt[k]),
               held_maxabs=float(max(abs(y1[j] - y0[j]) for j in range(3) if j != k)),
               u=dict(twist_deg=float(np.degrees(u[0])), varus_deg=float(np.degrees(u[1])), lengthening_mm=float(u[2])),
               newton_maxerr=hist, distal_bit_identical=bool(np.array_equal(V1[fr['w'] == 0], V[fr['w'] == 0])),
               n_moved=int((fr['w'] > 0).sum()), certificate=cert, cpu_s=time.process_time() - t0)
    return out, rep


# ------------------------------------------------------------------ (e) tissue layers
def vertex_normals(V, F):
    fn = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    N = np.zeros_like(V)
    for k in range(3):
        np.add.at(N, F[:, k], fn)
    return N / np.maximum(np.linalg.norm(N, axis=1), 1e-12)[:, None]


def _prism_volumes(Vo, Vi, F):
    """results/N7b/e_regions.py shell_tets + tet_vol (copied): exact prism volume per face (3 tetrahedra)."""
    a, b, c = F[:, 0], F[:, 1], F[:, 2]
    T = [np.stack([Vo[a], Vo[b], Vo[c], Vi[a]], 1), np.stack([Vo[b], Vo[c], Vi[a], Vi[b]], 1),
         np.stack([Vo[c], Vi[a], Vi[b], Vi[c]], 1)]
    vol = lambda X: np.abs(np.einsum('ij,ij->i', X[:, 1] - X[:, 0], np.cross(X[:, 2] - X[:, 0], X[:, 3] - X[:, 0]))) / 6  # noqa: E731
    return sum(vol(t) for t in T)


def op_tissue(inst, thickness_mm=None, source='UNKNOWN (no CT for this individual)'):
    """Attach tissue fields and compute the cortical layer shell. thickness_mm: per-vertex (NaN = not measurable)."""
    out = inst.copy()
    n = len(inst.V)
    t = np.full(n, np.nan) if thickness_mm is None else np.asarray(thickness_mm, float).copy()
    out.fields['cortical_thickness_mm'] = dict(values=t, kind='length', units=inst.units, source=source)
    cart = np.zeros(n)
    vals = inst.values()
    ha = vals.get(f'{inst.reg.part}/region/tissue.head_articular')
    if ha is not None:
        cart[ha['indices']] = 2.0
    out.fields['cartilage_thickness_mm'] = dict(values=cart, kind='length', units=inst.units,
                                                source='ASSUMED 2 mm on tissue.head_articular (N7b e; UNKNOWN)')
    I._log(out, 'tissue', cortical_source=source, n_measured=int(np.isfinite(t).sum()))
    rep = dict(materials=MATERIALS, cortical_source=source, measured_fraction=float(np.isfinite(t).mean()))
    if np.isfinite(t).any():
        N = vertex_normals(inst.V, inst.F)
        tt = np.where(np.isfinite(t), t, 0.0)
        Vi = inst.V - tt[:, None] * N
        pv = _prism_volumes(inst.V, Vi, inst.F)
        rep['cortical_shell_volume_mm3'] = float(pv.sum())
        per = {}
        for b, v in vals.items():
            if b.startswith(f'{inst.reg.part}/region/tissue.'):
                s = np.zeros(n, bool)
                s[v['indices']] = True
                fsel = s[inst.F].all(1)
                per[b.split('/')[-1]] = dict(median_thickness_mm=float(np.nanmedian(t[v['indices']])),
                                             measured_fraction=float(np.isfinite(t[v['indices']]).mean()),
                                             shell_volume_mm3=float(pv[fsel].sum()))
        rep['per_region'] = per
    return out, rep


# ------------------------------------------------------------------ (c) registration
def _umeyama(A, B, scale=True, w=None):
    return C.umeyama(A, B, scale=scale, w=w)


def _rigid_icp(X, dense, tree, back, iters=30):
    from scipy.spatial import cKDTree
    X = X.copy()
    for _ in range(iters):
        d, i = tree.query(X)
        dd, j = cKDTree(X).query(back)
        A = np.r_[X, X[j]]
        B = np.r_[dense[i], back]
        r = np.r_[d, dd]
        w = (r < 3 * np.median(r) + 1e-9).astype(float)
        s, R, t = _umeyama(A, B, scale=True, w=w)
        X = s * X @ R.T + t
    d, _ = tree.query(X)
    dd, _ = cKDTree(X).query(back)
    return X, float(np.sqrt(0.5 * (np.mean(d ** 2) + np.mean(dd ** 2))))


def _init_rigid(V, dense, tree, back):
    ct, cs = V.mean(0), dense.mean(0)
    _, _, At = np.linalg.svd(V - ct, full_matrices=False)
    _, _, Bt = np.linalg.svd(dense - cs, full_matrices=False)
    cands = [np.eye(3)]
    for sx in (1, -1):
        for sy in (1, -1):
            S = np.diag([sx, sy, 1.0])
            R = Bt.T @ S @ At
            if np.linalg.det(R) < 0:
                S[2, 2] = -1.0
                R = Bt.T @ S @ At
            cands.append(R)
    best = None
    for R in cands:
        X, err = _rigid_icp((V - ct) @ R.T + cs, dense, tree, back, iters=25)
        if best is None or err < best[1]:
            best = (X, err)
    return best


def _tps_icp(X, F, dense, dense_n, tree, back, seed=0,
             stages=((40, 3.0), (80, 1.0), (160, 0.3), (320, 0.1), (640, 0.03), (1200, 0.01)), iters=8):
    from scipy.interpolate import RBFInterpolator
    from scipy.spatial import cKDTree
    from .wrap import fps
    n = len(X)
    for ncp, lam in stages:
        for _ in range(iters):
            N = vertex_normals(X, F)
            d, i = tree.query(X)
            c = dense[i]
            ok = ((N * dense_n[i]).sum(1) > 0.0) & (d < max(5.0, 4 * np.median(d)))
            disp = np.zeros_like(X)
            cnt = np.zeros(n)
            disp[ok] += (c - X)[ok]
            cnt[ok] += 1
            dd, j = cKDTree(X).query(back)
            okb = dd < max(5.0, 4 * np.median(dd))
            np.add.at(disp, j[okb], (back - X[j])[okb])
            np.add.at(cnt, j[okb], 1)
            has = cnt > 0
            disp[has] /= cnt[has][:, None]
            cp = fps(X, min(ncp, n), seed=seed)
            k = max(4, int(3 * n / ncp))
            _, nb = cKDTree(X[has]).query(X[cp], k=min(k, int(has.sum())))
            dcp = disp[has][nb].mean(1)
            f = RBFInterpolator(X[cp] / 100.0, dcp, kernel='thin_plate_spline', smoothing=lam * ncp / 40.0)
            X = X + f(X / 100.0)
    return X


def op_register_surface(gm, dense, normals, subject='surface', stages=None, reproject='plane'):
    """(c) Register the template (population mean in correspondence) to a target surface sample (dense points +
    outward normals, mm). Returns (Instance in the target frame, report). Scale is removed afterwards so that the
    instance keeps the target's size (the rigid-ICP scale only initialises)."""
    from scipy.spatial import cKDTree
    t0 = time.process_time()
    dense = np.asarray(dense, float)
    back = dense[np.random.default_rng(0).choice(len(dense), min(10000, len(dense)), replace=False)]
    tree = cKDTree(dense)
    TV = gm.model.mu
    X, rig = _init_rigid(TV, dense, tree, back)
    kw = {} if stages is None else dict(stages=stages)
    Xpre = _tps_icp(X, gm.F, dense, np.asarray(normals, float), tree, back, **kw)
    d, i = tree.query(Xpre)
    nrm = np.asarray(normals, float)
    if reproject == 'nearest_sample':      # first version (N7c run 1): snaps vertices onto samples -> folds
        Xp = dense[i]
    else:                                  # point-to-plane projection onto the nearest sample's tangent plane
        Xp = Xpre - (np.einsum('ij,ij->i', Xpre - dense[i], nrm[i]))[:, None] * nrm[i]
    inst = I.Instance(gm.registry, Xp, gm.F, frame=f'observation@{subject}', subject=subject)
    I._log(inst, 'register_surface', rigid_rms_mm=rig, n_dense=len(dense), reproject=reproject)
    cert = C2.certificate_v2(TV, gm.F, Xp, gm.F, attachments=C2.instance_declarations(inst)[0],
                             joints=C2.instance_declarations(inst)[1], local=False, op='register_surface')
    dd, _ = cKDTree(Xp).query(back)
    return inst, dict(reproject=reproject, rigid_rms_mm=rig, reproj_shift_median_mm=float(np.median(d)),
                      surface_to_instance_median_mm=float(np.median(dd)), certificate=cert,
                      cpu_s=time.process_time() - t0)


# ------------------------------------------------------------------ (h) wraps + B35
def wraps_summary(inst):
    vals = inst.values()
    out = {}
    for b, e in inst.reg.entities.items():
        if e.kind == 'wrap':
            v = vals[b]
            out[e.id] = dict(R_mm=v['R'], fit_rms_mm=v['rms'], n_fit=v['n'], length_mm=v['L'])
    return out


def op_outlier(gm, landmarks, names=None):
    """BT-B35 leave-one-landmark-out Mahalanobis test (threshold chi2_3(0.999)); rows per landmark."""
    from .outlier import landmark_outlier_test
    rows = landmark_outlier_test(gm, landmarks, names)
    return dict(rows=rows, flagged=[r['landmark'] for r in rows if r['flag']], test='BT-B35 LOLO chi2_3(0.999)')
