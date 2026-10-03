"""geomgr.core: numpy/scipy-only numerics of the geometry manager (runs on the OVH box as well).

Population model (GPA without scale + PCA over all n-1 modes), feature catalog, Gaussian conditioning
on any mix of scalar measures and landmark coordinates, predictive covariance per vertex, calibration
metrics and the numeric part of the per-operation certificate.

Copied functions (source files are NOT imported so this runs without the workspace tree):
  * sphere_fit, unit, measures, knee_frame, umeyama, apply  <- results/P1/p1_common.py
    (sha256 8e04d424ee61ad2a...; identical definitions, comments shortened)
  * gpa, b2 formula (stage A)                               <- results/P1/p1_eval.py (sha256 ecce5d4b5c99c90d...)
  * fold_check / _patch_rotations                           <- results/X1b/geomgr_ll2/model.py
    (sha256 aa46af8a97106e0e...; rotation-invariant fold check, D1 raw variant dropped)
"""
from __future__ import annotations

import numpy as np
from scipy import linalg, stats
import scipy.sparse as sp

CHI2_3_90 = float(stats.chi2.ppf(0.90, 3))   # 6.2514
Z_90 = float(stats.norm.ppf(0.95))           # 1.6449 (two-sided 90 %)
KAPPA_GRID = np.round(np.arange(0.50, 4.0001, 0.01), 2)

LM_NAMES = ['SFH', 'AFH', 'LFH', 'PFH', 'MFH', 'IFH', 'SNI', 'ANI', 'INI', 'PNI', 'SGT', 'LT', 'PITC',
            'PLC', 'PMC', 'PPLC', 'PPMC', 'LEC', 'MEC', 'ICN', 'DLC', 'DMC']
HEAD6 = ('SFH', 'AFH', 'LFH', 'PFH', 'MFH', 'IFH')
MEAS_KEYS = ["L_tl", "W_epi", "L_mech", "D_head", "CCD", "AV"]
PT_NAMES = LM_NAMES + ['HC6']
PAIRS = [(i, j) for i in range(len(PT_NAMES)) for j in range(i + 1, len(PT_NAMES))]
PAIR_NAMES = [f"d({PT_NAMES[i]},{PT_NAMES[j]})" for i, j in PAIRS]
FEMUR_FEATURES = list(MEAS_KEYS) + PAIR_NAMES
X5 = ["L_mech", "W_epi", "D_head", "CCD", "AV"]
T5 = ["d(SGT,DMC)", "AV", "d(IFH,SGT)", "d(SGT,DLC)", "CCD"]      # results/H7/part2.json recommended_all70 order[:5]


def feature_units(name):
    if name in ('CCD', 'AV'):
        return 'deg'
    if name in ('height',):
        return 'cm'
    if name in ('weight',):
        return 'kg'
    if name in ('age',):
        return 'yr'
    return 'mm'


# ------------------------------------------------------------------ P1 copies
def sphere_fit(P):
    A = np.c_[2 * P, np.ones(len(P))]
    b = (P ** 2).sum(1)
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    ctr = c[:3]
    return ctr, float(np.sqrt(c[3] + ctr @ ctr))


def unit(v):
    return v / np.linalg.norm(v)


def measures(L, V):
    """P1 p1_common.measures (identical)."""
    head_c, head_r = sphere_fit(np.array([L[k] for k in HEAD6]))
    neck_c = np.mean([L[k] for k in ("SNI", "ANI", "INI", "PNI")], 0)
    knee_c = 0.5 * (L["LEC"] + L["MEC"])
    mech = unit(head_c - knee_c)
    proj = (V - knee_c) @ mech
    lo, hi = proj.min(), proj.max()
    cents = []
    for f in np.linspace(0.30, 0.70, 9):
        h = lo + f * (hi - lo)
        sl = V[np.abs(proj - h) < 3.0]
        if len(sl) > 10:
            cents.append(sl.mean(0))
    C = np.array(cents)
    _, _, vt = np.linalg.svd(C - C.mean(0))
    shaft = unit(vt[0] * np.sign(vt[0] @ mech))
    neck = unit(head_c - neck_c)
    ccd = float(np.degrees(np.arccos(np.clip(neck @ -shaft, -1, 1))))

    def pp(v):
        return v - (v @ mech) * mech
    med = unit(pp(L["PMC"] - L["PLC"]))
    ant = unit(pp(knee_c - 0.5 * (L["PMC"] + L["PLC"])))
    ant = unit(ant - (ant @ med) * med)
    nn = pp(neck)
    av = float(np.degrees(np.arctan2(nn @ ant, nn @ med)))
    return {"L_tl": float(np.linalg.norm(L["SGT"] - L["LEC"])), "W_epi": float(np.linalg.norm(L["MEC"] - L["LEC"])),
            "L_mech": float(np.linalg.norm(head_c - knee_c)), "D_head": 2 * head_r, "CCD": ccd, "AV": av,
            "_head_c": head_c, "_knee_c": knee_c, "_mech": mech, "_shaft": shaft}


def knee_frame(L):
    """P1 p1_common.knee_frame: origin epicondyle midpoint, x MEC->LEC, y anterior, z proximal (rows = axes)."""
    kc = 0.5 * (L["LEC"] + L["MEC"])
    x = unit(L["LEC"] - L["MEC"])
    y = kc - 0.5 * (L["PLC"] + L["PMC"])
    y = unit(y - (y @ x) * x)
    z = np.cross(x, y)
    return kc, np.stack([x, y, z])


def umeyama(A, B, scale=True, w=None):
    """s,R,t with s*R@A+t ~ B (rows = points)."""
    if w is None:
        w = np.ones(len(A))
    w = w / w.sum()
    ma, mb = w @ A, w @ B
    A0, B0 = A - ma, B - mb
    C = (B0 * w[:, None]).T @ A0
    U, S, Vt = np.linalg.svd(C)
    D = np.eye(3)
    D[2, 2] = np.sign(np.linalg.det(U @ Vt))
    R = U @ D @ Vt
    s = (S * np.diag(D)).sum() / (w @ (A0 ** 2).sum(1)) if scale else 1.0
    t = mb - s * R @ ma
    return s, R, t


def apply(s, R, t, X):
    return s * X @ R.T + t


def gpa(S, iters=6):
    """P1 p1_eval.gpa (identical): generalized Procrustes without scale."""
    ref = S[0].copy()
    A = S.copy()
    for _ in range(iters):
        for i in range(len(S)):
            s, Rm, t = umeyama(A[i], ref, scale=False)
            A[i] = apply(s, Rm, t, A[i])
        new = A.mean(0)
        s, Rm, t = umeyama(new, ref, scale=False)
        ref = apply(s, Rm, t, new)
    return A, ref


# ------------------------------------------------------------------ anchors / features
class Anchors:
    """Linear anchors as a sparse (m x N) matrix: point = A @ V. Names in order."""

    def __init__(self, names, indices, weights, n_vertices):
        self.names = list(names)
        rows = np.concatenate([[i] * len(ix) for i, ix in enumerate(indices)]).astype(int)
        cols = np.concatenate([np.asarray(ix, int) for ix in indices])
        vals = np.concatenate([np.asarray(w, float) for w in weights])
        self.M = sp.csr_matrix((vals, (rows, cols)), shape=(len(self.names), n_vertices))
        self.index = {n: i for i, n in enumerate(self.names)}

    def points(self, V):
        return np.asarray(self.M @ V)

    def dict(self, V):
        P = self.points(V)
        return {n: P[i] for i, n in enumerate(self.names)}


def femur_points(anchors, V):
    """22 landmarks + HC6 (sphere through the 6 head landmarks) as dict."""
    L = anchors.dict(V)
    L['HC6'] = sphere_fit(np.array([L[k] for k in HEAD6]))[0]
    return L


def femur_features_from_landmarks(L, V):
    """Feature vector in FEMUR_FEATURES order (H7 cand_from_landmarks definitions)."""
    m = measures(L, V)
    P = np.array([L[k] for k in LM_NAMES] + [m["_head_c"]])
    d = np.array([np.linalg.norm(P[i] - P[j]) for i, j in PAIRS])
    return np.r_[[m[k] for k in MEAS_KEYS], d]


def tibia_points(anchors, V):
    return anchors.dict(V)


def tibia_geom_features(L):
    return {'L_tib': float(np.linalg.norm(L['TPROX'] - L['MM'])),
            'W_plat': float(np.linalg.norm(L['MTP'] - L['LTP']))}


# ------------------------------------------------------------------ population model
class ShapeModel:
    """GPA (no scale) + PCA over all n-1 modes. Members' features M (n,k) with names; groups = person ids."""

    def __init__(self, S, M=None, feat_names=(), groups=None, keep_aligned=True):
        S = np.asarray(S, float)
        self.n, self.N = S.shape[0], S.shape[1]
        A, _ = gpa(S)
        self.mu = A.mean(0)
        X = (A - self.mu).reshape(self.n, -1)
        U, s, Vt = np.linalg.svd(X, full_matrices=False)
        r = self.n - 1
        self.V = Vt[:r].T                        # (3N, r) orthonormal
        self.B = U[:, :r] * s[:r]                # (n, r) scores
        self.lam = s[:r] ** 2 / (self.n - 1)
        self.r = r
        self.aligned = A if keep_aligned else None
        self.M = None if M is None else np.asarray(M, float)
        self.feat_names = list(feat_names)
        self.groups = None if groups is None else np.asarray(groups)
        self.sig2_oos = np.zeros(self.N)          # set by set_oos()

    # geometry
    def shape(self, b):
        return self.mu + (self.V @ b).reshape(self.N, 3)

    def G(self):
        return self.V.reshape(self.N, 3, self.r)

    def project(self, X):
        """Rigidly align X to mu, coefficients, residual per vertex (N,)."""
        s, R, t = umeyama(X, self.mu, scale=False)
        Xa = apply(s, R, t, X)
        b = self.V.T @ (Xa - self.mu).ravel()
        res = Xa - self.shape(b)
        return b, (res ** 2).sum(1)

    def set_oos(self, sig2):
        self.sig2_oos = np.asarray(sig2, float)


def oos_variance(S, groups, feat=None):
    """Model incompleteness per vertex: inner leave-one-group-out projection residual, isotropic per axis.
    Returns sig2 (N,) = mean |res_v|^2 / 3 over held-out shapes."""
    S = np.asarray(S, float)
    groups = np.asarray(groups)
    acc = np.zeros(S.shape[1])
    cnt = 0
    for g in np.unique(groups):
        tr = groups != g
        m = ShapeModel(S[tr], keep_aligned=False)
        for i in np.flatnonzero(~tr):
            _, r2 = m.project(S[i])
            acc += r2
            cnt += 1
    return acc / cnt / 3.0


# ------------------------------------------------------------------ conditioning
STEP_TRACE = None


class Posterior:
    pass


def _skew(v):
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])


def _expm_so3(w):
    th = np.linalg.norm(w)
    K = _skew(w)
    if th < 1e-12:
        return np.eye(3) + K
    return np.eye(3) + np.sin(th) / th * K + (1 - np.cos(th)) / th ** 2 * (K @ K)


def stage_a(model, feat_obs, R_feat):
    """P1 B2 (p1_eval.Model.b2) generalized to any member features incl. covariates.
    Returns b1, S1 (posterior over b), innovation d2, dof, per-feature standardized innovation."""
    r = model.r
    Sbb = model.B.T @ model.B / (model.n - 1)
    if not feat_obs:
        return np.zeros(r), Sbb, 0.0, 0, {}
    names = list(feat_obs)
    idx = [model.feat_names.index(k) for k in names]
    mu_m = model.M.mean(0)[idx]
    Mc = model.M[:, idx] - mu_m
    n1 = model.n - 1
    Sbm = model.B.T @ Mc / n1
    Smm = Mc.T @ Mc / n1 + R_feat
    m = np.array([feat_obs[k] for k in names], float)
    K = Sbm @ np.linalg.inv(Smm)
    b = K @ (m - mu_m)
    Sc = Sbb - K @ Sbm.T
    Sc = 0.5 * (Sc + Sc.T)
    innov = m - mu_m
    d2 = float(innov @ np.linalg.solve(Smm, innov))
    zf = {k: float(innov[i] / np.sqrt(Smm[i, i])) for i, k in enumerate(names)}
    return b, Sc, d2, len(names), zf


def _lm_fun(model, anchors, names, point_fn):
    """f(b) -> (L,3) model-frame points for names; J(b) -> (L,3,r)."""
    lin = [n for n in names if n in anchors.index]
    Msel = anchors.M[[anchors.index[n] for n in lin]] if lin else None
    G = model.G()
    Jlin = np.einsum('ln,nkr->lkr', Msel.toarray(), G) if lin else None

    def f(b):
        V = model.shape(b)
        P = point_fn(anchors, V)
        return np.array([P[n] for n in names])

    def J(b):
        out = np.zeros((len(names), 3, model.r))
        for i, n in enumerate(names):
            if n in anchors.index:
                out[i] = Jlin[lin.index(n)]
        nonlin = [i for i, n in enumerate(names) if n not in anchors.index]
        if nonlin:
            h = 1e-3 * np.sqrt(np.maximum(model.lam, 1e-12))
            for k in range(model.r):
                e = np.zeros(model.r)
                e[k] = h[k]
                fp, fm = f(b + e), f(b - e)
                for i in nonlin:
                    out[i, :, k] = (fp[i] - fm[i]) / (2 * h[k])
        return out
    return f, J


def condition(model, anchors, point_fn, feat_obs=None, R_feat=None, lm_obs=None, C_lm=None,
              frame_fn=None, max_iter=60, tol=1e-7):
    """Gaussian conditioning on any subset of features (stage A, empirical joint Gaussian) and landmark
    coordinates in an arbitrary observation frame (stage B, Gauss-Newton over b and rigid pose, flat pose prior).

    C_lm: name -> 3x3 noise covariance in the ANATOMICAL frame given by frame_fn(points)->(origin, rows);
    the model-incompleteness variance at each landmark is added isotropically.
    """
    feat_obs = dict(feat_obs or {})
    lm_obs = dict(lm_obs or {})
    post = Posterior()
    b1, S1, d2, k, zf = stage_a(model, feat_obs, R_feat if feat_obs else None)
    post.stage_a = dict(d2=d2, dof=k, p=float(stats.chi2.sf(d2, k)) if k else 1.0, z=zf)
    post.b1, post.S1 = b1, S1
    r = model.r
    if not lm_obs:
        post.b, post.Sbb, post.pose, post.Sigma, post.pose_identified = b1, S1, None, None, False
        post.stage_b = None
        return post
    names = list(lm_obs)
    Y = np.array([lm_obs[n] for n in names], float)
    f, J = _lm_fun(model, anchors, names, point_fn)
    # isotropic incompleteness at landmarks: anchor-weighted sig2 (HC6: mean of head anchors)
    s2 = model.sig2_oos
    s2_lm = []
    for n in names:
        if n in anchors.index:
            row = anchors.M[anchors.index[n]]
            s2_lm.append(float(row.data @ s2[row.indices]))
        else:
            s2_lm.append(float(np.mean([s2[anchors.M[anchors.index[h]].indices].mean() for h in HEAD6])))
    b = b1.copy()
    P = f(b)
    if len(names) >= 3:
        _, R, t = umeyama(P, Y, scale=False)
    else:
        R, t = np.eye(3), Y.mean(0) - P.mean(0)
    S1i = np.linalg.inv(S1)

    def noise_obs(R, t, b):
        Vp = model.shape(b) @ R.T + t
        o, Ax = frame_fn(point_fn(anchors, Vp))
        Cs = []
        for i, n in enumerate(names):
            Ca = C_lm.get(n, np.zeros((3, 3))) if C_lm else np.zeros((3, 3))
            Cs.append(Ax.T @ Ca @ Ax + s2_lm[i] * np.eye(3))
        return Cs

    Cs = noise_obs(R, t, b)
    Wb = linalg.block_diag(*[np.linalg.inv(C) for C in Cs])
    it = 0
    for it in range(max_iter):
        P = f(b)
        Jf = J(b)
        RP = P @ R.T
        res = (Y - (RP + t)).ravel()
        Jb = np.einsum('ij,ljr->lir', R, Jf).reshape(-1, r)
        Jw = np.concatenate([-_skew(p) for p in RP], 0)
        Jt = np.tile(np.eye(3), (len(names), 1))
        Jall = np.c_[Jb, Jw, Jt]
        H = Jall.T @ Wb @ Jall
        H[:r, :r] += S1i
        H[r:, r:] += 1e-9 * np.eye(6)
        g = Jall.T @ Wb @ res
        g[:r] -= S1i @ (b - b1)
        d = np.linalg.solve(H, g)
        b = b + d[:r]
        R = _expm_so3(d[r:r + 3]) @ R
        t = t + d[r + 3:]
        if it == 3:                     # re-evaluate the frame-dependent noise once the pose has settled
            Cs = noise_obs(R, t, b)
            Wb = linalg.block_diag(*[np.linalg.inv(C) for C in Cs])
        if STEP_TRACE is not None:
            STEP_TRACE.append(float(np.linalg.norm(d)))
        if np.linalg.norm(d) < tol and it > 3:
            break
    P = f(b)
    Jf = J(b)
    RP = P @ R.T
    res = (Y - (RP + t)).ravel()
    Jb = np.einsum('ij,ljr->lir', R, Jf).reshape(-1, r)
    Jall = np.c_[Jb, np.concatenate([-_skew(p) for p in RP], 0), np.tile(np.eye(3), (len(names), 1))]
    H = Jall.T @ Wb @ Jall
    H[:r, :r] += S1i
    # pose identifiability: Schur complement of the pose block
    Hpp = H[r:, r:] - H[r:, :r] @ np.linalg.solve(H[:r, :r], H[:r, r:])
    ev = np.linalg.eigvalsh(0.5 * (Hpp + Hpp.T))
    pose_ok = bool(ev.min() > 1e-8 * max(ev.max(), 1e-300))
    Hreg = H.copy()
    Hreg[r:, r:] += 1e-9 * np.eye(6)
    Sigma = np.linalg.inv(Hreg)
    Sigma = 0.5 * (Sigma + Sigma.T)
    ssr = float(res @ Wb @ res)
    dof = max(1, 3 * len(names) - (6 if pose_ok else 5))
    per_lm = {}
    for i, n in enumerate(names):
        ri = res[3 * i:3 * i + 3]
        per_lm[n] = float(ri @ np.linalg.solve(Cs[i], ri))
    post.b, post.Sbb, post.pose, post.Sigma = b, Sigma[:r, :r], (R, t), Sigma
    post.pose_identified = pose_ok
    post.stage_b = dict(ssr=ssr, dof=dof, p=float(stats.chi2.sf(ssr, dof)), iters=it + 1,
                        per_landmark_d2=per_lm, residual_rms_mm=float(np.sqrt(np.mean(res.reshape(-1, 3) ** 2) * 3)),
                        pose_min_eig=float(ev.min()))
    return post


def _psd_sqrt(S):
    w, Q = np.linalg.eigh(0.5 * (S + S.T))
    return Q * np.sqrt(np.clip(w, 0, None))


def predictive(model, post, placed=False, n_eff=None):
    """Mean vertices (N,3) and covariance (N,3,3). Shape frame = model frame; placed = observation frame
    (requires an identified pose). Includes (1 + 1/n) for the estimated mean and the incompleteness term."""
    n_eff = model.n if n_eff is None else n_eff
    infl = 1.0 + 1.0 / n_eff
    G = model.G()
    X = model.shape(post.b)
    if not placed:
        L = _psd_sqrt(post.Sbb)
        GL = np.einsum('nkr,rs->nks', G, L)
        C = np.einsum('nks,nls->nkl', GL, GL) * infl
        C += model.sig2_oos[:, None, None] * np.eye(3)
        return X, C
    R, t = post.pose
    Xp = X @ R.T + t
    RX = X @ R.T
    L = _psd_sqrt(post.Sigma)
    r = model.r
    Jb = np.einsum('ij,njr->nir', R, G)
    Jw = np.zeros((model.N, 3, 3))
    Jw[:, 0, 1], Jw[:, 0, 2] = RX[:, 2], -RX[:, 1]
    Jw[:, 1, 0], Jw[:, 1, 2] = -RX[:, 2], RX[:, 0]
    Jw[:, 2, 0], Jw[:, 2, 1] = RX[:, 1], -RX[:, 0]
    Jall = np.concatenate([Jb, Jw, np.broadcast_to(np.eye(3), (model.N, 3, 3))], 2)
    JL = np.einsum('nkq,qs->nks', Jall, L)
    C = np.einsum('nks,nls->nkl', JL, JL) * infl
    C += model.sig2_oos[:, None, None] * np.eye(3)
    return Xp, C


def landmark_predictive(model, anchors, point_fn, post, names, placed=False):
    """Mean and 3x3 covariance for named points (linear anchors exact, HC6 numerically linearized)."""
    f, J = _lm_fun(model, anchors, names, point_fn)
    b = post.b
    P = f(b)
    Jf = J(b)
    infl = 1.0 + 1.0 / model.n
    s2 = model.sig2_oos
    out = {}
    for i, n in enumerate(names):
        if n in anchors.index:
            row = anchors.M[anchors.index[n]]
            s2i = float(row.data @ s2[row.indices])
        else:
            s2i = float(np.mean([s2[anchors.M[anchors.index[h]].indices].mean() for h in HEAD6]))
        if not placed:
            C = Jf[i] @ post.Sbb @ Jf[i].T * infl + s2i * np.eye(3)
            out[n] = (P[i], C)
        else:
            R, t = post.pose
            p = R @ P[i] + t
            Ji = np.c_[R @ Jf[i], -_skew(R @ P[i]), np.eye(3)]
            C = Ji @ post.Sigma @ Ji.T * infl + s2i * np.eye(3)
            out[n] = (p, C)
    return out


# ------------------------------------------------------------------ metrics
def vertex_normals(V, F):
    t = V[F]
    fn = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    vn = np.zeros_like(V)
    for k in range(3):
        np.add.at(vn, F[:, k], fn)
    return vn / np.maximum(np.linalg.norm(vn, axis=1, keepdims=True), 1e-300)


def coverage_curves(e, C, normals=None):
    """Per-vertex 3D Mahalanobis coverage and 1D normal coverage on KAPPA_GRID (covariance scaled by kappa^2)."""
    Ci = np.linalg.inv(C)
    z2 = np.einsum('ni,nij,nj->n', e, Ci, e)
    k2 = KAPPA_GRID ** 2
    cov3 = np.array([(z2 <= CHI2_3_90 * k) .mean() for k in k2])
    out = dict(z2_median=float(np.median(z2)), cov3=cov3)
    if normals is not None:
        sn = np.sqrt(np.einsum('ni,nij,nj->n', normals, C, normals))
        zn = np.abs(np.einsum('ni,ni->n', normals, e)) / np.maximum(sn, 1e-12)
        out['cov1'] = np.array([(zn <= Z_90 * k).mean() for k in KAPPA_GRID])
        out['sig_n_median'] = float(np.median(sn))
    return out


def kappa_for(curves, target=0.90):
    """Smallest kappa on the grid whose pooled coverage reaches the target."""
    c = np.mean(np.asarray(curves), 0)
    i = int(np.searchsorted(c, target))
    return float(KAPPA_GRID[min(i, len(KAPPA_GRID) - 1)])


def at_kappa(curve, kappa):
    return float(np.asarray(curve)[int(np.argmin(np.abs(KAPPA_GRID - kappa)))])


# ------------------------------------------------------------------ fold check (X1b copy)
def _patch_rotations(V0, V1, F, rings=2):
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
    return np.einsum('fji,fj,fkj->fik', Vt, D, U)


def fold_check(V0, V1, F, rings=2):
    """X1b rotation-invariant fold check: folded = cos(normal change after removing patch rotation) <= 0."""
    V0 = np.asarray(V0, float)
    V1 = np.asarray(V1, float)
    F = np.asarray(F)

    def normals(V):
        t = V[F]
        nn = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
        return nn, np.linalg.norm(nn, axis=1)
    n0, a0 = normals(V0)
    n1, a1 = normals(V1)
    R = _patch_rotations(V0, V1, F, rings=rings)
    n0r = np.einsum('fij,fj->fi', R, n0)
    cosang = np.einsum('ij,ij->i', n0r, n1) / np.maximum(a0 * a1, 1e-300)
    t = V0[F]
    longest = np.max(np.stack([np.linalg.norm(t[:, 1] - t[:, 0], axis=1), np.linalg.norm(t[:, 2] - t[:, 1], axis=1),
                               np.linalg.norm(t[:, 0] - t[:, 2], axis=1)]), axis=0)
    degenerate = (a0 / longest) / longest < 1e-3
    flipped = cosang <= 0.0
    return dict(n_folded_faces=int(np.sum(flipped & ~degenerate)), n_degenerate_source_faces=int(degenerate.sum()),
                min_area_ratio=float(np.min(a1 / np.maximum(a0, 1e-300))),
                max_area_ratio=float(np.max(a1 / np.maximum(a0, 1e-300))),
                min_cos_normal_change=float(np.min(cosang)))


def signed_volume(V, F):
    t = V[F]
    return float(np.einsum('ij,ij->i', t[:, 0], np.cross(t[:, 1], t[:, 2])).sum() / 6.0)
