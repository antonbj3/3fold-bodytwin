"""N1 core: landmark morphs of TLEM (affine, affine-HJC, RBF φ=r, RBF φ=r²log r; two
landmark sources) and our population-conditioned morphs (PSM, PSM+R) with common error measurement.

Read-only import: results/P3 (p3_chain -> P1, X1), results/RM1 (rm1_core: kriging/REML for the RTPS reference).
No files outside results/N1 and external_media are written. See PREREG.md.
"""
from __future__ import annotations

import os
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_k, "2")
import sys
import json
from pathlib import Path

import numpy as np
import trimesh
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
RES = HERE.parent
sys.dont_write_bytecode = True
CLOUD = os.environ.get("N1_CLOUD")
if CLOUD:                                   # cloud run: vendored pure Python packages + dataset in the package (cloud_bundle.sh)
    _B = Path(CLOUD)
    sys.path.insert(0, str(_B / "vendor"))
sys.path[:0] = [str(RES / "P3"), str(RES / "RM1"), str(RES / "P1"), str(RES / "X1")]
if CLOUD:
    import p1_common as _pc
    from geomgr_ll import tlem as _tl, vsd as _vs
    _pc.ROOT = _B / "data"
    _pc.IMP_ZIP = _B / "data" / "femur_3D_surface_meshes.zip"
    _pc.VSD = _vs.VSD = _B / "data" / "VSD"
    _pc.TLEM = _tl.ROOT = _B / "data" / "TLEM2.0"
    _tl.BONES = _tl.ROOT / "TLEM 2.0 - Bones - Local Reference Frame"
    _tl.TABLES = _tl.ROOT / "TLEM 2.0 - Musculoskeletal Model Dataset"
import p3_chain as P  # noqa: E402
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_k, "2")
from geomgr_ll import morph, tlem, joints  # noqa: E402
from geomgr_ll.tps import TPS3, fit_affine, apply_affine  # noqa: E402
from p1_common import LM_NAMES, load_vsd_femur, load_imperial, measures  # noqa: E402

E = P.E
F = E.F
TMP = Path(os.environ.get("N1_TMP", "external_media"))
PRE = Path(os.environ.get("N1_PRE", "external_media"))
TMP.mkdir(parents=True, exist_ok=True)

HEAD6 = ("SFH", "AFH", "LFH", "PFH", "MFH", "IFH")
SETS = {
    "L9": ["HJC", "SGT", "LT", "MEC", "LEC", "PMC", "PLC", "DMC", "DLC"],
    "L17": ["HJC", "SNI", "ANI", "INI", "PNI", "SGT", "LT", "PITC", "PLC", "PMC", "PPLC", "PPMC", "LEC", "MEC",
            "ICN", "DLC", "DMC"],
}
SEED = 20260924

# ------------------------------------------------------------------ landmarks
def feat(Ld, names):
    """(m,3) landmarks in names order; HJC = geometric sphere through 6 femoral head points (as in morph.subject_targets)."""
    out = []
    for k in names:
        if k == "HJC":
            out.append(joints.sphere_fit(np.array([Ld[h] for h in HEAD6]))["center"])
        else:
            out.append(np.asarray(Ld[k], float))
    return np.array(out, float)


def feat_shape(V, names):
    return feat(E.landmarks(V), names)


def frame_idx(names):
    return names.index("HJC"), names.index("MEC"), names.index("LEC")


def frame_of(Y, names):
    i, m, l = frame_idx(names)
    return joints.femur_frame(Y[i], Y[m], Y[l])[1]      # columns = axes


def rot_block(Nloc, R):
    m = Nloc.shape[0] // 3
    B = np.kron(np.eye(m), R)
    return B @ Nloc @ B.T


def wkabsch(A, B, w=None):
    """R, t such that A @ R.T + t ≈ B (weighted, without scale)."""
    w = np.ones(len(A)) if w is None else np.asarray(w, float)
    w = w / w.sum()
    ma, mb = w @ A, w @ B
    H = ((B - mb) * w[:, None]).T @ (A - ma)
    U, _, Vt = np.linalg.svd(H)
    D = np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))])
    R = U @ D @ Vt
    return R, mb - R @ ma


# ------------------------------------------------------------------ VSD ratings
SUBJ, SIDE, ARR = P.load_vsd_raters()      # (19,22,20,3), mirrored to the right
N_RATER = 5


def ok_ratings(s):
    return [r for r in range(ARR.shape[2]) if not np.isnan(ARR[s, :, r]).any()]


def rating_dict(s, r):
    return {n: ARR[s, l, r] for l, n in enumerate(LM_NAMES)}


def mean_dict(s, rs=None):
    rs = ok_ratings(s) if rs is None else rs
    return {n: ARR[s, l, rs].mean(0) for l, n in enumerate(LM_NAMES)}


def m4_draws(s, n=10):
    """10 draws: 4 different raters, one random valid trial each (index r = rater*4 + trial)."""
    rng = np.random.default_rng(SEED + s)
    ok = set(ok_ratings(s))
    out = []
    for _ in range(n):
        raters = rng.choice(N_RATER, 4, replace=False)
        rs = []
        for q in raters:
            cand = [q * 4 + t for t in range(4) if q * 4 + t in ok]
            rs.append(int(rng.choice(cand)))
        out.append(sorted(rs))
    return out


def noise_cov(s, names, cond="S1"):
    """Second moment (in the local femur frame) of (input − mean-20) over the other 18 individuals.
    S1: each individual rating (corrected n/(n−1), as in RM1). M4: 10 M4 draws per individual. M20: 0."""
    d = len(names) * 3
    if cond == "M20":
        return np.zeros((d, d))
    dev = []
    for j in range(len(SUBJ)):
        if j == s:
            continue
        ok = ok_ratings(j)
        m = feat(mean_dict(j), names)
        R = frame_of(m, names)
        if cond == "S1":
            n = len(ok)
            for r in ok:
                dev.append(((feat(rating_dict(j, r), names) - m) @ R).ravel() * np.sqrt(n / (n - 1)))
        else:
            for rs in m4_draws(j):
                dev.append(((feat(mean_dict(j, rs), names) - m) @ R).ravel())
    dev = np.array(dev)
    return dev.T @ dev / len(dev)


_RG = P.RG
_RKEYS = P.RKEYS


def vsd_reg(sid):
    return _RG["Xp"][_RKEYS.index(f"vsd:{sid}")].astype(float)


def def_cov(s, names):
    """Definition deviation: mean-20 minus template landmarks on P1's registration of CT, other 18, block 3×3."""
    dev = []
    for j in range(len(SUBJ)):
        if j == s:
            continue
        m = feat(mean_dict(j), names)
        t = feat_shape(vsd_reg(SUBJ[j]), names)
        R = frame_of(m, names)
        dev.append((m - t) @ R)
    dev = np.array(dev)                      # (18,m,3)
    C = np.zeros((len(names) * 3,) * 2)
    for i in range(len(names)):
        C[3 * i:3 * i + 3, 3 * i:3 * i + 3] = dev[:, i].T @ dev[:, i] / len(dev)
    return C


# ------------------------------------------------------------------ TLEM sources
XT = P.XT.astype(float)                                   # template registered to TLEM (TLEM frame)
HIPJOINT = np.asarray(json.loads((RES / "geometry_inputs.json").read_text())["hip_joint_mm"], dtype=float)
_TAB = dict(morph.tlem_landmarks()[1])
_TMPL = E.landmarks(XT)
TLEM_MESH = tlem.load_mesh("femur")
TLEM_V = np.asarray(TLEM_MESH.vertices, float)
TLEM_F = np.asarray(TLEM_MESH.faces)


def tlem_src(source, names, hjc="hipjoint"):
    if source == "tab":
        L = dict(_TAB)
        if any(k not in L for k in names):
            return None
    else:
        L = dict(_TMPL)
    out = []
    for k in names:
        if k == "HJC":
            out.append(HIPJOINT if hjc == "hipjoint" else _TAB["HJC"])
        else:
            out.append(np.asarray(L[k], float))
    return np.array(out, float)


class RBF:
    """RBF + polynomial of degree 1. kernel 'r' (= TPS3, biharmonic 3D) or 'r2logr'."""

    def __init__(self, src, dst, kernel):
        self.c = np.asarray(src, float)
        self.k = kernel
        k = len(src)
        self.mu = self.c.mean(0)
        K = self._phi(cdist(self.c, self.c))
        Pm = np.c_[np.ones(k), self.c - self.mu]
        L = np.zeros((k + 4, k + 4))
        L[:k, :k] = K
        L[:k, k:] = Pm
        L[k:, :k] = Pm.T
        rhs = np.r_[np.asarray(dst, float) - self.c, np.zeros((4, 3))]
        self.w = np.linalg.solve(L, rhs)

    def _phi(self, r):
        if self.k == "r":
            return r
        rm = r / 1000.0                          # The reference formulation uses meters; r² log r is not scale invariant
        with np.errstate(divide="ignore", invalid="ignore"):
            v = rm ** 2 * np.log(rm)
        return np.nan_to_num(v)

    def __call__(self, X):
        X = np.atleast_2d(np.asarray(X, float))
        out = np.empty_like(X)
        for a in range(0, len(X), 20000):
            Xa = X[a:a + 20000]
            B = np.c_[self._phi(cdist(Xa, self.c)), np.ones(len(Xa)), Xa - self.mu]
            out[a:a + 20000] = Xa + B @ self.w
        return out


def ab_map(kind, src, dst, names):
    if kind == "AFF":
        M = fit_affine(src, dst)
        return lambda X: apply_affine(M, np.atleast_2d(X))
    if kind == "AFFH":
        h = names.index("HJC")
        A0, B0 = src - src[h], dst - dst[h]
        Lm, *_ = np.linalg.lstsq(A0, B0, rcond=None)       # (3,3): x_rel @ Lm
        return lambda X: dst[h] + (np.atleast_2d(X) - src[h]) @ Lm
    if kind == "RBF1":
        return RBF(src, dst, "r")
    if kind == "RBF2":
        return RBF(src, dst, "r2logr")
    raise ValueError(kind)


AB_KINDS = ["AFF", "AFFH", "RBF1", "RBF2"]
SRC_TLEM_KNEE = P.KNEE_T
TAB_DISTAL = np.array([_TAB[k] for k in P.DISTAL])
TMPL_DISTAL = np.array([_TMPL[k] for k in P.DISTAL])


def ab_arm(kind, source, y, names):
    """landmark-based morph of TLEM to the landmarks y. Returns (fem-dict for P.quantities, surface mesh vertices (TLEM topology),
    template correspondence vertices (XT through the same mapping, for rigid fitting in E_form))."""
    src = tlem_src(source, names)
    if src is None:
        return None
    M = ab_map(kind, src, y, names)
    Ls = _TAB if source == "tab" else _TMPL
    hjc = M(HIPJOINT)[0]
    mec, lec = M(np.asarray(Ls["MEC"]))[0], M(np.asarray(Ls["LEC"]))[0]
    distal = M(TAB_DISTAL if source == "tab" else TMPL_DISTAL)
    Xm = M(XT)
    shaft = measures(E.landmarks(Xm), Xm)["_shaft"]
    fem = dict(pts=dict(zip(P.FEM_IDS, M(P.FEM_XYZ))), hjc=hjc, mec=mec, lec=lec, distal=distal,
               knee=M(SRC_TLEM_KNEE), shaft=shaft, n_sphere=-1)
    return fem, M(TLEM_V), Xm


# ------------------------------------------------------------------ shape model (PSM)
class PSM:
    """Posterior shape model. Prior: PCA of GPA-fitted shapes (without scale). Landmark function linearized per mode."""

    def __init__(self, shapes, names, perm_seed=None):
        self.names = names
        A, _ = E.gpa(np.asarray(shapes, float))
        self.n = len(A)
        self.mu = A.mean(0)
        X = (A - self.mu).reshape(self.n, -1)
        U, s, Vt = np.linalg.svd(X, full_matrices=False)
        r = self.n - 1
        self.V = Vt[:r].T
        self.lam = s[:r] ** 2 / (self.n - 1)
        self.scores = U[:, :r] * s[:r]
        self.mu_f = feat_shape(self.mu, names).ravel()
        if perm_seed is None:
            W = np.empty((len(self.mu_f), r))
            for k in range(r):
                sg = np.sqrt(self.lam[k])
                vk = self.V[:, k].reshape(-1, 3)
                W[:, k] = (feat_shape(self.mu + sg * vk, names).ravel() - feat_shape(self.mu - sg * vk, names).ravel()) / (2 * sg)
            self.W = W
        else:                                   # countertest: regression of (permuted) landmarks on shape scores
            Fm = np.array([feat_shape(a, names).ravel() for a in A])
            Fm = Fm[np.random.default_rng(perm_seed).permutation(len(Fm))]
            Fc = Fm - Fm.mean(0)
            self.mu_f = Fm.mean(0)
            self.W = np.linalg.lstsq(self.scores, Fc, rcond=None)[0].T
        self.r = r

    def condition(self, y, Nloc, Soos_loc, K, iters=200, tol=0.01):
        """y (m,3) in the CT frame. Nloc: rater+definition noise, Soos_loc: model incompleteness (both in the local femur frame).
        Returns dict with the posterior of all modes (for MC) and K-truncated mean shape, residual field ĝ, pose."""
        W, lam, mu_f = self.W, self.lam, self.mu_f
        G_all = (W * lam) @ W.T
        m = len(self.names)
        w = np.ones(m)
        R, t = wkabsch(y, mu_f.reshape(-1, 3), w)
        fhat = mu_f.reshape(-1, 3)
        step = np.inf
        for it in range(iters):
            ym = y @ R.T + t
            Rf = frame_of(fhat, self.names)
            Nm = rot_block(Nloc, Rf)
            Om = rot_block(Soos_loc, Rf)
            S = G_all + Om + Nm
            alpha = np.linalg.solve(S, ym.ravel() - mu_f)
            fnew = (mu_f + (G_all + Om) @ alpha).reshape(-1, 3)
            Cd = np.array([np.trace((Om + Nm)[3 * i:3 * i + 3, 3 * i:3 * i + 3]) for i in range(m)])
            w = 1.0 / np.maximum(Cd, 1e-6)
            R, t = wkabsch(y, fnew, w)
            step = float(np.abs(fnew - fhat).max())
            fhat = fnew
            if step < tol and it > 2:
                break
        ym = y @ R.T + t
        Rf = frame_of(fhat, self.names)
        Nm, Om = rot_block(Nloc, Rf), rot_block(Soos_loc, Rf)
        S = G_all + Om + Nm
        Si = np.linalg.inv(S)
        alpha = Si @ (ym.ravel() - mu_f)
        b = lam * (W.T @ alpha)
        g = Om @ alpha
        # joint posterior for z = (b_all, g): P = diag(Λ, Om), H = [W, I]
        Pz = np.zeros((self.r + len(mu_f),) * 2)
        Pz[:self.r, :self.r] = np.diag(lam)
        Pz[self.r:, self.r:] = Om
        H = np.c_[W, np.eye(len(mu_f))]
        PH = Pz @ H.T
        Cz = Pz - PH @ Si @ PH.T
        Cz = 0.5 * (Cz + Cz.T)
        return dict(R=R, t=t, b=b, g=g, alpha=alpha, Cz=Cz, K=K, n_iter=it + 1, last_step=step,
                    resid_mm=float(np.sqrt(np.mean(np.sum((ym - fhat) ** 2, 1)))))

    def shape(self, post, K, resid=False, b=None, g=None):
        """Shape in the CT frame. K-truncated mean shape (+ residual-TPS if resid)."""
        b = post["b"] if b is None else b
        Sm = self.mu + (self.V[:, :K] @ b[:K]).reshape(-1, 3)
        if resid:
            g = post["g"] if g is None else g
            g = g + self.W[:, K:] @ b[K:]          # Σ_g = Σ_trunk + Σ_oos: landmark part of the truncated modes + incompleteness
            f_base = (self.mu_f + self.W[:, :K] @ b[:K]).reshape(-1, 3)
            Sm = TPS3(f_base).displace(Sm, g.reshape(-1, 3))
        return (Sm - post["t"]) @ post["R"]

    def draws(self, post, n, rng):
        L = np.linalg.cholesky(post["Cz"] + 1e-9 * np.eye(len(post["Cz"])))
        mean = np.r_[post["b"], post["g"]]
        Z = mean + rng.standard_normal((n, len(mean))) @ L.T
        return [(z[:self.r], z[self.r:]) for z in Z]


def psm_fem(Sct):
    """Femur side (P3 arm A rule) for a shape in template correspondence, CT frame."""
    return P.femur_from_template(Sct)


# ------------------------------------------------------------------ Σ_oos: model incompleteness (Imperial-LOSO)
def soos_imperial(names):
    f = TMP / f"soos_{'_'.join(names[:2])}_{len(names)}.npz"
    if f.exists():
        z = np.load(f)
        return z["C"], z["vvar"]
    Y = E.Y
    subj = E.SUBJ
    devs, vres = [], []
    for p in sorted(set(subj)):
        tr = [i for i, s in enumerate(subj) if s != p]
        te = [i for i, s in enumerate(subj) if s == p]
        A, _ = E.gpa(Y[tr])
        mu = A.mean(0)
        X = (A - mu).reshape(len(A), -1)
        _, _, Vt = np.linalg.svd(X, full_matrices=False)
        V = Vt[:len(A) - 1].T
        for i in te:
            s_, Rm, tt = E.umeyama(Y[i], mu, scale=False)
            Yi = E.apply(s_, Rm, tt, Y[i])
            for _ in range(3):                    # rigid fitting to the projection (iterate)
                proj = mu + (V @ (V.T @ (Yi - mu).ravel())).reshape(-1, 3)
                s_, Rm, tt = E.umeyama(Yi, proj, scale=False)
                Yi = E.apply(s_, Rm, tt, Yi)
            proj = mu + (V @ (V.T @ (Yi - mu).ravel())).reshape(-1, 3)
            ft, fp = feat_shape(Yi, names), feat_shape(proj, names)
            Rf = frame_of(ft, names)
            devs.append((ft - fp) @ Rf)
            nrm = trimesh.Trimesh(proj, F, process=False).vertex_normals
            vres.append(np.einsum("ij,ij->i", Yi - proj, nrm))
    devs = np.array(devs)
    C = np.zeros((len(names) * 3,) * 2)
    for i in range(len(names)):
        C[3 * i:3 * i + 3, 3 * i:3 * i + 3] = devs[:, i].T @ devs[:, i] / len(devs)
    vvar = np.mean(np.array(vres) ** 2, 0)
    np.savez(f, C=C, vvar=vvar)
    return C, vvar


# ------------------------------------------------------------------ error measurement
def samples(mesh_or_VF, n, seed=0):
    if isinstance(mesh_or_VF, tuple):
        m = trimesh.Trimesh(mesh_or_VF[0], mesh_or_VF[1], process=False)
    else:
        m = mesh_or_VF
    pts, _ = trimesh.sample.sample_surface(m, n, seed=seed)
    return pts


class Target:
    """Raw surface (CT or MR) with dense samples and KD trees."""

    def __init__(self, mesh, key, n=200000):
        f = TMP / f"tgt_{key}.npy"
        if f.exists():
            self.pts = np.load(f)
        else:
            self.pts = samples(mesh, n, seed=1)
            np.save(f, self.pts)
        self.tree = cKDTree(self.pts)

    def sym_rms_pts(self, ps):
        d1, _ = self.tree.query(ps)
        d2, _ = cKDTree(ps).query(self.pts[::2])
        return float(np.sqrt(0.5 * (np.mean(d1 ** 2) + np.mean(d2 ** 2))))


def eval_surface(tgt, V, Fc, Xcorr, Xtrue, n=100000):
    """E_yta (as placed) and E_form (after rigid fitting Xcorr -> Xtrue via template correspondence).
    The surface is sampled once; E_form moves the samples rigidly (same as sampling the moved surface)."""
    ps = samples((V, Fc), n, seed=2)
    e_placed = tgt.sym_rms_pts(ps)
    R, t = wkabsch(Xcorr, Xtrue)
    e_shape = tgt.sym_rms_pts(ps @ R.T + t)
    return e_placed, e_shape


def qsumm(q):
    return dict(r_ab=float(q["r_ab"]), R_BW=float(q["R_BW"]), gi_hjc=[float(v) for v in q["gi_hjc"]])


def att_err(fem, femT):
    d = [np.linalg.norm(np.asarray(fem["pts"][i]) - np.asarray(femT["pts"][i])) for i in P.FEM_IDS]
    return float(np.median(d))


# ------------------------------------------------------------------ precomputed truths (morph.build requires C extensions; locally)
import pickle  # noqa: E402


def truth(sid):
    """P3.truth(sid) precomputed locally by n1_precompute.py (pel, T1, T2, floor)."""
    with open(PRE / f"truth_{sid}.pkl", "rb") as f:
        return pickle.load(f)


def pel_identity():
    with open(PRE / "pel_identity.pkl", "rb") as f:
        return pickle.load(f)
