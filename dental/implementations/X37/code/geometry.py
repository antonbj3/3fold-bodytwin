"""Zip-streamed surface measurements. Matrices map column points; arrays store rows."""
from dental_release.paths import expand as _release_expand
import hashlib, io, re, zipfile
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
import trimesh
DATA = Path(_release_expand('@DENTAL_DATA_ROOT@/geometry/AlignerMovement_Zenodo11280343'))
CACHE = Path(_release_expand('@DENTAL_WORK_ROOT@/X37_ALIGNER_MOVEMENT'))
ROOT = Path(__file__).resolve().parents[1]
MANUAL = {('3485', 'OK'): [[-9.8421, -8.2143, 26.0685], [-27.851, -6.6035, -15.414], [25.8475, -3.8557, -2.5125]], ('3485', 'UK'): [[-7.3285, -7.6276, 25.0078], [26.3526, -9.18039, -1.72466], [-24.716, -6.17749, -13.3343]], ('6457', 'OK'): [[-1.9318, -9.1141, 26.3856], [-32.5699, -10.399, -13.6267], [30.3429, -10.2828, -14.3173]], ('6457', 'UK'): [[-1.4049, -10.1022, 20.1666], [29.6122, -6.189, -17.008], [-28.1837, -8.2212, -19.2044]]}

def digest(b):
    return hashlib.sha256(b).hexdigest()

def apply(H, p):
    return p @ H[:3, :3].T + H[:3, 3]

def mat(R=np.eye(3), t=np.zeros(3)):
    H = np.eye(4)
    H[:3, :3] = R
    H[:3, 3] = t
    return H

def frame(ps):
    (p1, p2, p3) = np.asarray(ps)
    o = (p1 + p2 + p3) / 3
    x = p3 - p2
    x /= np.linalg.norm(x)
    if x[0] < 0:
        x = -x
    y = np.cross(p2 - p1, p3 - p1)
    y /= np.linalg.norm(y)
    if y[1] < 0:
        y = -y
    z = np.cross(x, y)
    z /= np.linalg.norm(z)
    y = np.cross(z, x)
    R = np.array([x, y, z])
    return mat(R, -R @ o)

def member(p, src, t, fdi):
    seg = 'teethSeg' if src == 'Sirona' and t == 9 else 'teethSeg_noSnap_TolS'
    return f'{p}_Clinical_Trial/{src}/T{t}/{seg}/{p}_Z{fdi}.stl'

def initial_frame(z, p, src, t, jaw):
    if src == 'Sirona' and t == 9:
        return frame(MANUAL[p, jaw])
    name = f'{p}_Clinical_Trial/{src}/T{t}/teethSeg_noSnap_TolS/{p}_{jaw}_occPlane.mat'
    d = loadmat(io.BytesIO(z.read(name)), simplify_cells=True)['occPlaneState']
    return frame([d[k] for k in ['P1', 'P2', 'P3']])

def sample_mesh(blob, H, seed, n=16000):
    m = trimesh.load(io.BytesIO(blob), file_type='stl', process=False)
    tri = apply(H, m.vertices)[m.faces]
    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    norm = np.linalg.norm(cross, axis=1)
    keep = norm > 1e-12
    tri = tri[keep]
    cross = cross[keep]
    norm = norm[keep]
    area = norm / 2
    ns = cross / norm[:, None]
    rng = np.random.default_rng(seed)
    ids = rng.choice(len(tri), size=n, p=area / area.sum())
    ab = rng.random((n, 2))
    ab = np.sqrt(ab[:, 0])[:, None] * np.c_[1 - ab[:, 1], ab[:, 1]]
    pts = tri[ids, 0] + ab[:, 0, None] * (tri[ids, 1] - tri[ids, 0]) + ab[:, 1, None] * (tri[ids, 2] - tri[ids, 0])
    c = np.average(tri.mean(1), axis=0, weights=area)
    return (pts, ns[ids], c, area.sum())

def load_cloud(z, p, src, t, fdi):
    n = member(p, src, t, fdi)
    path = CACHE / f'{p}_{src}_T{t}_Z{fdi}.npz'
    if path.exists():
        d = np.load(path)
        return {k: d[k] for k in d.files}
    b = z.read(n)
    jaw = 'OK' if fdi // 10 in [1, 2] else 'UK'
    H = initial_frame(z, p, src, t, jaw)
    if src == 'Sirona' and t == 9:
        C = Rotation.from_euler('x', -90, degrees=True).as_matrix()
        if jaw == 'UK':
            C = np.diag([1, -1, 1]) @ C
        H = H @ mat(C)
    seed = int(digest(n.encode())[:8], 16)
    (pts, ns, c, area) = sample_mesh(b, H, seed)
    out = {'points': pts.astype('float32'), 'normals': ns.astype('float32'), 'centroid': c, 'area_mm2': np.array(area), 'source_sha256': np.array(digest(b)), 'frame': H}
    np.savez_compressed(path, **out)
    return out

def kabsch(a, b):
    ac = a.mean(0)
    bc = b.mean(0)
    (u, s, vt) = np.linalg.svd((a - ac).T @ (b - bc))
    R = vt.T @ u.T
    if np.linalg.det(R) < 0:
        vt[-1] *= -1
        R = vt.T @ u.T
    return mat(R, bc - R @ ac)

def icp(a, b, normals=None, initial=None, method='plane', iterations=60, trim=0.9):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    H = np.eye(4) if initial is None else initial.copy()
    tree = cKDTree(b)
    last = np.inf
    for i in range(iterations):
        v = apply(H, a)
        (d, idx) = tree.query(v, workers=1)
        mask = d <= np.quantile(d, trim)
        v0 = v[mask]
        q = b[idx[mask]]
        if method == 'point':
            inc = kabsch(v0, q)
        else:
            ns = np.asarray(normals)[idx[mask]]
            rr = np.einsum('ij,ij->i', ns, q - v0)
            center = v0.mean(0)
            A = np.c_[np.cross(v0 - center, ns), ns]
            sol = np.linalg.lstsq(A, rr, rcond=1e-09)[0]
            rot = sol[:3]
            theta = np.linalg.norm(rot)
            if theta > 0.15:
                rot *= 0.15 / theta
            tr = sol[3:]
            R = Rotation.from_rotvec(rot).as_matrix()
            inc = mat(R, center + tr - R @ center)
        H = inc @ H
        loss = float(np.mean(d[mask] ** 2))
        if abs(last - loss) < 1e-10:
            break
        last = loss
    v = apply(H, a)
    (d, idx) = tree.query(v, workers=1)
    mask = d <= np.quantile(d, trim)
    plane = float(np.sqrt(np.mean(np.einsum('ij,ij->i', np.asarray(normals)[idx[mask]], v[mask] - b[idx[mask]]) ** 2))) if normals is not None else None
    return (H, {'iterations': i + 1, 'trimmed_point_rms_mm': float(np.sqrt(np.mean(d[mask] ** 2))), 'trimmed_plane_rms_mm': plane, 'inlier_fraction': float(mask.mean())})

def angle(H):
    return float(np.degrees(Rotation.from_matrix(H[:3, :3]).magnitude()))

def union(clouds, fdis, n=1000):
    return (np.concatenate([clouds[k]['points'][:n] for k in fdis]), np.concatenate([clouds[k]['normals'][:n] for k in fdis]))
