import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from geometry import *

def raw_orientation(patient, t):
    xyz = [90, 0, 0]
    if t == 8:
        xyz = [180, 0, 0] if patient == '3485' else [0, 90, 180]
    if patient == '6457' and t == 4:
        xyz = [-100, 180, 40]
    Rx = Rotation.from_euler('x', xyz[0], degrees=True).as_matrix()
    Ry = Rotation.from_euler('y', xyz[1], degrees=True).as_matrix()
    Rz = Rotation.from_euler('z', xyz[2], degrees=True).as_matrix()
    return mat((Rz @ Ry @ Rx).T)

def load_full(z, p, t):
    cache = CACHE / f'{p}_full_upper_T{t}_100k.npz'
    name = f'{p}_Clinical_Trial/Sirona/T{t}/{p}_OnyxCeph3_Export_OK.stl'
    if not cache.exists():
        blob = z.read(name)
        H = initial_frame(z, p, 'Sirona', t, 'OK') @ raw_orientation(p, t)
        (pts, ns, c, area) = sample_mesh(blob, H, int(digest(name.encode())[:8], 16), 100000)
        np.savez_compressed(cache, points=pts.astype('float32'), normals=ns.astype('float32'), source_sha256=np.array(digest(blob)))
    d = np.load(cache)
    return (d['points'], d['normals'], str(d['source_sha256']))

def roi(pts, normals, crowns):
    (cp, cn) = union(crowns, sorted(crowns), 4000)
    (dd, _) = cKDTree(cp).query(pts, workers=4)
    zi = np.mean([c['centroid'][2] for (f, c) in crowns.items() if f % 10 == 1])
    zm = np.mean([c['centroid'][2] for (f, c) in crowns.items() if f % 10 == 6])
    lo = min(zm + 3, zi - 5)
    hi = max(zm + 3, zi - 5)
    m = (np.abs(pts[:, 0]) <= 8) & (pts[:, 2] >= lo) & (pts[:, 2] <= hi) & (dd >= 3)
    r = pts[m]
    n = normals[m]
    ext = np.ptp(r, axis=0) if len(r) else np.zeros(3)
    ok = len(r) >= 500 and ext[0] >= 8 and (ext[2] >= 8)
    return (r, n, {'points': len(r), 'extent_mm': ext.tolist(), 'z_bounds_mm': [lo, hi], 'coverage_pass': bool(ok), 'resolution': 'PER_SURFACE_REGION'})

def heldout(H, pts, target, normals):
    v = apply(H, pts)
    (d, i) = cKDTree(target).query(v, workers=1)
    mask = d <= np.quantile(d, 0.9)
    res = np.einsum('ij,ij->i', normals[i[mask]], v[mask] - target[i[mask]])
    return float(np.sqrt(np.mean(res ** 2)))
