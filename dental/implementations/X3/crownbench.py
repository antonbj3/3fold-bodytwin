"""Geometry primitives shared by public generation and private evaluation.
No function in this module loads a hidden surface. Unit is provisionally mm.
"""
from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import pathlib, json, hashlib, datetime, time
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import minimize, least_squares
from scipy.ndimage import distance_transform_edt
import trimesh
from skimage.measure import marching_cubes
P = pathlib.Path(__file__).resolve().parent
DATA = pathlib.Path(os.environ.get('X3_DATA_DIR', _release_expand('@DENTAL_WORK_ROOT@/X3_crown_benchmark')))
H = 0.2
LO = np.array([-9.0, -9.0, -10.0])
SHAPE = (91, 91, 101)
TARGETS = (2, 3, 5, 6)

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        while (b := f.read(1048576)):
            h.update(b)
    return h.hexdigest()

def put(p, x):
    p = pathlib.Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(x, indent=2, allow_nan=False, ensure_ascii=False) + '\n')
    tmp.replace(p)

def state(phase, last, nxt, **extra):
    j = json.load(open(P / 'CURRENT_WORK_STATE.json'))
    j.update(phase=phase, last_gate=last, next_operation=nxt, updated_utc=now(), **extra)
    put(P / 'CURRENT_WORK_STATE.json', j)

def disk_guard():
    total = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    if total > 2900000000:
        raise RuntimeError('INTERMEDIATE_LIMIT_3GB imminent')
    return total

def mesh(v, f):
    return trimesh.Trimesh(vertices=v, faces=f, process=False)

def load_obj(p):
    vs = []
    fs = []
    with open(p) as stream:
        for line in stream:
            if line.startswith('v '):
                vs.append([float(x) for x in line.split()[1:4]])
            elif line.startswith('f '):
                ids = [int(x.split('/')[0]) - 1 for x in line.split()[1:]]
                if len(ids) != 3:
                    raise ValueError('NON_TRIANGULAR_FACE')
                fs.append(ids)
    return (np.asarray(vs, dtype=np.float64), np.asarray(fs, dtype=np.int32))

def select(v, f, labels, k):
    keep = labels == k
    inds = np.flatnonzero(keep)
    lut = np.full(len(v), -1, dtype=np.int32)
    lut[inds] = np.arange(len(inds))
    faces = f[np.all(keep[f], axis=1)]
    return (v[inds], lut[faces])

def mask(v, f, labels, k):
    keep = labels != k
    inds = np.flatnonzero(keep)
    lut = np.full(len(v), -1, dtype=np.int32)
    lut[inds] = np.arange(len(inds))
    faces = f[np.all(keep[f], axis=1)]
    return (v[inds], lut[faces], labels[inds])

def swap(k):
    return {1: 2, 2: 1, 3: 4, 4: 3}[k // 10] * 10 + k % 10

def unit(v):
    n = np.linalg.norm(v)
    if n < 1e-10:
        raise ValueError('DEGENERATE_DIRECTION')
    return v / n

def boundary(v, f):
    e = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    (_, idx, cnt) = np.unique(np.sort(e, axis=1), axis=0, return_index=True, return_counts=True)
    edges = e[idx[cnt == 1]]
    return (edges, np.unique(edges))

def close(v, f):
    (e, ids) = boundary(v, f)
    if len(e) == 0:
        return (v, f, ids)
    adj = {}
    for (a, b) in e:
        adj.setdefault(int(a), set()).add(int(b))
        adj.setdefault(int(b), set()).add(int(a))
    seen = set()
    vs = [*v]
    fs = [*f]
    for root in adj:
        if root in seen:
            continue
        stack = [root]
        component = set()
        while stack:
            u = stack.pop()
            if u in component:
                continue
            component.add(u)
            stack.extend(adj[u] - component)
        seen |= component
        center = np.mean(v[list(component)], axis=0)
        ci = len(vs)
        vs.append(center)
        for (a, b) in e:
            if int(a) in component:
                fs.append([b, a, ci])
    return (np.asarray(vs), np.asarray(fs, dtype=np.int32), ids)

def tooth_dict(v, f, labels):
    return {int(k): select(v, f, labels, k) for k in np.unique(labels) if k > 0 and (labels == k).sum() >= 300}

def axis(teeth):
    cents = np.array([a.mean(0) for (a, b) in teeth.values()])
    (_, _, vh) = np.linalg.svd(cents - cents.mean(0), full_matrices=False)
    n = vh[-1]
    diff = []
    for (v, f) in teeth.values():
        (_, ids) = boundary(v, f)
        if len(ids):
            diff.append((v.mean(0) - v[ids].mean(0)) @ n)
    if diff and np.median(diff) < 0:
        n = -n
    return n

def reflection(teeth, q):
    pairs = [(q * 10 + i, (q + 1) * 10 + i) for i in range(1, 9) if q * 10 + i in teeth and (q + 1) * 10 + i in teeth]
    if len(pairs) < 3:
        raise ValueError('FEWER_THAN3_VISIBLE_HOMOLOG_PAIRS')
    L = np.array([teeth[a][0].mean(0) for (a, b) in pairs])
    R = np.array([teeth[b][0].mean(0) for (a, b) in pairs])
    c = ((L + R) / 2).mean(0)
    n0 = unit((L - R).mean(0))

    def fun(n):
        n = unit(n)
        pred = R - 2 * ((R - c) @ n)[:, None] * n
        return float(np.mean(np.sum((pred - L) ** 2, axis=1)))
    res = minimize(fun, n0, method='SLSQP', constraints={'type': 'eq', 'fun': lambda n: n @ n - 1}, options={'maxiter': 200, 'ftol': 1e-10})
    if not res.success:
        raise ValueError('REFLECTION_FIT_FAILED:' + res.message)
    n = unit(res.x)
    return (n, c, {'pairs': pairs, 'rms_mm': float(np.sqrt(fun(n))), 'normal': n.tolist(), 'origin': c.tolist()})

def reflect(v, n, c):
    return v - 2 * ((v - c) @ n)[:, None] * n

def local_frame(teeth, k, n):
    i = k % 10
    q = k // 10
    if q * 10 + i - 1 in teeth and q * 10 + i + 1 in teeth:
        d = teeth[q * 10 + i + 1][0].mean(0) - teeth[q * 10 + i - 1][0].mean(0)
    elif i == 1 and q * 10 + 2 in teeth and (swap(k) in teeth):
        d = teeth[q * 10 + 2][0].mean(0) - teeth[swap(k)][0].mean(0)
    elif q * 10 + i - 1 in teeth:
        d = teeth[k][0].mean(0) - teeth[q * 10 + i - 1][0].mean(0)
    else:
        raise ValueError('NO_TANGENT_NEIGHBOR')
    x = unit(d - n * (d @ n))
    y = unit(np.cross(n, x))
    return np.column_stack([x, y, n])

def sampled(v, f, count, seed):
    rng = np.random.default_rng(seed)
    t = v[f]
    areas = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2
    if not (len(f) and areas.sum() > 0):
        raise ValueError('EMPTY_SURFACE')
    ii = rng.choice(len(f), size=count, p=areas / areas.sum())
    r = np.sqrt(rng.random(count))
    s = rng.random(count)
    return (1 - r)[:, None] * t[ii, 0] + (r * (1 - s))[:, None] * t[ii, 1] + (r * s)[:, None] * t[ii, 2]

def seed(key):
    return int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)

def distances(points, v, f):
    m = mesh(v, f)
    parts = []
    for i in range(0, len(points), 256):
        parts.append(trimesh.proximity.closest_point(m, points[i:i + 256])[1])
    return np.concatenate(parts)

def signed_grid(v, f):
    (cv, cf, _) = close(v, f)
    vox = mesh(cv, cf).voxelized(pitch=H).fill()
    xyz = vox.points
    ids = np.rint((xyz - LO) / H).astype(np.int32)
    if np.any(ids < 0) or np.any(ids >= np.array(SHAPE)):
        raise ValueError('SDF_GRID_OVERFLOW')
    occ = np.zeros(SHAPE, dtype=bool)
    occ[tuple(ids.T)] = True
    phi = (distance_transform_edt(~occ) - distance_transform_edt(occ)) * H
    return phi.astype(np.float32)

def zero_mesh(phi):
    (v, f, _, _) = marching_cubes(phi, level=0, spacing=(H, H, H))
    return (v + LO, f.astype(np.int32))

def sdf_roundtrip(v, f, F, center):
    local = (v - center) @ F
    phi = signed_grid(local, f)
    (lv, lf) = zero_mesh(phi)
    centers = lv[lf].mean(1)
    keep = distances(centers, local, f) <= 1.5 * H
    lf = lf[keep]
    return (lv @ F.T + center, lf, phi, {'cap_triangles_removed': int((~keep).sum()), 'grid_h_mm': H, 'occupancy_convention': 'negative inside; voxel center +/-h/2 biased zero surface'})

def neighbor_fit(v, f, teeth, k, F):
    c = v.mean(0)
    loc = (v - c) @ F
    sp = sampled(v, f, 768, 6103)
    sl = (sp - c) @ F
    nn = [k - 1, k + 1]
    trees = [cKDTree(teeth[z][0]) for z in nn]
    source_nn = [swap(z) for z in nn]
    target_height = np.mean([np.quantile(teeth[z][0] @ F[:, 2], 0.85) for z in nn])
    if all((z in teeth for z in source_nn)):
        source_height = np.mean([np.quantile(teeth[z][0] @ F[:, 2], 0.85) for z in source_nn])
        height = np.quantile(v @ F[:, 2], 0.85) + target_height - source_height
    else:
        height = np.quantile(v @ F[:, 2], 0.85)

    def transformed(p):
        z = sl.copy()
        z[:, 0] *= p[3]
        return z @ F.T + c + F @ p[:3]

    def residual(p):
        z = transformed(p)
        gap = [(np.quantile(t.query(z, workers=1)[0], 0.02) - 0.1) / 0.2 for t in trees]
        h = (np.quantile(z @ F[:, 2], 0.85) - height) / 0.5
        return np.r_[gap, h, 0.35 * np.asarray(p[:3]), 0.35 * (p[3] - 1) / 0.2]
    res = least_squares(residual, [0, 0, 0, 1], bounds=([-2, -2, -2, 0.75], [2, 2, 2, 1.25]), max_nfev=60, diff_step=0.001, ftol=1e-06, xtol=1e-06, gtol=1e-06)
    out = loc.copy()
    out[:, 0] *= res.x[3]
    out = out @ F.T + c + F @ res.x[:3]
    return (out, f, {'params': res.x.tolist(), 'cost': float(res.cost), 'nfev': int(res.nfev), 'success': bool(res.success), 'height_objective_mm': float(height), 'residual': residual(res.x).tolist(), 'fit_seconds': None})

def write_surface(p, v, f):
    p = pathlib.Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(p, vertices=v.astype(np.float64), faces=f.astype(np.int32))

def read_surface(p):
    with np.load(p) as d:
        return (d['vertices'], d['faces'])
