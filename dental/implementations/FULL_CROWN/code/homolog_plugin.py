"""Public-only full-surface generator. No original target or scores are accepted."""
import os, sys, json, time, resource
from pathlib import Path
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from skimage.measure import marching_cubes
import trimesh
sys.path.insert(0, '/deps')
import mesh2sdf

def sample(tri, n=600):
    a = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    if not len(tri) or a.sum() <= 0:
        raise ValueError('empty observed surface')
    i = np.searchsorted(np.cumsum(a), (np.arange(n) + 0.5) * a.sum() / n)
    return tri[i].mean(1)

def main_component(tri):
    m = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=True)
    parts = m.split(only_watertight=False)
    if not len(parts):
        raise ValueError('No connected homolog')
    keep = max(parts, key=lambda x: x.area)
    return (keep, float(1 - keep.area / m.area), len(parts))

def cap(m):
    f = np.asarray(m.faces)
    v = np.asarray(m.vertices)
    ed = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    (_, inv, ct) = np.unique(np.sort(ed, axis=1), axis=0, return_inverse=True, return_counts=True)
    b = ed[ct[inv] == 1]
    pending = {tuple(e) for e in b.tolist()}
    new = []
    vertices = v.tolist()
    maxz = -1000000000.0
    while pending:
        (a, b) = pending.pop()
        loop = [a, b]
        start = a
        while loop[-1] != start:
            choices = [e for e in pending if e[0] == loop[-1]]
            if not choices:
                break
            e = min(choices)
            pending.remove(e)
            loop.append(e[1])
        if loop[-1] != start:
            raise ValueError('Nonclosed mesh boundary; no fabricated surface accepted')
        pp = v[loop[:-1]]
        center = pp.mean(0)
        maxz = max(maxz, float(pp[:, 2].max()))
        center[2] = min(center[2], v[:, 2].min() - 0.5)
        k = len(vertices)
        vertices.append(center.tolist())
        for (x, y) in zip(loop[:-1], loop[1:]):
            new.append([y, x, k])
    mesh = trimesh.Trimesh(np.array(vertices), np.vstack([f, np.array(new, int).reshape(-1, 3)]), process=False)
    trimesh.repair.fix_normals(mesh, multibody=True)
    return (mesh, dict(cap_faces=len(new), cap_boundary_max_z_mm=maxz, watertight=bool(mesh.is_watertight), winding=bool(mesh.is_winding_consistent)))

def fit(p, method):
    (donor, drop, ncomp) = main_component(p['homolog_triangles'])
    dv = donor.vertices.copy()
    dm = float(np.quantile(dv[donor.faces].reshape(-1, 3)[:, 2], 0.15))
    dv[:, 0] *= -1
    donor.faces = donor.faces[:, ::-1]
    dtri = dv[donor.faces]
    mask = (dtri[:, :, 2].max(1) <= dm + 0.5) & (dtri[:, :, 2].min(1) >= dm - 1.0)
    if mask.sum() < 30 or len(p['collar_triangles']) < 30:
        raise ValueError('insufficient observed collar')
    source = sample(dtri[mask])
    target = sample(p['collar_triangles'])
    center = np.array([np.median(source[:, 0]), np.median(source[:, 1]), dm])
    src = source - center
    v = dv - center
    m = float(p['margin_z'])
    t0 = np.array([np.median(target[:, 0]), np.median(target[:, 1]), m])
    A0 = np.eye(3)
    (pv, pf, _, _) = marching_cubes(p['preparation'], 0, spacing=(float(p['step']),) * 3)
    pv += p['origin']
    (pm, _drop, _ncomp) = main_component(pv[pf])
    top = float(pm.vertices[:, 2].max() + 1.28)
    if method == 'homolog_affine':
        A0[2, 2] = np.clip((top - m) / max(dv[:, 2].max() - dm, 0.1), 0.8, 1.25)
    A = A0.copy()
    t = t0.copy()
    tree = cKDTree(target)
    hist = []
    for iteration in range(35):
        curr = src @ A.T + t
        (d, idx) = tree.query(curr, workers=1)
        (dr, ir) = cKDTree(curr).query(target, workers=1)
        sa = np.vstack([src, src[ir]])
        ta = np.vstack([target[idx], target])
        dis = np.r_[d, dr]
        w = np.minimum(1.0, 1.0 / np.maximum(dis, 1e-12))
        hist.append(float(np.sqrt(np.mean(dis ** 2))))
        if method == 'homolog_rigid':
            w /= w.sum()
            cs = (sa * w[:, None]).sum(0)
            ct = (ta * w[:, None]).sum(0)
            H = (sa - cs).T @ ((ta - ct) * w[:, None])
            (U, S, Vt) = np.linalg.svd(H)
            B = Vt.T @ U.T
            if np.linalg.det(B) < 0:
                Vt[-1] *= -1
                B = Vt.T @ U.T
            rv = Rotation.from_matrix(B).as_rotvec()
            ang = np.linalg.norm(rv)
            if ang > np.deg2rad(25):
                B = Rotation.from_rotvec(rv * np.deg2rad(25) / ang).as_matrix()
            A = B
            t = ct - cs @ A.T
        else:
            X = np.c_[sa, np.ones(len(sa))]
            W = w[:, None]
            prior = np.vstack([A0.T, t0])
            ridge = np.diag([0.1 * len(sa)] * 3 + [0.0])
            coef = np.linalg.solve(X.T @ (W * X) + ridge, X.T @ (W * ta) + ridge @ prior)
            A = coef[:3].T
            t = coef[3]
            (U, S, Vt) = np.linalg.svd(A)
            A = U @ np.diag(np.clip(S, 0.8, 1.25)) @ Vt
        if not np.isfinite(A).all():
            raise ValueError('nonfinite fit')
    donor.vertices = v @ A.T + t
    (closed, capinfo) = cap(donor)
    return (closed, dict(A=A.tolist(), translation_mm=t.tolist(), center_mm=center.tolist(), collar_fit_rmse_mm=hist[-1], fit_history_mm=hist, donor_rejected_area_fraction=drop, donor_components=ncomp, **capinfo))

def extract(p, m, kind):
    low = np.minimum(m.bounds[0], p['origin'])
    high = np.maximum(m.bounds[1], p['origin'] + float(p['step']) * (np.array(p['preparation'].shape) - 1))
    center = (low + high) / 2
    side = float((high - low).max() + 0.6)
    n = int(np.ceil(side / 0.15))
    h = side / n
    origin = center - side / 2
    outer = mesh2sdf.compute((m.vertices - center) / (side / 2), m.faces, size=n, fix=False) * (side / 2)
    if outer.min() >= 0:
        raise ValueError('no signed interior')
    axes = [origin[k] + h * np.arange(n) for k in range(3)]
    X = np.meshgrid(*axes, indexing='ij')
    ix = np.array([(X[k] - p['origin'][k]) / float(p['step']) for k in range(3)])
    cavity = ndimage.map_coordinates(p['cavity'], ix, order=1, mode='constant', cval=100)
    cut = np.broadcast_to(float(p['margin_z']) - X[2], outer.shape)
    field = np.maximum(outer, cut) if kind == 'outer' else np.maximum.reduce([outer, -cavity, cut])
    (v, f, _, _) = marching_cubes(field, 0, spacing=(h,) * 3, allow_degenerate=False)
    v += origin
    pts = v[f].mean(1)
    q = ((pts - origin) / h).T
    if kind == 'outer':
        roles = np.where(ndimage.map_coordinates(outer, q, order=1) >= float(p['margin_z']) - pts[:, 2], 0, 2).astype(np.int8)
    else:
        roles = np.argmax(np.stack([ndimage.map_coordinates(outer, q, order=1), ndimage.map_coordinates(-cavity, q, order=1), float(p['margin_z']) - pts[:, 2]]), axis=0).astype(np.int8)
    mesh = trimesh.Trimesh(v, f, process=False)
    if mesh.volume < 0:
        f = f[:, ::-1]
    return (v, f, roles, dict(grid_step_mm=h, grid_n=n, signed_mesh_watertight=bool(m.is_watertight)))

class Plugin:

    def generate(self, task):
        allowed = {'public', 'method', 'kind'}
        if set(task) != allowed:
            raise ValueError('Unknown or missing public task keys')
        (m, info) = fit(task['public'], task['method'])
        (v, f, r, gi) = extract(task['public'], m, task['kind'])
        return dict(vertices=v, faces=f, face_roles=r, info={**info, **gi})

def run():
    records = json.loads(Path('/inputs/RECORDS.json').read_text())
    rows = []
    st = time.perf_counter()
    for rec in records:
        with np.load(Path('/inputs') / (rec['key'] + '.npz')) as a:
            p = dict(a)
        for method in ['homolog_rigid', 'homolog_affine']:
            t = time.perf_counter()
            try:
                (m, info) = fit(p, method)
            except Exception as e:
                for kind in ['outer', 'shell']:
                    rows.append(dict(rec, participant=method + '_' + kind, status='ABSTAIN', reason=str(e)))
                continue
            for kind in ['outer', 'shell']:
                name = method + '_' + kind
                row = dict(key=rec['key'], family=rec['family'], case_key=rec['case_key'], source_fdi=rec['source_fdi'], participant=name, status='FAILED', kind=kind)
                try:
                    (v, f, roles, grid) = extract(p, m, kind)
                    dest = Path('/output') / name / rec['key']
                    dest.mkdir(exist_ok=True, parents=True)
                    np.savez_compressed(dest / 'mesh.npz', vertices=v, faces=f, face_roles=roles)
                    trimesh.Trimesh(v @ p['source_R'].T + p['source_base'], f, process=False).export(dest / 'crown.stl')
                    (dest / 'FIT.json').write_text(json.dumps(dict(info, **grid), indent=2))
                    row.update(status='EXPORTED', fit=info, **grid)
                except Exception as e:
                    row['reason'] = str(e)
                row['seconds'] = time.perf_counter() - t
                rows.append(row)
            print(rec['key'], method, rows[-1]['status'], 'fit', round(info['collar_fit_rmse_mm'], 3), flush=True)
        Path('/output/RECORDS.json').write_text(json.dumps(rows, indent=2))
    Path('/output/COST.json').write_text(json.dumps(dict(seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024), indent=2))
if __name__ == '__main__':
    run()
