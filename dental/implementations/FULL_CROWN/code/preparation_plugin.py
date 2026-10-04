import os, sys, json, time, resource
from pathlib import Path
import numpy as np
from scipy import ndimage, sparse
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree
from skimage.measure import marching_cubes
import trimesh
from homolog_plugin import fit, cap, extract, sample, main_component

def construct(p):
    (mesh, info) = fit(p, 'homolog_rigid')
    margin = float(p['margin_z'])
    h = float(p['step'])
    (v, f, _, _) = marching_cubes(p['preparation'], 1.28, spacing=(h,) * 3)
    v += p['origin']
    tri = v[f]
    norm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    norm /= np.maximum(np.linalg.norm(norm, axis=1, keepdims=True), 1e-12)
    center = tri.mean(1)
    iz = ((center - p['origin']) / h).T
    grad = np.gradient(p['preparation'], h)
    nz = ndimage.map_coordinates(grad[2], iz, order=1)
    gn = np.sqrt(sum((ndimage.map_coordinates(g, iz, order=1) ** 2 for g in grad)))
    mask = (nz / np.maximum(gn, 1e-12) > 0.55) & (center[:, 2] > margin + 2.0)
    if mask.sum() < 10:
        raise ValueError('no supported inverse occlusal patch')
    tc = sample(p['collar_triangles'], 1200)
    to = sample(tri[mask], 2400)
    trees = [cKDTree(tc), cKDTree(to)]
    pts = mesh.vertices.copy()
    hist = []

    def assign(q):
        mesh.vertices = q
        norms = mesh.vertex_normals
        mc = (q[:, 2] >= margin - 1.0) & (q[:, 2] <= margin + 0.5)
        mo = (norms[:, 2] > 0.55) & (q[:, 2] > margin + 2.0)
        mc[mo] = False
        dest = np.zeros_like(q)
        w = np.zeros(len(q))
        ds = []
        for (m, t) in [(mc, trees[0]), (mo, trees[1])]:
            if not m.any():
                continue
            (dist, ix) = t.query(q[m], workers=1)
            dest[m] = t.data[ix]
            w[m] = np.minimum(1.0, 1.0 / np.maximum(dist, 1e-09))
            ds.extend(dist.tolist())
        return (dest, w, float(np.sqrt(np.mean(np.square(ds)))) if ds else 100.0)
    for it in range(25):
        (dest, w, err) = assign(pts)
        hist.append(err)
        good = w > 0
        if good.sum() < 30:
            raise ValueError('insufficient patch correspondences')
        a = pts[good]
        b = dest[good]
        ww = w[good]
        ww /= ww.sum()
        c0 = np.sum(a * ww[:, None], 0)
        c1 = np.sum(b * ww[:, None], 0)
        (U, S, Vt) = np.linalg.svd((a - c0).T @ ((b - c1) * ww[:, None]))
        R = Vt.T @ U.T
        if np.linalg.det(R) < 0:
            Vt[-1] *= -1
            R = Vt.T @ U.T
        pts = (pts - c0) @ R.T + c1
    initial = pts.copy()
    edges = mesh.edges_unique
    ii = np.r_[edges[:, 0], edges[:, 1]]
    jj = np.r_[edges[:, 1], edges[:, 0]]
    adj = sparse.coo_matrix((np.ones(len(ii)), (ii, jj)), shape=(len(pts), len(pts))).tocsr()
    deg = np.asarray(adj.sum(1)).ravel()
    L = sparse.eye(len(pts)) - sparse.diags(1 / np.maximum(deg, 1)) @ adj
    pen = 5 * (L.T @ L) + 0.02 * sparse.eye(len(pts))
    hist2 = []
    for it in range(6):
        (dest, w, err) = assign(pts)
        hist2.append(err)
        M = (pen + sparse.diags(w)).tocsc()
        d = splu(M).solve(w[:, None] * (dest - initial))
        norm = np.linalg.norm(d, axis=1)
        d *= np.minimum(1, 2 / np.maximum(norm, 1e-12))[:, None]
        pts = initial + d
    mesh.vertices = pts
    info.update(inverse_occlusal_faces=int(mask.sum()), joint_rigid_fit_history_mm=hist, joint_deform_fit_history_mm=hist2, max_displacement_mm=float(np.linalg.norm(pts - initial, axis=1).max()), shape_prior='same-case homolog; joint collar plus conditional inverse-prep observations, no hidden original')
    return (mesh, info)

def run():
    records = json.loads(Path('/inputs/RECORDS.json').read_text())
    out = Path('/output')
    rows = []
    start = time.perf_counter()
    for rec in records:
        with np.load(Path('/inputs') / (rec['key'] + '.npz')) as a:
            p = dict(a)
        st = time.perf_counter()
        try:
            (m, info) = construct(p)
        except Exception as e:
            for kind in ['outer', 'shell']:
                rows.append(dict(key=rec['key'], family=rec['family'], case_key=rec['case_key'], participant='homolog_preparation_' + kind, status='ABSTAIN', reason=str(e), kind=kind))
            (out / 'RECORDS.json').write_text(json.dumps(rows, indent=2))
            continue
        for kind in ['outer', 'shell']:
            name = 'homolog_preparation_' + kind
            row = dict(key=rec['key'], family=rec['family'], case_key=rec['case_key'], source_fdi=rec['source_fdi'], participant=name, status='FAILED', kind=kind)
            try:
                (v, f, roles, grid) = extract(p, m, kind)
                dest = out / name / rec['key']
                dest.mkdir(exist_ok=True, parents=True)
                np.savez_compressed(dest / 'mesh.npz', vertices=v, faces=f, face_roles=roles)
                trimesh.Trimesh(v @ p['source_R'].T + p['source_base'], f, process=False).export(dest / 'crown.stl')
                (dest / 'FIT.json').write_text(json.dumps(dict(info, **grid), indent=2))
                row.update(status='EXPORTED', fit=info, **grid)
            except Exception as e:
                row['reason'] = str(e)
            row['seconds'] = time.perf_counter() - st
            rows.append(row)
        (out / 'RECORDS.json').write_text(json.dumps(rows, indent=2))
        print(rec['key'], rows[-1]['status'], 'fit', round(info['joint_deform_fit_history_mm'][-1], 4), flush=True)
    (out / 'COST.json').write_text(json.dumps(dict(seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024), indent=2))
if __name__ == '__main__':
    run()
