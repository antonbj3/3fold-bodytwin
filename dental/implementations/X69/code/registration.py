from dental_release.paths import expand as _release_expand
import pathlib, json, time, itertools, hashlib, resource
import numpy as np
import nibabel as nib
from scipy.spatial import cKDTree
from skimage.measure import marching_cubes
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = pathlib.Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/dental_sol_night/X69'))

def dump(path, obj):
    pathlib.Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False))

def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()

def stl(path):
    b = pathlib.Path(path).read_bytes()
    n = int.from_bytes(b[80:84], 'little')
    if len(b) != 84 + 50 * n:
        raise ValueError('STL format')
    dt = np.dtype([('normal', '<f4', (3,)), ('vertices', '<f4', (3, 3)), ('attribute', '<u2')])
    return np.frombuffer(b, dt, count=n, offset=84)['vertices'].astype(float)

def transform(p, T):
    return p @ T[:3, :3].T + T[:3, 3]

def sample(tri, n, rng):
    a = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    ids = rng.choice(len(tri), n, p=a / a.sum())
    uv = rng.random((n, 2))
    uv[uv.sum(1) > 1] = 1 - uv[uv.sum(1) > 1]
    return tri[ids, 0] + uv[:, 0, None] * (tri[ids, 1] - tri[ids, 0]) + uv[:, 1, None] * (tri[ids, 2] - tri[ids, 0])

def procrustes(p, q):
    pc = p.mean(0)
    qc = q.mean(0)
    (u, s, vt) = np.linalg.svd((p - pc).T @ (q - qc))
    R = vt.T @ u.T
    if np.linalg.det(R) < 0:
        vt[-1] *= -1
        R = vt.T @ u.T
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = qc - R @ pc
    return T

def icp(p, q, T, trim=0.75, iters=80):
    tree = cKDTree(q)
    T = T.copy()
    last = np.inf
    for i in range(iters):
        x = transform(p, T)
        (d, idx) = tree.query(x, workers=4)
        keep = np.argsort(d, kind='stable')[:int(len(d) * trim)]
        loss = float(np.mean(d[keep] ** 2))
        D = procrustes(x[keep], q[idx[keep]])
        T = D @ T
        if abs(last - loss) < 1e-10:
            break
        last = loss
    (d, _) = tree.query(transform(p, T), workers=4)
    loss = float(np.mean(np.sort(d)[:int(len(d) * trim)] ** 2))
    return (T, loss, i + 1)

def initializations(p, q):
    (_, Ep) = np.linalg.eigh(np.cov(p.T))
    (_, Eq) = np.linalg.eigh(np.cov(q.T))
    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product([-1, 1], repeat=3):
            P = np.eye(3)[:, perm] * np.array(signs)[None, :]
            R = Eq @ P @ Ep.T
            if np.linalg.det(R) < 0:
                continue
            T = np.eye(4)
            T[:3, :3] = R
            T[:3, 3] = q.mean(0) - R @ p.mean(0)
            out.append(T)
    return out

def stats(d):
    return {'n': len(d), 'mean_mm': float(np.mean(d)), 'rms_mm': float(np.sqrt(np.mean(d * d))), 'p95_mm': float(np.percentile(d, 95)), 'max_mm': float(np.max(d)), 'resolution': 'PER_SURFACE_REGION'}

def extract():
    cache = DATA / 'r1_geometry.npz'
    if cache.exists():
        return np.load(cache)
    im = nib.load(DATA / '001_cbct.nii.gz')
    a = np.asarray(im.dataobj)
    sh = np.array(a.shape)
    low = (sh * np.array([0.2, 0.2, 0])).astype(int)
    high = (sh * np.array([0.8, 0.8, 1])).astype(int)
    coarse = a[low[0]:high[0]:2, low[1]:high[1]:2, ::2]
    idx = np.argwhere(coarse > 2200) * 2 + low
    margin = np.ceil(8 / np.array(im.header.get_zooms())).astype(int)
    lo = np.maximum(low, idx.min(0) - margin)
    hi = np.minimum(high, idx.max(0) + margin + 1)
    crop = a[tuple((slice(l, h) for (l, h) in zip(lo, hi)))].astype(np.float32)
    (v, f, n, _) = marching_cubes(crop, level=1800, step_size=2, allow_degenerate=False)
    target = nib.affines.apply_affine(im.affine, v + lo)
    rng = np.random.default_rng(69)
    points = []
    jaws = []
    regions = []
    supports = []
    for jaw in ['lower', 'upper']:
        tri = stl(DATA / f'001_{jaw.title()}JawScan.stl')
        cent = tri.mean(1)
        mid = np.median(cent[:, 1])
        crown = tri[cent[:, 1] >= mid] if jaw == 'lower' else tri[cent[:, 1] <= mid]
        p = sample(crown, 9000, rng)
        ter = np.quantile(p[:, 0], [1 / 3, 2 / 3])
        reg = np.digitize(p[:, 0], ter)
        points.append(p)
        jaws.extend([jaw] * len(p))
        regions.extend(reg.tolist())
        supports.append({'jaw': jaw, 'native_y_boundary': float(mid), 'x_tercile_boundaries': ter.tolist(), 'faces_total': len(tri), 'faces_crown_closure': len(crown), 'faces_rejected': len(tri) - len(crown), 'rejected_fraction': 1 - len(crown) / len(tri)})
    np.savez_compressed(cache, points=np.concatenate(points), jaw=np.array(jaws), region=np.array(regions), target=target, faces=f, affine=im.affine, crop_lo=lo)
    dump(ROOT / 'raw/EXTRACTION.json', {'CBCT_shape': im.shape, 'voxel_mm': [float(x) for x in im.header.get_zooms()], 'unit': im.header.get_xyzt_units()[0], 'NIfTI_affine': im.affine.tolist(), 'crop_lo_ijk': lo.tolist(), 'crop_hi_ijk': hi.tolist(), 'surface_vertices': len(v), 'surface_triangles': len(f), 'isovalue': 1800, 'isovalue_is_HU': False, 'supports': supports, 'semantic_labels': 'UNKNOWN: isosurface is not a tooth segmentation', 'IOS_metric_unit': 'mm assumption, uncalibrated STL', 'data_array_sha256': sha(cache)})
    return np.load(cache)

def nearest_surface(points, vertices, faces):
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    pts = vtk.vtkPoints()
    pts.SetData(numpy_to_vtk(vertices, deep=True))
    cells = vtk.vtkCellArray()
    conn = np.column_stack((np.full(len(faces), 3, dtype=np.int64), faces)).ravel().astype(np.int64)
    cells.SetCells(len(faces), numpy_to_vtkIdTypeArray(conn, deep=True))
    mesh = vtk.vtkPolyData()
    mesh.SetPoints(pts)
    mesh.SetPolys(cells)
    loc = vtk.vtkStaticCellLocator()
    loc.SetDataSet(mesh)
    loc.BuildLocator()
    ds = []
    qs = []
    ids = []
    for p in points:
        closest = [0.0, 0.0, 0.0]
        cell = vtk.reference(0)
        sub = vtk.reference(0)
        dist = vtk.reference(0.0)
        loc.FindClosestPoint(p, closest, cell, sub, dist)
        ds.append(float(dist) ** 0.5)
        qs.append(closest)
        ids.append(int(cell))
    return (np.array(ds), np.array(qs), np.array(ids))

def run_r1():
    start = time.monotonic()
    g = extract()
    p = g['points']
    target = g['target']
    reg = g['region']
    fit = reg != 1
    rng = np.random.default_rng(69)
    q = target[rng.choice(len(target), min(18000, len(target)), replace=False)]
    rows = []
    best = None
    for (i, T) in enumerate(initializations(p[fit], q)):
        (T, loss, n) = icp(p[fit], q, T)
        rows.append({'start': i, 'fit_loss_mm2': loss, 'iterations': n, 'transform': T.tolist()})
        if best is None or loss < best[0]:
            best = (loss, T)
        print('R1 start', i, 'fit_loss', loss, flush=True)
    T = best[1]
    x = transform(p, T)
    (d, witness, cell) = nearest_surface(x, target, g['faces'])
    regions = []
    for jaw in ['lower', 'upper']:
        for k in range(3):
            sel = (g['jaw'] == jaw) & (reg == k)
            regions.append({'jaw': jaw, 'native_x_tercile': k, 'heldout_from_fit': k == 1, **stats(d[sel]), 'gate_p95_mm': 0.5, 'gate_pass': float(np.percentile(d[sel], 95)) <= 0.5})
    np.savez_compressed(DATA / 'r1_witnesses.npz', source_points=p, transformed_points=x, nearest_surface=witness, surface_cell=cell, residual_mm=d, region=reg, jaw=g['jaw'], source_to_target=T)
    result = {'construction': 'R1', 'claim_type': ['capability', 'information_link'], 'actual_pairs_acquired': 1, 'unit': 'mm conditional on STL source-unit assumption', 'source_to_target': T.tolist(), 'fit_summary': stats(d[fit]), 'heldout_summary': stats(d[~fit]), 'regions': regions, 'heldout_gate_pass': all((r['gate_pass'] for r in regions if r['heldout_from_fit'])), 'independent_landmark_TRE': 'UNKNOWN', 'reference_pose_error': 'UNKNOWN_NO_PUBLISHER_TRANSFORM', 'semantic_tooth_ownership': 'UNKNOWN', 'seconds': time.monotonic() - start, 'maxrss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'witnesses_path': str(DATA / 'r1_witnesses.npz'), 'witnesses_sha256': sha(DATA / 'r1_witnesses.npz'), 'starts': rows, 'uniform_error_bound_mm': None, 'physical_accuracy': 'UNKNOWN', 'external_referent': {'kind': 'published_dataset', 'locator': 'https://doi.org/10.6084/m9.figshare.26965903.v3', 'compared_quantity': 'source-declared paired acquired crown geometry; surface consistency with released CBCT intensity isosurface', 'refutes_us': True}}
    dump(ROOT / 'raw/R1_RESULTS.json', result)
    np.savez_compressed(DATA / 'figure_r1.npz', target=target[::max(1, len(target) // 9000)], source=x, residual=d, jaw=g['jaw'])
    print(json.dumps({k: result[k] for k in ['heldout_summary', 'heldout_gate_pass', 'seconds', 'maxrss_MiB']}, indent=2))
if __name__ == '__main__':
    run_r1()
