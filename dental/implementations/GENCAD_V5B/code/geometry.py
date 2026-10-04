from common import *
import trimesh
from skimage.measure import marching_cubes
from scipy import ndimage
MET = load_module('x55_metrology_readonly', D / 'LANE_X55_METROLOGY/code/metrology.py')
PR = read(P / 'PREREG_R1.json')
Q = PR['parameters']
CACHE = {}

def prep_mesh(p):
    h = float(p['step'])
    (v, f, _, _) = marching_cubes(p['preparation'], 0, spacing=(h,) * 3, allow_degenerate=False)
    v += p['origin']
    m = trimesh.Trimesh(v, f, process=False)
    if m.volume < 0:
        m.invert()
    return m

def preparation(rec):
    key = rec['key']
    if key not in CACHE:
        p = load_npz(rec['prep_path'])
        m = prep_mesh(p)
        T = np.eye(4)
        T[:3, :3] = p['source_R']
        T[:3, 3] = p['source_base']
        MET.validate_transform(T)
        CACHE[key] = (p, m)
    return CACHE[key]

def targets(points, margin):
    return np.where(points[:, 2] < margin + Q['marginal_band_mm'], Q['marginal_spacer_mm'], Q['internal_spacer_mm'])

def field_values(pm, points):
    (c, d, f) = MET.closest(pm, points)
    n = pm.face_normals[f]
    signed = np.sum((points - c) * n, axis=1)
    return (signed, c, n)

def sample_faces(v, f, n=512):
    tri = v[f]
    areas = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    if not len(areas) or areas.sum() <= 0:
        raise ValueError('Empty positive-area surface')
    ix = np.searchsorted(np.cumsum(areas), (np.arange(n) + 0.5) / n * areas.sum())
    return (tri[ix].mean(1), ix, np.full(n, areas.sum() / n))

def topo(v, f):
    m = trimesh.Trimesh(v, f, process=False)
    return dict(watertight=bool(m.is_watertight), winding_consistent=bool(m.is_winding_consistent), positive_volume=bool(m.volume > 0), min_triangle_area_mm2=float(m.area_faces.min()), volume_mm3=float(m.volume))

def adapt(rec):
    z = load_npz(rec['mesh_path'])
    v0 = z['vertices']
    v = v0.copy().astype(float)
    f = z['faces']
    roles = z['face_roles']
    (p, pm) = preparation(rec)
    margin = float(p['margin_z'])
    inner = np.unique(f[roles == 1])
    marginal = np.unique(f[roles == 2]) if np.any(roles == 2) else np.array([], int)
    if not len(marginal):
        marginal = np.flatnonzero(v[:, 2] <= margin + float(p['step']) / 2)
    moved = np.union1d(inner, marginal)
    outer_only = np.setdiff1d(np.arange(len(v)), moved)
    if not len(inner):
        raise ValueError('NO_INTAGLIO')
    for k in range(Q['projection_iterations']):
        (gap, c, n) = field_values(pm, v[inner])
        v[inner] = c + targets(c, margin)[:, None] * n
        if len(marginal):
            v[marginal, 2] = margin
    (ip, ix, w) = sample_faces(v, f[roles == 1], Q['probes'])
    (g, cp, n) = field_values(pm, ip)
    target = targets(cp, margin)
    extfaces = f[roles == 0]
    em = trimesh.Trimesh(v, extfaces, process=False)
    (ep, t, ef) = MET.closest(em, ip)
    sign = np.sum((ep - ip) * n, axis=1)
    signed_t = np.where(sign >= 0, t, -t)
    (before, _, _) = field_values(pm, sample_faces(v0, f[roles == 1], Q['probes'])[0])
    tpg = topo(v, f)
    identity = float(np.max(abs(v[outer_only] - v0[outer_only]))) if len(outer_only) else 0.0
    gate = dict(closed_positive_oriented=tpg['watertight'] and tpg['winding_consistent'] and tpg['positive_volume'] and (tpg['min_triangle_area_mm2'] > 1e-14), intaglio_gap=bool(np.quantile(abs(g - target), 0.95) <= 0.025), marginal_plane=bool(len(marginal) and np.max(abs(v[marginal, 2] - margin)) <= 1e-09), outer_only_identity=identity == 0, nonnegative_gap=bool(g.min() >= 0), positive_thickness=bool(signed_t.min() > 0), displacement=bool(np.quantile(np.linalg.norm(v - v0, axis=1), 0.95) <= 0.5))
    out = dict(uid=rec['uid'], key=rec['key'], participant=rec['participant'], track=rec['track'], family=rec['family'], status='ADAPTED', geometry_pass=all(gate.values()), gates=gate, topology=tpg, margin_kind='PUBLIC_VIRTUAL_FINISH_PLANE', outer_only_identity_error_mm=identity, displacement_p95_mm=float(np.quantile(np.linalg.norm(v - v0, axis=1), 0.95)), gap_p05_um=float(np.quantile(g, 0.05) * 1000), gap_mean_um=float(g.mean() * 1000), gap_p95_um=float(np.quantile(g, 0.95) * 1000), gap_min_um=float(g.min() * 1000), gap_error_p95_um=float(np.quantile(abs(g - target), 0.95) * 1000), original_gap_p95_um=float(np.quantile(before, 0.95) * 1000), thickness_min_mm=float(signed_t.min()), thickness_p05_mm=float(np.quantile(signed_t, 0.05)), thickness_median_mm=float(np.median(signed_t)), thickness_mean_mm=float(np.mean(signed_t)), thickness_p95_mm=float(np.quantile(signed_t, 0.95)), original_shape_p95_mm=rec['reconstruction_p95_mm'], resolution='PER_TOOTH', field_resolution='PER_POINT', sample_count=len(g), continuous_geometry_enclosure='MISSING', seated_film_status='UNKNOWN_UNMEASURED', source_mesh_sha256=sha(rec['mesh_path']), prep_sha256=sha(rec['prep_path']))
    folder = DATA / 'R1' / rec['uid']
    folder.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(folder / 'mesh.npz', vertices=v, faces=f, face_roles=roles)
    np.savez_compressed(folder / 'fields.npz', sites_mm=ip, normals=n, area_weights_mm2=w, gap_mm=g, target_mm=target, thickness_mm=signed_t, source_face_ix=ix)
    trimesh.Trimesh(v @ p['source_R'].T + p['source_base'], f, process=False).export(folder / 'crown.stl')
    out.update(mesh_path=str(folder / 'mesh.npz'), field_path=str(folder / 'fields.npz'), mesh_sha256=sha(folder / 'mesh.npz'), field_sha256=sha(folder / 'fields.npz'), stl_sha256=sha(folder / 'crown.stl'))
    return out

def run(limit=None):
    rows = []
    t0 = time.perf_counter()
    for (i, r) in enumerate(cohort() if limit is None else cohort()[:limit]):
        if r['status'] != 'SCORED':
            rows.append(dict(uid=r['uid'], participant=r['participant'], key=r['key'], status='SOURCE_UNSCORED', reason=r.get('reason', r['status'])))
            continue
        st = time.perf_counter()
        try:
            out = adapt(r)
        except (ValueError, RuntimeError, IndexError) as e:
            out = dict(uid=r['uid'], participant=r['participant'], key=r['key'], status='ADAPTATION_FAILED', reason=str(e))
        out['seconds'] = time.perf_counter() - st
        rows.append(out)
        dump(P / 'raw/R1_CHECKPOINT.json', rows)
        print(i + 1, r['participant'], out['status'], out.get('geometry_pass'), round(out['seconds'], 2), flush=True)
    dump(P / ('raw/R1_PILOT.json' if limit else 'raw/R1_GEOMETRY.json'), dict(rows=rows, seconds=time.perf_counter() - t0))
    return rows
if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else None)
