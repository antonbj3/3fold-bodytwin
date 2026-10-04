from common import *
from geometry import MET, preparation, field_values, sample_faces, targets, topo
import trimesh
from scipy import ndimage
from skimage.measure import marching_cubes
try:
    from trimesh.ray.ray_pyembree import RayMeshIntersector
except ImportError:
    from trimesh.ray.ray_triangle import RayMeshIntersector
Q = read(P / 'PREREG_R3.json')['parameters']
FINE = {}

def fine_prep(rec):
    if rec['key'] not in FINE:
        (p0, pm) = preparation(rec)
        p = dict(p0)
        p['origin'] = p0['origin'].copy()
        h = Q['remesh_step_mm']
        p['origin'][2] = float(p['margin_z']) - np.floor((float(p['margin_z']) - p0['origin'][2]) / h) * h
        offset = (p['origin'] - p0['origin']) / float(p0['step'])
        step = float(p['step'])
        a = p['preparation']
        n = np.ceil((np.array(a.shape) - 1) * step / h).astype(int) + 1
        pp = np.empty(n, dtype=np.float32)
        yz = np.stack(np.meshgrid(np.arange(n[1]) * h / step + offset[1], np.arange(n[2]) * h / step + offset[2], indexing='ij'))
        for i in range(n[0]):
            pp[i] = ndimage.map_coordinates(a, np.vstack([np.full(yz[0].size, i * h / step + offset[0]), yz.reshape(2, -1)]), order=1, mode='nearest').reshape(n[1:])
        grad = tuple(np.gradient(pp, h))
        norm = np.sqrt(sum((g * g for g in grad)))
        zs = p['origin'][2] + h * np.arange(n[2])
        g = np.where(zs < float(p['margin_z']) + 0.5, 0.03, 0.08).astype(np.float32)
        cavity = pp - g[None, None, :] * np.maximum(norm, Q['gradient_floor'])
        collar = pp - (g[None, None, :] + Q['collar_width_mm']) * np.maximum(norm, Q['gradient_floor'])
        cut = (float(p['margin_z']) - zs).astype(np.float32)
        collar = np.maximum(collar, (zs - float(p['margin_z']) - Q['collar_height_mm'])[None, None, :])
        FINE[rec['key']] = (p, pm, cavity, collar, cut)
    return FINE[rec['key']]

def sampled_shape(v, f, roles, rec):
    a = load_npz(V4 / 'payload/whole_private' / rec['key'] / 'reference.npz')
    (p, _) = preparation(rec)
    tri = a['source_triangles']
    tri = tri[tri.mean(1)[:, 2] >= float(p['margin_z'])]
    source = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=False)
    ext = trimesh.Trimesh(v, f[roles == 0], process=False)
    ps = sample_faces(source.vertices, source.faces, Q['source_shape_probes'])[0]
    pe = sample_faces(v, f[roles == 0], Q['source_shape_probes'])[0]
    return max(float(np.quantile(MET.closest(ext, ps)[1], 0.95)), float(np.quantile(MET.closest(source, pe)[1], 0.95)))

def construct(rec):
    (p, pm, cavity, collar, cut) = fine_prep(rec)
    h = Q['remesh_step_mm']
    origin = p['origin']
    z = load_npz(rec['mesh_path'])
    old = trimesh.Trimesh(z['vertices'], z['faces'], process=False)
    vox = old.voxelized(h).fill()
    ix = np.rint((vox.points - origin) / h).astype(int)
    inside = np.all((ix >= 0) & (ix < np.array(cavity.shape)), axis=1)
    occ = np.zeros(cavity.shape, bool)
    ix = ix[inside]
    occ[tuple(ix.T)] = True
    occ |= cavity <= 0
    out = (ndimage.distance_transform_edt(~occ, sampling=h) - ndimage.distance_transform_edt(occ, sampling=h)).astype(np.float32)
    outer = np.minimum(out, collar)
    field = np.maximum(np.maximum(outer, -cavity), cut[None, None, :])
    (v, f, _, _) = marching_cubes(field, 0, spacing=(h,) * 3, allow_degenerate=False)
    v += origin
    m = trimesh.Trimesh(v, f, process=False)
    if m.volume < 0:
        f = f[:, ::-1]
    cent = v[f].mean(1)
    coords = ((cent - origin) / h).T
    vals = np.stack([ndimage.map_coordinates(outer, coords, order=1), ndimage.map_coordinates(-cavity, coords, order=1), float(p['margin_z']) - cent[:, 2]])
    roles = np.argmax(vals, axis=0).astype(np.int8)
    it = f[roles == 1]
    (sites, fi, w) = sample_faces(v, it, Q['probes'])
    (g, cp, n) = field_values(pm, sites)
    tg = targets(cp, float(p['margin_z']))
    ext = trimesh.Trimesh(v, f[roles == 0], process=False)
    ray = RayMeshIntersector(ext)
    (hit, ri, _) = ray.intersects_location(sites + n * 1e-05, n, multiple_hits=False)
    t = np.full(len(sites), np.nan)
    distances = np.sum((hit - sites[ri]) * n[ri], axis=1)
    good = (distances > 0) & (distances <= Q['ray_limit_mm'])
    t[ri[good]] = distances[good]
    tp = topo(v, f)
    margverts = np.unique(f[roles == 2])
    marginerr = float(np.max(abs(v[margverts, 2] - float(p['margin_z'])))) if len(margverts) else float('inf')
    shape = sampled_shape(v, f, roles, rec)
    valid = np.isfinite(t)
    gate = dict(closed_positive_oriented=tp['watertight'] and tp['winding_consistent'] and tp['positive_volume'], intaglio_gap=bool(np.quantile(abs(g - tg), 0.95) <= 0.025), nonnegative_gap=bool(g.min() >= 0), marginal_plane=bool(marginerr <= 1e-06), ray_coverage=bool(valid.mean() >= 0.95), postadapt_shape=bool(shape - rec['reconstruction_p95_mm'] <= 0.35))
    folder = DATA / 'R3' / rec['uid']
    folder.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(folder / 'mesh.npz', vertices=v, faces=f, face_roles=roles)
    np.savez_compressed(folder / 'fields.npz', sites_mm=sites, normals=n, area_weights_mm2=w, gap_mm=g, target_mm=tg, thickness_mm=t, source_face_ix=fi)
    out = dict(uid=rec['uid'], key=rec['key'], participant=rec['participant'], track=rec['track'], family=rec['family'], status='ADAPTED', gates=gate, geometry_pass=all(gate.values()), topology=tp, margin_error_mm=marginerr, original_shape_p95_mm=rec['reconstruction_p95_mm'], adapted_shape_p95_mm=shape, shape_worsening_mm=shape - rec['reconstruction_p95_mm'], gap_error_p95_um=float(np.quantile(abs(g - tg), 0.95) * 1000), gap_p05_um=float(np.quantile(g, 0.05) * 1000), gap_mean_um=float(g.mean() * 1000), gap_p95_um=float(np.quantile(g, 0.95) * 1000), gap_min_um=float(g.min() * 1000), ray_coverage=float(valid.mean()), thickness_median_mm=float(np.nanmedian(t)), thickness_p05_mm=float(np.nanquantile(t, 0.05)), thickness_p95_mm=float(np.nanquantile(t, 0.95)), resolution='PER_TOOTH', field_resolution='PER_POINT', continuous_geometry_enclosure='MISSING', physical_seated_film='UNKNOWN_UNMEASURED', mesh_path=str(folder / 'mesh.npz'), field_path=str(folder / 'fields.npz'), mesh_sha256=sha(folder / 'mesh.npz'), field_sha256=sha(folder / 'fields.npz'), source_mesh_sha256=sha(rec['mesh_path']), prep_sha256=sha(rec['prep_path']))
    return out

def run(limit=None):
    rows = []
    start = time.perf_counter()
    items = cohort() if limit is None else cohort()[:limit]
    for (i, r) in enumerate(items):
        if r['status'] != 'SCORED':
            rows.append(dict(uid=r['uid'], key=r['key'], participant=r['participant'], status='SOURCE_UNSCORED', reason=r.get('reason', r['status'])))
            continue
        st = time.perf_counter()
        try:
            out = construct(r)
        except (ValueError, RuntimeError, IndexError) as e:
            out = dict(uid=r['uid'], key=r['key'], participant=r['participant'], status='ADAPTATION_FAILED', reason=str(e))
        out['seconds'] = time.perf_counter() - st
        rows.append(out)
        dump(P / 'raw/R3_CHECKPOINT.json', rows)
        print(i + 1, r['participant'], out['status'], out.get('geometry_pass'), round(out['seconds'], 2), flush=True)
    dump(P / ('raw/R3_PILOT.json' if limit else 'raw/R3_GEOMETRY.json'), dict(rows=rows, seconds=time.perf_counter() - start))
    return rows
if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else None)
