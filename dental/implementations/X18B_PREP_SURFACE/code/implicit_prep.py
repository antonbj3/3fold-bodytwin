"""Euclidean erosion of an explicitly closed virtual insertion envelope."""
from common import *
from scipy import ndimage
from scipy.spatial import cKDTree
from skimage.measure import marching_cubes
import trimesh, resource

def mesh_field(field, origin, step, level):
    (v, f, _, _) = marching_cubes(field.astype(np.float32), level=level, spacing=(step,) * 3, allow_degenerate=False)
    return (v + origin, f)

def predict():
    pr = json.loads((H / 'PREREG_R4.json').read_text())
    s = pr['cases'][0]
    start = time.perf_counter()
    (data, man) = pair(s['case'])
    a = data[s['jaw']]
    tri = a['tri'][a['owner'] == s['fdi']]
    step = pr['grid_step_mm']
    lo = tri.min((0, 1)) - 3 * step
    hi = tri.max((0, 1)) + 3 * step
    axes = [np.arange(lo[k], hi[k] + step * 0.1, step) for k in range(3)]
    (xx, yy) = np.meshgrid(axes[0], axes[1], indexing='ij')
    xy = np.c_[xx.ravel(), yy.ravel()]
    (height, _) = query_height(tri, xy, False)
    height = height.reshape(xx.shape)
    base = float(tri[:, :, 2].min())
    zz = axes[2][None, None, :]
    occupied = np.isfinite(height)[:, :, None] & (zz >= base) & (zz <= height[:, :, None])
    if occupied.sum() == 0:
        raise RuntimeError('Empty virtual source envelope')
    inside = ndimage.distance_transform_edt(occupied, sampling=(step,) * 3)
    outside = ndimage.distance_transform_edt(~occupied, sampling=(step,) * 3)
    phi = outside - inside
    eroded = inside >= pr['depth_mm']
    if eroded.sum() == 0:
        raise RuntimeError('Empty implicit preparation after erosion')
    origin = np.array([q[0] for q in axes])
    (v, f) = mesh_field(phi, origin, step, -pr['depth_mm'])
    (v0, f0) = mesh_field(phi, origin, step, 0)
    path = D / '127_36_implicit_preparation.stl'
    bad = D / '127_36_injected_no_erosion.stl'
    write_stl(path, v, f)
    write_stl(bad, v0, f0)
    raw = D / '127_36_implicit_fields.npz'
    np.savez_compressed(raw, phi=phi.astype(np.float32), occupied=occupied, eroded=eroded, origin=origin, grid_step_mm=step, source_height=height, basal_z_mm=base, vertices=v, faces=f, unprepared_vertices=v0, unprepared_faces=f0)
    pred = dict(case=s['case'], fdi=s['fdi'], resolution='PER_TOOTH', source=man, grid_shape=list(occupied.shape), occupied_voxels=int(occupied.sum()), eroded_voxels=int(eroded.sum()), virtual_volume_mm3=float(occupied.sum() * step ** 3), virtual_preparation_volume_mm3=float(eroded.sum() * step ** 3), virtual_removed_volume_mm3=float((occupied.sum() - eroded.sum()) * step ** 3), erosion_depth_mm=pr['depth_mm'], basal_z_mm=base, preparation_mesh=artifact(path), fault_mesh=artifact(bad), field_arrays=artifact(raw), runtime_wall_s=time.perf_counter() - start, source_margin='UNKNOWN; virtual basal plane, no cervical annotation', physical_preparation_validation='NOT_RUN')
    dump(H / 'raw/PREDICTIONS_R4.json', pred)
    files = [pred['preparation_mesh'], pred['fault_mesh'], pred['field_arrays'], artifact(H / 'raw/PREDICTIONS_R4.json')]
    dump(H / 'FROZEN_PREDICTIONS_R4.json', dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(H / 'PREREG_R4.json'), files=files, code_sha256={str(p): sha(p) for p in [H / 'code/implicit_prep.py', H / 'code/common.py']}, physical_measurement_status='NOT_RUN', external_carrier='R1 exact exterior retained; implicit representation applies only to virtual preparation'))
    state('R4_PREDICTIONS_FROZEN', 'Implicit preparation and zero-depth fault frozen', 'Independent nearest-background distance, ray monotonicity and topology checks')
    print(json.dumps(pred, indent=2))

def evaluate():
    verify_freeze('R4')
    pr = json.loads((H / 'PREREG_R4.json').read_text())
    p = json.loads((H / 'raw/PREDICTIONS_R4.json').read_text())
    a = np.load(p['field_arrays']['path'])
    occupied = a['occupied']
    eroded = a['eroded']
    step = pr['grid_step_mm']
    v = a['vertices']
    f = a['faces']
    origin = a['origin']
    start = time.perf_counter()
    exterior = ndimage.binary_dilation(occupied, structure=np.ones((3, 3, 3))) & ~occupied
    background = np.argwhere(exterior) * step + origin
    tree = cKDTree(background)
    distance = tree.query(v, workers=1)[0]
    bad_distance = tree.query(a['unprepared_vertices'], workers=1)[0]
    err = float(np.max(np.abs(distance - pr['depth_mm'])))
    starts = eroded & ~np.concatenate([np.zeros(eroded.shape[:2] + (1,), bool), eroded[:, :, :-1]], axis=2)
    raymax = int(starts.sum(2).max())
    edges = np.sort(np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1)
    counts = np.unique(edges, axis=0, return_counts=True)[1]
    nonmanifold = int(np.sum(counts != 2))
    mesh = trimesh.Trimesh(vertices=v, faces=f, process=False)
    carrier = np.load(D / '127_36_carrier.npz')
    ext = carrier['external_vertices'][carrier['external_faces']]
    (data, _) = pair(127)
    src = data['lower']['tri'][data['lower']['owner'] == 36]
    exterror = float(np.max(np.abs(src - ext)))
    m = trimesh.Trimesh(vertices=carrier['external_vertices'], faces=carrier['external_faces'], process=False)
    sample = np.linspace(0, len(v) - 1, min(2048, len(v)), dtype=int)
    dd = []
    for chunk in np.array_split(sample, max(1, int(np.ceil(len(sample) / 256)))):
        dd.extend(trimesh.proximity.closest_point(m, v[chunk])[1].tolist())
    gates = dict(manifold=nonmanifold == 0 and bool(mesh.is_watertight), depth=err <= pr['metrics']['erosion_distance_tolerance_mm'], external=exterror <= pr['metrics']['external_unchanged_max_mm'], insertion_rays=raymax <= pr['metrics']['ray_occupancy_interval_max'], nonempty=bool(eroded.any()), injected_zero_depth=float(np.max(np.abs(bad_distance - pr['depth_mm']))) > pr['metrics']['erosion_distance_tolerance_mm'])
    out = dict(round='R4', claim_type='capability', external_referent=pr['external_referent'], additional_external_referent=pr['additional_external_referent'], decision='IMPLICIT_SPECIMEN_PASS' if all(gates.values()) else 'IMPLICIT_SPECIMEN_FAIL', gates=gates, prediction=p, voxel_distance_error_mm=err, voxel_distance_range_mm=[float(distance.min()), float(distance.max())], injected_no_erosion_distance_range_mm=[float(bad_distance.min()), float(bad_distance.max())], prep_mesh_watertight=bool(mesh.is_watertight), nonmanifold_edges=nonmanifold, vertical_ray_max_intervals=raymax, external_retention_error_mm=exterror, sampled_source_triangle_distance_mm=dict(n=len(dd), min=float(min(dd)), p50=float(np.median(dd)), p95=float(np.quantile(dd, 0.95)), resolution='PER_POINT', status='Geometric diagnostic against raw predicted tooth surface, includes extraction boundaries; not an actual preparation thickness measurement'), uncertainty=dict(voxel_spacing_mm=step, lattice_discretization_debt_mm=pr['metrics']['voxel_sdf_error_debt_mm'], continuous_boundary_error_bound='UNKNOWN; height envelope interpolation/source-boundary aliasing and scan/FDI error not included in lattice debt', cervical_margin='UNKNOWN', pulp_remaining_dentin='UNKNOWN', clinical_or_material_risk='UNKNOWN'), cost=dict(validation_wall_s=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    dump(H / 'rounds/R4.json', out)
    (H / 'HANDOFF_R4.md').write_text(f"R4 {out['decision']}: Euclidean-eroded insertion-envelope preparation specimen, watertight={out['prep_mesh_watertight']}, ray intervals={raymax}, lattice depth error={err:.6g}mm. Exact R1 exterior retained. Virtual basal closure and sampled source envelope are declared closures; actual cervical margin/pulp/preparation and physical pressure unknown. Next: independent cervical/FDI annotation, registered pre/post preparation scans and same-time pressure map/resultant measurement.\n")
    state('R4_DECIDED', out['decision'], 'Package executable demo and exact unresolved measurement ports')
    print(json.dumps(dict(decision=out['decision'], gates=gates, nonmanifold_edges=nonmanifold, depth_error_mm=err), indent=2))
if __name__ == '__main__':
    {'predict': predict, 'evaluate': evaluate}[sys.argv[1]]()
