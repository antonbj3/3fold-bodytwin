from dental_release.paths import expand as _release_expand
import datetime, hashlib, json, resource, time
from pathlib import Path
import numpy as np
from full_geometry import write_stl, project_cylinder
from printable_geometry import mesh_from_mask, topology, triangle_cylinder_bracket, inside_mesh
ROOT = Path(__file__).resolve().parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X8-guide-nerve-risk'))

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2) + '\n')
start = time.perf_counter()
pr = json.loads((ROOT / 'PREREG_R6.json').read_text())
assert sha(ROOT / 'PREREG_R6.json') == (ROOT / 'PREREG_R6.sha256').read_text().strip()
assert sha(ROOT / 'PHANTOM_SCENE.json') == pr['source_hashes']['scene']
assert sha(DATA / 'P394_44_CROPPED_SDF.npz') == pr['source_hashes']['sdf_npz']
with np.load(DATA / 'P394_44_CROPPED_SDF.npz') as arr:
    mask = arr['canal_mask']
    sp = arr['spacing_zyx_mm']
scene = json.loads((ROOT / 'PHANTOM_SCENE.json').read_text())
E = np.array(scene['entry_local_xyz_mm'])
a = np.array(scene['axis_local_xyz'])
L = scene['length_mm']
radius = scene['radius_mm']
dest = ROOT / 'phantom_manifold'
dest.mkdir(exist_ok=True)
meshes = {}
stats = {}
for (name, occupancy) in [('CANAL_LUMEN.stl', mask), ('PHANTOM_BLOCK_WITH_LUMEN.stl', ~mask)]:
    tri = mesh_from_mask(occupancy, sp)
    write_stl(dest / name, tri)
    stats[name] = topology(tri)
    if name == 'CANAL_LUMEN.stl':
        canal = tri
    print(name, stats[name], flush=True)
gap = triangle_cylinder_bracket(canal, E, a, L, radius)
mid = E + 0.5 * L * a
contained = inside_mesh(mid, canal)
topopass = all((all((s[k] == 0 for k in ['boundary_edges', 'edges_incidence_not2', 'orientation_inconsistent_edges', 'bad_vertex_links', 'duplicate_triangles', 'zero_area_triangles'])) for s in stats.values()))
old = np.array(scene['whole_annotated_canal_clearance_mm'])
change = float(max(abs(gap['lower_mm'] - old[0]), abs(gap['upper_mm'] - old[1])))
bad = topology(np.concatenate([canal, canal[:1]]))
bad_topology_rejected = bad['edges_incidence_not2'] > 0 and bad['duplicate_triangles'] > 0
w = min(gap['records'], key=lambda r: r['upper_mm'])
point = np.array(w['point_xyz_mm'])
d = float(np.linalg.norm(point - project_cylinder(point, E, a, L, radius)))
faults = {'duplicated_facet_rejected': bad_topology_rejected, 'wrong_lower_plus_.1_rejected': gap['lower_mm'] + 0.1 > d, 'wrong_upper_minus_.1_rejected': gap['upper_mm'] - 0.1 < gap['lower_mm'], 'wrong_fidelity_delta_plus_1_rejected': change + 1 > 0.3}
assert all(faults.values())
out = {'construction': 'R6', 'topology': stats, 'gates': {'closed_edge_manifold': topopass, 'vertex_link_manifold': all((s['bad_vertex_links'] == 0 for s in stats.values())), 'clearance_numeric_bracket': gap['gap_mm'] <= 0.0001, 'annotation_clearance_preservation': change <= 0.3, 'manufacture_metrology': 'UNKNOWN'}, 'signed_field': 'Signed grid-centre EDT with outside padding; piecewise affine Freudenthal6-tetrahedra extraction', 'continuous_cylinder_to_exported_canal_surface_distance': gap, 'cylinder_midpoint_inside_canal': contained, 'interpretation': 'Surface gap equals solid gap here only with positive surface separation and midpoint outside. Not a general inside-solid distance API.', 'source_voxel_gap_mm': old.tolist(), 'maximum_gap_change_mm': change, 'practice2mm_class': 'ABOVE' if gap['lower_mm'] >= 2 else 'BELOW' if gap['upper_mm'] < 2 else 'UNKNOWN', 'controls': {'method_outcome': 'TIE: conventional tetrahedral isosurface and convex triangle minimization; new consumer is audited printable representation', 'injected_faults': faults}, 'external_referent': pr['external_referent'], 'cost': {'wall_seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1, 'gpu': False, 'physical_measurements': 0}, 'limits': ['Manifold/closed topology does not certify fabrication', 'Signed field smoothing changes the source voxel geometry', 'Cylinder STL polygonization separately adds<=0.002410mm ideal geometric discrepancy', 'Actual anatomy, metrology, drill-path, clinical risk UNKNOWN']}
assert not contained
dump('RAW_R6.json', out)
pred = {'frozen_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prereg_sha256': sha(ROOT / 'PREREG_R6.json'), 'expected_ideal_cylinder_to_manifold_canal_gap_mm': [gap['lower_mm'], gap['upper_mm']], 'files': {p.name: sha(p) for p in dest.glob('*.stl')}, 'source_voxel_gap_mm': old.tolist(), 'source_gap_change_mm': change, 'physical_measurements_performed': False, 'manufactured_prediction': 'UNKNOWN until process and metrology are fixed; frozen values refer to CAD surface'}
if (ROOT / 'FROZEN_MANIFOLD_PHANTOM_PREDICTIONS.json').exists():
    previous = json.loads((ROOT / 'FROZEN_MANIFOLD_PHANTOM_PREDICTIONS.json').read_text())
    assert previous['files'] == pred['files']
    assert previous['expected_ideal_cylinder_to_manifold_canal_gap_mm'] == pred['expected_ideal_cylinder_to_manifold_canal_gap_mm']
else:
    dump('FROZEN_MANIFOLD_PHANTOM_PREDICTIONS.json', pred)
dump('CURRENT_WORK_STATE.json', {'lane': 'X8-guide-nerve-risk', 'status': 'R6_COMPLETE', 'latest_gate': out['gates'], 'next_operation': 'Independent slicer/material/process review, fabrication and signed gap metrology; no clinical risk calibration'})
(ROOT / 'HANDOFF_R6.md').write_text('R6 completed; RAW_R6.json records all fixed topology/fidelity gates and CAD prediction. Originals and earlier freezes retained. Manifold files under phantom_manifold/, frozen before any manufacture. The smoothed reference has a different clearance than the source voxel volume; do not compare a printed smoothed part to the original ideal-voxel prediction. Next: material/process/registration and independent metrology. Clinical risk UNKNOWN.\n')
print(json.dumps({k: v for (k, v) in out.items() if k != 'continuous_cylinder_to_exported_canal_surface_distance'}, indent=2))
