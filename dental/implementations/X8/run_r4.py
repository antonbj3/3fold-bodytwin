from dental_release.paths import expand as _release_expand
import datetime, hashlib, json, resource, sys, time, zipfile
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from full_geometry import voxel_cylinder_bracket, voxel_triangles, cylinder_triangles, write_stl
from risk_operator import body_bound
ROOT = Path(__file__).resolve().parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X8-guide-nerve-risk'))
sys.path.insert(0, str(Path(__import__('os').environ.get('DENTAL_PROJECT_ROOT', str(ROOT.parent.parent))) / 'cells/geometry'))
sys.path.insert(0, str(Path(__import__('os').environ.get('DENTAL_PROJECT_ROOT', str(ROOT.parent.parent))) / 'cells/procedure'))
import tf2_io, canal_risk

def dump(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
start = time.perf_counter()
assert sha(ROOT / 'PREREG_R4.json') == (ROOT / 'PREREG_R4.sha256').read_text().strip()
pr = json.loads((ROOT / 'PREREG_R4.json').read_text())
case = pr['case']
fdi = pr['fdi']
member = tf2_io.ROOT + '/labelsTr/' + case + '.mha'
with zipfile.ZipFile(tf2_io.ZIP) as z:
    b = z.read(member)
label_hash = hashlib.sha256(b).hexdigest()
(lab, spacing, header) = tf2_io.read_mha_bytes(b)
sp = np.array(spacing)
nv = Path(__import__('os').environ['DENTAL_PROJECT_ROOT']) / 'results/NV1_canals/per_case' / f'{case}.json'
rec = json.loads(nv.read_text())['teeth'][str(fdi)]
E = np.array(rec['crest_entry'])
a = -np.array(rec['axis'])
a /= np.linalg.norm(a)
L = pr['length_mm']
radius = pr['radius_mm']
mask = np.isin(lab, [3, 4])
centers = np.argwhere(mask) * sp
result = voxel_cylinder_bracket(centers, sp, E, a, L, radius)
rawsite = next((r for r in json.loads((DATA / 'K3_SITES.json').read_text())['sites'] if r['case'] == case and r['fdi'] == fdi))
P1 = E + (L + 4) * a
lo = np.maximum(np.floor((np.minimum(E, P1) - 14) / sp).astype(int), 0)
hi = np.minimum(np.ceil((np.maximum(E, P1) + 14) / sp).astype(int), lab.shape)
sub = mask[tuple((slice(x, y) for (x, y) in zip(lo, hi)))]
edt = ndi.distance_transform_edt(~sub, sampling=sp).astype(np.float32)

def field(P):
    return ndi.map_coordinates(edt, (P / sp - lo).T, order=1, mode='nearest') - 0.15
(e1, e2) = canal_risk.frame(a)
old = float(canal_risk.min_dist(field, E, a, e1, e2, L, np.zeros(1), radius, {'d_ent': np.zeros((1, 3)), 'tilt': np.zeros((1, 2))}, canal_risk.drill_template(radius, flat=True))[0])
replay = abs(old - rawsite['d_nom_implant_body'])
delta_min = max(result['lower_mm'] - old, old - result['upper_mm'], 0.0)
ph = ROOT / 'phantom'
ph.mkdir(exist_ok=True)
end = E + L * a
lo = np.maximum(np.floor((np.minimum(E, end) - pr['roi_margin_mm']) / sp).astype(int), 0)
hi = np.minimum(np.ceil((np.maximum(E, end) + pr['roi_margin_mm']) / sp).astype(int), lab.shape)
sl = tuple((slice(x, y) for (x, y) in zip(lo, hi)))
cm = mask[sl]
origin_zyx = lo * sp - sp / 2
entry_xyz = (E - origin_zyx)[::-1]
axis_xyz = a[::-1]
counts = {}
counts['CANAL_LUMEN.stl'] = write_stl(ph / 'CANAL_LUMEN.stl', voxel_triangles(cm, sp))
counts['PHANTOM_BLOCK_WITH_LUMEN.stl'] = write_stl(ph / 'PHANTOM_BLOCK_WITH_LUMEN.stl', voxel_triangles(~cm, sp))
counts['PLANNED_IMPLANT.stl'] = write_stl(ph / 'PLANNED_IMPLANT.stl', cylinder_triangles(entry_xyz, axis_xyz, L, radius, pr['cylinder_facets']))
sdf = (ndi.distance_transform_edt(~cm, sampling=sp) - ndi.distance_transform_edt(cm, sampling=sp)).astype(np.float32)
np.savez_compressed(DATA / 'P394_44_CROPPED_SDF.npz', canal_mask=cm, canal_grid_centre_sdf_mm=sdf, spacing_zyx_mm=sp, origin_zyx_mm=origin_zyx)
scene = {'case': case, 'fdi': fdi, 'source_zip': tf2_io.ZIP, 'source_member': member, 'source_member_sha256': label_hash, 'plan_source': str(nv), 'plan_source_sha256': sha(nv), 'units': 'mm', 'stl_axes': 'x,y,z; original computational z,y,x coordinates reversed', 'stl_origin_global_zyx_mm': origin_zyx.tolist(), 'voxel_centers_local_xyz_mm': '(index_xyz+0.5)*spacing_xyz', 'entry_local_xyz_mm': entry_xyz.tolist(), 'axis_local_xyz': axis_xyz.tolist(), 'length_mm': L, 'radius_mm': radius, 'block_size_xyz_mm': ((hi - lo) * sp)[::-1].tolist(), 'triangles': counts, 'whole_annotated_canal_clearance_mm': [result['lower_mm'], result['upper_mm']], 'sdf_note': 'Grid-centre signed EDT is an approximate SDF; exact voxel-union surfaces are exported. Do not interpret as an anatomical error certificate.', 'manufacturing_status': 'CANDIDATE_NOT_FABRICATED; voxel staircases and diagonal contacts require slicer/manufacturing review', 'clinical_use': 'NONE; simulated virtual tooth-axis plan; no patient recommendation'}
dump('PHANTOM_SCENE.json', scene)
profiles = json.loads((ROOT / 'GUIDE_PROFILES.json').read_text())['profiles']
geo_risks = [{'guide': p['id'], 'whole_label_body_probability_upper_if_moments_and_B_valid': float(body_bound(result['lower_mm'], p, 0.3)), 'true_anatomy_probability': 'UNKNOWN', 'clinical_injury_probability': 'UNKNOWN'} for p in profiles]
w = result['witness']
point = np.array(w['witness_zyx_mm'])
n = np.array(w['normal_zyx'])
center = np.array(w['center_zyx_mm'])
distance = lambda pp: np.linalg.norm(pp - __import__('full_geometry').project_cylinder(pp, E, a, L, radius))
check_lower = max(0.0, float(distance(point) + n @ (np.where(n >= 0, center - sp / 2, center + sp / 2) - point)))
upper = float(distance(point))
fault = {'claimed_lower_plus_.1_rejected_by_feasible_point': result['lower_mm'] + 0.1 > upper, 'claimed_upper_minus_.1_rejected_by_supporting_plane': result['upper_mm'] - 0.1 < check_lower, 'K3_distance_plus_.1_rejected': abs(old - (rawsite['d_nom_implant_body'] + 0.1)) > 1e-05, 'manufactured_status_without_measurement_rejected': len(list(ph.glob('*measurement*'))) == 0}
assert all(fault.values())
out = {'construction': 'R4', 'case': case, 'fdi': fdi, 'full_voxel_union_distance': result, 'old_K3_mm': old, 'archived_K3_mm': rawsite['d_nom_implant_body'], 'K3_replay_error_mm': replay, 'minimum_difference_from_full_bracket_mm': delta_min, 'gates': {'whole_geometry_gap': result['gap_mm'] <= 0.0001, 'K3_replay': replay <= 1e-05, 'sampled_full_equivalence': delta_min <= 0.1, 'manufacture_metrology_validation': 'UNKNOWN'}, 'controls': {'same_information_method': 'TIE to conventional convex box-distance optimization; no algorithm novelty claim', 'injected_faults': fault}, 'guide_conditions': geo_risks, 'external_referent': pr['external_referent'], 'phantom_scene': 'PHANTOM_SCENE.json', 'cost': {'wall_seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'gpu': False, 'threads': 1, 'physical_measurements': 0}}
dump('RAW_R4.json', out)
frozen = {'frozen_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prereg_sha256': sha(ROOT / 'PREREG_R4.json'), 'scene_sha256': sha(ROOT / 'PHANTOM_SCENE.json'), 'expected_annotated_voxel_union_gap_mm': [result['lower_mm'], result['upper_mm']], 'measurement_tolerance': 'Must be fixed by independent lab metrology before measurement; no empirical acceptance gate invented here', 'files': {p.name: sha(p) for p in ph.glob('*.stl')}, 'physical_measurements_performed': False}
if (ROOT / 'FROZEN_PHANTOM_PREDICTIONS.json').exists():
    previous = json.loads((ROOT / 'FROZEN_PHANTOM_PREDICTIONS.json').read_text())
    assert previous['files'] == frozen['files']
else:
    dump('FROZEN_PHANTOM_PREDICTIONS.json', frozen)
dump('CURRENT_WORK_STATE.json', {'lane': 'X8-guide-nerve-risk', 'status': 'R4_COMPLETE_PACKAGING', 'latest_gate': out['gates'], 'next_operation': 'Independent review of raw data and lab acquisition of signed clearance-loss; no measured clinical hazard on disk'})
print(json.dumps({k: v for (k, v) in out.items() if k != 'full_voxel_union_distance'}, indent=2))
