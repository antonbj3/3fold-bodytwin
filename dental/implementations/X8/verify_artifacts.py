"""Verify numerical certificates and force each critical guard to reject a fault."""
import csv, hashlib, json
from pathlib import Path
import numpy as np
from full_geometry import project_cylinder
from risk_operator import body_bound, scalar_control
ROOT = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
for i in range(1, 7):
    assert sha(ROOT / f'PREREG_R{i}.json') == (ROOT / f'PREREG_R{i}.sha256').read_text().strip()
r5 = json.loads((ROOT / 'RAW_R5.json').read_text())
raw = Path(r5['raw_certificates']['path'])
assert sha(raw) == r5['raw_certificates']['sha256']
records = [json.loads(x) for x in raw.read_text().splitlines()]
checks = 0
failed_faults = 0
for row in records:
    E = np.array(row['entry_zyx_mm'])
    a = np.array(row['axis_zyx'])
    L = row['length_mm']
    r = row['radius_mm']
    point = np.array(row['witness']['witness_zyx_mm'])
    feasible = float(np.linalg.norm(point - project_cylinder(point, E, a, L, r)))
    assert row['lower_mm'] <= feasible + 1e-08
    assert abs(feasible - row['upper_mm']) <= 1e-08
    assert row['gap_mm'] <= 0.0001
    assert row['lower_mm'] + 0.1 > feasible
    assert row['upper_mm'] - 0.1 < row['lower_mm']
    for box in row['records']:
        pp = np.array(box['witness_zyx_mm'])
        delta = pp - project_cylinder(pp, E, a, L, r)
        d = float(np.linalg.norm(delta))
        assert abs(d - box['upper_mm']) <= 1e-08
        if d > 0:
            assert np.linalg.norm(delta / d - np.array(box['normal_zyx'])) <= 1e-08
        assert box['lower_mm'] <= box['upper_mm'] + 1e-08
        checks += 1
profiles = json.loads((ROOT / 'GUIDE_PROFILES.json').read_text())['profiles']
bound_errors = [abs(float(body_bound(8, p)) - scalar_control(8, p)) for p in profiles]
assert max(bound_errors) < 1e-10
assert all((abs(float(body_bound(8, p)) + 0.1 - scalar_control(8, p)) > 1e-10 for p in profiles))

def clinical_output_guard(row):
    if row.get('clinical_risk') != 'UNKNOWN' and (not all((row.get(k) == 'True' for k in ['true_anatomy_validated', 'tail_population_validated', 'paired_sensory_outcome_validated']))):
        raise ValueError('Clinical numeric risk lacks measured prerequisites')
with (ROOT / 'PER_SITE_FULL_GEOMETRY_RISK.csv').open() as f:
    rr = list(csv.DictReader(f))
for row in rr:
    clinical_output_guard(row)
bad = {**rr[0], 'clinical_risk': '0.0001'}
try:
    clinical_output_guard(bad)
except ValueError:
    clinical_fault = True
else:
    clinical_fault = False
assert clinical_fault
frozen = [('FROZEN_PREDICTIONS.json', 'PER_SITE_RISK.csv'), ('FROZEN_FULL_GEOMETRY_PREDICTIONS.json', 'PER_SITE_FULL_GEOMETRY_RISK.csv')]
for (name, target) in frozen:
    j = json.loads((ROOT / name).read_text())
    assert sha(ROOT / target) == j['per_site_prediction_sha256']
pf = json.loads((ROOT / 'FROZEN_PHANTOM_PREDICTIONS.json').read_text())
assert sha(ROOT / 'PHANTOM_SCENE.json') == pf['scene_sha256']
for (name, h) in pf['files'].items():
    assert sha(ROOT / 'phantom' / name) == h
pm = json.loads((ROOT / 'FROZEN_MANIFOLD_PHANTOM_PREDICTIONS.json').read_text())
for (name, h) in pm['files'].items():
    assert sha(ROOT / 'phantom_manifold' / name) == h
r6 = json.loads((ROOT / 'RAW_R6.json').read_text())
scene = json.loads((ROOT / 'PHANTOM_SCENE.json').read_text())
E = np.array(scene['entry_local_xyz_mm'])
a = np.array(scene['axis_local_xyz'])
L = scene['length_mm']
radius = scene['radius_mm']
for box in r6['continuous_cylinder_to_exported_canal_surface_distance']['records']:
    verts = np.array(box['vertices_xyz_mm'])
    uv = np.array(box['barycentric_uv'])
    assert np.all(uv >= -1e-10) and uv.sum() <= 1 + 1e-10
    p = verts[0] + uv @ (verts[1:] - verts[0])
    assert np.linalg.norm(p - np.array(box['point_xyz_mm'])) < 1e-09
    delta = p - project_cylinder(p, E, a, L, radius)
    d = float(np.linalg.norm(delta))
    normal = delta / d if d else np.zeros(3)
    lower = max(0.0, d + float(((verts - p) @ normal).min()))
    assert abs(lower - box['lower_mm']) < 1e-08 and abs(d - box['upper_mm']) < 1e-08
assert r6['controls']['injected_faults']['duplicated_facet_rejected']
result = {'status': 'PASS', 'frozen_preregs_checked': 6, 'sites_checked': len(records), 'convex_box_witnesses_checked': checks, 'voxel_and_manifold_phantom_freezes_checked': True, 'manifold_triangle_witnesses_checked': len(r6['continuous_cylinder_to_exported_canal_surface_distance']['records']), 'bound_scalar_control_max_abs_error': max(bound_errors), 'injected_wrong_distance_rejected_at_every_site': True, 'injected_wrong_probability_bound_rejected': True, 'injected_clinical_numeric_risk_rejected': clinical_fault, 'scope': 'Producer numerical and artifact verification only; no independent scientific review, anatomy or clinical validity'}
(ROOT / 'VERIFICATION.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
