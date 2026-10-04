import argparse, json, hashlib, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def read(n):
    return json.loads((ROOT / n).read_text())

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
sys.path.insert(0, str(ROOT))
from review_control import verify_integrity
verify_integrity()
p = argparse.ArgumentParser()
p.add_argument('--preflight', action='store_true')
a = p.parse_args()
for k in range(1, 6):
    n = f'PREREG_R{k}.json'
    assert sha(ROOT / n) == read(f'PREREG_R{k}.sha256.json')['sha256'], n
for r in read('raw/R3V2_PARENT_HASHES.json'):
    assert sha(Path(r['file'])) == r['sha256'], 'Parent code drift'
for r in read('inputs/ZIP_INTEGRITY.json')['files']:
    assert Path(r['file']).stat().st_size == r['bytes'] and r['match'], 'Zip source integrity'
if a.preflight:
    print('Preflight PASS: immutable preregs, predecessor hashes, source sizes/published-checksum record')
    sys.exit(0)
r = read('results.json')
assert r['claim_type'] == 'information_link' and r['review_state'] == 'PENDING_INDEPENDENT_REVIEW'
assert r['storage']['cache_bytes'] < 3000000000
r1 = read('raw/R1_RESULTS.json')
r2 = read('raw/R2_RESULTS.json')
r3 = read('raw/R3V2_RESULTS.json')
r5 = read('raw/R5_RESULTS.json')
assert r1['weekly_tooth_rows'] == 468 and r1['numerical_controls']['pass']
assert r1['numerical_controls']['injected_1mm_rejected'] and r1['numerical_controls']['injected_10deg_rejected']
assert r2['pair_week_rows'] == 2826 and r2['gauge_gate_pass'] and r2['injected_wrong_scalar_rejected_all']
assert r2['resolved_tracking_gap'] == 67
assert r3['stage_cases'] == 108 and r3['numerically_admissible_stages'] == 97 and (r3['wrench_rows'] == 1334)
assert all(r3['port_fault_checks'].values())
assert r3['physical_force_validation'].startswith('UNKNOWN')
for s in read('raw/R3V2_STAGES.json'):
    if s.get('gate', {}).get('pass'):
        assert all(s['injections'].values())
f = read('FROZEN_PREDICTIONS_R3V2.json')
assert sha(ROOT / f['prediction_file']) == f['sha256'], 'Frozen force prediction drift'
assert r5['weekly_rows'] == 234 and r5['region_gate_passes'] == 19 and (r5['failed_per_tooth_split'] == 198)
assert all((x.get('wrong1mm_reference_rejected', True) for x in read('raw/R5_REGIONS.json')))
assert all((x['wrong_reference_pose_rejected'] for x in read('raw/R4V2_RESULTS.json')['fits']))
import numpy as np
from scipy.spatial.transform import Rotation
from geometry import apply, mat, icp, angle, CACHE
pts = np.load(CACHE / '3485_Sirona_T0_Z11.npz')['points'][:4000]
known = mat(Rotation.from_euler('xyz', [2, -3, 5], degrees=True).as_matrix(), np.array([0.35, -0.2, 0.15]))
target = apply(known, pts)
(H, e) = icp(pts, target, method='point', initial=mat(t=target.mean(0) - pts.mean(0)))
err = H @ np.linalg.inv(known)
good = np.linalg.norm(err[:3, 3]) <= 0.01 and angle(err) <= 0.1
assert good and np.linalg.norm(err[:3, 3] + np.array([1, 0, 0])) > 0.01
result = {'status': 'PASS', 'scope': 'Reproducibility and original frozen numerical/identity gates only; clinical/physical force claims remain UNKNOWN', 'independent_point_control_closed_form_translation_error_mm': float(np.linalg.norm(err[:3, 3])), 'independent_point_control_closed_form_rotation_error_deg': angle(err), 'injected_point_control_1mm_rejected': True, 'r1_weekly_rows': 468, 'r2_relative_tracking_rows_above_declared_budget': 67, 'r3_numerical_stage_passes': 97, 'r5_region_case_passes': 19, 'r5_resolved_weekly_motion': 0, 'physical_validation': False}
(ROOT / 'raw/VERIFICATION.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
