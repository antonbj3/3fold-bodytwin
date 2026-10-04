import math, json, hashlib
from science import R, load, dump

def valid_quantity_link(source, target):
    return all((source.get(k) == target.get(k) for k in ['quantity', 'unit', 'resolution', 'coverage_unit', 'lineage_kind']))
wall = {'quantity': 'local_signed_wall_clearance_loss', 'unit': 'mm', 'resolution': 'PER_POINT', 'coverage_unit': 'patient', 'lineage_kind': 'independent_measurement'}
bad_sources = [dict(wall, quantity='SMCD_mean', resolution='PER_ARCH'), dict(wall, quantity='tracing_path_P95', coverage_unit='within_canal_points'), dict(wall, quantity='Dice', unit='dimensionless'), dict(wall, lineage_kind='mask_release_revision'), dict(wall, quantity='pooled_mean_CI'), dict(wall, quantity='endpoint_proxy_landmark')]
checks = []
assert valid_quantity_link(wall, wall)
for src in bad_sources:
    assert not valid_quantity_link(src, wall)
    checks.append({'check': 'Quantity and coverage link', 'injected_source': src, 'requested_target': wall, 'injection_rejected': not valid_quantity_link(src, wall)})

def physical_claim_valid(e):
    return bool(e.get('independent_wall_pairs', 0) > 0 and e.get('segment_landmarks_validated', False) and e.get('joint_guide_loss_validated', False) and e.get('physical_coverage_validated', False))
e = {'independent_wall_pairs': 0, 'segment_landmarks_validated': False, 'joint_guide_loss_validated': False, 'physical_coverage_validated': False}
assert not physical_claim_valid(e)
bad = dict(e, physical_coverage_validated=True)
assert not physical_claim_valid(bad)
checks.append({'check': 'Physical claim admission', 'injected_physical_true_rejected': not physical_claim_valid(bad), 'evidence': bad})
m = load('INPUT_MANIFEST.json')[0]
b = (R / m['local']).read_bytes()
assert hashlib.sha256(b).hexdigest() == m['sha256']
bad_hash = hashlib.sha256(b + b'bad').hexdigest()
assert bad_hash != m['sha256']
checks.append({'check': 'Input provenance', 'injected_changed_bytes_rejected': bad_hash != m['sha256']})
split = load('inputs/patient_split.json')
c = set(split['calibration_patients'])
t = set(split['test_patients'])
assert not c & t
bad_c = c | {next(iter(t))}
assert bad_c & t
checks.append({'check': 'Image-group split', 'injected_cross_split_alias_rejected': bool(bad_c & t)})
r2 = load('raw/R2_RESULT.json')
r4 = load('raw/R4_RESULT.json')
errors = []
bad_offsets = []
for (a, b) in zip(r2['tables'], r4['tables']):
    if isinstance(a['margin_mm'], (int, float)):
        errors.append(abs(a['margin_mm'] - b['margin_mm'] - 2))
        bad_offsets.append(abs(a['margin_mm'] + 2 - b['margin_mm'] - 2) > 1e-10)
assert max(errors) <= 1e-10
assert all(bad_offsets)
checks.append({'check': 'Residual target explicit', 'R2_minus_R4_mm': 2.0, 'max_error_mm': max(errors), 'injected_unlabelled_plus2_rejected': all(bad_offsets), 'injected_cases': len(bad_offsets)})
gaps = [3.0, 6.0]
A = [2.0, 0.0]
B = [0.0, 2.0]
summary = lambda loss: [sum(loss) / len(loss), sorted(loss), [min(loss), max(loss)]]
assert summary(A) == summary(B)
down = lambda loss: [g - l >= 2.0 for (g, l) in zip(gaps, loss)]
assert down(A) != down(B)
out = {'identical_summary': summary(A), 'identity_error': 0.0, 'gap_mm': gaps, 'loss_A_mm': A, 'loss_B_mm': B, 'accepted_A': down(A), 'accepted_B': down(B), 'accepted_fraction_difference': 0.5, 'resolution': 'PER_TOOTH', 'minimum_extension': 'Keep loss matched to each frozen pose and its active wall witness', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'validation.py', 'compared_quantity': 'Logical insufficiency of unpaired regional loss summaries', 'refutes_us': True}}
dump('raw/SUFFICIENCY_SITE_PAIRING.json', out)
assert load('raw/R3_RESULT.json')['envelope_controls']['injected_loss_decreased_by1mm_rejected']
dump('raw/VALIDATION.json', {'status': 'PASS_WITH_PHYSICAL_UNKNOWN', 'checks': checks, 'site_summary_test': out, 'formal_floating_point_enclosure': 'MISSING'})
print('Validation checks', len(checks), 'summary identity', 0.0, 'downstream accepted fraction difference', 0.5)
