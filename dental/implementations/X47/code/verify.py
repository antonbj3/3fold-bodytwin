import copy, json, math, hashlib, csv
import numpy as np
from pathlib import Path
from extract import extract
from lab_port import ROOT, sha, allocation, verify_allocation
from compare_lab import score

def source_check(rows):
    fresh = {r['row_id']: r for r in extract()}
    if len(fresh) != len(rows):
        raise ValueError('GROUP_COUNT_MISMATCH')
    for r in rows:
        e = fresh[r['row_id']]
        for k in ['mean_force_N', 'sd_force_N']:
            if abs(r[k] - e[k]) > 1e-09:
                raise ValueError('SOURCE_TABLE_VALUE_MISMATCH')
        for k in ['unit', 'endpoint', 'n', 'source_sha256', 'source_path', 'doi', 'table_id', 'height_mm', 'angle_kind', 'total_convergence_deg', 'thermal_cycles', 'cement', 'substrate']:
            if r[k] != e[k]:
                raise ValueError('SOURCE_METADATA_MISMATCH:' + k)
        if sha(Path(r['source_path'])) != r['source_sha256']:
            raise ValueError('SOURCE_HASH_CHANGED')
    return True

def fixture():
    rows = allocation()
    for (i, r) in enumerate(rows):
        r.update(force_N=100 + i % 4, force_unit='N', failure_mode='cement_debonding', material_batch='TEST_ONLY', cement_product='TEST_ONLY', cement_batch='TEST_ONLY', surface_treatment='TEST_ONLY', die_material_batch='TEST_ONLY', actual_TC_cycles=r['thermal_cycles'], TC_low_C=5, TC_high_C=55, dwell_s=30, transfer_s=5, post_storage_hours=200, pull_speed_mm_min=1, measured_pull_axis_deg=0, force_trace_sha256='a' * 64, image_sha256='b' * 64)
    return rows

def run(rows):
    tests = {}
    source_check(rows)
    for (name, key, value) in [('double_force', 'mean_force_N', rows[0]['mean_force_N'] * 2), ('wrong_unit', 'unit', 'MPa'), ('wrong_endpoint', 'endpoint', 'planar_shear_bond'), ('missing_locator', 'table_id', 'MISSING'), ('wrong_TOC_convention', 'angle_kind', 'half_angle'), ('wrong_source_hash', 'source_sha256', '0' * 64)]:
        bad = copy.deepcopy(rows)
        bad[0][key] = value
        try:
            source_check(bad)
        except (ValueError, KeyError):
            tests[name] = True
        else:
            tests[name] = False
    for (name, mut) in [('reuse_ultimate', lambda x: x.append(copy.deepcopy(x[0]))), ('duplicate_specimen', lambda x: x[1].update(specimen_id=x[0]['specimen_id'])), ('wrong_endpoint_allocation', lambda x: x[0].update(ultimate_endpoint='two_destructive_endpoints'))]:
        bad = allocation()
        mut(bad)
        try:
            verify_allocation(bad)
        except ValueError:
            tests[name] = True
        else:
            tests[name] = False
    good = fixture()
    tests['valid_synthetic_score'] = score(good)['status'] == 'SCORED'
    retention_idx = next((i for (i, r) in enumerate(good) if r['ultimate_endpoint'] == 'retention'))
    aged_idx = next((i for (i, r) in enumerate(good) if r['age_state'] == 'TC10000'))
    for (name, idx, key, value) in [('wrong_lab_force_unit', retention_idx, 'force_unit', 'MPa'), ('nonaxial_pull', retention_idx, 'measured_pull_axis_deg', 30), ('wrong_cement_batch', retention_idx, 'cement_batch', 'WRONG'), ('wrong_thermal_dose', aged_idx, 'actual_TC_cycles', 9999), ('wrong_thermal_temperature', aged_idx, 'TC_high_C', 45), ('nonzero_premature_loss', aged_idx, 'failure_mode', 'premature_decementation'), ('missing_trace', retention_idx, 'force_trace_sha256', '')]:
        bad = copy.deepcopy(good)
        bad[idx][key] = value
        try:
            score(bad)
        except ValueError:
            tests[name] = True
        else:
            tests[name] = False
    bad = copy.deepcopy(good)
    bad[retention_idx]['failure_mode'] = 'grip_failure'
    tests['competing_grip_failure_not_debonding'] = any((r['status'].startswith('UNKNOWN') for r in score(bad)['designs']))
    bad = copy.deepcopy(good)
    for r in bad:
        if r['age_state'] == 'TC10000':
            r['force_N'] *= 0.4
    tests['aging_null_rejects_injected_loss'] = all((r['status'] == 'REJECT_AGING_NULL' for r in score(bad)['designs']))
    bad = copy.deepcopy(good)
    bad[aged_idx].update(failure_mode='premature_decementation', force_N=0)
    tests['premature_zero_not_dropped'] = any((r['status'] == 'UNKNOWN_INTERVAL_ZERO_OR_CENSORED' and sum(r['premature_losses']) > 0 for r in score(bad)['designs']))
    bad = copy.deepcopy(good)
    for r in bad:
        r['force_N'] = 100
    tests['zero_variance_not_certificate'] = all((r['status'] == 'UNKNOWN_MEASUREMENT_VARIANCE' for r in score(bad)['designs']))
    with (ROOT / 'LAB_MEASUREMENTS_TEMPLATE.csv').open() as f:
        template = list(csv.DictReader(f))
    tests['empty_template_pending'] = score(template)['status'] == 'PENDING_MEASUREMENT'
    r2 = json.loads((ROOT / 'rounds/RESULTS_R2.json').read_text())
    tests['GLS_parity'] = r2['model']['gls_control_max_abs'] < 1e-08
    from model import design, covariance
    eligible = [r for r in rows if r['independent_first_state'] and r['substrate'] == 'titanium' and (r['height_mm'] is not None)]
    X = design(eligible, 'reduced')
    y = np.log([r['mean_force_N'] for r in eligible])
    V = covariance(eligible, r2['model']['tau2'])
    control_beta = np.linalg.solve(X.T @ np.linalg.solve(V, X), X.T @ np.linalg.solve(V, y))

    def gls_matches(candidate):
        candidate = np.asarray(candidate, dtype=float)
        return bool(candidate.shape == control_beta.shape and np.isfinite(candidate).all() and (np.max(np.abs(candidate - control_beta)) < 1e-08))
    tests['GLS_parity'] = gls_matches(r2['model']['beta'])
    bad_beta = np.array(r2['model']['beta'])
    bad_beta[0] += 1
    tests['GLS_rejects_wrong_beta'] = not gls_matches(bad_beta)
    tests['summary_exact_identity'] = r2['sufficiency_test']['numeric_identity_error'] == 0 and r2['sufficiency_test']['summary_byte_identity']
    tests['summary_fails_downstream'] = r2['sufficiency_test']['relative_difference'] > 0.1
    for name in ['PREREG_RETENTION_R1', 'PREREG_RETENTION_R2', 'PREREG_LAB_PORT_R3', 'FROZEN_PREDICTIONS']:
        tests[name + '_frozen_hash'] = sha(ROOT / (name + '.json')) == (ROOT / (name + '.sha256')).read_text().split()[0]
    r3 = json.loads((ROOT / 'rounds/RESULTS_R3.json').read_text())
    for r in r3['geometry']:
        tests[r['design_id'] + '_area_partition'] = r['area_partition_identity_error_mm2'] <= 1e-10
        tests[r['design_id'] + '_array_hash'] = sha(ROOT / r['local_face_artifact']) == r['local_face_sha256']
    for r in r3['source_manifest']:
        tests[Path(r['path']).parent.name + '_' + Path(r['path']).name + '_source_hash'] = sha(Path(r['path'])) == r['sha256']
    result = {'pass': all(tests.values()), 'checks': tests, 'source_replayed_groups': len(rows), 'force_fixture_kind': 'our_own_fixture', 'fixture_scope': 'Fault detection only, no physical validation or external reference', 'physical_measurements_performed': False}
    (ROOT / 'VERIFICATION.json').write_text(json.dumps(result, indent=2) + '\n')
    if not result['pass']:
        raise ValueError('VERIFICATION_FAIL:' + str([k for (k, v) in tests.items() if not v]))
    return result
if __name__ == '__main__':
    r = run(json.loads((ROOT / 'raw/RETENTION_GROUPS.json').read_text()))
    print(len(r['checks']), 'checks', r['pass'])
