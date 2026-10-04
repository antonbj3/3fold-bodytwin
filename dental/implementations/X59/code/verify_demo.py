"""Executed falsifiers for source/observation contracts, not scientific admission."""
import copy, csv, json, tempfile
from pathlib import Path
import numpy as np
from extract_sources import ROOT, sha, save, extract_cunali
from analyze_r2 import contrast, ecological_check, witness, query_x13
from analyze_r3 import ols_enclosure
from analyze_r4 import rank
import lab_bridge

def run():
    checks = []

    def record(name, ok, details=''):
        checks.append(dict(name=name, pass_=bool(ok), details=details))
        if not ok:
            raise AssertionError(name + ': ' + str(details))
    for rnd in [1, 2, 3, 4]:
        p = ROOT / f'PREREG_R{rnd}.json'
        record(f'R{rnd}_freeze_hash', sha(p) == p.with_suffix('.sha256').read_text().strip())
    rows = extract_cunali()
    base = ecological_check(rows)
    corrupt = copy.deepcopy(rows)
    corrupt[0]['ct_mean_um'] = str(float(corrupt[0]['ct_mean_um']) + 1000)
    record('external_OLS_control_rejects_1000um_source_error', ecological_check(corrupt)[0]['external_coefficient_gate'] == 'FAIL', 'Cunali Figure4 independent published coefficients are falsifier; R2 Dentsply failure is retained.')
    q = contrast(rows[0])
    (lo, hi) = map(float, q['sample_difference_sd_identified_set_um'])
    actual_sds = []
    z = np.array([-1.0, -1.0, 1.0, 1.0])
    z /= z.std(ddof=1)
    for paired_ct in [z, z[::-1]]:
        actual_sds.append(float(np.std(float(rows[0]['replica_sd_um']) * z - float(rows[0]['ct_sd_um']) * paired_ct, ddof=1)))
    record('sharp_covariance_endpoints_reproduced', abs(min(actual_sds) - lo) < 1e-12 and abs(max(actual_sds) - hi) < 1e-12)
    record('covariance_enclosure_rejects_injected_2x_upper_SD', not lo <= hi * 2 <= hi)
    w = witness()
    record('same_summary_sufficiency_exact', w['identity_max_abs_error_um'] == 0 and w['downstream_classification_error_difference'] == 4)
    xx = json.loads((ROOT / 'raw/X13_INPUT_SNAPSHOT.json').read_text())
    donors = [contrast(r) for r in rows]
    target = copy.deepcopy(next((r for r in xx if r['method'] == 'silicone_replica' and r['material'] == 'zirconia' and (r['restoration'] == 'crown'))))
    original_target = copy.deepcopy(target)
    target['material'] = 'resin'
    control = query_x13([original_target], donors, persist=False)
    rejected = query_x13([target], donors, persist=False)
    record('material_scope_rejects_injected_wrong_material', control['conditional_eligible_cells'] == 1 and rejected['conditional_eligible_cells'] == 0)
    target = copy.deepcopy(original_target)
    target['method'] = 'direct_microscopy'
    rejected = query_x13([target], donors, persist=False)
    record('method_scope_rejects_wrong_operator', rejected['conditional_eligible_cells'] == 0)
    r3 = json.loads((ROOT / 'RESULTS_R3.json').read_text())
    record('source_rounding_enclosure_checked', all((r['independent_corner_checks'] == 256 for r in r3['rounding_ols'])))
    record('rounding_enclosure_rejects_1000um_false_intercept', all((not r['intercept_rounding_enclosure_um'][0] <= 1000 <= r['intercept_rounding_enclosure_um'][1] for r in r3['rounding_ols'])))
    obs = [[1, 0, 0, 1, 0], [1, 1, 0, 1, 0], [1, 1, 1, 1, 0], [1, 1, 1, 0, 1]]
    record('rank_control_refutes_false_absolute_gap_identification', rank(obs + [[1, 0, 0, 0, 0]]) > rank(obs))

    def must_reject(name, fn):
        try:
            fn()
        except ValueError as e:
            record(name, True, str(e))
        else:
            record(name, False, 'Corruption was accepted')
    must_reject('lab_no_actual_data_refused', lambda : lab_bridge.fit([]))
    fixtures = []
    for i in range(10):
        fixtures.append(dict(specimen_id=f'CONTRACT_FIXTURE_{i}', protocol_id='CONTRACT_TEST', region='marginal', point_id='p1', batch_id='CALIBRATION_FIXTURE_BATCH', point_pair_verified='true', ct_state_bias_calibrated='true', single_observation_bound_um='0', dry_ct_um=str(80 + i), in_situ_replica_ct_um=str(90 + i), sectioned_replica_ct_um=str(92 + i), sectioned_replica_optical_um=str(95 + i)))
    bad = copy.deepcopy(fixtures)
    bad[0]['point_pair_verified'] = 'false'
    must_reject('lab_unregistered_point_refused', lambda : lab_bridge.fit(bad))
    bad = copy.deepcopy(fixtures)
    bad[0]['ct_state_bias_calibrated'] = 'false'
    must_reject('lab_CT_closure_refused', lambda : lab_bridge.fit(bad))
    bad = copy.deepcopy(fixtures)
    bad[0]['dry_ct_um'] = 'nan'
    must_reject('lab_nan_refused', lambda : lab_bridge.fit(bad))
    bad = copy.deepcopy(fixtures)
    bad[0]['single_observation_bound_um'] = '-1'
    must_reject('lab_negative_bound_refused', lambda : lab_bridge.fit(bad))
    bad = copy.deepcopy(fixtures)
    bad[0]['batch_id'] = ''
    must_reject('lab_missing_batch_identity_refused', lambda : lab_bridge.fit(bad))
    must_reject('lab_point_pseudoreplication_refused', lambda : lab_bridge.fit([dict(r, specimen_id='ONE_COPING', point_id=str(i)) for (i, r) in enumerate(fixtures)]))
    with tempfile.TemporaryDirectory(prefix='contract_fixture_', dir=ROOT) as td:
        td = Path(td)

        def write(name, rows):
            p = td / name
            with p.open('w') as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            return p
        cal = write('cal.csv', fixtures)
        design = write('design.csv', [dict(specimen_id='HELD_CONTRACT_FIXTURE', protocol_id='CONTRACT_TEST', region='marginal', batch_id='HELD_FIXTURE_BATCH', replica_mean_um='120')])
        frozen = td / 'frozen.json'
        lab_bridge.freeze(cal, design, frozen)
        must_reject('lab_cannot_overwrite_predictions', lambda : lab_bridge.freeze(cal, design, frozen))
        leakage = write('leakage.csv', [dict(specimen_id='HELD2', protocol_id='CONTRACT_TEST', region='marginal', batch_id='HELD_FIXTURE_BATCH', replica_mean_um='120', ct_mean_um='105')])
        must_reject('lab_validation_leakage_refused', lambda : lab_bridge.freeze(cal, leakage, td / 'leaked.json'))
        validation = write('val.csv', [dict(fixtures[0], specimen_id='HELD_CONTRACT_FIXTURE', batch_id='HELD_FIXTURE_BATCH', dry_ct_um='105', in_situ_replica_ct_um='115', sectioned_replica_ct_um='117', sectioned_replica_optical_um='120')])
        result = lab_bridge.score(frozen, validation, td / 'score.json')
        record('lab_freeze_then_score_software_contract', result['overall_gate'] == 'PASS', 'Own fixture; not external or physical evidence.')
        validation_bad = write('val_bad.csv', [dict(fixtures[0], specimen_id='HELD_CONTRACT_FIXTURE', batch_id='HELD_FIXTURE_BATCH', dry_ct_um='1105', in_situ_replica_ct_um='115', sectioned_replica_ct_um='117', sectioned_replica_optical_um='120')])
        result = lab_bridge.score(frozen, validation_bad, td / 'bad_score.json')
        record('lab_validation_control_fails_injected_1000um_CT', result['overall_gate'] == 'FAIL')
        wrong_batch = write('wrong_batch.csv', [dict(fixtures[0], specimen_id='HELD_CONTRACT_FIXTURE', batch_id='WRONG_BATCH', dry_ct_um='105', in_situ_replica_ct_um='115', sectioned_replica_ct_um='117', sectioned_replica_optical_um='120')])
        must_reject('lab_wrong_validation_batch_refused', lambda : lab_bridge.score(frozen, wrong_batch, td / 'wrong_batch_score.json'))
        frozen.write_text(frozen.read_text() + ' ')
        must_reject('lab_frozen_hash_drift_refused', lambda : lab_bridge.score(frozen, validation, td / 'drift.json'))
    save('DELIVERY_CHECK.json', dict(status='PASS', checks=checks, interpretation='Software/contract and exact-algebra checks only. Every active comparator has an executed adverse input; no independent scientific review.'))
    print(f'{len(checks)} falsifier/contract checks PASS; no scientific admission')
    return checks
if __name__ == '__main__':
    run()
