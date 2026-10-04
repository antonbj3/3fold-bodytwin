"""Scientific checks and deliberately corrupted inputs; no physical fixture as truth."""
import hashlib, json, math, pathlib, copy
import numpy as np
from collections import Counter
from model import feature
from calibration import monotone_inverse, lagrange
from enclosures import certified_derivative
from query import query
ROOT = pathlib.Path(__file__).resolve().parents[1]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
read = lambda n: json.loads((ROOT / n).read_text())

def source_row_valid(row, external_value, tolerance=1e-12):
    return bool(row['locator'] and row['torque_Ncm'] > 0 and (abs(row['torque_Ncm'] - external_value) <= tolerance))

def main():
    checks = {}
    faults = {}
    rows = read('raw/arms.json')
    primary = {r['id']: r for r in rows if r['primary']}
    allrows = {r['id']: r for r in rows}
    checks['prereg_and_prediction_receipts'] = all((sha(ROOT / f'PREREG_{r}.json') == (ROOT / f'PREREG_{r}.sha256').read_text().strip() and sha(ROOT / f'FROZEN_PREDICTIONS_{r}.json') == read(f'FROZEN_PREDICTIONS_{r}.sha256.json')['sha256'] and (sha(ROOT / 'raw/arms.json') == read(f'FROZEN_PREDICTIONS_{r}.json')['data_sha256']) for r in ['R1', 'R2', 'R3', 'R4']))
    sources = read('raw/source_manifest.json')
    checks['source_hashes'] = all((sha(pathlib.Path(s['path'])) == s['sha256'] for s in sources.values()))
    example = next(iter(sources.values()))
    b = pathlib.Path(example['path']).read_bytes()
    faults['tampered_source_hash'] = hashlib.sha256(b + b'corruption').hexdigest() != example['sha256']
    anchors = {'B_T2_2.7': 39.2, 'B_T2_3': 35.4, 'B_T2_3.3': 26.0, 'B_T2_3.6': 16.2, 'B_T5_3': 36.4, 'B_T5_3.3': 32.8, 'B_T5_3.6': 24.8, 'B_T5_3.8': 16.2}
    checks['publisher_table_anchors'] = all((source_row_valid(allrows[i], t) for (i, t) in anchors.items()))
    wrong = copy.deepcopy(allrows['B_T2_3.3'])
    wrong['torque_Ncm'] *= 10
    faults['wrong_external_torque'] = not source_row_valid(wrong, 26.0)
    isq = read('raw/ISQ_SOURCE_CHECK.json')['rows']
    actual = read('results_R3.json')['summary_sufficiency']
    checks['published_equal_torque_ISQ_anchor'] = actual['total_torque_Ncm'] == [v['torque_mean_Ncm'] for v in isq] and actual['observed_ISQ'] == [v['ISQ_mean'] for v in isq]
    faults['wrong_ISQ_anchor'] = actual['observed_ISQ'] != [10 * v['ISQ_mean'] for v in isq]
    frozen = read('FROZEN_PREDICTIONS_R1.json')
    maxqr = 0.0
    leakfree = True
    for entry in frozen['models']:
        (split, heldout) = (entry['split'], entry['heldout'])
        train = [primary[i] for i in entry['train_ids']]
        test = [r for r in primary.values() if r[split] == heldout]
        leakfree &= all((r[split] != heldout for r in train))
        if entry['method'] != 'equally_informed_WLS':
            continue
        m = entry['model']
        raw = feature(train, m['kind'])
        X = np.column_stack([np.ones(len(train)), (raw - np.array(m['mean'])) / np.array(m['std'])])
        count = Counter((r['study'] for r in train))
        w = np.array([1 / count[r['study']] for r in train])
        w = w / w.sum() * len(count)
        y = np.log([r['torque_Ncm'] / (r['D_mm'] * r['L_mm']) for r in train])
        beta = np.linalg.lstsq(np.sqrt(w)[:, None] * X, np.sqrt(w) * y, rcond=None)[0]
        Xt = np.column_stack([np.ones(len(test)), (feature(test, m['kind']) - m['mean']) / m['std']])
        pred = np.exp(Xt @ beta) * [r['D_mm'] * r['L_mm'] for r in test]
        old = {v['id']: v['predicted_Ncm'] for v in frozen['predictions'] if v['split'] == split and v['heldout'] == heldout and (v['method'] == entry['method']) and (v['kind'] == m['kind'])}
        maxqr = max(maxqr, max((abs(float(p) - old[r['id']]) / old[r['id']] for (r, p) in zip(test, pred))))
    checks['no_fold_training_leakage'] = bool(leakfree)
    checks['independent_weighted_design_WLS_agreement'] = maxqr <= 1e-08
    faults['wrong_equal_information_WLS_output'] = max((abs(2 * float(p) - old[r['id']]) / old[r['id']] for (r, p) in zip(test, pred))) > 1e-08
    r1 = read('results_R1.json')

    def study_metric(split, kind, mult=1.0):
        errors = {}
        for v in frozen['predictions']:
            if v['split'] == split and v['kind'] == kind and (v['method'] == 'candidate_ridge'):
                errors.setdefault(v['study'], []).append(abs(math.log(mult * v['predicted_Ncm'] / primary[v['id']]['torque_Ncm'])))
        return float(np.mean([np.median(e) for e in errors.values()]))
    checks['R1_metric_reproduction'] = all((abs(study_metric(s, k) - r1['scores'][s]['candidate_ridge'][k]) <= 1e-12 for s in ['study', 'family'] for k in ['practice', 'cortex', 'drilling', 'both']))
    faults['corrupted_R1_gain_and_absolute_gate'] = all((study_metric(s, 'both', 10) > math.log(1.5) for s in ['study', 'family']))
    r2 = read('FROZEN_PREDICTIONS_R2.json')

    def max_error(field, mult=1.0):
        return max((abs(math.log(mult * v[field] / allrows[v['id']]['torque_Ncm'])) for v in r2['predictions']))
    checks['R2_candidate_frozen_absolute_gate'] = max_error('predicted_Ncm') <= math.log(1.2)
    checks['R2_frozen_information_gain_gate'] = np.mean([abs(math.log(v['predicted_Ncm'] / allrows[v['id']]['torque_Ncm'])) for v in r2['predictions']]) <= 0.5 * np.mean([abs(math.log(v['practice_Ncm'] / allrows[v['id']]['torque_Ncm'])) for v in r2['predictions']])
    faults['corrupted_R2_candidate'] = max_error('predicted_Ncm', 2) > math.log(1.2)
    checks['R2_ordinary_control_absolute_check'] = max_error('ordinary_linear_Ncm') <= math.log(1.2)
    faults['corrupted_ordinary_control'] = max_error('ordinary_linear_Ncm', 2) > math.log(1.2)
    r3 = read('FROZEN_PREDICTIONS_R3.json')
    cert = []
    eq = []
    for v in r3['predictions']:
        rs = [allrows[i] for i in v['calibration_ids']]
        D = rs[0]['D_mm']
        df = allrows[v['test_id']]['Df_mm']
        x = [D * D - r['Df_mm'] ** 2 for r in rs]
        y = [r['torque_Ncm'] for r in rs]
        expected = lagrange(x, y, D * D - df * df)
        eq.append(abs(expected - v['predicted_Ncm']) / abs(expected))
        cert.append(certified_derivative(v['poly_c_b_a'], D, min((r['Df_mm'] for r in rs)), max((r['Df_mm'] for r in rs))))
    checks['R3_Lagrange_check'] = max(eq) <= 1e-10
    checks['R3_frozen_heldout_torque_gate'] = all((abs(math.log(v['predicted_Ncm'] / allrows[v['test_id']]['torque_Ncm'])) <= math.log(1.2) for v in r3['predictions']))
    faults['corrupted_Lagrange_control'] = abs(2 * v['predicted_Ncm'] - expected) / abs(expected) > 1e-10
    checks['R3_certified_derivative_gate'] = all((v['certified_nonincreasing_fitted_torque'] for v in cert))
    v = r3['predictions'][0]
    rs = [allrows[i] for i in v['calibration_ids']]
    reverse = certified_derivative([-x for x in v['poly_c_b_a']], rs[0]['D_mm'], 2.7, 3.6)
    faults['reversed_polynomial_derivative'] = not reverse['certified_nonincreasing_fitted_torque']
    bad = monotone_inverse([(2.7, 39.2), (3.0, 35.4), (3.6, 392.0)], 25.0)
    faults['nonmonotone_calibration'] = bad['status'] == 'INCONSISTENT'
    points = [(3.6, 16.2), (2.7, 39.2), (3.0, 35.4)]
    bounds = [(15.0, 17.0), (38.0, 40.0), (34.0, 36.0)]
    interval_query = monotone_inverse(points, 25.0, bounds)
    checks['unsorted_interval_bounds_keep_their_measurements'] = interval_query['Df_outer_interval_mm'] == [3.0, 3.6]
    checks['pressure_stays_unidentified'] = interval_query['radial_pressure_MPa'] is None
    r4 = read('FROZEN_PREDICTIONS_R4.json')
    checks['R4_finite_menu_query_gate'] = sum((20 <= allrows[v['selected']['id']]['torque_Ncm'] <= 30 for v in r4['predictions'])) / len(r4['predictions']) >= 0.8
    checks['R4_uses_only_calibration_response_for_choice'] = all((v['selected']['id'] == min(v['menu'], key=lambda x: (abs(math.log(x['predicted_Ncm'] / 25.0)), -x['Df_mm']))['id'] and v['selected']['id'] not in v['endpoint_ids'] for v in r4['predictions']))
    faults['corrupted_R4_selected_torque'] = all((not 20 <= 10 * allrows[v['selected']['id']]['torque_Ncm'] <= 30 for v in r4['predictions']))
    inputs = read('raw/LAB_QUERY_INPUT.json')
    queries = [query(v) for v in inputs]
    checks['lab_query_reproduces_frozen_discrete_choices'] = [q['selected']['Df_mm'] for q in queries] == [v['selected']['Df_mm'] for v in r4['predictions']]
    wrong = copy.deepcopy(inputs[0])
    wrong['same_regime_confirmed'] = False
    faults['lab_query_unknown_regime_rejected'] = query(wrong)['status'] == 'UNKNOWN_REGIME_COMPATIBILITY'
    checks['all_injected_faults_rejected'] = all(faults.values())
    out = dict(status='PASS' if all(checks.values()) else 'FAIL', checks=checks, injected_faults_rejected=faults, independent_WLS_max_relative_difference=maxqr, independent_Lagrange_max_relative_difference=max(eq), numerical_physics_scope='Checks certify arithmetic, provenance, holdouts and frozen criteria; physical/model uncertainty is not certified', external_referent=dict(kind='independent_measurement', locator='https://jpis.org/pdf/10.5051/jpis.2013.43.1.30 pp.33-34', compared_quantity='Published intermediate and endpoint insertion torque means Ncm', refutes_us=True))
    (ROOT / 'VERIFICATION.json').write_text(json.dumps(out, indent=2, allow_nan=False, default=lambda v: v.item()) + '\n')
    print(json.dumps(out, indent=2, default=lambda v: v.item()))
    assert all(checks.values()), 'Verification failed; inspect VERIFICATION.json'
if __name__ == '__main__':
    main()
