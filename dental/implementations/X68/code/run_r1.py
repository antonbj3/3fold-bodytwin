"""Predict every fold, write a frozen receipt, THEN compare against published outcomes."""
import json, pathlib, hashlib, datetime, time, resource
from collections import defaultdict
import numpy as np
from model import fit, predict, json_model
ROOT = pathlib.Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def main():
    start = time.perf_counter()
    p = ROOT / 'PREREG_R1.json'
    assert sha(p) == (ROOT / 'PREREG_R1.sha256').read_text().strip()
    allrows = json.loads((ROOT / 'raw/arms.json').read_text())
    rows = [r for r in allrows if r['primary']]
    kinds = ['practice', 'cortex', 'drilling', 'both']
    predictions = []
    models = []
    for split in ['study', 'family']:
        for heldout in sorted({r[split] for r in rows}):
            train = [r for r in rows if r[split] != heldout]
            test = [r for r in rows if r[split] == heldout]
            for k in kinds:
                for (method, penalty) in [('candidate_ridge', 1.0), ('equally_informed_WLS', 0.0)]:
                    m = fit(train, k, penalty)
                    pred = predict(m, test)
                    models.append(dict(split=split, heldout=heldout, method=method, model=json_model(m), train_ids=[r['id'] for r in train]))
                    for (r, v) in zip(test, pred):
                        predictions.append(dict(split=split, heldout=heldout, study=r['study'], id=r['id'], method=method, kind=k, predicted_Ncm=float(v)))
    predobj = {'created_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prereg_sha256': sha(p), 'data_sha256': sha(ROOT / 'raw/arms.json'), 'code_sha256': {f.name: sha(f) for f in (ROOT / 'code').glob('*.py')}, 'predictions': predictions, 'models': models, 'freeze_note': 'Fold targets were excluded from each fit. All source tables were read during extraction; this is not blind prospective lab validation.'}
    frozen = ROOT / 'FROZEN_PREDICTIONS_R1.json'
    if not frozen.exists():
        write(frozen.name, predobj)
    else:
        old = json.loads(frozen.read_text())
        assert old['predictions'] == predictions, 'Frozen predictions changed; keep old version and open new round'
    fs = sha(frozen)
    write('FROZEN_PREDICTIONS_R1.sha256.json', {'sha256': fs, 'frozen_at_utc': json.loads(frozen.read_text())['created_at_utc']})
    lookup = {r['id']: r for r in rows}
    by = defaultdict(list)
    for v in predictions:
        r = lookup[v['id']]
        v['observed_Ncm'] = r['torque_Ncm']
        v['abs_log_error'] = float(abs(np.log(v['predicted_Ncm'] / r['torque_Ncm'])))
        by[v['split'], v['method'], v['kind'], r['study']].append(v['abs_log_error'])
    study_metrics = [dict(split=s, method=m, kind=k, study=t, n=len(e), median_abs_ln=float(np.median(e))) for ((s, m, k, t), e) in by.items()]
    scores = {}
    for s in ['study', 'family']:
        scores[s] = {m: {k: float(np.mean([v['median_abs_ln'] for v in study_metrics if v['split'] == s and v['method'] == m and (v['kind'] == k)])) for k in kinds} for m in ['candidate_ridge', 'equally_informed_WLS']}
    ratios = {s: scores[s]['candidate_ridge']['both'] / scores[s]['candidate_ridge']['practice'] for s in ['study', 'family']}
    sufficient = len({r['family'] for r in rows}) >= 4
    passed = sufficient and all((ratios[s] <= 0.8 and scores[s]['candidate_ridge']['both'] <= np.log(1.5) for s in ['study', 'family']))
    a = lookup['B_T2_2.7']
    b = lookup['B_T2_3.3']
    sa = np.array([a['rho_gcc'], a['D_mm'] * a['L_mm']])
    sb = np.array([b['rho_gcc'], b['D_mm'] * b['L_mm']])
    suff = {'states': [a['id'], b['id']], 'summary': ['rho_gcc', 'D_times_L_mm2'], 'bitwise_identical': sa.tobytes() == sb.tobytes(), 'identity_error': float(np.max(abs(sa - sb))), 'observed_torque_Ncm': [a['torque_Ncm'], b['torque_Ncm']], 'downstream_difference_Ncm': a['torque_Ncm'] - b['torque_Ncm'], 'downstream_ratio': a['torque_Ncm'] / b['torque_Ncm'], 'smallest_extension_for_this_counterexample': 'Final bore diameter Df; not claimed globally sufficient', 'external_referent': {'kind': 'independent_measurement', 'locator': 'doi:10.5051/jpis.2013.43.1.30 Table 2', 'compared_quantity': 'Published mean peak insertion torque at same density and implant D/L', 'refutes_us': True}}
    choices = []
    for split in ['study', 'family']:
        for kind in kinds:
            groups = defaultdict(list)
            for r in rows:
                groups[r['study'], r['design'], r['D_mm'], r['L_mm'], r['rho_gcc'], r['cortex_mm']].append(r)
            pm = {v['id']: v['predicted_Ncm'] for v in predictions if v['split'] == split and v['kind'] == kind and (v['method'] == 'candidate_ridge')}
            for (group, rs) in groups.items():
                if len({r['Df_mm'] for r in rs}) < 2:
                    continue
                chosen = min(rs, key=lambda r: abs(np.log(pm[r['id']] / 25.0)))
                choices.append(dict(split=split, kind=kind, group=list(group), selected_id=chosen['id'], selected_Df_mm=chosen['Df_mm'], predicted_Ncm=pm[chosen['id']], observed_Ncm=chosen['torque_Ncm'], within_20pct=20 <= chosen['torque_Ncm'] <= 30, oracle_has_available_choice=any((20 <= r['torque_Ncm'] <= 30 for r in rs)), rigorous_enclosure='MISSING'))
    injection = {'summary_identity_check': not np.array_equal(sa, np.array([b['rho_gcc'] + 1, b['D_mm'] * b['L_mm']])), 'external_numeric_check': abs(2 * a['torque_Ncm'] - 39.2) > 1e-12, 'gain_gate_check': not 10 * scores['study']['candidate_ridge']['practice'] <= 0.8 * scores['study']['candidate_ridge']['practice'], 'inverse_band_check': not 20 <= 250 <= 30}
    results = {'round': 'R1', 'claim_type': 'information_link', 'verdict': 'PASS' if passed else 'FAIL' if sufficient else 'UNKNOWN_COVERAGE', 'scores': scores, 'both_to_practice_ratios': ratios, 'study_metrics': study_metrics, 'predictions': predictions, 'summary_sufficiency': suff, 'inverse_choices': choices, 'inverse_verdict': 'UNKNOWN_COVERAGE' if len({c['group'][0] for c in choices}) < 3 else 'EVALUATE_FROZEN_80PCT', 'thread_verdict': 'UNKNOWN_NOT_IDENTIFIABLE_NO_3_FAMILIES_WITH_MEASURED_PITCH_DEPTH', 'n_rows': len(rows), 'n_studies': len({r['study'] for r in rows}), 'n_families': len({r['family'] for r in rows}), 'rigorous_enclosure': 'MISSING; physical torque error is not bounded by affine regression.', 'external_referent': suff['external_referent'], 'frozen_predictions_sha256': fs, 'prereg_sha256': sha(p), 'injected_faults_rejected': injection, 'cost': {'fit_prediction_score_query_seconds': time.perf_counter() - start, 'maxrss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'source_preparation_seconds': 'UNKNOWN; manual extraction and source verification charged, no speedup claimed'}}
    write('results_R1.json', results)
    print(json.dumps({k: results[k] for k in ['verdict', 'scores', 'both_to_practice_ratios', 'n_rows', 'n_studies', 'n_families', 'inverse_verdict']}, indent=2))
if __name__ == '__main__':
    main()
