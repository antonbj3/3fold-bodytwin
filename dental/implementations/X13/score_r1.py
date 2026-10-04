from pathlib import Path
import json, hashlib, math, datetime, copy, time
from calibrate_gap import REGIONS, fit
P = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p, x):
    Path(p).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def balanced(values, rows):
    return sum((sum((v for (v, d) in zip(values, rows) if d['study'] == s)) / sum((d['study'] == s for d in rows)) for s in sorted(set((d['study'] for d in rows))))) / len(set((d['study'] for d in rows)))

def rmse(pred, rows):
    return math.sqrt(balanced([(p - d['measured_mean_um']) ** 2 for (p, d) in zip(pred, rows)], rows))

def basic_gate(pred, rows):
    return rmse(pred, rows) <= 20

def family_rmse(pred, rows):
    fams = sorted(set((d['source_family'] for d in rows)))
    return math.sqrt(sum((sum(((p - d['measured_mean_um']) ** 2 for (p, d) in zip(pred, rows) if d['source_family'] == s)) / sum((d['source_family'] == s for d in rows)) for s in fams)) / len(fams))

def verify():
    for tag in ['R1', 'R1_LINEAGE']:
        assert sha(P / f'PREREG_{tag}.json') == (P / f'PREREG_{tag}.sha256').read_text().strip()
    f = P / 'FROZEN_PREDICTIONS_R1.json'
    assert sha(f) == f.with_suffix('.sha256').read_text().strip()
    fr = json.loads(f.read_text())
    assert sha(P / 'measurements.json') == fr['source_measurements_sha256']
    for s in json.loads((P / 'SOURCE_MANIFEST.json').read_text())['sources']:
        assert sha(s['path']) == s['sha256']
    return fr

def dose_support(rows, region):
    q = [d for d in rows if d['region'] == region]
    contexts = {}
    for d in q:
        key = tuple((d[k] for k in ['source_family', 'cad', 'machine', 'material', 'method', 'state', 'restoration']))
        contexts.setdefault(key, set()).add(d['internal_spacer_um'])
    varied = {k[0] for (k, v) in contexts.items() if len(v) >= 2}
    return dict(independent_families_with_multiple_settings=sorted(varied), count=len(varied), required=3, inverse_supported=len(varied) >= 3)

def run():
    started = time.perf_counter()
    fr = verify()
    rows = json.loads((P / 'measurements.json').read_text())
    lookup = {d['row_id']: d for d in rows}
    metrics = {}
    tests = []
    for region in REGIONS:
        p = [d for d in fr['predictions'] if d['region'] == region]
        q = [lookup[d['row_id']] for d in p]
        mu = [d['prediction_um'] for d in p]
        pred_rmse = rmse(mu, q)
        gls = rmse([d['gls_um'] for d in p], q)
        ctrl = {k: rmse([d['controls_um'][k] for d in p], q) for k in ['mean', 'offset', 'affine']}
        coverage = balanced([d['lower_um'] <= t['measured_mean_um'] <= d['upper_um'] for (d, t) in zip(p, q)], q)
        eligible = [i for (i, d) in enumerate(p) if d['practice_um'] is not None]
        pe = [p[i] for i in eligible]
        qe = [q[i] for i in eligible]
        prac = rmse([d['practice_um'] for d in pe], qe)
        matched_pred = rmse([d['prediction_um'] for d in pe], qe)
        practice_ratio = matched_pred / prac if prac else None
        maxglsdiff = max((abs(d['gls_um'] - d['prediction_um']) for d in p))
        passed = pred_rmse <= 20 and practice_ratio is not None and (practice_ratio <= 0.8) and (coverage >= 0.9) and (len(set((d['source_family'] for d in q))) >= 3)
        metrics[region] = dict(rows=len(p), studies=len(set((d['study'] for d in q))), independent_families=len(set((d['source_family'] for d in q))), candidate_rmse_um=pred_rmse, gls_rmse_um=gls, simple_controls_rmse_um=ctrl, internal_identity_rmse_um=rmse([d['internal_identity_um'] for d in p], q), practice_rmse_um=prac, practice_candidate_same_support_rmse_um=matched_pred, practice_ratio=practice_ratio, practice_support_rows=len(eligible), interval_study_balanced_coverage=coverage, maximum_gls_prediction_difference_um=maxglsdiff, equivalence_pass=maxglsdiff <= 0.001, prediction_gate='PASS' if passed else 'FAIL', dose_support=dose_support(rows, region))
        metrics[region]['candidate_source_family_balanced_rmse_um'] = family_rmse(mu, q)
        for (method, v) in [('candidate', mu), ('GLS', [d['gls_um'] for d in p])] + [(k, [d['controls_um'][k] for d in p]) for k in ctrl]:
            bad = [x + 1000 for x in v]
            assert not basic_gate(bad, q)
            tests.append(dict(region=region, control=method, injection='+1000 um prediction to every held cell', mutated_rmse_um=rmse(bad, q), error_gate_rejects=True))
        for (name, values, truth) in [('CAD_explicit', [d['practice_um'] for d in pe], qe), ('internal_identity', [d['internal_identity_um'] for d in p], q)]:
            bad = [v + 1000 for v in values]
            assert not basic_gate(bad, truth)
            tests.append(dict(region=region, control=name, injection='+1000 um nominal prediction', mutated_rmse_um=rmse(bad, truth), error_gate_rejects=True))
        assert abs(p[0]['gls_um'] + 1000 - p[0]['prediction_um']) > 0.001
    mutated = copy.deepcopy(rows)
    mutated[0]['measured_mean_um'] += 1000
    mutated_bytes = (json.dumps(mutated, indent=2, ensure_ascii=False) + '\n').encode()
    assert hashlib.sha256(mutated_bytes).hexdigest() != fr['source_measurements_sha256']
    tests.append(dict(control='source_integrity', injection='one extracted measured mean +1000 um', expected_sha256=fr['source_measurements_sha256'], tampered_sha256=hashlib.sha256(mutated_bytes).hexdigest(), rejected=True))
    tests.append(dict(control='GLS_equivalence', injection='one GLS prediction +1000 um', parity_gate_rejects=True))
    models = {region: fit([d for d in rows if d['region'] == region]) for region in REGIONS}
    dump(P / 'REGIONAL_MODELS_R1.json', dict(models=models, input_sha256=sha(P / 'measurements.json'), use='research prediction only; inverse unsupported unless dose gate passes'))
    result = dict(round='R1', review_state='PENDING_INDEPENDENT_REVIEW', prediction_outcome='PASS' if all((d['prediction_gate'] == 'PASS' for d in metrics.values())) else 'FAIL', inverse_outcome='SUPPORTED' if all((d['dose_support']['inverse_supported'] for d in metrics.values())) else 'UNKNOWN', method_outcome='CLASSICAL_GLS_TIE', metrics=metrics, corruption_tests=tests, scoring_seconds=time.perf_counter() - started, prediction_sha256=sha(P / 'FROZEN_PREDICTIONS_R1.json'), external_referent={'kind': 'independent_measurement', 'locator': str(P / 'SOURCE_MANIFEST.json'), 'compared_quantity': 'Held-source-family regional reported group-mean gap in um', 'refutes_us': True}, limits=['Retrospective extracted group means, not blinded physical measurements', 'No calibration at an unseen machine batch', 'No specimen-level probabilities', 'Pulpl/internal-overall/absolute discrepancy are not pooled into occlusal or marginal gap', 'Predictive intervals rely on unverified Gaussian/linear closures'])
    dump(P / 'RESULTS_R1.json', result)
    dump(P / 'CURRENT_WORK_STATE.json', dict(lane='X13-cement-gap', phase='R1_SCORED', latest_gate=result['prediction_outcome'], next_operation='change information: measure a regional anchor in held study, predict other arms without seeing their outcomes', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    for (k, v) in metrics.items():
        print(k, v['candidate_rmse_um'], v['practice_ratio'], v['prediction_gate'], v['dose_support']['count'])
if __name__ == '__main__':
    run()
