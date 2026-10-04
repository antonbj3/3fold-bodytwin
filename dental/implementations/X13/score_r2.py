from pathlib import Path
import json, datetime, time, copy, hashlib
from score_r1 import rmse, balanced, basic_gate, sha, dump
from calibrate_gap import REGIONS
P = Path(__file__).resolve().parent

def run():
    started = time.perf_counter()
    f = P / 'FROZEN_PREDICTIONS_R2.json'
    assert sha(f) == f.with_suffix('.sha256').read_text().strip()
    fr = json.loads(f.read_text())
    assert sha(P / 'measurements.json') == fr['source_measurements_sha256']
    assert sha(P / 'PREREG_R2.json') == fr['prereg_sha256']
    assert sha(P / 'FROZEN_PREDICTIONS_R1.json') == fr['r1_prediction_sha256']
    rows = {d['row_id']: d for d in json.loads((P / 'measurements.json').read_text())}
    metrics = {}
    tests = []
    for region in REGIONS:
        p = [d for d in fr['predictions'] if d['region'] == region]
        q = [rows[d['row_id']] for d in p]
        er = rmse([d['prediction_um'] for d in p], q)
        base = rmse([d['R1_prediction_um'] for d in p], q)
        ratio = er / base
        cov = balanced([d['lower_um'] <= r['measured_mean_um'] <= d['upper_um'] for (d, r) in zip(p, q)], q)
        controls = {k: rmse([d[k] for d in p], q) for k in ['gaussian_control_um', 'constant_anchor_control_um', 'additive_control_um']}
        parity = max((abs(d['prediction_um'] - d['gaussian_control_um']) for d in p))
        passed = er <= 20 and ratio <= 0.8 and (cov >= 0.9) and (len(set((d['source_family'] for d in q))) >= 3)
        metrics[region] = dict(candidate_rmse_um=er, R1_same_support_rmse_um=base, ratio=ratio, interval_coverage=cov, control_rmse_um=controls, maximum_gaussian_control_difference_um=parity, evaluation_rows=len(p), studies=len(set((d['study'] for d in q))), source_families=len(set((d['source_family'] for d in q))), prediction_gate='PASS' if passed else 'FAIL')
        for name in ['prediction_um', 'gaussian_control_um', 'constant_anchor_control_um', 'additive_control_um']:
            bad = [d[name] + 1000 for d in p]
            assert not basic_gate(bad, q)
            tests.append(dict(region=region, control=name, injection='+1000 um to held prediction', mutated_rmse_um=rmse(bad, q), rejected=True))
    anchorgroups = {a['row_id'].rsplit(':', 1)[0] for a in fr['anchors']}
    families = set((d['source_family'] for d in fr['predictions']))
    result = dict(round='R2', outcome='PASS' if all((m['prediction_gate'] == 'PASS' for m in metrics.values())) else 'FAIL', metrics=metrics, inverse_outcome='UNKNOWN_UNCHANGED', method_outcome='CLASSICAL_GAUSSIAN_CONDITIONING_TIE', anchor_values=len(fr['anchors']), distinct_group_region_prefixes=len(anchorgroups), studies_with_evaluation=len(set((d['study'] for d in fr['predictions']))), calibration_cost='Published group means from n specimens; actual lab acquisition time/cost UNKNOWN', corruption_tests=tests, prediction_sha256=sha(f), score_seconds=time.perf_counter() - started, external_referent={'kind': 'independent_measurement', 'locator': str(P / 'SOURCE_MANIFEST.json'), 'compared_quantity': 'held other-arm regional group-mean gap after deterministic regional anchor', 'refutes_us': True}, review_state='PENDING_INDEPENDENT_REVIEW')
    dump(P / 'RESULTS_R2.json', result)
    dump(P / 'CURRENT_WORK_STATE.json', dict(lane='X13-cement-gap', phase='R2_SCORED', latest_gate=result['outcome'], next_operation='separate spacer-dose identification from cross-arm bias; paired local dose calibration only where external doses exist', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    for (region, m) in metrics.items():
        print(region, m['candidate_rmse_um'], m['ratio'], m['prediction_gate'])
if __name__ == '__main__':
    run()
