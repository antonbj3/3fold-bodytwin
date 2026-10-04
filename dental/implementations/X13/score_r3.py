from pathlib import Path
import json, time, datetime
from score_r1 import rmse, sha, dump, basic_gate
P = Path(__file__).resolve().parent

def run():
    started = time.perf_counter()
    f = P / 'FROZEN_PREDICTIONS_R3.json'
    assert sha(f) == f.with_suffix('.sha256').read_text().strip()
    fr = json.loads(f.read_text())
    assert sha(P / 'measurements.json') == fr['source_measurements_sha256']
    assert sha(P / 'PREREG_R3.json') == fr['prereg_sha256']
    lookup = {d['row_id']: d for d in json.loads((P / 'measurements.json').read_text())}
    p = fr['predictions']
    q = [lookup[d['row_id']] for d in p]
    errors = [d['prediction_um'] - r['measured_mean_um'] for (d, r) in zip(p, q)]
    er = rmse([d['prediction_um'] for d in p], q)
    ctrl = {k: rmse([d[k] for d in p], q) for k in ['reciprocal_control_um', 'affine_control_um', 'constant_control_um']}
    ratio = er / ctrl['affine_control_um']
    maxer = max(map(abs, errors))
    passed = er <= 20 and maxer <= 20 and (ratio <= 0.8) and (len(set((d['study'] for d in q))) >= 2)
    parity = max((abs(d['prediction_um'] - d['reciprocal_control_um']) for d in p))
    assert parity <= 0.001
    tests = []
    for method in ['prediction_um', 'reciprocal_control_um', 'affine_control_um', 'constant_control_um']:
        bad = [d[method] + 1000 for d in p]
        assert not basic_gate(bad, q)
        tests.append(dict(control=method, injection='+1000 um to prediction', mutated_rmse_um=rmse(bad, q), rejected=True))
    cells = [dict(row_id=d['row_id'], spacer_um=d['spacer_um'], prediction_um=d['prediction_um'], observed_um=r['measured_mean_um'], error_um=e, affine_um=d['affine_control_um'], extrapolation=True) for (d, r, e) in zip(p, q, errors)]
    result = dict(round='R3', pilot_outcome='PASS' if passed else 'FAIL', full_regional_inverse_outcome='UNKNOWN', candidate_rmse_um=er, maximum_absolute_error_um=maxer, control_rmse_um=ctrl, ratio_to_affine=ratio, maximum_reciprocal_control_difference_um=parity, cells=cells, corruption_tests=tests, prediction_sha256=sha(f), score_seconds=time.perf_counter() - started, method_outcome='CLASSICAL_RECIPROCAL_INTERPOLATION_TIE', external_referent=json.loads((P / 'PREREG_R3.json').read_text())['external_referent'], review_state='PENDING_INDEPENDENT_REVIEW', limitations=['Only 2 studies and 3 held higher-dose cells', 'Retrospective, chosen after earlier negative construction', 'Same-information reciprocal control is exactly equivalent', 'No measured internal regional dose response', 'Anchor SD sensitivity is not a predictive confidence interval', 'Observed successful extrapolation does not authorize extrapolation in inverse queries'])
    dump(P / 'RESULTS_R3.json', result)
    dump(P / 'CURRENT_WORK_STATE.json', dict(lane='X13-cement-gap', phase='R3_SCORED', latest_gate=result['pilot_outcome'], full_goal='UNKNOWN', next_operation='couple observation semantics and hydraulic conductance; preserve local marginal pilot and specify same-machine regional perturbation experiment', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    print(er, maxer, ratio, result['pilot_outcome'])
if __name__ == '__main__':
    run()
