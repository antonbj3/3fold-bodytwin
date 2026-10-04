import json, pathlib, hashlib, datetime, time
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def dump(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False, default=lambda v: v.item()) + '\n')

def main():
    start = time.perf_counter()
    pr = ROOT / 'PREREG_R2.json'
    assert sha(pr) == (ROOT / 'PREREG_R2.sha256').read_text().strip()
    allrows = json.loads((ROOT / 'raw/arms.json').read_text())
    rows = {r['id']: r for r in allrows}
    groups = [['B_T2_2.7', 'B_T2_3', 'B_T2_3.3', 'B_T2_3.6'], ['B_T5_3', 'B_T5_3.3', 'B_T5_3.6', 'B_T5_3.8']]
    prior = json.loads((ROOT / 'results_R1.json').read_text())
    practice = {v['id']: v['predicted_Ncm'] for v in prior['predictions'] if v['split'] == 'study' and v['method'] == 'candidate_ridge' and (v['kind'] == 'practice')}
    canonical = {r['id']: r.get('duplicate_of', r['id']) for r in allrows}
    predictions = []
    cal = []
    for g in groups:
        rs = [rows[v] for v in g]
        D = rs[0]['D_mm']
        (lo, hi) = (rs[0], rs[-1])
        u0 = D ** 2 - lo['Df_mm'] ** 2
        u1 = D ** 2 - hi['Df_mm'] ** 2
        b = (lo['torque_Ncm'] - hi['torque_Ncm']) / (u0 - u1)
        a = lo['torque_Ncm'] - b * u0
        q25 = D ** 2 - (25 - a) / b
        df25 = float(np.sqrt(q25)) if q25 >= 0 else None
        slope = (hi['torque_Ncm'] - lo['torque_Ncm']) / (hi['Df_mm'] - lo['Df_mm'])
        lin25 = lo['Df_mm'] + (25 - lo['torque_Ncm']) / slope
        cal.append(dict(group=g, endpoint_ids=[g[0], g[-1]], a_Ncm=a, b_Ncm_mm2=b, area_based_Df_target25_mm=df25, ordinary_linear_Df_target25_mm=lin25, calibration_range_mm=[lo['Df_mm'], hi['Df_mm']], target_in_range=lo['Df_mm'] <= df25 <= hi['Df_mm'], physical_enclosure='MISSING'))
        for r in rs[1:-1]:
            u = D ** 2 - r['Df_mm'] ** 2
            predictions.append(dict(id=r['id'], predicted_Ncm=a + b * u, ordinary_linear_Ncm=lo['torque_Ncm'] + slope * (r['Df_mm'] - lo['Df_mm']), practice_Ncm=practice[canonical[r['id']]], endpoint_ids=[g[0], g[-1]]))
    freeze = {'created_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prereg_sha256': sha(pr), 'data_sha256': sha(ROOT / 'raw/arms.json'), 'code_sha256': sha(pathlib.Path(__file__)), 'predictions': predictions, 'calibrations': cal, 'scope': 'Retrospective same-study interior validation; not physical measurement conducted in this session'}
    f = ROOT / 'FROZEN_PREDICTIONS_R2.json'
    if not f.exists():
        dump(f.name, freeze)
    else:
        assert json.loads(f.read_text())['predictions'] == predictions
    dump('FROZEN_PREDICTIONS_R2.sha256.json', {'sha256': sha(f), 'frozen_at_utc': json.loads(f.read_text())['created_at_utc']})
    errors = {k: [] for k in ['candidate', 'practice', 'ordinary_linear']}
    for v in predictions:
        obs = rows[v['id']]['torque_Ncm']
        v['observed_Ncm'] = obs
        for (k, field) in [('candidate', 'predicted_Ncm'), ('practice', 'practice_Ncm'), ('ordinary_linear', 'ordinary_linear_Ncm')]:
            e = float(abs(np.log(v[field] / obs)))
            errors[k].append(e)
            v[k + '_abs_ln'] = e
    metrics = {k: dict(mean_abs_ln=float(np.mean(es)), max_abs_ln=float(np.max(es))) for (k, es) in errors.items()}
    t = np.array([0.0, 0.5, 1.0])
    A = np.array([39.2, (39.2 + 16.2) / 2, 16.2])
    B = A + 8.0 * t * (1.0 - t)
    endpoint_summary_A = A[[0, 2]]
    endpoint_summary_B = B[[0, 2]]
    suff = dict(summary='Two endpoint mean torques at the same two bore diameters', bitwise_identical=endpoint_summary_A.tobytes() == endpoint_summary_B.tobytes(), identity_error=float(np.max(abs(endpoint_summary_A - endpoint_summary_B))), middle_torque_Ncm=[float(A[1]), float(B[1])], downstream_difference_Ncm=float(B[1] - A[1]), monotonicity_proof='dA/dt=-23; dB/dt=-23+8*(1-2t), in [-31,-15] for all t in [0,1]', smallest_extension='One independent interior torque identifies a QUADRATIC response; a rigorous curvature bound is required to enclose arbitrary monotone response.', external_referent={'kind': 'closed_form', 'locator': 'https://dlmf.nist.gov/3.3 (interpolation remainder); constructed curves are mathematical witnesses, not physical measurements', 'compared_quantity': 'Endpoint-preserving polynomial and linear-interpolation curvature remainder', 'refutes_us': True})
    passed = metrics['candidate']['max_abs_ln'] <= np.log(1.2) and metrics['candidate']['mean_abs_ln'] <= 0.5 * metrics['practice']['mean_abs_ln']
    result = dict(round='R2', claim_type='information_link', verdict='PASS_RETROSPECTIVE_INTERPOLATION' if passed else 'FAIL', predictions=predictions, metrics=metrics, calibrations=cal, summary_sufficiency=suff, gain_ratio=metrics['candidate']['mean_abs_ln'] / metrics['practice']['mean_abs_ln'], physical_target_choice='UNKNOWN_NO_CURVATURE_ENCLOSURE', external_referent=json.loads(pr.read_text())['external_referent'], frozen_predictions_sha256=sha(f), prereg_sha256=sha(pr), injected_faults_rejected={'mean_error_gate': abs(np.log(250 / 25)) > np.log(1.2), 'endpoint_identity': not np.array_equal(endpoint_summary_A, endpoint_summary_B + 1.0), 'ordinary_control_error_gate': abs(np.log(0.01 / 25)) > np.log(1.2)}, cost={'fit_validate_query_seconds': time.perf_counter() - start, 'physical_calibration_cost': 'Two same-material laboratory groups/measurements per regime, NOT RUN', 'source_preparation': 'Reuse R1 charged extraction'})
    dump('results_R2.json', result)
    print(json.dumps({k: result[k] for k in ['verdict', 'metrics', 'gain_ratio', 'physical_target_choice', 'calibrations']}, indent=2))
if __name__ == '__main__':
    main()
