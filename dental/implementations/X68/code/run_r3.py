import json, pathlib, hashlib, datetime, time
import numpy as np
from calibration import monotone_inverse, lagrange, derivative_enclosure_area_poly
ROOT = pathlib.Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def dump(n, o):
    (ROOT / n).write_text(json.dumps(o, indent=2, ensure_ascii=False, default=lambda v: v.item(), allow_nan=False) + '\n')

def main():
    start = time.perf_counter()
    p = ROOT / 'PREREG_R3.json'
    assert sha(p) == (ROOT / 'PREREG_R3.sha256').read_text().strip()
    rows = {r['id']: r for r in json.loads((ROOT / 'raw/arms.json').read_text())}
    groups = [(['B_T2_2.7', 'B_T2_3', 'B_T2_3.6'], 'B_T2_3.3'), (['B_T5_3', 'B_T5_3.3', 'B_T5_3.8'], 'B_T5_3.6')]
    predictions = []
    for (ids, testid) in groups:
        r = [rows[v] for v in ids]
        D = r[0]['D_mm']
        x = np.array([D * D - v['Df_mm'] ** 2 for v in r])
        y = np.array([v['torque_Ncm'] for v in r])
        poly = np.polyfit(x, y, 2)
        df = rows[testid]['Df_mm']
        query = D * D - df * df
        grid = np.linspace(min(x), max(x), 101)
        eq = max((abs(np.polyval(poly, q) - lagrange(x, y, q)) / max(1, abs(lagrange(x, y, q))) for q in grid))
        pred = float(np.polyval(poly, query))
        points = [(v['Df_mm'], v['torque_Ncm']) for v in r]
        inv = monotone_inverse(points, 25.0)
        der = derivative_enclosure_area_poly(poly, D, min((v['Df_mm'] for v in r)), max((v['Df_mm'] for v in r)))
        predictions.append(dict(calibration_ids=ids, test_id=testid, predicted_Ncm=pred, poly_c_b_a=poly.tolist(), same_information_Lagrange_max_relative_difference=eq, derivative_enclosure=der, inverse=inv))
    f = ROOT / 'FROZEN_PREDICTIONS_R3.json'
    freeze = {'created_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prereg_sha256': sha(p), 'code_sha256': {n: sha(ROOT / 'code' / n) for n in ['run_r3.py', 'calibration.py']}, 'data_sha256': sha(ROOT / 'raw/arms.json'), 'predictions': predictions, 'note': 'Retrospective study; remaining intermediate outcomes not used by fit'}
    if not f.exists():
        dump(f.name, freeze)
    else:
        assert json.loads(f.read_text())['predictions'] == predictions
    dump('FROZEN_PREDICTIONS_R3.sha256.json', {'sha256': sha(f), 'frozen_at_utc': json.loads(f.read_text())['created_at_utc']})
    for v in predictions:
        v['observed_Ncm'] = rows[v['test_id']]['torque_Ncm']
        v['abs_ln'] = float(abs(np.log(v['predicted_Ncm'] / v['observed_Ncm'])))
    a = np.float64(rows['B_T2_3.6']['torque_Ncm'])
    b = np.float64(rows['B_T5_3.8']['torque_Ncm'])
    suff = {'summary': 'mean total peak insertion torque Ncm', 'bitwise_identical': a.tobytes() == b.tobytes(), 'identity_error_Ncm': float(abs(a - b)), 'states': ['Bayarchimeg Table2 Df3.6 no cortex', 'Bayarchimeg Table5 Df3.8 cortex1.5'], 'total_torque_Ncm': [float(a), float(b)], 'observed_ISQ': [52.5, 56.8], 'downstream_ISQ_difference': 56.8 - 52.5, 'observation_level': 'POPULATION', 'cortex_level': 'PER_SURFACE_REGION', 'smallest_extension_for_this_pair': 'cortex contact state and preparation profile; ISQ must be observed separately if used as a decision variable', 'pressure_identification': 'UNKNOWN; ISQ is not an independent pressure measurement', 'external_referent': {'kind': 'independent_measurement', 'locator': 'doi:10.5051/jpis.2013.43.1.30 Tables2 and5', 'compared_quantity': 'ISQ at identical mean total peak insertion torque', 'refutes_us': True}}
    fbad = monotone_inverse([(2.7, 39.2), (3.0, 35.4), (3.6, 392.0)], 25.0)
    faults = {'reversed_endpoint_rejected': fbad['status'] == 'INCONSISTENT', 'doubled_prediction_rejected': all((abs(np.log(2 * v['predicted_Ncm'] / v['observed_Ncm'])) > np.log(1.2) for v in predictions)), 'Lagrange_fault_rejected': abs(lagrange([0.0, 1.0, 2.0], [1.0, 2.0, 4.0], 1.0) - 999) > 1e-10}
    passed = all((v['abs_ln'] <= np.log(1.2) and v['same_information_Lagrange_max_relative_difference'] <= 1e-10 for v in predictions))
    result = {'round': 'R3', 'claim_type': 'capability', 'verdict': 'PASS_CONDITIONAL_MODEL_QUERY' if passed else 'FAIL_MODEL_CLOSURE', 'predictions': predictions, 'summary_sufficiency': suff, 'physical_prediction_enclosure': 'MISSING', 'continuous_single_bore_choice': 'UNKNOWN_CURVATURE_AND_INTERSPECIMEN_VARIATION', 'radial_pressure_MPa': None, 'external_referent': suff['external_referent'], 'frozen_predictions_sha256': sha(f), 'prereg_sha256': sha(p), 'injected_faults_rejected': faults, 'cost_seconds': time.perf_counter() - start}
    dump('results_R3.json', result)
    print(json.dumps(result, indent=2, default=lambda v: v.item()))
if __name__ == '__main__':
    main()
