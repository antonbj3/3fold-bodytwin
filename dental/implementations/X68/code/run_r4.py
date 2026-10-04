"""Freeze a finite available-bore choice before comparing its published outcome."""
import json, pathlib, hashlib, datetime, math, time
from fractions import Fraction as F
import numpy as np
from enclosures import outward
ROOT = pathlib.Path(__file__).resolve().parents[1]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
read = lambda n: json.loads((ROOT / n).read_text())

def dump(n, o):
    (ROOT / n).write_text(json.dumps(o, indent=2, allow_nan=False) + '\n')

def choose(menu, field):
    return min(menu, key=lambda v: (abs(math.log(v[field] / 25.0)), -v['Df_mm']))

def main():
    start = time.perf_counter()
    pr = ROOT / 'PREREG_R4.json'
    assert sha(pr) == (ROOT / 'PREREG_R4.sha256').read_text().strip()
    rows = {r['id']: r for r in read('raw/arms.json')}
    r2 = read('FROZEN_PREDICTIONS_R2.json')
    r1 = read('results_R1.json')
    predictions = []
    for c in r2['calibrations']:
        endpoints = [rows[i] for i in c['endpoint_ids']]
        (a, b) = endpoints
        D = a['D_mm']
        slope = (b['torque_Ncm'] - a['torque_Ncm']) / (b['Df_mm'] - a['Df_mm'])
        menu = []
        for i in c['group']:
            d = rows[i]['Df_mm']
            menu.append(dict(id=i, Df_mm=d, predicted_Ncm=c['a_Ncm'] + c['b_Ncm_mm2'] * (D * D - d * d), ordinary_Ncm=a['torque_Ncm'] + slope * (d - a['Df_mm'])))
        selected = choose(menu, 'predicted_Ncm')
        ordinary = choose(menu, 'ordinary_Ncm')
        old = next((v for v in r1['inverse_choices'] if v['split'] == 'study' and v['kind'] == 'practice' and (v['group'][-1] == a['cortex_mm'])))
        predictions.append(dict(endpoint_ids=c['endpoint_ids'], menu=menu, selected=selected, ordinary_selected=ordinary, prior_practice_selected_id=old['selected_id'], rigorous_physical_enclosure='MISSING', guaranteed_fixed_Df_mm=None))
    f = ROOT / 'FROZEN_PREDICTIONS_R4.json'
    freeze = dict(created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(pr), data_sha256=sha(ROOT / 'raw/arms.json'), code_sha256=sha(pathlib.Path(__file__)), predictions=predictions, scope='Retrospective finite-menu query; selected interior torques excluded from selection. All source means previously exposed.')
    if not f.exists():
        dump(f.name, freeze)
    else:
        assert read(f.name)['predictions'] == predictions, 'Preserve freeze and open new round'
    dump('FROZEN_PREDICTIONS_R4.sha256.json', dict(sha256=sha(f), frozen_at_utc=read(f.name)['created_at_utc']))
    for v in predictions:
        for field in ['selected', 'ordinary_selected']:
            v[field]['observed_Ncm'] = rows[v[field]['id']]['torque_Ncm']
            v[field]['within_band'] = 20 <= v[field]['observed_Ncm'] <= 30
        v['practice_observed_Ncm'] = rows[v['prior_practice_selected_id']]['torque_Ncm']
        v['practice_within_band'] = 20 <= v['practice_observed_Ncm'] <= 30
    (d0, d1, T0, T1) = map(lambda x: F(float(x)), [2.7, 3.6, 39.2, 16.2])
    h = d1 - d0
    s = (T0 - T1) / (d1 * d1 - d0 * d0)

    def A(d):
        q = F(float(d))
        w = (q * q - d0 * d0) / (d1 * d1 - d0 * d0)
        return float((1 - w) * T0 + w * T1)

    def B(d):
        q = F(float(d))
        t = (q - d0) / h
        w = (q * q - d0 * d0) / (d1 * d1 - d0 * d0)
        return float((1 - w) * T0 + w * T1 - 20 * t * (1 - t))
    ds = [2.7, 3.0, 3.3, 3.6]
    sa = np.array([A(ds[0]), A(ds[-1])])
    sb = np.array([B(ds[0]), B(ds[-1])])
    am = [dict(Df_mm=d, predicted_Ncm=A(d)) for d in ds]
    bm = [dict(Df_mm=d, predicted_Ncm=B(d)) for d in ds]
    da = choose(am, 'predicted_Ncm')['Df_mm']
    db = choose(bm, 'predicted_Ncm')['Df_mm']
    ga = [-2 * s * d for d in [d0, d1]]
    gb = [-2 * s * d - 20 / h * (1 - 2 * (d - d0) / h) for d in [d0, d1]]
    suff = dict(summary='Two endpoint torque means and same finite nominal-bore menu', bitwise_identical=sa.tobytes() == sb.tobytes(), identity_error_Ncm=float(np.max(abs(sa - sb))), endpoint_means_Ncm=sa.tolist(), selected_Df_mm=[da, db], downstream_choice_difference_mm=abs(da - db), curve_A_menu=am, curve_B_menu=bm, rigorous_derivative_bounds_Ncm_per_mm={'A': [outward(min(ga), True), outward(max(ga), False)], 'B': [outward(min(gb), True), outward(max(gb), False)]}, exact_rational_derivative_bounds={'A': [str(min(ga)), str(max(ga))], 'B': [str(min(gb)), str(max(gb))]}, monotonicity_proof='Affine derivatives attain extrema at endpoints; exact rational maxima <0. Both responses positive by endpoint positivity and monotonicity.', smallest_extension='One measured interior response separates these two quadratic curves; arbitrary physical response still requires an independently justified error enclosure.', external_referent=dict(kind='our_own_fixture', locator='code/run_r4.py exact-rational endpoint-preserving construction; not physical truth', compared_quantity='Finite-menu choice ambiguity for identical calibration summary', refutes_us=True))
    hits = sum((v['selected']['within_band'] for v in predictions))
    ordinary = sum((v['ordinary_selected']['within_band'] for v in predictions))
    passed = len(predictions) >= 2 and hits / len(predictions) >= 0.8

    def band(x):
        return 20 <= x <= 30
    faults = dict(corrupted_selected_measurement=all((not band(10 * v['selected']['observed_Ncm']) for v in predictions)), corrupted_ordinary_measurement=all((not band(10 * v['ordinary_selected']['observed_Ncm']) for v in predictions)), corrupted_endpoint_identity=not np.array_equal(sa, sb + 1.0), reversed_derivative_certificate=not max([-g for g in ga]) <= 0)
    result = dict(round='R4', claim_type='capability', verdict='PASS_RETROSPECTIVE_DISCRETE_QUERY' if passed else 'FAIL', selected_hits=hits, evaluable_groups=len(predictions), hit_fraction=hits / len(predictions), ordinary_control_hits=ordinary, prior_practice_hits=sum((v['practice_within_band'] for v in predictions)), independent_study_families=1, general_transfer_verdict='UNKNOWN_COVERAGE', physical_fixed_bore_guarantee=None, radial_pressure_MPa=None, predictions=predictions, summary_sufficiency=suff, external_referent=read(pr.name)['external_referent'], injected_faults_rejected=faults, prereg_sha256=sha(pr), frozen_predictions_sha256=sha(f), cost_seconds=time.perf_counter() - start)
    dump('results_R4.json', result)
    assert all(faults.values())
    print(json.dumps({k: result[k] for k in ['verdict', 'selected_hits', 'evaluable_groups', 'ordinary_control_hits', 'prior_practice_hits', 'general_transfer_verdict', 'summary_sufficiency']}, indent=2))
if __name__ == '__main__':
    main()
