"""New controls preserve the original weak source-mutation outcome."""
import datetime as dt, hashlib, json
from fractions import Fraction as F
from pathlib import Path
import ports
from sources import load
from make_demo import write
ROOT = Path(__file__).resolve().parent

def main():
    reg = ROOT / 'PREREG_CONTROL_REPAIR.json'
    h = hashlib.sha256(reg.read_bytes()).hexdigest()
    receipt = ROOT / 'PREREG_CONTROL_REPAIR_FREEZE.json'
    if not receipt.exists():
        write(receipt.name, {'frozen_utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'prereg_sha256': h})
    if json.loads(receipt.read_text())['prereg_sha256'] != h:
        raise ValueError('control prereg drift')
    data = load()
    controls = []
    for (g, rows) in data['ISQ'].items():
        for p in rows:
            bad = p['mean_ISQ'] + 10
            controls.append({'name': g + '_day' + str(p['day']) + '_returned_ISQ+10', 'detected': abs(bad - p['mean_ISQ']) > 1})
    for (g, rows) in data['WEAR'].items():
        for p in rows:
            tol = max(0.1, 1.96 * p['SEM_mg'])
            bad = p['mean_mg'] + 10
            controls.append({'name': g + '_N' + str(p['cycles']) + '_returned_mass+10mg', 'detected': abs(bad - p['mean_mg']) > tol})
    for (grade, states) in data['LTD'].items():
        v = states['aged']
        correct = ports.weibull_quantile(v['m'], v['sigma0_MPa'], 0.05)
        bad = ports.weibull_quantile(states['unaged']['m'], v['sigma0_MPa'], 0.05)
        controls.append({'name': grade + '_aged_shape_wrongly_unaged', 'detected': abs(bad - correct) / correct > 1e-12})
    calls = [('reverse_temporal_knots', lambda : ports.interpolate_enclosed([90, 28, 14, 0], [70, 68, 66, 71], 21)), ('reverse_wear_inventory', lambda : ports.wear_envelope([5, 20, 120], [2, 1, 0.5], 40)), ('clinical_ISQ_to_BIC', lambda : ports.promote_isq('BIC')), ('mass_to_volume_without_rho', lambda : ports.mass_to_volume(1))]
    for (name, fn) in calls:
        try:
            fn()
            caught = False
        except ports.PortError:
            caught = True
        controls.append({'name': name, 'detected': caught})
    enclosure_checks = []
    for (g, rows) in data['WEAR'].items():
        train = rows[:3]
        xs = [p['cycles'] for p in train]
        ys = [F(str(p['source_mean_1e-4g'])) / 10 for p in train]
        for intercept in [False, True]:
            ak = ports.fit_rate(xs, ys, intercept)
            for n in [40000, 80000, 120000]:
                exact = ak[0] + n * ak[1]
                bounds = ports.rate_query(ak, n)['arithmetic_enclosure_mg']
                enclosure_checks.append(F.from_float(bounds[0]) <= exact <= F.from_float(bounds[1]))
    write('raw/CONTROL_REPAIR.json', {'controls': controls, 'count': len(controls), 'detected_count': sum((c['detected'] for c in controls)), 'all_detected': all((c['detected'] for c in controls)), 'exact_rational_enclosure_checks': len(enclosure_checks), 'all_enclosed': all(enclosure_checks), 'preserved_original_weak_mutation': 'raw/RESULT_K34_WEAR_R1.json / Grandio120000;1/15 source+10mg mutations not detected', 'scope': 'Specified faults only. Physical validity, forged new source data and statistical coverage are not certified.'})
    print('Controls repaired', len(controls), 'detected', sum((c['detected'] for c in controls)))
if __name__ == '__main__':
    main()
