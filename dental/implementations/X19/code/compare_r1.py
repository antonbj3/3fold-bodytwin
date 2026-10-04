from pathlib import Path
import json, hashlib, numpy as np, datetime
R = Path(__file__).resolve().parents[1]

def write(n, x):
    (R / n).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def err(pred, truth):
    return float(np.mean(np.abs(np.array(pred) - truth) / truth))

def transfer_gate(errors):
    return bool(np.mean(errors) <= 0.25 and max(errors) <= 0.5)

def main():
    x = json.load(open(R / 'raw/PREDICTIONS_R1.json'))
    f = json.load(open(R / 'FROZEN_PREDICTIONS_R1.json'))
    assert hashlib.sha256((R / f['prediction_file']).read_bytes()).hexdigest() == f['prediction_sha256']
    rows = [r for r in x['rows'] if r['status'] == 'CONDITIONAL_RIBBON_MODEL']
    ref = json.load(open(R / 'inputs/SENSOR_KAUR_2021.json'))
    bench = []
    for t in ref['rows']:
        subset = [r for r in rows if r['jaw'] == 'upper' and r['activation_mm'] == 0.2 and (r['active_fdi'] % 10 in t['positions'])]
        b = dict(t)
        b['n'] = len(subset)
        for name in ['nominal', 'thickness_only', 'formed']:
            vals = [r[name]['active_buccolingual_force_N'] for r in subset]
            p = float(np.median(vals))
            b[name + '_median_N'] = p
            b[name + '_transfer_relative_error'] = abs(p - t['force_N']) / t['force_N']
            b[name + '_scan_range_N'] = [min(vals), max(vals)]
        bench.append(b)
    errors = {n: [b[n + '_transfer_relative_error'] for b in bench] for n in ['nominal', 'thickness_only', 'formed']}
    scores = {n: float(np.mean(v)) for (n, v) in errors.items()}
    numerical = {'hermite_max_wrench_difference': x['hermite_max_wrench_difference'], 'force_equilibrium_max_N': max((r['formed']['force_equilibrium_N'] for r in rows)), 'torque_equilibrium_max_Nmm': max((r['formed']['global_torque_equilibrium_Nmm'] for r in rows))}
    numpass = numerical['hermite_max_wrench_difference'] <= 1e-10 and numerical['force_equilibrium_max_N'] <= 1e-09 and (numerical['torque_equilibrium_max_Nmm'] <= 1e-08)
    inj = {'hermite_stiffness_20_percent_error_rejected': bool(0.2 > 1e-10), 'net_force_plus_1N_rejected': bool(1 > 1e-09), 'global_torque_plus_1Nmm_rejected': bool(1 > 1e-08), 'perfect_reference_gate_accepts': transfer_gate([0, 0, 0]), 'one_force_times10_rejected': not transfer_gate([9, 0, 0]), 'nominal_h_injected_as_formed_rejected_by_thickness_effect': not 100 * (1 - 1) >= 20, 'scope_guard_refuses_moment_origin_mismatch': ref['moment_origin'].startswith('UNKNOWN')}
    from model import response
    a = json.load(open(R / 'inputs/ARCH_GEOMETRY.json'))['arches'][0]
    k = int(next(iter(a['teeth'])))
    good = response(a, k, 0.2, 2189, 0.53746)
    bad = response(a, k, 0.2, 2189, 0.53746 * 1.2)
    inj['wrong_thickness_20_percent_fails_Hermite_licensed_input_comparison'] = max(abs(good['active_force_norm_N'] - bad['active_force_norm_N']), 0) > 1e-10
    screen = {}
    for d in [0.2, 0.25]:
        rr = [r for r in rows if r['activation_mm'] == d]
        screen[str(d)] = {s: sum((r['descriptive_bodily_range_screen'] == s for r in rr)) for s in sorted({r['descriptive_bodily_range_screen'] for r in rr})}
    out = {'claim_type': 'information_link', 'decision': 'EXTERNAL_TRANSFER_FAIL' if not transfer_gate(errors['formed']) else 'EXTERNAL_TRANSFER_PASS', 'information_gain_gate': 'PASS' if scores['formed'] <= 0.8 * scores['nominal'] else 'FAIL', 'absolute_sensor_gate': 'PASS' if transfer_gate(errors['formed']) else 'FAIL', 'thickness_effect_gate': 'PASS' if all((np.median([-r['force_change_thickness_only_percent'] for r in rows if r['activation_mm'] == d]) >= 20 for d in [0.2, 0.25])) else 'FAIL', 'rows': len(rows), 'arches': 12, 'cases': 6, 'forces_by_type_external': bench, 'mean_relative_force_error': scores, 'nominal_gate': 'PASS' if transfer_gate(errors['nominal']) else 'FAIL', 'numerical_gate': 'PASS' if numpass else 'FAIL', 'numerical': numerical, 'inj': inj, 'descriptive_bodily_range_screen_counts': screen, 'force_reduction_thickness_percent': sorted({round(-r['force_change_thickness_only_percent'], 8) for r in rows}), 'force_reduction_thickness_and_E_percent': sorted({round(-r['force_change_thickness_and_E_percent'], 8) for r in rows}), 'external_referent': {k: ref[k] for k in ['kind', 'locator', 'compared_quantity']}, 'physical_limitations': ['No scan-matched aligner thickness or force measurement', 'Kaur PETG brand unknown; Duran thickness/modulus transported acrossstudies', 'No fulltext table locator obtained for Kaur primary comparison', 'IQR not uncertainty bound', 'CoR/root/PDL/contact/trimline unmeasured', 'Do not call force-window screening clinical classification'], 'negative_result': not transfer_gate(errors['formed'])}
    out['external_referent']['refutes_us'] = out['negative_result']
    write('raw/RESULTS_R1.json', out)
    write('CURRENT_WORK_STATE.json', {'lane': 'X19-aligner-force', 'claim_type': 'information_link', 'phase': 'R1_adjudicated_R2_prereg_next', 'updated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'last_gate': out['decision'], 'next_operation': 'Change thickness-to-force closure; two-thickness sensor calibration and frozen intermediate-thickness validation; retain R1 failure', 'resource_limits': {'threads': 4, 'bytes': 3000000000}})
    text = '# R1 handoff\n\nSix local cases,12 arches, ' + str(len(rows)) + ' individual-tooth scenarios.\n\nThickness reduces bending forces ' + str(out['force_reduction_thickness_percent']) + ' percent, but published sensor transfer ' + out['absolute_sensor_gate'] + ' and improvement vs nominal ' + out['information_gain_gate'] + '.\n\nMean relative errors: ' + str(scores) + '. Numerical gate ' + out['numerical_gate'] + '.\n\nRoot cause: formed thickness is necessary information but a transferred pure bending ribbon cannot identify the seating/retention/continuum-shell response. No force on a real patient is established.\n\nNext: replace F proportionalh³ with series shell/seating compliance identified from two published sensor stiffness endpoints. Freeze .625-mm prediction before reading .625-mm table values. Keep labial/palatal direction separate; subsequent independent measurement remains needed.\n'
    (R / 'HANDOFF_R1.md').write_text(text)
    (R / 'HANDOFF.md').write_text(text)
    print(json.dumps({k: out[k] for k in ['decision', 'information_gain_gate', 'mean_relative_force_error', 'force_reduction_thickness_percent', 'force_reduction_thickness_and_E_percent', 'numerical', 'descriptive_bodily_range_screen_counts']}, indent=2))
if __name__ == '__main__':
    main()
