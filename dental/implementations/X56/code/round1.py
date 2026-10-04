"""LOSO on source-compatible observation operators; numerical model ≠ attribution."""
from dental_release.paths import expand as _release_expand
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '4'
from pathlib import Path
import json, math, time, hashlib
import numpy as np
from motion_port import local_threshold
LANE = Path(__file__).resolve().parents[1]

def features(d):
    return np.array([1, math.log(d['rho_g_cm3'] / 0.32), d['cortex_mm'], math.sin(math.radians(d['angle_deg'])) ** 2, d['implant_tilt_deg'] / 30])

def fit(train):
    x = np.array([features(d) for d in train])
    y = np.log([d['motion_um'] / d['force_N'] for d in train])
    w = np.array([1 / sum((q['study'] == d['study'] for q in train)) for d in train])
    a = x.T @ (w[:, None] * x) + np.diag([0, 1, 1, 1, 1])
    return (np.linalg.solve(a, x.T @ (w * y)), float(np.average(y, weights=w)))

def main():
    start = time.monotonic()
    data = json.loads((LANE / 'raw/measurements.json').read_text())
    folds = []
    for sid in sorted({d['study'] for d in data}):
        test = [d for d in data if d['study'] == sid]
        train = [d for d in data if d['study'] != sid and d['observation_operator'] == 'abutment_crosshead']
        if test[0]['observation_operator'] != 'abutment_crosshead' or not train:
            folds.append({'study': sid, 'verdict': 'UNKNOWN_UNSEEN_OBSERVATION_OPERATOR', 'n_test': len(test), 'error_ln': None, 'predictions': []})
            continue
        (beta, base) = fit(train)
        pr = []
        for d in test:
            y = d['force_N'] * math.exp(float(features(d) @ beta))
            b = d['force_N'] * math.exp(base)
            pr.append({'id': d['id'], 'observed_um': d['motion_um'], 'predicted_um': y, 'same_information_intercept_um': b, 'abs_log_error': abs(math.log(y / d['motion_um'])), 'check_abs_log_error': abs(math.log(b / d['motion_um'])), 'resolution': 'POPULATION'})
        folds.append({'study': sid, 'verdict': 'EXPLORATORY_ONLY_FEWER_THAN_3_OPERATOR_COMPATIBLE_STUDIES', 'n_test': len(test), 'n_train_studies': len({d['study'] for d in train}), 'error_ln': float(np.median([p['abs_log_error'] for p in pr])), 'check_error_ln': float(np.median([p['check_abs_log_error'] for p in pr])), 'coefficients': beta.tolist(), 'predictions': pr})
    compatible = [p for f in folds for p in f['predictions']]
    rawbase = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/datasets/dental_3fold_extra/scratch_MICRO/runs'))
    runs = {}
    for f in [] if os.environ.get('X56_PORTABLE') == '1' else sorted(rawbase.glob('*.jsonl')):
        for line in f.read_text().splitlines():
            d = json.loads(line)
            if d.get('ok'):
                runs[d['id']] = {**d, 'raw_path': str(f), 'raw_sha256': hashlib.sha256(f.read_bytes()).hexdigest()}
    if not runs:
        runs = {d['id']: d for d in json.loads((LANE / 'raw/historical_runs.json').read_text())}
    perturb = []
    for (original, modified, obs) in [('HO1_00', 'HO1_00_E66posthoc', 358.7), ('HO1_45', 'HO1_45_E66posthoc', 550)]:
        (b, m) = (runs[original], runs[modified])
        perturb.append({'base_id': original, 'modified_id': modified, 'measured_system_um': obs, 'baseline_predicted_system_um': b['u_load_along_F_um'], 'modified_predicted_system_um': m['u_load_along_F_um'], 'baseline_predicted_slip_um': b['slip_max_um'], 'modified_predicted_slip_um': m['slip_max_um'], 'system_gain': m['u_load_along_F_um'] / b['u_load_along_F_um'], 'slip_gain': m['slip_max_um'] / b['slip_max_um'], 'remaining_measured_over_predicted': obs / m['u_load_along_F_um'], 'closed_log_gap_fraction': math.log(m['u_load_along_F_um'] / b['u_load_along_F_um']) / math.log(obs / b['u_load_along_F_um']), 'changed_material_moduli_MPa': {'foam': [b['spec']['E_s'], m['spec']['E_s']], 'cortex': [b['spec']['E_c'], m['spec']['E_c']]}, 'causal_attribution': 'NOT_IDENTIFIED: BOTH foam and cortex E change, also penalty and recalibrated pressure; 66 MPa is specimen-system pretest fit', 'resolution': 'PER_POINT (FE system); PER_SURFACE_REGION (FE peak slip); POPULATION (measured system)', 'raw_paths': [b['raw_path'], m['raw_path']]})
    controls = []
    for (var, key) in [('contact penalty', 'CV_K4'), ('cortex modulus', 'HO1_45_E78'), ('prestress thermal parameter', 'HO1_45_Tc0')]:
        (b, m) = (runs['HO1_45'], runs[key])
        controls.append({'parameter': var, 'run': key, 'system_gain': m['u_load_along_F_um'] / b['u_load_along_F_um'], 'slip_gain': m['slip_max_um'] / b['slip_max_um'], 'measured_cause': False, 'rigorous_enclosure': 'MISSING', 'warning': 'single actual historical simulated perturbation, NOT measurement or worst-case bound'})
    ratio = np.exp(np.mean([math.log(d['remaining_measured_over_predicted']) for d in perturb]))
    injected = [{**p, 'predicted_um': 1000 * p['predicted_um']} for p in compatible]
    negative = {'x1000_injection_median_ln_error': float(np.median([abs(math.log(p['predicted_um'] / p['observed_um'])) for p in injected])), 'x1000_fails_factor2': all((abs(math.log(p['predicted_um'] / p['observed_um'])) > math.log(2) for p in injected)), 'wrong_observable_interface_gate_rejects_all': all((local_threshold([d['motion_um'], d['motion_um']], d['interface_relative_direct']).startswith('UNKNOWN') for d in data))}
    out = {'round': 'R1', 'claim_type': 'information_link', 'verdict': 'TRANSFERABILITY_AND_INTERFACE_GATE_NOT_ESTABLISHED', 'external_referent': {'kind': 'independent_measurement', 'locator': 'PMC5577443#T2; PMC7150554#T1; PMC8705369#materials-14-07886-t001; PMC10792929#Tab1', 'compared_quantity': 'source-defined loaded system or optical-marker displacement (not paired local interface slip)', 'refutes_us': True}, 'counts': {'rows': len(data), 'primary_studies': 4, 'operator_compatible_studies': 2, 'direct_interface_studies': 0, 'LOSO_scored_rows': len(compatible)}, 'LOSO_folds': folds, 'LOSO_median_abs_log_error': float(np.median([p['abs_log_error'] for p in compatible])), 'LOSO_study_balanced_median_abs_log_error': float(np.median([f['error_ln'] for f in folds if f['error_ln'] is not None])), 'primary_gate': {'coverage_pass': False, 'required_compatible_studies': 3, 'threshold': math.log(2), 'error_pass_exploratory': float(np.median([p['abs_log_error'] for p in compatible])) <= math.log(2), 'current_practice_comparison': 'Historical K27 predictions exist only for matched 2017/2019 conditions; no fabricated prediction for new implants', 'direct_interface_threshold_decisions': 0}, 'historical_measured_material_link': perturb, 'remaining_geometric_system_underprediction_factor': float(ratio), 'other_actual_FE_perturbations': controls, 'causal_bisection': 'Magnitude gap partially reduced by source-linked PAIRED foam AND cortex moduli; friction, p0, geometry, crush and fixture are unresolved, not individually identified', 'negative_controls': negative, 'full_cost': {'runtime_s': time.monotonic() - start, 'training': '2 whole-study fits (one training study each)', 'source_extraction': '1414 candidate paths, local XML only', 'new_FE_runs': 0, 'past_FE_runtime_charged': 'raw t_total per input run in source snapshot; not a warm-query claim', 'historical_FE_runtime_s': sum((runs[k].get('t_total', 0) for k in ['HO1_00', 'HO1_45', 'HO1_00_E66posthoc', 'HO1_45_E66posthoc', 'CV_K4', 'HO1_45_E78', 'HO1_45_Tc0'])), 'historical_cost_status': 'reused past physical-model preparation/solve; full experimental acquisition cost UNKNOWN', 'sensitivity_rigorous_enclosure': 'MISSING; perturbation ratios are not linear uncertainty propagation'}, 'resolution': 'POPULATION for table means; distinct FE spatial output levels retained'}
    (LANE / 'rounds/R1').mkdir(parents=True, exist_ok=True)
    (LANE / 'rounds/R1/results.json').write_text(json.dumps(out, indent=2) + '\n')
    (LANE / 'raw/historical_runs.json').write_text(json.dumps([runs[k] for k in ['HO1_00', 'HO1_45', 'HO1_00_E66posthoc', 'HO1_45_E66posthoc', 'CV_K4', 'HO1_45_E78', 'HO1_45_Tc0']], indent=2) + '\n')
    print(json.dumps({k: out[k] for k in ['verdict', 'counts', 'LOSO_median_abs_log_error', 'remaining_geometric_system_underprediction_factor', 'negative_controls']}, indent=2))
if __name__ == '__main__':
    main()
