"""Verify frozen scientific decisions, units/origins and real corruptions."""
import copy, json, hashlib
from pathlib import Path
import numpy as np
from calibration_port import fit
from whole_arch import build
R = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    checks = {}
    records = []
    for i in range(1, 10):
        p = R / f'PREREG_R{i}.json'
        checks[f'prereg_R{i}_immutable'] = sha(p) == json.loads((R / f'PREREG_R{i}.sha256.json').read_text())['sha256']
    for record in json.loads((R / 'SOURCE_MANIFEST.json').read_text()):
        records.append(dict(path=record['path'], expected=record['sha256'], observed=sha(Path(record['path']))))
    checks['source_hashes_match'] = all((r['expected'] == r['observed'] for r in records))
    frozen = json.loads((R / 'FROZEN_PREDICTIONS_R1.json').read_text())
    checks['R1_prediction_immutable'] = sha(R / 'raw/PREDICTIONS_R1.json') == frozen['predictions_sha256']
    checks['R1_numeric_gates_pass'] = all((json.loads((R / 'round1/results.json').read_text())['gates'][k] for k in ['numerical_parity', 'nonpenetration', 'complementarity', 'stationarity', 'force_balance', 'torque_balance']))
    p = json.loads((R / 'raw/PREDICTIONS_R1.json').read_text())
    ge = json.loads((R / 'PREREG_R1.json').read_text())['gates']
    row = p['rows'][0]['regional']
    regional = row['regional']
    force = np.sum([r['force_N'] for r in regional if r['fdi'] == row['active_fdi']], axis=0)
    reported = np.array(row['wrenches'][str(row['active_fdi'])][:3])
    checks['regional_force_reconstructs_tooth'] = np.linalg.norm(force - reported) <= ge['numerical_wrench_parity_abs']
    corrupt = reported.copy()
    corrupt[0] += 0.01
    checks['injected_tooth_force_refuted_by_region_balance'] = np.linalg.norm(force - corrupt) > ge['force_balance_N']
    ref = json.loads((R / 'round1/results.json').read_text())
    altered = np.array([r['regional_N'] + 10 for r in ref['summary']])
    truth = np.array([r['external_N'] for r in ref['summary']])
    checks['injected_plus10N_external_refuted'] = bool(np.all(np.abs(altered - truth) / truth > ge['force_relative_error_max']))
    checks['failed_x10_injection_retained'] = ref['injection_checks']['force_x10_rejected'] is False
    checks['inherited_gap_corruption_refuted'] = json.loads((R / 'raw/INHERITED_GAP_CHECK.json').read_text())['corruption_rejected']
    r2 = json.loads((R / 'round2/results.json').read_text())
    checks['R2_corruptions_refuted'] = r2['gates']['injected_errors_rejected']
    r3 = json.loads((R / 'round3/results.json').read_text())
    checks['R3_calibration_gates'] = bool(r3['all_fixture_gates_pass'])
    r4 = json.loads((R / 'round4/results.json').read_text())
    checks['R4_envelope_gates'] = all(r4['gates'].values())
    r5 = json.loads((R / 'round5/results.json').read_text())
    checks['R5_unfit_gates'] = all(r5['gates'].values())
    r6 = json.loads((R / 'round6/results.json').read_text())
    checks['R6_stronger_refutation_retained'] = not r6['gates']['neighbour_transfer'] and r6['neighbour_groups_pass'] == 1
    checks['R6_corruptions_refuted'] = r6['gates']['corruptions_refuted']
    r7 = json.loads((R / 'round7/results.json').read_text())
    checks['R7_numerical_failure_retained'] = not r7['gates']['numerical_parity'] and (not r7['gates']['nonpenetration'])
    r8 = json.loads((R / 'round8/results.json').read_text())
    checks['R8_numerical_repair_passes'] = all((r8['gates'][k] for k in ['numerical_parity', 'force_balance', 'torque_balance', 'nonpenetration', 'null_modes']))
    checks['R8_neighbour_refutation_retained'] = sum((r['gate'] for r in r8['neighbour_rows'])) == 2 and (not r8['gates']['neighbours'])
    r9 = json.loads((R / 'round9/results.json').read_text())
    checks['R9_membrane_refutation_and_dropout_retained'] = not r9['gates']['neighbours'] and len(r9['failures']) == 2
    for i in [6, 7, 8, 9]:
        f = json.loads((R / f'FROZEN_PREDICTIONS_R{i}.json').read_text())
        if 'predictions' in f and 'predictions_sha256' in f:
            checks[f'R{i}_prediction_payload_immutable'] = hashlib.sha256(json.dumps(f['predictions'], sort_keys=True).encode()).hexdigest() == f['predictions_sha256']
    checks['whole_arch_corrupted_support_refuted'] = all((r['gates']['injected_support_force_rejected'] for r in [r7, r8, r9]))
    checks['whole_arch_plus10N_external_refuted'] = all((abs(r['predicted_N'] + 10 - r['force_N']) / abs(r['force_N']) > 0.25 for result in [r8, r9] for r in result['neighbour_rows']))
    arches = {a['case']: a for a in json.loads((R / 'inputs/ARCHES.json').read_text())}
    regions = json.loads((R / 'inputs/PARK_REGIONS.json').read_text())
    geometry = {case: build(arch, regions)['patch'] for (case, arch) in arches.items()}
    moment_residuals = []
    for row in json.loads((R / 'FROZEN_PREDICTIONS_R8.json').read_text())['predictions']['rows']:
        reactions = row['candidate']['regional']
        for (fdi, wrench) in row['candidate']['wrenches'].items():
            reconstructed = sum((np.cross(p['lever'], r['force_N']) for (p, r) in zip(geometry[row['case']], reactions) if str(p['fdi']) == fdi), np.zeros(3))
            moment_residuals.append((reconstructed, np.array(wrench[3:])))
    checks['whole_arch_regional_moment_reconstructs'] = all((np.max(np.abs(a - b)) <= 1e-05 for (a, b) in moment_residuals))
    checks['whole_arch_injected_moment_rejected'] = all((np.max(np.abs(a - (b + np.array([0.01, 0, 0])))) > 1e-05 for (a, b) in moment_residuals))
    observations = json.loads((R / 'raw/R3_FIXTURE_OBSERVATIONS.json').read_text())
    o = copy.deepcopy(observations[0])
    o['optical_housing_rigidity_residual_mm'] = 0.06
    checks['nonrigid_housing_refuted'] = fit([o], 6.0)['status'] == 'RIGID_HOUSING_CLOSURE_REFUTED'
    template = json.loads((R / 'inputs/MATCHED_LAB_MEASUREMENT_TEMPLATE.json').read_text())
    checks['missing_lab_measurements_refused'] = fit(template, 6.0)['status'] == 'NEEDS_MATCHED_MEASUREMENT'
    checks = {k: bool(v) for (k, v) in checks.items()}
    out = dict(checks=checks, all_checks_pass=all(checks.values()), source_verification=records, retained_failed_gate='R1 prereg force x10 corruption does NOT reject; preserve as failure. New +10N corruption does reject same force threshold.', corrections=[dict(field='raw/PREDICTIONS_R1.json relaxation[].wrench_difference', reason='Legacy Euclidean norm mixes N and Nmm; not a physical scalar', replacement='Package exports force difference N and moment difference Nmm separately. Frozen raw output is preserved.'), dict(field='held_30N simulation', reason='Housing translations reach ~3 mm, beyond small-motion closure', replacement='Conditional stress test only; no quantitative physical 30N seating/relaxation claim')])
    (R / 'raw/VERIFICATION.json').write_text(json.dumps(out, indent=2) + '\n')
    print('Verification:', sum(checks.values()), '/', len(checks), 'including retained failed x10 gate')
    if not out['all_checks_pass']:
        raise RuntimeError('Failed verification: ' + ','.join((k for (k, v) in checks.items() if not v)))
if __name__ == '__main__':
    main()
