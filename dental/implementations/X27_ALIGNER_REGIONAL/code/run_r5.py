import json, hashlib, time
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from contact_model import make_system, solve
R = Path(__file__).resolve().parents[1]

def put(n, x):
    (R / n).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def main():
    start = time.perf_counter()
    pr = json.loads((R / 'PREREG_R5.json').read_text())
    ge = pr['gates']
    assert hashlib.sha256((R / 'PREREG_R5.json').read_bytes()).hexdigest() == json.loads((R / 'PREREG_R5.sha256.json').read_text())['sha256']
    r1 = json.loads((R / 'round1/results.json').read_text())
    cal = next((r for r in r1['summary'] if r['fdi'] == 11))
    gamma = cal['external_N'] / cal['regional_N']
    gamma_nominal = cal['external_N'] / cal['practice_N']
    pred = [dict(fdi=r['fdi'], force_N=gamma * r['regional_N'], nominal_calibrated_force_N=gamma_nominal * r['practice_N']) for r in r1['summary'] if r['fdi'] != 11]
    frozen = dict(frozen_utc=datetime.now(timezone.utc).isoformat(), prereg_sha256=hashlib.sha256((R / 'PREREG_R5.json').read_bytes()).hexdigest(), calibration_fdi=11, gain=gamma, nominal_ablation_gain=gamma_nominal, predictions=pred, source_visibility='All three published endpoint means seen before prereg; two types are unfit, not blinded')
    h = hashlib.sha256(json.dumps(dict(gain=gamma, predictions=pred), sort_keys=True).encode()).hexdigest()
    frozen['predictions_sha256'] = h
    if (R / 'FROZEN_PREDICTIONS_R5.json').exists():
        if json.loads((R / 'FROZEN_PREDICTIONS_R5.json').read_text())['predictions_sha256'] != h:
            raise RuntimeError('R5 prediction drift')
    else:
        put('FROZEN_PREDICTIONS_R5.json', frozen)
    scores = []
    controls = []
    injections = []
    for row in pred:
        ref = next((r for r in r1['summary'] if r['fdi'] == row['fdi']))
        truth = ref['external_N']
        error = abs(row['force_N'] - truth) / truth
        ab = abs(row['nominal_calibrated_force_N'] - truth) / truth
        scores.append(dict(**row, external_N=truth, external_SD_N=ref['external_SD_N'], relative_error=error, nominal_calibrated_relative_error=ab, uncalibrated_practice_relative_error=ref['practice_relative_error'], unfit_gate=error <= ge['unfit_relative_error_max'], resolution_level='PER_TOOTH', uncertainty='Mean-force type transfer only; external SD not an accuracy confidence interval'))
        corrupt = row['force_N'] + 10
        injections.append(dict(fdi=row['fdi'], injected_force_N=corrupt, relative_error=abs(corrupt - truth) / truth, rejected=abs(corrupt - truth) / truth > ge['unfit_relative_error_max']))
    regions = json.loads((R / 'inputs/PARK_REGIONS.json').read_text())
    arches = json.loads((R / 'inputs/ARCHES.json').read_text())
    raw = json.loads((R / 'raw/PREDICTIONS_R1.json').read_text())['rows']
    for row in raw:
        arch = next((a for a in arches if a['case'] == row['case']))
        system = make_system(arch, row['active_fdi'], regions)
        control = solve(system, method='primal', relaxation_factor=gamma)
        parity = max((float(np.max(np.abs(np.array(control['wrenches'][k]) - gamma * np.array(row['regional']['wrenches'][k])))) for k in control['wrenches']))
        controls.append(dict(case=row['case'], fdi=row['active_fdi'], max_wrench_difference=parity))
    gates = dict(unfit_types=all((r['unfit_gate'] for r in scores)), improvement=float(np.median([r['relative_error'] for r in scores])) <= (1 - ge['median_error_improvement_fraction']) * float(np.median([r['uncalibrated_practice_relative_error'] for r in scores])), scaled_primal_control=max((c['max_wrench_difference'] for c in controls)) <= ge['numerical_scaled_wrench_parity_abs'], wrong_force_rejected=all((r['rejected'] for r in injections)))
    result = dict(round='R5', claim_type='information_link', outcome='CALIBRATED_UNFIT_TYPE_TRANSFER_PASS' if all(gates.values()) else 'CALIBRATED_UNFIT_TYPE_TRANSFER_FAIL', external_referent=pr['external_referent'], rows=scores, gain=gamma, forbidden_interpretation_as_modulus_MPa=gamma * 2746, modulus_interpretation='REJECTED: source-conditioned omitted stiffness gain; not measured PETG E', gates=gates, controls=controls, injections=injections, source_visibility=frozen['source_visibility'], retained_prereg_x10_force_corruption_rejected=r1['injection_checks']['force_x10_rejected'], limitation='No table locator or paired specimen geometry in Kaur abstract; no independently held study; no physical individual-force/moment validation', resolution_level='PER_TOOTH', physical_validation='UNKNOWN', wall_s=time.perf_counter() - start, next_operation='Measure matched finished shell and same-arch baseline/three independent perturbations/holdout; test whether effective omitted stiffness is membrane, retained prestrain, housing or contact geometry')
    put('round5/results.json', result)
    put('CURRENT_WORK_STATE.json', dict(lane='X27-aligner-regional', state='R5_COMPLETE_PACKAGING', latest_gate=result['outcome'], next_operation='One-command replay, source/corruption verification, graph coverage proposal and final handoff'))
    (R / 'round5/HANDOFF.md').write_text('# R5 force-calibrated type handoff\n\n' + result['outcome'] + '\n\n' + json.dumps(scores, indent=2) + '\n\nOne measured central-incisor mean fitted an effective gain. Canine/premolar were unfit but known before prereg. No independent-study or specimen force validation. The gain cannot be called a material modulus. R1 x10 corruption rejection failed and is retained.\n\n' + result['next_operation'] + '\n')
    print(json.dumps(dict(outcome=result['outcome'], gain=gamma, rows=scores, gates=gates), indent=2))
if __name__ == '__main__':
    main()
