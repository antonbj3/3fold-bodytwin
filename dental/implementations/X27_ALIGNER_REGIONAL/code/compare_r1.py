import json, hashlib
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parents[1]

def put(n, x):
    (R / n).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def main():
    pr = json.loads((R / 'PREREG_R1.json').read_text())
    fr = json.loads((R / 'FROZEN_PREDICTIONS_R1.json').read_text())
    p = fr['predictions']
    truth = json.loads((R / 'inputs/KAUR_REFERENT.json').read_text())
    summ = []
    for (fdi, ref) in zip([11, 13, 15], truth['rows']):
        rs = [r for r in p['rows'] if r['active_fdi'] == fdi]
        if not rs:
            raise RuntimeError('No predictions for external tooth type')
        row = dict(fdi=fdi, external_N=ref['force_N'], external_SD_N=ref['force_SD_N'], resolution_level='PER_TOOTH', uncertainty='Median across three unmatched source arches; external SD is specimen scatter, not validation CI')
        for k in ['regional', 'nominal_contact', 'tooth_type_only']:
            v = [r[k]['active_buccolingual_N'] for r in rs]
            row[k + '_N'] = float(np.median(v))
            row[k + '_relative_error'] = float(abs(np.median(v) - ref['force_N']) / ref['force_N'])
            row[k + '_range_N'] = [min(v), max(v)]
        v = [r['practice']['active_buccolingual_force_N'] for r in rs]
        row['practice_N'] = float(np.median(v))
        row['practice_relative_error'] = float(abs(np.median(v) - ref['force_N']) / ref['force_N'])
        row['held_30N_N'] = float(np.median([r['held_30N']['active_buccolingual_N'] for r in rs]))
        summ.append(row)
    ge = pr['gates']
    candidate = np.array([r['regional_relative_error'] for r in summ])
    base = np.array([r['practice_relative_error'] for r in summ])
    maxres = {k: max((r['regional']['residuals'][k] for r in p['rows'])) for k in p['rows'][0]['regional']['residuals']}
    gates = dict(external_absolute=bool(np.all(candidate <= ge['force_relative_error_max'])), improvement=bool(np.median(candidate) <= (1 - ge['median_error_improvement_fraction']) * np.median(base)), numerical_parity=bool(max((c['max_wrench_difference'] for c in p['controls'])) <= ge['numerical_wrench_parity_abs']), nonpenetration=maxres['penetration_mm'] <= ge['nonpenetration_mm'], complementarity=maxres['complementarity_Nmm'] <= ge['complementarity_Nmm'], stationarity=maxres['stationarity_N'] <= ge['stationarity_N'], force_balance=maxres['global_force_balance_N'] <= ge['force_balance_N'], torque_balance=maxres['global_torque_balance_Nmm'] <= ge['torque_balance_Nmm'])
    inj = dict(force_x10_rejected=bool(np.any(np.abs(10 * np.array([r['regional_N'] for r in summ]) - np.array([r['external_N'] for r in summ])) / np.array([r['external_N'] for r in summ]) > ge['force_relative_error_max'])), negative_gap_rejected=0.01 > ge['nonpenetration_mm'], wrong_wrench_control_rejected=0.01 > ge['numerical_wrench_parity_abs'], unbalanced_force_rejected=0.01 > ge['force_balance_N'], unbalanced_torque_rejected=0.01 > ge['torque_balance_Nmm'])
    result = dict(round='R1', claim_type='information_link', outcome='REGIONAL_CONTACT_MODEL_EXTERNAL_TRANSFER_PASS' if gates['external_absolute'] and gates['improvement'] else 'REGIONAL_CONTACT_MODEL_EXTERNAL_TRANSFER_FAIL', external_referent=pr['external_referent'], summary=summ, gates=gates, injection_checks=inj, max_residuals=maxres, median_candidate_relative_error=float(np.median(candidate)), median_practice_relative_error=float(np.median(base)), insertion_endpoint_path_difference=p['insertion_endpoint_path_difference'], insertion_result='Unique convex frictionless model has no retained path/branch; cannot represent field reduced-cap branch result', relaxation_max_failure_of_scalar_scaling=max((r['wrench_difference'] for r in p['relaxation'])), resolution_level=['PER_SURFACE_REGION', 'PER_TOOTH'], physical_moments_validated=False, prediction_sha256=fr['predictions_sha256'], source_visibility=pr['source_visibility'], negative_result=not (gates['external_absolute'] and gates['improvement']), next_construction='Source-conditioned measured tangent wrench port, calibrated on Cho .3 mm groups and tested on unfit .6 mm groups; preserve separate tooth/material and apex origin')
    put('round1/results.json', result)
    (R / 'round1/HANDOFF.md').write_text(f"# R1 regional contact handoff\n\nOutcome: {result['outcome']}.\nGates: {json.dumps(gates)}.\nMedian external force error: {result['median_candidate_relative_error']:.4f}; practice proxy: {result['median_practice_relative_error']:.4f}.\n\nPhysical transfer remains unvalidated. Park marginals are not a matched tooth-region map. Convex, fixed-normal contact has a unique final state and cannot retain the field-cap starting branch.\n\nNext construction: {result['next_construction']}\n")
    put('CURRENT_WORK_STATE.json', dict(lane='X27-aligner-regional', state='R1_COMPLETE_R2_STARTING', latest_gate=result['outcome'], next_operation=result['next_construction']))
    print(json.dumps(dict(outcome=result['outcome'], summary=summ, gates=gates), indent=2))
if __name__ == '__main__':
    main()
