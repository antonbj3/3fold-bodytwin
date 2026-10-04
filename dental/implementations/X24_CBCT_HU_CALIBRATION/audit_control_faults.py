"""Delivery checks with altered values; does not change historical gates."""
from pathlib import Path
import hashlib
import json
import os
import numpy as np
from calibration import apply, bone_branch_b
from run_r2 import lagrange
P = Path(__file__).resolve().parent

def audit():
    r1 = json.loads((P / 'results_R1.json').read_text())
    r2 = json.loads((P / 'results_R2.json').read_text())
    r3 = json.loads((P / 'results_R3.json').read_text())
    frozen = json.loads((P / 'FROZEN_RESPONSE_R2.json').read_text())
    rows = []

    def record(name, clean_value, fault_value, threshold, mode='upper'):
        if mode == 'upper':
            (clean_pass, fault_pass) = (clean_value <= threshold, fault_value <= threshold)
        else:
            (clean_pass, fault_pass) = (clean_value > threshold, fault_value > threshold)
        rows.append({'control': name, 'clean_value': float(clean_value), 'injected_value': float(fault_value), 'threshold': threshold, 'predicate': mode, 'clean_pass': bool(clean_pass), 'injected_pass': bool(fault_pass), 'fault_rejected': bool(not fault_pass)})
    row = next((x for x in r1['rows'] if x['device'] != 'Varian'))
    coeff = np.array(next((x for x in r1['scanner_calibration'] if x['device'] == row['device']))['a_inverse_Href'])
    altered = coeff.copy()
    altered[0] += 1
    record('R1 equal-information arithmetic, intercept +1 HU', row['control_disagreement_HU'], abs(float(apply(altered, row['g'])) - row['H_calibrated']), 1e-08)
    record('Legacy air measurement, +10 gray', max((abs(x['air_replay_error']) for x in r1['scanner_calibration'])), 10, 1)
    row = next((x for x in r2['rows'] if x['variant'] == 'quadratic3_candidate' and x['device'] != 'Varian'))
    model = frozen['models'][row['device']]
    reference = float(lagrange(model['x'], model['y'], row['gray']))
    record('R2 equal-information Lagrange, predicted offset +1 HU', abs(row['H_calibrated'] - reference), abs(row['H_calibrated'] + 1 - reference), 1e-08)
    record('R2 increasing response, derivative sign corrupted', min((x['min_derivative'] for x in r2['monotonicity'])), -1, 0, 'lower')
    row = next((x for x in r3['answers'] if x['device'] != 'Varian'))
    H = row['H_reference_closure']
    wrong = float(bone_branch_b(H + 1)['E_MPa']) / row['E_reference_closure_MPa'] - 1
    record('R3 affine inverse, recovered H offset +1 HU', r3['max_affine_inverse_relative_error'], abs(wrong), 1e-09)
    record('R3 shared-anchor invariance, one anchor +100 gray', max((x['normalization_relative_error'] for x in r3['fault_injections'])), min((x['anchor_plus100gray_relative_error'] for x in r3['fault_injections'])), 1e-09)
    status_errors = sum((c['status'] != 'ABSTAIN_UNKNOWN' for c in r3['K42_certificates']))
    record('K42 UNKNOWN required, one output changed to ACCEPT_CONDITIONAL', status_errors, status_errors + 1, 0)
    blob = (P / 'PREREG_R1.json').read_bytes()
    expected = (P / 'PREREG_R1.sha256').read_text().strip()
    record('Frozen artifact identity, one extra byte', int(hashlib.sha256(blob).hexdigest() != expected), int(hashlib.sha256(blob + b' ').hexdigest() != expected), 0)
    delta = row['HU_symmetric_budget_5percent_E']
    E = row['E_reference_closure_MPa']
    errors = [max(abs(float(bone_branch_b(H + k * delta)['E_MPa']) / E - 1), abs(float(bone_branch_b(H - k * delta)['E_MPa']) / E - 1)) for k in [1, 2]]
    record('Consumer 5% E budget, double permitted HU error', errors[0], errors[1], 0.05 + 1e-12)
    d = r3['dropout']
    conservation = abs(d['source_sites'] - d['retained_sites'] - d['rejected_sites'])
    record('Dropout conservation, one extra retained site', conservation, conservation + 1, 0)
    record('Regional accuracy already fails; target +1000 HU also rejected', max((abs(x['error_HU']) for x in r2['rows'])), min((abs(x['error_HU'] + 1000) for x in r2['rows'])), 40)
    report = {'status': 'PASS' if all((x['fault_rejected'] for x in rows)) else 'FAIL', 'historical_R2_plus100_HU_fault_gate': r2['gates']['error_injection'], 'scope': 'Delivery fault checks; historical PREREG/gates preserved, no accuracy improvement claimed', 'controls': rows}
    out = Path(os.environ.get('X24_OUTPUT_DIR', str(P)))
    (out / 'CONTROL_FAULT_AUDIT.json').write_text(json.dumps(report, indent=2) + '\n')
    assert report['status'] == 'PASS'
    print('PASS', len(rows), 'delivery control fault checks; R2 original fault gate remains FALSE')
    return report
if __name__ == '__main__':
    audit()
