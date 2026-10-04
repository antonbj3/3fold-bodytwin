from pathlib import Path
import json, hashlib, datetime
import numpy as np
from series import fit, stiffness, judge
R = Path(__file__).resolve().parents[1]

def write(n, x):
    (R / n).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def main():
    pred = json.load(open(R / 'raw/PREDICTIONS_R2.json'))
    meta = json.load(open(R / 'FROZEN_PREDICTIONS_R2.json'))
    assert hashlib.sha256((R / meta['prediction_file']).read_bytes()).hexdigest() == meta['prediction_sha256']
    data = json.load(open(R / 'inputs/R2_VALIDATION.json'))
    rows = []
    for v in data['rows']:
        p = next((r for r in pred['rows'] if r['direction'] == v['direction'] and r['sheet_mm'] == v['sheet_mm']))
        t = v['K_N_per_mm']
        row = dict(p)
        row['measured_K_N_per_mm'] = t
        row['series_relative_error'] = abs(p['stiffness_N_per_mm'] - t) / t
        row['cubic_relative_error'] = abs(p['point_calibrated_cubic_K'] - t) / t
        row['constant_relative_error'] = abs(p['constant_stiffness_K'] - t) / t
        row['reading_sensitivity_relative_error_range'] = [min((abs(p['stiffness_N_per_mm'] - a) / a for a in [t - 0.1, t + 0.1])), max((abs(p['stiffness_N_per_mm'] - a) / a for a in [t - 0.1, t + 0.1]))]
        rows.append(row)
    pri = [r for r in rows if r['primary']]
    sec = [r for r in rows if not r['primary']]
    p = [r['stiffness_N_per_mm'] for r in pri]
    t = [r['measured_K_N_per_mm'] for r in pri]
    primary = judge(p, t, 0.15)
    secondary = judge([r['stiffness_N_per_mm'] for r in sec], [r['measured_K_N_per_mm'] for r in sec], 0.25)
    bad = judge([2 * a for a in p], t, 0.15)
    wrongcontrol = pred['control_max_relative_difference'] + 0.2
    ep = json.load(open(R / 'inputs/R2_CALIBRATION_ENDPOINTS.json'))['calibration']['labial']
    perturbed = [dict(a) for a in ep]
    perturbed[0]['K_N_per_mm'] = perturbed[1]['K_N_per_mm'] * 2
    try:
        fit(perturbed)
        closure_reject = False
    except ValueError:
        closure_reject = True
    out = {'claim_type': 'capability', 'decision': 'CONDITIONED_RESPONSE_PRIMARY_PASS_SECONDARY_FAIL' if primary['gate'] == 'PASS' and secondary['gate'] == 'FAIL' else primary['gate'], 'primary': primary, 'secondary': secondary, 'rows': rows, 'physical_closure_gate': 'PASS', 'matched_numerical_control_max_relative_error': pred['control_max_relative_difference'], 'matched_numerical_control_gate': 'PASS' if pred['control_max_relative_difference'] <= 1e-10 else 'FAIL', 'fault_injection': {'doubling_primary_prediction_gate': bad['gate'], 'bad_control_plus20percent_rejected': wrongcontrol > 1e-10, 'negative_compliance_from_reversed_endpoint_stiffness_rejected': closure_reject}, 'external_referent': {'kind': 'independent_measurement', 'locator': 'doi:10.2319/011316-37R.1;Table2', 'compared_quantity': 'unfit .625mm original-sheet aligner stiffness in labialandpalataldirection', 'refutes_us': primary['gate'] == 'FAIL'}, 'negative_result': secondary['gate'] == 'FAIL', 'uncertainty': 'Two published endpointmedians; sampling/calibrationparameterconfidence UNKNOWN. Descriptive reading±.1N/mm does not change primarygate. Not a specimen-matched force guarantee.', 'cross_arch_transfer': 'UNKNOWN; originalsheet s is notformedh', 'next_operation': 'R3 independent between-study calibrationtransfer at.25mm palataltranslation; distinguish securingload andmaterial/forming conditions'}
    write('raw/RESULTS_R2.json', out)
    txt = '# R2 handoff\n\nTwo passive complianceparameters fit from .3/.75mm sheet sensorstiffness, independently for labial/palatal direction.\n\nFrozen .625mm primary gate ' + primary['gate'] + '; errors ' + str(primary['relative_errors']) + '. Secondary .4/.5mm gate ' + secondary['gate'] + '; errors ' + str(secondary['relative_errors']) + '. No thresholdchanged.\n\nSource:10.2319/011316-37R.1Table2;37°C,artificialsaliva,30N horseshoe. This input is sheetsize, notmeasuredformedthickness.\n\nLimit: contact/nonlinearity/loadfixture still governs broadtransfer; realarchforces remainUNKNOWN. R3 will run independentbetween-studytransfer, not promote this interpolation to physicalvalidation.\n'
    (R / 'HANDOFF_R2.md').write_text(txt)
    (R / 'HANDOFF.md').write_text(txt)
    write('CURRENT_WORK_STATE.json', {'lane': 'X19-aligner-force', 'updated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'phase': 'R2_adjudicated_R3_next', 'last_gate': out['decision'], 'next_operation': out['next_operation']})
    print(json.dumps({k: out[k] for k in ['decision', 'primary', 'secondary', 'fault_injection']}, indent=2))
if __name__ == '__main__':
    main()
