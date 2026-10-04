"""Compare independently measured regional vertical forces to a frozen query.
No parameter fit. Synthetic interface checks are explicitly not physical reference.
"""
from dental_release.paths import expand as _release_expand
import argparse, json
from pathlib import Path
import numpy as np
from contact_model import P, sha, write, local_path

def compare(record, prediction, geometry_sha):
    required = ['kind', 'locator', 'case', 'geometry_sha256', 'force_target_N', 'quantity', 'observation_matrix', 'values_N', 'applicability_evidence_locator']
    for k in required:
        if k not in record:
            raise ValueError('MISSING_' + k)
    if record['kind'] != 'independent_measurement':
        raise ValueError('NOT_INDEPENDENT_MEASUREMENT')
    if not record['locator'] or not record['applicability_evidence_locator']:
        raise ValueError('MISSING_LOCATOR')
    if record['geometry_sha256'] != geometry_sha:
        raise ValueError('GEOMETRY_HASH_MISMATCH')
    if record['case'] != prediction['case']:
        raise ValueError('CASE_MISMATCH')
    if record['quantity'] != 'regional_vertical_force_N':
        raise ValueError('WRONG_QUANTITY_OR_UNIT')
    if abs(record['force_target_N'] - prediction['force_control']['force_target_N']) > 1e-05:
        raise ValueError('LOAD_MISMATCH')
    matrix = np.asarray(record['observation_matrix'], float)
    values = np.asarray(record['values_N'], float)
    if matrix.ndim != 2 or matrix.shape != (len(values), len(prediction['rows'])):
        raise ValueError('OBSERVATION_SHAPE_MISMATCH')
    uncertainty_provided = 'uncertainty_N' in record
    uncertainty = np.asarray(record.get('uncertainty_N', np.zeros(len(values))), float)
    if uncertainty.shape != values.shape or np.any(uncertainty < 0):
        raise ValueError('UNCERTAINTY_INVALID')
    if not all((np.isfinite(a).all() for a in [matrix, values, uncertainty])):
        raise ValueError('NONFINITE_MEASUREMENT')
    predicted = matrix @ np.array([r['vector_on_upper_N'][2] for r in prediction['rows']])
    error = np.abs(predicted - values)
    accepted = error <= 5.0 + uncertainty
    referent_kind = 'our_own_fixture' if record['locator'].startswith('our_own_fixture:') else 'independent_measurement'
    return dict(claim_type='capability', predicted_N=predicted, measured_N=values, error_N=error, uncertainty_N=uncertainty if uncertainty_provided else None, uncertainty_status='PROVIDED' if uncertainty_provided else 'UNKNOWN', gate=bool(accepted.all()), criterion='each absolute residual<=5N+reported uncertainty; absent uncertainty checks only5N point residual, not physical validation', external_referent=dict(kind=referent_kind, locator=record['locator'], compared_quantity='regional_vertical_force_N', refutes_us=True), physical_validation='UNKNOWN' if not uncertainty_provided or referent_kind == 'our_own_fixture' else 'CONDITIONAL_ON_OBSERVATION_APPLICABILITY', physical_scope='Only supplied observation/preload/geometry applicability; no parameters refitted', applicability_evidence_locator=record['applicability_evidence_locator'])

def main():
    a = argparse.ArgumentParser()
    a.add_argument('measurement')
    a.add_argument('--prediction', default=str(P / 'rounds/R2/case001_force100.json'))
    a.add_argument('--out', default=str(P / 'INDEPENDENT_MEASUREMENT_COMPARISON.json'))
    args = a.parse_args()
    freeze = P / 'FROZEN_PREDICTIONS_R2.json'
    if sha(freeze) != (P / 'FROZEN_PREDICTIONS_R2.json.sha256').read_text().strip():
        raise ValueError('FREEZE_HASH_MISMATCH')
    requested = Path(args.prediction).resolve()
    entries = json.loads(freeze.read_text())['manifest']
    match = [x for x in entries if local_path(x['path']).resolve() == requested or (str(x['path']).startswith(_release_expand('@DENTAL_IMPLEMENTATIONS@/X54/')) and (P / Path(x['path']).relative_to(_release_expand('@DENTAL_IMPLEMENTATIONS@/X54'))).resolve() == requested)]
    if len(match) != 1 or sha(requested) != match[0]['sha256']:
        raise ValueError('PREDICTION_NOT_FROZEN_OR_HASH_CHANGED')
    rec = json.loads(Path(args.measurement).read_text())
    pred = json.loads(requested.read_text())
    original = json.loads((P / 'rounds/R1' / f"case{pred['case']:03d}_d0.05.json").read_text())
    g = sha(P / 'raw' / f"geometry_{pred['case']:03d}_h02.json")
    if g != original['geometry_sha256']:
        raise ValueError('ORIGINAL_GEOMETRY_HASH_CHANGED')
    result = compare(rec, pred, g)
    result['prediction_sha256'] = sha(args.prediction)
    write(args.out, result)
    print('PASS' if result['gate'] else 'FAIL')
if __name__ == '__main__':
    main()
