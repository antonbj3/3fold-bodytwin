"""Complete preregistered contrast/uncertainty scores without changing predictions."""
import json, pathlib, itertools, math
import numpy as np
from enclosures import certified_derivative
ROOT = pathlib.Path(__file__).resolve().parents[1]

def main():
    rows = [r for r in json.loads((ROOT / 'raw/arms.json').read_text()) if r['primary']]
    frozen = json.loads((ROOT / 'FROZEN_PREDICTIONS_R1.json').read_text())
    pairs = []
    keys = ['study', 'design', 'D_mm', 'L_mm', 'rho_gcc']
    for (a, b) in itertools.combinations(rows, 2):
        if any((a[k] != b[k] for k in keys)):
            continue
        variable = None
        if a['cortex_mm'] == b['cortex_mm'] and a['Df_mm'] != b['Df_mm']:
            variable = 'drilling'
        if a['Df_mm'] == b['Df_mm'] and a['cortex_mm'] != b['cortex_mm']:
            variable = 'cortex'
        if variable:
            pairs.append((variable, a, b))
    records = []
    for split in ['study', 'family']:
        for kind in ['practice', 'cortex', 'drilling', 'both']:
            lookup = {v['id']: v['predicted_Ncm'] for v in frozen['predictions'] if v['split'] == split and v['method'] == 'candidate_ridge' and (v['kind'] == kind)}
            for (variable, a, b) in pairs:
                observed = math.log(b['torque_Ncm'] / a['torque_Ncm'])
                predicted = math.log(lookup[b['id']] / lookup[a['id']])
                error = abs(predicted - observed)
                sign = np.sign(observed) == np.sign(predicted)
                records.append(dict(split=split, kind=kind, variable=variable, ids=[a['id'], b['id']], observed_log_ratio=observed, predicted_log_ratio=predicted, abs_log_effect_error=error, direction_correct=bool(sign), within_effect_tolerance=error <= math.log(1.5), both_criteria_met=bool(sign and error <= math.log(1.5)), level='POPULATION', locators=[a['locator'], b['locator']]))
    summary = []
    for (split, kind, variable) in itertools.product(['study', 'family'], ['practice', 'cortex', 'drilling', 'both'], ['drilling', 'cortex']):
        r = [x for x in records if (x['split'], x['kind'], x['variable']) == (split, kind, variable)]
        summary.append(dict(split=split, kind=kind, variable=variable, pairs=len(r), direction_correct=sum((x['direction_correct'] for x in r)), both_criteria_met=sum((x['both_criteria_met'] for x in r)), median_abs_log_effect_error=float(np.median([x['abs_log_effect_error'] for x in r]))))
    empirical = []
    lookup = {r['id']: r for r in rows}
    for split in ['study', 'family']:
        for kind in ['practice', 'both']:
            residuals = [abs(math.log(v['predicted_Ncm'] / lookup[v['id']]['torque_Ncm'])) for v in frozen['predictions'] if v['split'] == split and v['kind'] == kind and (v['method'] == 'candidate_ridge')]
            quantile = float(np.quantile(residuals, 0.95))
            empirical.append(dict(split=split, kind=kind, n=len(residuals), pooled_arm_95pct_abs_log_residual=quantile, retrospective_multiplicative_error_factor=math.exp(quantile), method='Descriptive pooled-arm residual quantile; not calibrated future coverage', rigorous_physical_enclosure='MISSING', level='POPULATION'))
    r3 = json.loads((ROOT / 'FROZEN_PREDICTIONS_R3.json').read_text())
    allrows = {r['id']: r for r in json.loads((ROOT / 'raw/arms.json').read_text())}
    enclosures = []
    for v in r3['predictions']:
        r = [allrows[i] for i in v['calibration_ids']]
        e = certified_derivative(v['poly_c_b_a'], r[0]['D_mm'], min((x['Df_mm'] for x in r)), max((x['Df_mm'] for x in r)))
        e['test_id'] = v['test_id']
        enclosures.append(e)
    out = dict(round='R1/R3 preregistered scoring completion', note='No refit, changed metric, post-hoc gate or altered frozen prediction', contrast_summary=summary, contrasts=records, empirical_uncertainty=empirical, certified_fitted_polynomial_derivative_bounds=enclosures, all_fitted_polynomials_certified_nonincreasing=all((e['certified_nonincreasing_fitted_torque'] for e in enclosures)))
    (ROOT / 'SCORES_SUPPLEMENT.json').write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(contrast_summary=[x for x in summary if x['kind'] == 'both'], empirical_uncertainty=empirical, derivative_bounds=enclosures), indent=2))
if __name__ == '__main__':
    main()
