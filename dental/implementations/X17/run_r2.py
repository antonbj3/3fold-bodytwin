import json, time, hashlib
import numpy as np
from run_r1 import ROOT, write, freeze, state, gate

def main():
    start = time.monotonic()
    p = ROOT / 'PREREG_R2.json'
    assert hashlib.sha256(p.read_bytes()).hexdigest() == (ROOT / (p.name + '.sha256')).read_text().split()[0]
    src = json.loads((ROOT / 'raw/round_R1.json').read_text())['external_comparisons']
    pred = []
    for x in src:
        det = x['ap1_mm'] / x['ap0_mm'] * (x['la1_mm'] / x['la0_mm'])
        a = x['area0_mm2'] * det
        strong = x['area0_mm2'] * float(np.linalg.det(np.diag([x['ap1_mm'] / x['ap0_mm'], x['la1_mm'] / x['la0_mm']])))
        pred.append({'group': x['group'], 'area1_predicted_mm2': a, 'strong_control_mm2': strong, 'dose_mm': x['dose_mm']})
    digest = freeze('FROZEN_PREDICTIONS_R2.json', {'predictions': pred, 'input_reference': 'raw/round_R1.json external AP/lateral observations; area1 target not used in formula'})
    rows = []
    for (x, y) in zip(src, pred):
        target = x['area1_observed_mm2']
        cand = y['area1_predicted_mm2']
        err = abs(cand - target)
        wrong = cand * 2
        rows.append(dict(y, observed_area1_mm2=target, absolute_error_mm2=err, relative_error=err / target, predicted_change_per_mm_mm2=(cand - x['area0_mm2']) / x['dose_mm'], volume_control_absolute_error_mm2=abs(x['proportional_volume_predicted_area1_mm2'] - target), strong_control_difference_mm2=abs(cand - y['strong_control_mm2']), injected_lateral_scale=2, injected_relative_error=abs(wrong - target) / target, injected_error_rejected=not gate(abs(wrong - target) / target, 0.2)))
    ratio = sum((r['absolute_error_mm2'] for r in rows)) / sum((r['volume_control_absolute_error_mm2'] for r in rows))
    controls = {'external_area': all((gate(r['relative_error'], 0.2) for r in rows)), 'advantage_vs_volume': gate(ratio, 0.5), 'strong_affine_match': all((gate(r['strong_control_difference_mm2'], 1e-10) for r in rows))}
    out = {'round': 'R2', 'outcome': 'INFORMATION_GAIN_WITH_METHOD_TIE' if all(controls.values()) else 'FAIL', 'external_comparisons': rows, 'controls': controls, 'injected_error_rejected': {'external_area': all((r['injected_error_rejected'] for r in rows)), 'advantage_vs_volume': not gate(2.0, 0.5), 'strong_affine_match': not gate(10.0, 1e-10)}, 'candidate_to_volume_absolute_error_ratio': ratio, 'prediction_sha256': digest, 'cost': {'seconds': time.monotonic() - start, 'fit_parameters': 0, 'post_measurements_required_per_cross_section': ['AP diameter mm', 'lateral diameter mm']}, 'limitations': ['Published group means, two groups in one cohort: not two independent studies.', 'Post-diameters are measured input, not predicted by baseline geometry.', 'Minimum plane can change; determinant of mean diameters does not equal mean individual determinant.', 'No toothfairy case has a paired MAD scan; no patient validation.'], 'external_referent': json.loads(p.read_text())['external_referent']}
    write('raw/round_R2.json', out)
    state('R2 adjudicated', out['outcome'], 'R3 replace fixed transfer coefficient with spatial feasible set constrained by volume; test whether minimum is identifiable')
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    main()
