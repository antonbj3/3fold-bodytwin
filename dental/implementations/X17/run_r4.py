import json, time, hashlib
import numpy as np
from run_r1 import ROOT, write, freeze, state, gate

def main():
    start = time.monotonic()
    p = ROOT / 'PREREG_R4.json'
    assert hashlib.sha256(p.read_bytes()).hexdigest() == (ROOT / (p.name + '.sha256')).read_text().split()[0]
    groups = json.loads((ROOT / 'raw/round_R1.json').read_text())['external_comparisons']
    pred = []
    for (j, tgt) in enumerate(groups):
        src = groups[1 - j]
        kap = (src['ap1_mm'] - src['ap0_mm']) / src['dose_mm']
        kl = (src['la1_mm'] - src['la0_mm']) / src['dose_mm']
        kv = (src['v1_cm3'] / src['v0_cm3'] - 1) / src['dose_mm']
        beta = np.linalg.lstsq([[src['dose_mm']]], [[src['ap1_mm'] - src['ap0_mm'], src['la1_mm'] - src['la0_mm']]], rcond=None)[0][0]
        area = tgt['area0_mm2'] * (1 + kap * tgt['dose_mm'] / tgt['ap0_mm']) * (1 + kl * tgt['dose_mm'] / tgt['la0_mm'])
        baseline = tgt['area0_mm2'] * (1 + beta[0] * tgt['dose_mm'] / tgt['ap0_mm']) * (1 + beta[1] * tgt['dose_mm'] / tgt['la0_mm'])
        vol = tgt['area0_mm2'] * (1 + kv * tgt['dose_mm'])
        pred.append(dict(training_group=src['group'], target_group=tgt['group'], k_AP=kap, k_lateral=kl, k_volume_relative=kv, predicted_area_mm2=area, volume_control_area_mm2=vol, strong_control_area_mm2=baseline))
    digest = freeze('FROZEN_PREDICTIONS_R4.json', {'predictions': pred, 'target_post_fields_used_in_prediction': []})
    rows = []
    for (x, tgt) in zip(pred, groups):
        y = tgt['area1_observed_mm2']
        err = abs(x['predicted_area_mm2'] - y)
        rows.append(dict(x, observed_area_mm2=y, relative_error=err / y, absolute_error_mm2=err, volume_control_absolute_error_mm2=abs(x['volume_control_area_mm2'] - y), strong_match_error_mm2=abs(x['strong_control_area_mm2'] - x['predicted_area_mm2'])))
    ratio = sum((x['absolute_error_mm2'] for x in rows)) / sum((x['volume_control_absolute_error_mm2'] for x in rows))
    controls = {'held_out_area': all((gate(x['relative_error'], 0.2) for x in rows)), 'advantage_vs_volume': gate(ratio, 0.5), 'conventional_regression_match': all((gate(x['strong_match_error_mm2'], 1e-10) for x in rows))}
    out = {'round': 'R4', 'outcome': 'FAIL_SHARED_RESPONSE_TRANSFER' if not controls['held_out_area'] else 'AGGREGATE_ONLY', 'comparisons': rows, 'controls': controls, 'injected_error_rejected': {'held_out_area': not gate(1.0, 0.2), 'advantage_vs_volume': not gate(1.0, 0.5), 'conventional_regression_match': not gate(1.0, 1e-10)}, 'candidate_to_volume_error_ratio': ratio, 'prediction_sha256': digest, 'cost': {'seconds': time.monotonic() - start, 'source_groups': 1, 'held_out_groups': 1, 'independent_studies': 1}, 'external_referent': json.loads(p.read_text())['external_referent'], 'implication': 'AP/lateral post measurement improves reconstruction, but a common transfer law inferred from one group does not provide a validated baseline-to-response curve.'}
    write('raw/round_R4.json', out)
    state('R4 adjudicated', out['outcome'], 'Package lab demo and next paired-movement acquisition construction')
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    main()
