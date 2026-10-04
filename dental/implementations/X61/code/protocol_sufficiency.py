"""Exact observation fibers: no fitting, numerical exploration or invented measurement."""
import json, math, hashlib
from pathlib import Path
from scipy.stats import norm
R = Path(__file__).resolve().parents[1]

def run():
    p = R / 'PREREG_R1.json'
    assert hashlib.sha256(p.read_bytes()).hexdigest() == (R / 'PREREG_R1.json.sha256').read_text().strip()
    a = {'cement_effective_offset_um': 64.0, 'sinter_factor': 0.75}
    b = {'cement_effective_offset_um': 32.75, 'sinter_factor': 0.78125}
    obs = lambda x: x['cement_effective_offset_um'] + 1000 * (x['sinter_factor'] - 0.75)
    states = [dict(summary=obs(x), downstream_post_length_mm=10 * x['sinter_factor'], **x) for x in [a, b]]
    mean_a = sum([4.0, 4.0, 4.0, 4.0]) / 4
    mean_b = sum([1.0, 1.0, 7.0, 7.0]) / 4
    tail_a = sum((x < 2 for x in [4.0, 4.0, 4.0, 4.0])) / 4
    tail_b = sum((x < 2 for x in [1.0, 1.0, 7.0, 7.0])) / 4
    pplus = norm.sf(2 / math.sqrt(3.5))
    pminus = norm.sf(2 / math.sqrt(0.5))
    out = dict(round='R1', claim_type='capability', outcome='PROTOCOL_ONLY_CAUSAL_LOCALIZATION_AND_PHYSICAL_PARAMETERS_NOT_IDENTIFIED', identical_summary_states=states, summary_identity_error=abs(obs(a) - obs(b)), downstream_length_difference_mm=abs(10 * a['sinter_factor'] - 10 * b['sinter_factor']), minimum_extension='Direct paired pre/post reference length plus independently assembled cement-film measurement', mean_force_sufficiency=dict(summary_mean_N=[mean_a, mean_b], identity_error=abs(mean_a - mean_b), downstream_lower_tail_probability=[tail_a, tail_b], difference=abs(tail_a - tail_b), minimum_extension='Keep specimen distribution, at least a tail statistic; mean alone cannot calibrate Weibull m'), posterior_mean_and_marginals_sufficiency=dict(summary_mean=[0.0, 0.0], summary_variances=[1.0, 1.0], identity_error=0.0, covariances=[0.75, -0.75], downstream_P_sum_above_2=[float(pplus), float(pminus)], difference=float(pplus - pminus), minimum_extension='Keep shared covariance or joint posterior, never independent marginal certificates'), physical_sigma0='UNKNOWN: force law identifies a force scale, not a material stress scale without calibrated specimen hazard', mechanism_vs_measurement='Same low recorded force can be true weak fracture or transcription error; trace consistency required', affine_enclosure='Missing for physical dry-gap map; collision is a synthetic mathematical witness', external_referent=json.loads(p.read_text())['external_referent'])
    (R / 'raw/R1_SUFFICIENCY.json').write_text(json.dumps(out, indent=2) + '\n')
    (R / 'HANDOFF_R1.md').write_text('R1: exact summary collision (0 identity error) gives 0.3125mm post-length difference. Mean force alone hides a 0.5 lower-tail probability difference. Same Gaussian marginal means/variances hide different downstream joint tails. Original CSV cannot identify physical sinter or cement film, nor separate weak specimen from bad force record. Next R2 changes acquisition: paired reference-length, assembled-film and trace-consistency sidecar, plus finite joint model posterior; preserve UNKNOWN without channels.\n')
    (R / 'CURRENT_WORK_STATE.json').write_text(json.dumps(dict(lane='X61-lab-alarm', status='R1_DECIDED', latest_gate=out['outcome'], next_operation='Freeze R2 profile and categorical sequential gates before generation'), indent=2) + '\n')
    print(out['outcome'])
if __name__ == '__main__':
    run()
