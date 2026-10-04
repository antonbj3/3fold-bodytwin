import math
import numpy as np
from scipy.stats import t
from common import read, save, state, verify_freeze, ROOT

def main():
    verify_freeze('PREREG_R1.json')
    verify_freeze('FROZEN_PREDICTIONS_R1.json')
    base = read('cbct_baseline')
    held = read('cbct_held')
    mats = [v['material'] for v in base[0]['inserts']]
    ref = next((r for r in base if r['device'] == 'Varian'))
    h = np.array([v['gray_median'] for v in ref['inserts']])
    (i, j) = (mats.index('AIR'), mats.index('Polystyrene'))
    cb = []
    for q in base:
        if q['device'] == 'Varian':
            continue
        g = np.array([x['gray_median'] for x in q['inserts']])
        slope = (h[j] - h[i]) / (g[j] - g[i])
        intercept = h[i] - slope * g[i]
        selected = [x for x in held if x['device'] == q['device'] and 'test-retest2_' in x['zip_member'] or (x['device'] == q['device'] and 'test-retest3_' in x['zip_member'])]
        greps = np.vstack([g] + [[v['gray_median'] for v in x['inserts']] for x in selected])
        residuals = intercept + slope * greps - h
        drift = np.ptp(greps, axis=0) * abs(slope)
        avg = intercept + slope * greps.mean(axis=0) - h
        cb.append({'device': q['device'], 'materials': mats, 'repeat_scans': len(greps), 'repeat_gray_medians': greps.tolist(), 'repeat_range_HU_ref': drift.tolist(), 'baseline_signed_error_HU_ref': residuals[0].tolist(), 'mean_of_repeats_signed_error_HU_ref': avg.tolist(), 'baseline_max_abs_error_HU_ref': float(abs(residuals[0]).max()), 'repeat_mean_max_abs_error_HU_ref': float(abs(avg).max()), 'independent_control_max_error': float(abs(h[i] + (greps - g[i]) / (g[j] - g[i]) * (h[j] - h[i]) - (intercept + slope * greps)).max()), 'resolution': 'PER_SURFACE_REGION', 'timescale': 'SIMULTANEOUS', 'scope': 'Three reported repeat scans; finite observed range, no unseen-case enclosure'})
    a = np.array([-1.0, 1.0])
    b = a[::-1].copy()
    truth = np.array([0.5, 0.5])
    identity = max(abs(abs(a).mean() - abs(b).mean()), abs(abs(a).max() - abs(b).max()))
    decisions_a = (truth + a <= 0).tolist()
    decisions_b = (truth + b <= 0).tolist()
    suff = {'summary_A': {'MAE': float(abs(a).mean()), 'max_abs': float(abs(a).max())}, 'summary_B': {'MAE': float(abs(b).mean()), 'max_abs': float(abs(b).max())}, 'identity_error': float(identity), 'residual_A': a.tolist(), 'residual_B': b.tolist(), 'truth': truth.tolist(), 'threshold': 0, 'decisions_A': decisions_a, 'decisions_B': decisions_b, 'downstream_decision_difference_count': sum((x != y for (x, y) in zip(decisions_a, decisions_b))), 'minimum_extension': 'Signed error interval keyed by decision location and state; no claim that this alone identifies unobserved physical errors', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'code/run_r1.py exact residual permutation', 'compared_quantity': 'Summary sufficiency only', 'refutes_us': True}}
    gap = read('net')['edges']
    gap = next((e for e in gap if e['id'] == 'D-E-K23'))['gap']['value']
    old_mae_sd = 0.71 * math.sqrt(2 / math.pi)
    correction = {'reported_proxy_MAE_degC': gap, 'swarm_other_cohort_sd_degC': 0.71, 'swarm_gaussian_expected_abs_degC': old_mae_sd, 'swarm_quadrature_number_degC': math.sqrt(gap ** 2 - old_mae_sd ** 2), 'identified_model_error_degC': None, 'identified_sensor_noise_degC': None, 'reason': 'Unpaired rotary cohort SD mixes biology/process/sensor; laser proxy comparison is another observable. MAE has no additive variance identity.', 'resolution': 'POPULATION', 'source': 'BT-DW48-NET-K23 RESULTS.md checked against R3/X33 observation contract', 'claims_rejected': ['94 percent is sensor noise', 'physics error is 0.20degC', 'proxy MAE is an individual error bound']}
    p = read('preload')
    n = p['n']
    delta = p['first']['mean_N'] - p['tenth']['mean_N']
    hw = t.ppf(0.975, n - 1) * (p['first']['sd_N'] + p['tenth']['sd_N']) / math.sqrt(n)
    preload = {'difference_N': delta, 'worst_covariance_pointwise95_interval_N': [delta - hw, delta + hw], 'sampling_closure': 'Normal iid specimen pairs; rho in [-1,1], worst SD=s_first+s_tenth; not individual/support guarantee', 'published_adjusted_contrast_intervals': p['rows'], 'resolution': 'POPULATION', 'decision': 'Same-protocol mean preload declines; mechanical clamp threshold UNKNOWN', 'locator': 'doi:10.1186/s40729-023-00473-3 Table1 and Results'}
    out = {'claim_type': 'capability', 'cbct_repeat_test': cb, 'sufficiency': suff, 'noise_claim_correction': correction, 'preload_group_answer': preload, 'full_physical_chain_closed': 0, 'outcome': 'LOCATED_DECISION_STATE_REQUIRED; EXISTING_REPEAT_DATA_DO_NOT_REMOVE_CALIBRATION_BIAS', 'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    save('rounds/R1_results.json', out)
    (ROOT / 'rounds/HANDOFF_R1.md').write_text('R1: scalar MAE/max fail exact sufficiency (identity0; two located decisions change). Other-cohort SD does not identify laser sensor error. RadRepro three-repeat mean retains systematic calibration error. Next construction: decision-specific finite anchor selection and covariance-robust experiment allocation, with genuine held protocols.\n')
    state('R1_COMPLETE', out['outcome'], 'Freeze R2: finite anchor holdout design and minimum population-contrast specimen allocation')
    print('R1', [(x['device'], x['baseline_max_abs_error_HU_ref'], x['repeat_mean_max_abs_error_HU_ref']) for x in cb], preload['worst_covariance_pointwise95_interval_N'])
if __name__ == '__main__':
    main()
