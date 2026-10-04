import json, time, resource
from pathlib import Path
import numpy as np
from scipy.stats import fisher_exact
from r1 import assert_freeze, dump
from r4 import exact_interval

def aggregate_loo_brier(good, bad):
    n = good + bad
    prediction_for_bad = (bad - 0.5) / n
    prediction_for_good = (bad + 0.5) / n
    return float(np.sum(bad * (1 - prediction_for_bad) ** 2 + good * prediction_for_good ** 2) / np.sum(n))

def literal_count_deletion(good, bad, pooled):
    loss = []
    for category in range(len(good)):
        for (label, num) in [(0, good[category]), (1, bad[category])]:
            g = good.copy()
            b = bad.copy()
            if label:
                b[category] -= 1
            else:
                g[category] -= 1
            prediction = (b.sum() + 0.5) / (g.sum() + b.sum() + 1) if pooled else (b[category] + 0.5) / (g[category] + b[category] + 1)
            loss.append(num * (prediction - label) ** 2)
    return float(sum(loss) / (good.sum() + bad.sum()))

def main():
    assert_freeze('R6')
    start = time.perf_counter()
    s = json.loads(Path('raw/R5_TRIAL_INPUTS.json').read_text())
    cats = ['normal', 'widening', 'lesion']
    counts = s['CBCT_completed_counts']
    good = np.array([counts['P'][c] for c in cats], int)
    bad = np.array([counts['fallback_RCT'][c] for c in cats], int)
    observed = []
    for (c, g, b) in zip(cats, good, bad):
        observed.append({'CBCT': c, 'completed_P': int(g), 'fallback_RCT': int(b), 'n': int(g + b), 'sample_risk': b / (g + b), 'IID_simultaneous95_interval': exact_interval(int(b), int(g + b), 0.05 / 3)})
    allocations = []
    controls = []
    for j in range(3):
        bf = bad.copy()
        bf[j] += 1
        abnormal_good = int(good[1:].sum())
        abnormal_bad = int(bf[1:].sum())
        table = [[abnormal_bad, abnormal_good], [int(bf[0]), int(good[0])]]
        f = fisher_exact(table, alternative='greater')
        diff = abnormal_bad / (abnormal_bad + abnormal_good) - bf[0] / (bf[0] + good[0])
        detailed = aggregate_loo_brier(good, bf)
        pooled = aggregate_loo_brier(np.array([good.sum()]), np.array([bf.sum()]))
        check = literal_count_deletion(good, bf, False)
        check_p = literal_count_deletion(good, bf, True)
        controls.append(max(abs(detailed - check), abs(pooled - check_p)))
        allocations.append({'missing_CBCT_assigned_to': cats[j], 'observed_counts_if_known': {'good': good.tolist(), 'bad': bf.tolist()}, 'risk_difference_abnormal_minus_normal': float(diff), 'one_sided_Fisher_p': float(f.pvalue), 'odds_ratio': float(f.statistic), 'LOO_Brier_CBCT': detailed, 'LOO_Brier_pooled': pooled, 'Brier_relative_gain': 1 - detailed / pooled, 'association_gate': bool(f.pvalue <= 0.05 and diff >= 0.15), 'forecast_gate': bool(1 - detailed / pooled >= 0.05)})
    observed_brier = aggregate_loo_brier(good, bad)
    pooled_observed = aggregate_loo_brier(np.array([good.sum()]), np.array([bad.sum()]))
    r = {'round': 'R6', 'claim_type': 'information_link', 'observed_categories': observed, 'missing_CBCT_cases': 1, 'fallback_assignment_denominator': 86, 'observed_only_LOO_gain': 1 - observed_brier / pooled_observed, 'all_missing_allocations': allocations, 'control_max_error': max(controls), 'gates': {'source_cohort_association': 'PASS' if all((r['association_gate'] for r in allocations)) else 'FAIL', 'all_missing_forecast_gain_ge_5percent': 'PASS' if all((r['forecast_gate'] for r in allocations)) else 'FAIL', 'literal_LOO_control': 'PASS' if max(controls) <= 1e-08 else 'FAIL'}, 'external_forecast_validation': 'UNKNOWN', 'clinical_scan_benefit': 'UNKNOWN', 'individual_treatment_choice': 'UNKNOWN', 'within_category_timing': 'UNKNOWN', 'cost': {'seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'prior_cost_inherited': ['R1', 'R2', 'R3', 'R4', 'R5'], 'clinical_acquisition': 'Published primary counts; no scan performed or recommended.'}}
    dump('raw/R6_RESULTS.json', r)
    dump('CURRENT_WORK_STATE.json', {'tag': 'XBREAK-hunt-4', 'stage': 'R6_COMPLETE', 'latest_gate': r['gates'], 'next_operation': 'Freeze R7: preserve correlation of preop category and within-procedure timing. Construct opposite timing worlds with same published means,SD and path counts.'})
    print(json.dumps(r, indent=2))
if __name__ == '__main__':
    main()
