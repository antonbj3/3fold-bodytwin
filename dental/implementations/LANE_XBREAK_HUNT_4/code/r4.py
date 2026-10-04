import json, time, resource
from pathlib import Path
import numpy as np
from scipy.stats import beta, binom
from scipy.optimize import brentq
from r1 import assert_freeze, dump

def exact_interval(k, n, alpha):
    return [0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1)), 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))]

def inverse_binomial(k, n, alpha):
    return [0.0 if k == 0 else brentq(lambda p: binom.sf(k - 1, n, p) - alpha / 2, 0, 1, xtol=1e-14), 1.0 if k == n else brentq(lambda p: binom.cdf(k, n, p) - alpha / 2, 0, 1, xtol=1e-14)]

def main():
    assert_freeze('R4')
    start = time.perf_counter()
    tables = np.array(json.loads(Path('raw/R3_TABLES.json').read_text()), int)
    alpha = 0.05 / 4
    checks = []
    cases = []
    for c in range(2):
        for p in range(2):
            pairs = sorted(set(((int(q[1, c, p]), int(q[:, c, p].sum())) for q in tables)))
            intervals = []
            for (k, n) in pairs:
                interval = exact_interval(k, n, alpha)
                oracle = inverse_binomial(k, n, alpha)
                checks.append({'C': c, 'P': p, 'k': k, 'n': n, 'interval': interval, 'oracle': oracle, 'error': max((abs(interval[i] - oracle[i]) for i in range(2)))})
                intervals.append(interval)
            interval = [min((r[0] for r in intervals)), max((r[1] for r in intervals))]
            cases.append({'C': c, 'P': p, 'conditional_population_interval': interval, 'width': interval[1] - interval[0], 'decision_at_0_8': 'D1' if interval[0] > 0.8 else 'D0' if interval[1] < 0.8 else 'UNKNOWN', 'transport_to_unselected_clinic': 'UNKNOWN'})
    error = max((r['error'] for r in checks))
    r = {'round': 'R4', 'claim_type': 'capability', 'alpha_per_pattern': alpha, 'cases': cases, 'raw_interval_checks': checks, 'gates': {'control_agreement': 'PASS' if error <= 1e-08 else 'FAIL', 'all_width_le_0_10': 'PASS' if all((c['width'] <= 0.1 for c in cases)) else 'FAIL', 'both_C1_decide_0_8': 'PASS' if all((c['decision_at_0_8'] != 'UNKNOWN' for c in cases if c['C'] == 1)) else 'FAIL'}, 'cluster_correct_coverage': 'UNKNOWN', 'unselected_clinic_and_treatment': 'UNKNOWN', 'cost': {'seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'inherited_cost': ['R1', 'R2', 'R3'], 'data_acquisition': 'No new patient data'}}
    dump('raw/R4_RESULTS.json', r)
    dump('CURRENT_WORK_STATE.json', {'tag': 'XBREAK-hunt-4', 'stage': 'R4_COMPLETE', 'latest_gate': r['gates'], 'next_operation': 'Change decision target: full assigned treatment pathway with intraoperative fallback, not absence of chamber bleeding.'})
    print(json.dumps({'cases': cases, 'gates': r['gates'], 'max_control_error': error, 'cost': r['cost']}, indent=2))
if __name__ == '__main__':
    main()
