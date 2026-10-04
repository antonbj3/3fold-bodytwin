import itertools, json, time, resource
from pathlib import Path
import numpy as np
from r1 import assert_freeze, constraints, ratio_lp, dump

def joint_rows(atoms):
    return (np.array([[int(c == 1 and p == j) for (d, c, p) in atoms] for j in [0, 1, 2]], float), np.array([116, 264, 0], float))

def main():
    assert_freeze('R2')
    t = time.perf_counter()
    s = json.loads(Path('raw/R1_SOURCE_COUNTS.json').read_text())
    (atoms, A, b) = constraints(s)
    (extra, totals) = joint_rows(atoms)
    cases = []
    for (c, p) in itertools.product(range(2), repeat=2):
        oracle = ratio_lp(A, b, atoms, c, p, extra, totals)
        interval = [r['probability'] for r in oracle]
        cases.append({'C': c, 'P': p, 'posterior_interval': interval, 'width': interval[1] - interval[0], 'oracle': oracle, 'decisions': [{'threshold_scenario': k / 10, 'decision': 'D1' if interval[0] > k / 10 + 1e-08 else 'D0' if interval[1] < k / 10 - 1e-08 else 'UNKNOWN'} for k in range(1, 10)]})
    low = max(0, s['cold_margins'][1][1] + s['percussion_margins'][1][1] - s['class_totals'][1], 264 - min(s['cold_margins'][0][1], s['percussion_margins'][0][1]))
    high = min(s['cold_margins'][1][1], s['percussion_margins'][1][1], 264)
    check11 = [low / 264, high / 264]
    check10 = [(307 - high) / 116, (307 - low) / 116]
    err = max((abs(cases[3]['posterior_interval'][i] - check11[i]) for i in range(2)))
    err = max(err, max((abs(cases[2]['posterior_interval'][i] - check10[i]) for i in range(2))))
    r = {'round': 'R2', 'claim_type': 'information_link', 'cases': cases, 'gates': {'all_width_le_0_10': 'PASS' if all((c['width'] <= 0.1 for c in cases)) else 'FAIL', 'fixed_denominator_control': 'PASS' if err <= 1e-08 else 'FAIL'}, 'count_control_max_error': err, 'outcome': 'Unlabeled joint tally alone does not decide all patterns. No new reference-conditioned observation was supplied.', 'cost': {'total_seconds': time.perf_counter() - t, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'prospective_acquisition': 'NOT_RUN', 'prior_R1_cost_inherited': True}}
    dump('raw/R2_RESULTS.json', r)
    dump('CURRENT_WORK_STATE.json', {'tag': 'XBREAK-hunt-4', 'stage': 'R2_COMPLETE', 'latest_gate': r['gates'], 'next_operation': 'Freeze R3; add rounded joint reference PPV with source precision and enumerate every admissible integer table.'})
    print(json.dumps({'gates': r['gates'], 'cases': [{k: c[k] for k in ['C', 'P', 'posterior_interval', 'width']} for c in cases], 'cost': r['cost']}, indent=2))
if __name__ == '__main__':
    main()
