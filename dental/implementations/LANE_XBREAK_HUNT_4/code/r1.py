import json, time, itertools, resource, hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')

def assert_freeze(name):
    p = Path(f'PREREG_{name}.json')
    assert hashlib.sha256(p.read_bytes()).hexdigest() == Path(f'PREREG_{name}.sha256').read_text().strip()

def constraints(source):
    atoms = list(itertools.product(range(2), range(3), range(3)))
    A = []
    b = []
    for d in range(2):
        for c in range(3):
            A.append([int(i == d and j == c) for (i, j, k) in atoms])
            b.append(source['cold_margins'][d][c])
        for p in range(3):
            A.append([int(i == d and k == p) for (i, j, k) in atoms])
            b.append(source['percussion_margins'][d][p])
    return (atoms, np.asarray(A, float), np.asarray(b, float))

def ratio_lp(A, b, atoms, c, p, extra_A=None, extra_b=None):
    if extra_A is not None:
        A = np.vstack([A, extra_A])
        b = np.r_[b, extra_b]
    event = np.array([int(j == c and k == p) for (i, j, k) in atoms], float)
    disease = np.array([int(i == 1 and j == c and (k == p)) for (i, j, k) in atoms], float)
    ae = np.vstack([np.c_[A, -b], np.r_[event, 0]])
    be = np.r_[np.zeros(len(b)), 1]
    results = []
    for sign in [1, -1]:
        r = linprog(np.r_[sign * disease, 0], A_eq=ae, b_eq=be, bounds=(0, None), method='highs')
        if not r.success:
            raise RuntimeError(r.message)
        x = r.x[:-1] / r.x[-1]
        results.append({'probability': float(disease @ r.x[:-1]), 'counts': x.tolist(), 'max_count_residual': float(np.max(np.abs(A @ x - b))), 'fractional_residual': float(np.max(np.abs(ae @ r.x - be)))})
    return results

def main():
    assert_freeze('R1v2')
    start = time.perf_counter()
    source = json.loads(Path('raw/R1_SOURCE_COUNTS.json').read_text())
    (atoms, A, b) = constraints(source)
    rows = []
    query_times = []
    for (c, p) in itertools.product(range(2), repeat=2):
        t = time.perf_counter()
        bounds = []
        for d in range(2):
            nd = source['class_totals'][d]
            mc = source['cold_margins'][d][c]
            mp = source['percussion_margins'][d][p]
            bounds.append([max(0, mc + mp - nd), min(mc, mp)])
        (v, n) = bounds
        lo = n[0] / (n[0] + v[1])
        hi = n[1] / (n[1] + v[0])
        query_times.append(time.perf_counter() - t)
        oracle = ratio_lp(A, b, atoms, c, p)
        expected = [source['cold_margins'][d][c] * source['percussion_margins'][d][p] / source['class_totals'][d] for d in range(2)]
        prod = expected[1] / sum(expected)
        decisions = []
        for threshold in np.arange(1, 10) / 10:
            decision = 'D1' if lo > threshold + 1e-08 else 'D0' if hi < threshold - 1e-08 else 'UNKNOWN'
            decisions.append({'threshold_scenario': float(threshold), 'robust_decision': decision, 'product_decision': 'D1' if prod > threshold else 'D0'})
        rows.append({'C': c, 'P': p, 'count_bounds_D0_D1': bounds, 'posterior_interval': [lo, hi], 'width': hi - lo, 'product_posterior': prod, 'oracle': oracle, 'decisions': decisions})
    disagreement = max((abs(row['posterior_interval'][i] - row['oracle'][i]['probability']) for row in rows for i in range(2)))
    residual = max((r['max_count_residual'] for row in rows for r in row['oracle']))
    valid = disagreement <= 1e-08 and residual <= 1e-08
    result = {'round': 'R1v2', 'claim_type': 'capability', 'atoms': atoms, 'cases': rows, 'gates': {'validity': 'PASS' if valid else 'FAIL', 'all_width_le_0_10': 'PASS' if all((row['width'] <= 0.1 for row in rows)) else 'FAIL'}, 'lp_agreement_max': disagreement, 'lp_count_residual_max': residual, 'algorithm_comparison': 'TIE with standard equally informed linear-fractional LP. No algorithm superiority.', 'cost': {'total_seconds': time.perf_counter() - start, 'analytic_query_seconds': query_times, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'prep_reasoning': 'UNKNOWN', 'external_acquisition': 'Published source reused; prospective paired measurement NOT_RUN', 'fit_seconds': 0, 'gpu_seconds': 0}}
    dump('raw/R1_RESULTS.json', result)
    dump('CURRENT_WORK_STATE.json', {'tag': 'XBREAK-hunt-4', 'stage': 'R1_COMPLETE', 'latest_gate': result['gates'], 'next_operation': 'Freeze R2: acquire joint event totals and rounded joint reference counts from Table2; retain all compatible integer tables.'})
    print(json.dumps({k: v for (k, v) in result.items() if k not in ['atoms', 'cases']}, indent=2))
if __name__ == '__main__':
    main()
