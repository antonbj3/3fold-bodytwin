"""One-sided gap reduction, then joint force-simplex reconstruction."""
import json, time, resource, warnings
import numpy as np
from scipy.optimize import linprog, OptimizeWarning
from run_r1 import ROOT, sha, write

def caps(baseline, minus, h, T, work_error):
    upper = np.minimum(T, np.maximum(0.0, (baseline - np.array(minus) + 2 * work_error) / h))
    lower = np.maximum(0.0, T - upper.sum() + upper)
    return (lower, upper)

def main():
    start = time.perf_counter()
    cpu = time.process_time()
    p = ROOT / 'PREREG_R3.json'
    cfg = json.loads(p.read_text())
    frozen = sha(p)
    write('raw/PREREG_R3_FROZEN.json', {'sha256': frozen, 'file': str(p)})
    worlds = json.loads((ROOT / 'raw/R1_WORLDS.json').read_text())
    obs = json.loads((ROOT / 'raw/R1_OBSERVATIONS.json').read_text())
    targets = list(range(11, 19)) + list(range(21, 29))
    T = 205
    records = []
    maxdiff = 0.0
    lp_calls = 0
    lp_s = 0.0
    warnings.filterwarnings('ignore', category=OptimizeWarning, message='Unrecognized options detected.*')
    for w in worlds:
        ww = [o for o in obs if o['case'] == w['case'] and o['scenario'] == w['scenario']]
        top = set(sorted(targets, key=lambda t: (-w['native_force_N'].get(str(t), 0), t))[:2])
        for h in cfg['h_mm']:
            minus = [next((o['minus_work_Nmm'] for o in ww if o['FDI'] == t and o['h_mm'] == h)) for t in targets]
            for eta in [0, 0.0001, 0.001]:
                (lo, hi) = caps(w['native_V_Nmm'], minus, h, T, T * eta)
                if hi.sum() < T - 1e-05:
                    raise ArithmeticError('Infeasible upper-force simplex')
                for (i, tooth) in enumerate(targets):
                    tick = time.perf_counter()
                    objective = np.eye(len(targets))[i]
                    low = linprog(objective, A_eq=np.ones((1, len(targets))), b_eq=[T], bounds=list(zip(np.zeros(len(targets)), hi)), method='highs', options={'threads': 4})
                    high = linprog(-objective, A_eq=np.ones((1, len(targets))), b_eq=[T], bounds=list(zip(np.zeros(len(targets)), hi)), method='highs', options={'threads': 4})
                    lp_s += time.perf_counter() - tick
                    lp_calls += 2
                    if not low.success or not high.success:
                        raise ArithmeticError('LP failure')
                    diff = max(abs(low.fun - lo[i]), abs(-high.fun - hi[i]))
                    maxdiff = max(maxdiff, diff)
                    q = w['native_force_N'].get(str(tooth), 0.0)
                    records.append({'case': w['case'], 'scenario': w['scenario'], 'FDI': tooth, 'h_mm': h, 'equivalent_displacement_error_mm': eta, 'facit_N': q, 'lower_N': float(lo[i]), 'upper_N': float(hi[i]), 'width_pp': float(100 * (hi[i] - lo[i]) / T), 'covered': bool(lo[i] - 1e-05 <= q <= hi[i] + 1e-05), 'top_two': tooth in top, 'LP_difference_N': float(diff)})
    summary = []
    for eta in [0, 0.0001, 0.001]:
        for h in cfg['h_mm']:
            rr = [r for r in records if r['equivalent_displacement_error_mm'] == eta and r['h_mm'] == h]
            top = [r for r in rr if r['top_two']]
            summary.append({'error_mm': eta, 'h_mm': h, 'n': len(rr), 'covered': sum((r['covered'] for r in rr)), 'fraction_all_le5pp': float(np.mean([r['width_pp'] <= 5 for r in rr])), 'fraction_top_two_le5pp': float(np.mean([r['width_pp'] <= 5 for r in top])), 'median_width_all_pp': float(np.median([r['width_pp'] for r in rr])), 'max_width_all_pp': max((r['width_pp'] for r in rr)), 'median_width_top_two_pp': float(np.median([r['width_pp'] for r in top]))})

    def gate(eta):
        return any((s['fraction_all_le5pp'] >= 0.8 and s['fraction_top_two_le5pp'] >= 0.8 for s in summary if s['error_mm'] == eta))
    out = {'round': 'R3', 'claim_type': cfg['claim_type'], 'status': 'PENDING_INDEPENDENT_REVIEW', 'prereg_sha256': frozen, 'coverage_gate': all((r['covered'] for r in records)), 'ideal_useful_gate': gate(0), 'noisy_useful_gate': gate(0.001), 'LP_parity_gate': bool(maxdiff <= 1e-07), 'max_LP_difference_N': float(maxdiff), 'summaries': summary, 'records': records, 'cost': {'wall_seconds': time.perf_counter() - start, 'CPU_seconds': time.process_time() - cpu, 'max_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'LP_calls': lp_calls, 'LP_seconds': lp_s, 'new_equilibrium_calls': 0, 'R1_upstream_QPs': 873, 'R1_upstream_wall_seconds': 4.156158004887402, 'physical_acquisition': 'NOT_RUN; 441 work reads for all h from scratch'}, 'physical_validation': 'UNKNOWN; ideal work revisited to isolate the sign-access intervention; R2 integration failure remains open'}
    write('rounds/R3/results.json', out)
    assert sha(p) == frozen
    print(json.dumps({k: out[k] for k in ['coverage_gate', 'ideal_useful_gate', 'noisy_useful_gate', 'LP_parity_gate', 'summaries', 'cost']}, indent=2))
if __name__ == '__main__':
    main()
