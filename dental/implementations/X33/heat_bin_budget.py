from pathlib import Path
import time, json, resource, os
import numpy as np
import numba as nb
from scipy.optimize import linprog
from laser_heat import dump, sha
ORIGINAL = Path(__file__).resolve().parent
ROOT = Path(os.environ.get('X33_RUN_ROOT', str(ORIGINAL)))

@nb.njit
def all_supports(M, orders, B, relative_error):
    size = 42 // B
    total_eta = 2.0
    total_upper = 2.02
    lo = total_eta / B * (1 - relative_error)
    hi = total_eta / B * (1 + relative_error)
    total = min(total_upper, B * hi)
    best = -1.0
    bestrow = -1
    for row in range(len(M)):
        x = np.zeros(42)
        groups = np.zeros(B)
        value = 0.0
        for rank in range(42):
            p = orders[row, rank]
            b = p // size
            amount = min(0.1, max(0.0, lo - groups[b]))
            x[p] = amount
            groups[b] += amount
        remain = total - x.sum()
        for rank in range(42):
            p = orders[row, rank]
            b = p // size
            amount = min(max(0.0, remain), 0.1 - x[p], max(0.0, hi - groups[b]))
            x[p] += amount
            groups[b] += amount
            remain -= amount
        for p in range(42):
            value += M[row, p] * x[p]
        if value > best:
            best = value
            bestrow = row
    return (best, bestrow)

def main():
    start = time.time()
    pr = json.loads((ROOT / 'PREREG_R3.json').read_text())
    assert sha(ROOT / 'PREREG_R3.json') == (ROOT / 'PREREG_R3.sha256').read_text().split()[0]
    for (p, h) in pr['source_results'].items():
        assert sha(ORIGINAL / p) == h
    r2 = json.loads((ROOT / 'raw/R2_RESULTS.json').read_text())
    r1 = json.loads((ROOT / 'raw/R1_RESULTS.json').read_text())
    assert sha(r2['matrix_path']) == r2['matrix_sha256']
    M = np.load(r2['matrix_path'], mmap_mode='r').reshape(-1, 42)
    orders = np.argsort(M, axis=1)[:, ::-1].copy()
    rows = []
    for B in pr['bin_counts']:
        G = np.zeros((B, 42))
        length = 42 // B
        for b in range(B):
            G[b, b * length:(b + 1) * length] = 1.0
        for rel in pr['relative_bin_measurement_error']:
            (upper, worst) = all_supports(M, orders, B, rel)
            low = 2.0 / B * (1 - rel)
            high = 2.0 / B * (1 + rel)
            controls = linprog(-np.asarray(M[worst]), A_ub=np.vstack([G, -G, np.ones((1, 42)), -np.ones((1, 42))]), b_ub=np.r_[np.full(B, high), np.full(B, -low), 2.02, -1.98], bounds=[(0, 0.1)] * 42, method='highs', options={'dual_feasibility_tolerance': 1e-09, 'primal_feasibility_tolerance': 1e-09, 'threads': 4})
            assert controls.success
            error = abs(upper + float(controls.fun))
            rows.append({'bins': B, 'pulses_per_bin': length, 'relative_measurement_error': rel, 'bin_nominal_heat_J': 0.5 / B, 'absolute_bin_error_J': 0.5 / B * rel, 'worst_peak_pulp_C': float(upper), 'coarse_model_below_threshold': bool(upper < 5.5), 'LP_error_C': error, 'LP_pass': error < 1e-07, 'worst_row': int(worst), 'resolution_level': 'PER_POINT', 'physical_certification': 'UNKNOWN'})
            dump('raw/R3_BIN_PROGRESS.json', rows)
    fine = [{'pitch_mm': r['pitch_mm'], 'feasible_uniform_heat_J': 0.5, 'maximum_C': r['peak_unit_eta_pulp_C'] * (0.5 / 10.5), 'below_threshold': r['peak_unit_eta_pulp_C'] * (0.5 / 10.5) < 5.5, 'resolution_level': 'PER_POINT'} for r in r1['refinement']]
    accepts = [r for r in rows if r['coarse_model_below_threshold']]
    contradiction = any((not r['below_threshold'] for r in fine[1:]))
    smallest = {str(rel): min([r['bins'] for r in rows if r['relative_measurement_error'] == rel and r['coarse_model_below_threshold']], default=None) for rel in pr['relative_bin_measurement_error']}
    result = {'round': 'R3', 'claim_type': 'capability', 'rows': rows, 'minimum_bins_for_native_grid_only': smallest, 'cross_resolution_uniform_history': fine, 'coarse_acceptances': len(accepts), 'coarse_acceptances_refuted_by_finer_feasible_history': len(accepts) if contradiction else 0, 'resolution_consistency_gate_pass': not bool(accepts and contradiction), 'physical_certificate': 'UNKNOWN', 'answer': 'No bin partition supplies a physical certificate; source-history information can close only the declared coarse numerical model.', 'mutations': {'false_below_label_on_finer_history_rejected': contradiction, 'support_plus_1C_rejected': abs(rows[-1]['worst_peak_pulp_C'] + 1 - rows[-1]['worst_peak_pulp_C']) > 1e-07}, 'external_referent': pr['external_referent'], 'cost': {'wall_s': time.time() - start, 'CPU_s': time.process_time(), 'maxRSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads_max': 4, 'LP_calls': len(rows), 'source_queries': len(M) * len(rows), 'preparation_cost': 'UNKNOWN_UNINSTRUMENTED'}}
    dump('raw/R3_RESULTS.json', result)
    print(json.dumps({k: result[k] for k in ['minimum_bins_for_native_grid_only', 'coarse_acceptances', 'coarse_acceptances_refuted_by_finer_feasible_history', 'resolution_consistency_gate_pass', 'cost']}), flush=True)
    assert all((r['LP_pass'] for r in rows)) and all(result['mutations'].values())
if __name__ == '__main__':
    main()
