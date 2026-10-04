"""Query-conditioned reaction projection; rate law is deliberately hidden."""
import json, time, resource, warnings
import numpy as np
from scipy.optimize import linprog, OptimizeWarning
from run_r1 import ROOT, sha, write

def project(y0, yv, h, d, de, fe, matched_velocity=False, matched_history=False, position_closure=False):
    if not (matched_velocity and matched_history and position_closure):
        return {'status': 'UNKNOWN', 'force_interval_N': None, 'reason': 'Unverified query-state/velocity/history or position closure'}
    if h <= 0 or d - de <= h or y0 - fe < 0:
        return {'status': 'UNKNOWN', 'force_interval_N': None, 'reason': 'Unsupported projection regime'}
    lo = ((d - de) / h - 1) * (y0 - fe) + (yv - fe)
    hi = ((d + de) / h - 1) * (y0 + fe) + (yv + fe)
    return {'status': 'CONDITIONAL_POINT_BOUND', 'force_interval_N': [lo, hi]}

def main():
    start = time.perf_counter()
    cpu = time.process_time()
    p = ROOT / 'PREREG_R6.json'
    cfg = json.loads(p.read_text())
    frozen = sha(p)
    write('raw/PREREG_R6_FROZEN.json', {'sha256': frozen, 'file': str(p)})
    source = json.loads((ROOT / 'rounds/R5/results.json').read_text())
    records = []
    guards = []
    maxdiff = 0.0
    h = 0.005
    de = cfg['native_displacement_error_mm']
    warnings.filterwarnings('ignore', category=OptimizeWarning, message='Unrecognized options detected.*')
    for w in source['worlds']:
        d = w['native_d_mm']
        v = w['native_velocity_mm_per_s']
        for law in ['linear', 'nonlinear']:
            for i in range(2):
                k = w['k_N_per_mm'][i]
                eta = w['eta_Ns_per_mm'][i]
                R = eta * v + (0 if law == 'linear' else (1 if i == 0 else -1) * w['nonlinear_native_delta_N'])
                y0_true = k * h
                yv_true = k * h + R
                facit = w['native_force_N' if law == 'linear' else 'nonlinear_force_N'][i]
                for fe in cfg['reaction_error_bounds_N']:
                    y0 = y0_true + fe * 0.5 * np.sin(w['world'] + i + 1)
                    yv = yv_true + fe * 0.6 * np.cos(w['world'] - i + 1)
                    ans = project(y0, yv, h, d, de, fe, matched_velocity=True, matched_history=True, position_closure=True)
                    (lo, hi) = ans['force_interval_N']
                    objective = [(d - de) / h - 1, 1]
                    omin = linprog(objective, bounds=[(y0 - fe, y0 + fe), (yv - fe, yv + fe)], method='highs', options={'threads': 4})
                    objective = [-((d + de) / h - 1), -1]
                    omax = linprog(objective, bounds=[(y0 - fe, y0 + fe), (yv - fe, yv + fe)], method='highs', options={'threads': 4})
                    if not omin.success or not omax.success:
                        raise ArithmeticError('LP failure')
                    diff = max(abs(lo - omin.fun), abs(hi + omax.fun))
                    maxdiff = max(maxdiff, diff)
                    records.append({'world': w['world'], 'law': law, 'support': i, 'reaction_error_N': fe, 'y0_true_N': y0_true, 'yv_true_N': yv_true, 'y0_observed_N': y0, 'yv_observed_N': yv, 'native_d_mm': d, 'native_velocity_point_mm_per_s': v, 'lower_N': lo, 'upper_N': hi, 'facit_N': facit, 'width_pp': 100 * (hi - lo) / 410, 'covered': bool(lo - 1e-05 <= facit <= hi + 1e-05), 'LP_difference_N': float(diff)})
                    for flag in ['matched_velocity', 'matched_history', 'position_closure']:
                        kwargs = {'matched_velocity': True, 'matched_history': True, 'position_closure': True}
                        kwargs[flag] = False
                        guard = project(y0, yv, h, d, de, fe, **kwargs)
                        guards.append({'world': w['world'], 'law': law, 'support': i, 'failed_contract': flag, **guard})
                    guard = project(y0, yv, h, d, de, fe)
                    guards.append({'world': w['world'], 'law': law, 'support': i, 'failed_contract': 'all_contracts_unspecified', **guard})
    summary = []
    for fe in cfg['reaction_error_bounds_N']:
        rr = [r for r in records if r['reaction_error_N'] == fe]
        summary.append({'reaction_error_N': fe, 'n': len(rr), 'covered': sum((r['covered'] for r in rr)), 'n_width_le5pp': int(sum((r['width_pp'] <= 5 for r in rr))), 'max_width_pp': float(max((r['width_pp'] for r in rr)))})
    out = {'round': 'R6', 'claim_type': 'information_link', 'status': 'PENDING_INDEPENDENT_REVIEW', 'prereg_sha256': frozen, 'force_coverage_gate': all((r['covered'] for r in records)), 'useful_force_gate': all((r['width_pp'] <= 5 for r in records if r['reaction_error_N'] == 0.005)), 'same_information_control_gate': bool(maxdiff <= 1e-07), 'regime_guard_gate': all((g['status'] == 'UNKNOWN' and g['force_interval_N'] is None for g in guards)), 'max_LP_difference_N': float(maxdiff), 'summaries': summary, 'records': records, 'actual_guard_calls': guards, 'cost': {'wall_seconds': time.perf_counter() - start, 'CPU_seconds': time.process_time() - cpu, 'max_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'LP_calls': 2 * len(records), 'new_reaction_observation_pairs': len(records), 'dynamic_probe_travel_seconds': h / source['worlds'][0]['native_velocity_mm_per_s'], 'physical_acquisition': 'NOT_RUN; static setup, matching and history hold costs UNKNOWN', 'fit': 'NONE'}, 'physical_native_force_certificate': 'UNKNOWN; position/rate/history matching and metrology not measured', 'validity': 'PER_POINT exact velocity under separable memoryless force law with linear position closure. No validated full-jaw, nonlinear-position, friction or irreversible tissue transfer.'}
    write('rounds/R6/results.json', out)
    assert sha(p) == frozen
    print(json.dumps({k: out[k] for k in ['force_coverage_gate', 'useful_force_gate', 'same_information_control_gate', 'regime_guard_gate', 'summaries', 'cost']}, indent=2))
if __name__ == '__main__':
    main()
