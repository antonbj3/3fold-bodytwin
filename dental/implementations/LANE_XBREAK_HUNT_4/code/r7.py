import json, time, resource
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from r1 import assert_freeze, dump

def group_coefficients(counts, weights):
    cats = ['normal', 'widening', 'lesion', 'missing']
    labels = [c for c in cats for _ in range(counts[c])]
    return (np.array([weights.get(c, 0.0) for c in labels]), labels)

def support(a, mu, s):
    n = len(a)
    ac = a - a.mean()
    norm = np.linalg.norm(ac)
    radius = s * np.sqrt(n - 1)
    middle = mu * a.sum()
    delta = radius * norm
    if norm == 0:
        v = np.zeros(n)
        v[0] = 1
        v[1] = -1
        v /= np.linalg.norm(v)
    else:
        v = ac / norm
    witnesses = [np.full(n, mu) - radius * v, np.full(n, mu) + radius * v]
    return ([middle - delta, middle + delta], witnesses)

def original_timing_oracle(a, mu, s, maximize):
    n = len(a)
    variance = (n - 1) * s * s
    sign = -1 if maximize else 1

    def variance_slack(x):
        return variance - np.sum((x - mu) ** 2)
    constraints = [{'type': 'eq', 'fun': lambda x: np.sum(x) - n * mu, 'jac': lambda x: np.ones(n)}, {'type': 'ineq', 'fun': variance_slack, 'jac': lambda x: -2 * (x - mu)}]
    solve = minimize(lambda x: sign * (a @ x), np.full(n, mu), jac=lambda x: sign * a, bounds=[(0, None)] * n, constraints=constraints, method='SLSQP', options={'ftol': 1e-10, 'maxiter': 2000})
    return {'value': float(a @ solve.x), 'success': bool(solve.success), 'message': solve.message, 'equality_residual': float(abs(solve.x.sum() - n * mu)), 'variance_slack': float(variance_slack(solve.x)), 'minimum_time': float(min(solve.x)), 'iterations': int(solve.nit)}

def main():
    assert_freeze('R7')
    start = time.perf_counter()
    s = json.loads(Path('raw/R5_TRIAL_INPUTS.json').read_text())
    groups = ['P', 'fallback_RCT']
    counts = s['CBCT_completed_counts']
    mom = s['time_first_visit']
    cats = ['normal', 'widening', 'lesion']
    queries = []
    for c in cats:
        den = counts['P'][c] + counts['fallback_RCT'][c]
        queries.append((c, {c: 1 / den}))
    queries.append(('normal_minus_lesion', {'normal': 1 / (34 + 5), 'lesion': -1 / (10 + 9)}))
    rows = []
    worlds = {}
    for (name, weights) in queries:
        intervals = []
        checks = []
        worlds[name] = {'type': 'CONSTRUCTED_COMPATIBLE_FIXTURE_NOT_MEASURED_PATIENTS', 'groups': {}}
        for g in groups:
            (a, labels) = group_coefficients(counts[g], weights)
            mu = mom[g]['mean']
            sd = mom[g]['sd']
            (interval, witness) = support(a, mu, sd)
            oracles = [original_timing_oracle(a, mu, sd, k) for k in [False, True]]
            err = max((abs(interval[i] - oracles[i]['value']) for i in range(2)))
            res = []
            for x in witness:
                res.append({'mean_residual': float(abs(x.mean() - mu)), 'sample_SD_residual': float(abs(x.std(ddof=1) - sd)), 'minimum_min': float(min(x))})
            intervals.append(interval)
            checks.append({'group': g, 'interval': interval, 'oracles': oracles, 'oracle_max_error': err, 'witness_checks': res})
            worlds[name]['groups'][g] = {'CBCT_labels': labels, 'timing_world_lower': witness[0].tolist(), 'timing_world_upper': witness[1].tolist(), 'mean': mu, 'sample_sd': sd, 'n': len(a)}
        total = np.array(intervals).sum(axis=0).tolist()
        rows.append({'query': name, 'interval_min': total, 'width_min': total[1] - total[0], 'groups': checks})
    category_rows = rows[:3]
    contrast = rows[3]['interval_min']
    errors = [x['oracle_max_error'] for r in rows for x in r['groups']]
    moment_error = max((max(q['mean_residual'], q['sample_SD_residual']) for r in rows for g in r['groups'] for q in g['witness_checks']))
    minimum = min((q['minimum_min'] for r in rows for g in r['groups'] for q in g['witness_checks']))
    interval_gate = all((r['interval_min'][1] <= 100 and r['width_min'] <= 20 for r in category_rows))
    valid = max(errors) <= 1e-05 and moment_error <= 1e-08 and (minimum >= 0)
    result = {'round': 'R7', 'claim_type': 'capability', 'query_intervals': rows, 'contrast_order': 'normal_faster' if contrast[1] < 0 else 'normal_slower' if contrast[0] > 0 else 'UNKNOWN: both signs attained', 'gates': {'attained_compatible_worlds': 'PASS' if moment_error <= 1e-08 and minimum >= 0 else 'FAIL', 'independent_control': 'PASS' if max(errors) <= 1e-05 else 'FAIL', '100min_budget_and_20min_width': 'PASS' if interval_gate else 'FAIL', 'category_order_identified': 'PASS' if contrast[0] * contrast[1] > 0 else 'FAIL'}, 'max_oracle_error_min': max(errors), 'max_mean_SD_residual_min': moment_error, 'minimum_witness_time_min': minimum, 'full_rounding_population_patient_prediction': 'UNKNOWN: nominal compatible worlds refute identification; not a sampling interval or per-patient guarantee.', 'cost': {'seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'SLSQP_fits': 16, 'prior_all_rounds_inherited': True, 'prospective_timing_panel': 'NOT_RUN'}}
    dump('raw/R7_TIMING_WORLDS.json', worlds)
    dump('raw/R7_RESULTS.json', result)
    dump('CURRENT_WORK_STATE.json', {'tag': 'XBREAK-hunt-4', 'stage': 'R7_COMPLETE', 'latest_gate': result['gates'], 'next_operation': 'Acquire/import joint assignment + preop observations + guard/fallback + complete timing panel. Build consumer that refuses fabricated conditional-time independence.'})
    print(json.dumps({'gates': result['gates'], 'queries': [{k: r[k] for k in ['query', 'interval_min', 'width_min']} for r in rows], 'contrast_order': result['contrast_order'], 'oracle_error': max(errors), 'moments': moment_error, 'cost': result['cost']}, indent=2))
if __name__ == '__main__':
    main()
