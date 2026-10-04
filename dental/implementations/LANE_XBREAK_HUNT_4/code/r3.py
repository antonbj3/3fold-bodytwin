import itertools, json, time, resource, math
from pathlib import Path
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from r1 import assert_freeze, constraints, dump
from r2 import joint_rows

def possible_rows(row_total, col):
    for a in range(min(row_total, col[0]) + 1):
        for b in range(min(row_total - a, col[1]) + 1):
            c = row_total - a - b
            if 0 <= c <= col[2]:
                yield np.array([a, b, c], int)

def table_candidates(source):
    result = []
    for k10 in range(116 + 1):
        if not 81.5 <= 100 * k10 / 116 < 82.5:
            continue
        for k11 in range(264 + 1):
            if not 79.5 <= 100 * k11 / 264 < 80.5:
                continue
            fixed = [np.array([116 - k10, 264 - k11, 0]), np.array([k10, k11, 0])]
            choices = []
            for d in range(2):
                if sum(fixed[d]) != source['cold_margins'][d][1]:
                    break
                remaining = np.array(source['percussion_margins'][d]) - fixed[d]
                if min(remaining) < 0:
                    break
                group = []
                for miss in possible_rows(source['cold_margins'][d][2], remaining):
                    first = remaining - miss
                    if sum(first) == source['cold_margins'][d][0]:
                        group.append(np.vstack([first, fixed[d], miss]))
                choices.append(group)
            if len(choices) == 2:
                result.extend((np.stack(pair) for pair in itertools.product(*choices)))
    return result

def milp_constraints(source, atoms, A, b):
    (extra, totals) = joint_rows(atoms)
    a = np.vstack([A, extra])
    lo = list(b) + list(totals)
    hi = lo.copy()
    for (p, n, r) in [(0, 116, 82), (1, 264, 80)]:
        mask = np.array([int(d == 1 and c == 1 and (k == p)) for (d, c, k) in atoms], float)
        a = np.vstack([a, mask])
        lo.append(math.ceil(n * (r - 0.5) / 100))
        hi.append(math.ceil(n * (r + 0.5) / 100) - 1)
    return LinearConstraint(a, np.array(lo), np.array(hi))

def integer_ratio_oracle(atoms, con, c, p, maximize):
    e = np.array([int(j == c and k == p) for (d, j, k) in atoms], float)
    n = np.array([int(d == 1 and j == c and (k == p)) for (d, j, k) in atoms], float)
    ratio = 0.5
    iterations = 0
    while iterations < 30:
        objective = (n - ratio * e) * (-1 if maximize else 1)
        solve = milp(objective, integrality=np.ones(len(atoms)), bounds=Bounds(np.zeros(len(atoms)), np.full(len(atoms), 708)), constraints=con, options={'mip_rel_gap': 0})
        if not solve.success:
            raise RuntimeError(solve.message)
        den = e @ solve.x
        if den <= 0:
            raise RuntimeError('Absent query event')
        updated = float(n @ solve.x / den)
        iterations += 1
        if abs(updated - ratio) < 1e-12:
            return (updated, iterations, solve.x.tolist())
        ratio = updated
    raise RuntimeError('Dinkelbach did not converge')

def main():
    assert_freeze('R3')
    start = time.perf_counter()
    s = json.loads(Path('raw/R1_SOURCE_COUNTS.json').read_text())
    (atoms, A, b) = constraints(s)
    tables = table_candidates(s)
    assert tables
    con = milp_constraints(s, atoms, A, b)
    cases = []
    for (c, p) in itertools.product(range(2), repeat=2):
        vals = np.array([q[1, c, p] / q[:, c, p].sum() for q in tables])
        oracle = [integer_ratio_oracle(atoms, con, c, p, m) for m in [False, True]]
        interval = [float(min(vals)), float(max(vals))]
        cases.append({'C': c, 'P': p, 'posterior_interval': interval, 'width': interval[1] - interval[0], 'count_pairs_D0_D1': sorted(set((tuple((int(x) for x in q[:, c, p])) for q in tables))), 'MILP_oracle': oracle, 'decisions': [{'threshold_scenario': k / 10, 'decision': 'D1' if interval[0] > k / 10 + 1e-08 else 'D0' if interval[1] < k / 10 - 1e-08 else 'UNKNOWN'} for k in range(1, 10)]})
    errors = [abs(c['posterior_interval'][i] - c['MILP_oracle'][i][0]) for c in cases for i in range(2)]
    report = json.loads(Path('raw/R3_PAIRED_REPORTS.json').read_text())
    crosschecks = []
    for row in report['joint']:
        (c, p) = (row['C'], row['P'])
        pos = int(tables[0][1, c, p])
        neg = int(tables[0][0, c, p])
        nd = sum(s['cold_margins'][1][:2])
        nv = sum(s['cold_margins'][0][:2])
        pred = {'PPV_percent': 100 * pos / (pos + neg), 'SN_percent': 100 * pos / nd, 'SP_percent': 100 * (nv - neg) / nv, 'TA_percent': 100 * (pos + nv - neg) / (nd + nv), 'NPV_percent': 100 * (nv - neg) / (nv - neg + nd - pos)}
        for (key, value) in pred.items():
            crosschecks.append({'C': c, 'P': p, 'metric': key, 'computed': value, 'source_display': row[key], 'within_reporting_band': row[key] - 0.5 <= value < row[key] + 0.5})
    results = {'round': 'R3', 'claim_type': 'information_link', 'admissible_integer_tables': len(tables), 'cases': cases, 'source_crosschecks': crosschecks, 'MILP_agreement_max': max(errors), 'gates': {'all_width_le_0_10': 'PASS' if all((c['width'] <= 0.1 for c in cases)) else 'FAIL', 'independent_integer_control': 'PASS' if max(errors) <= 1e-08 else 'FAIL', 'source_faithfulness': 'PASS' if all((r['within_reporting_band'] for r in crosschecks)) else 'FAIL'}, 'cost': {'total_seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'prior_cost_inherited': ['R1', 'R2'], 'prospective_paired_acquisition': 'NOT_RUN', 'population_validation': 'UNKNOWN'}}
    dump('raw/R3_TABLES.json', [q.tolist() for q in tables])
    dump('raw/R3_RESULTS.json', results)
    dump('CURRENT_WORK_STATE.json', {'tag': 'XBREAK-hunt-4', 'stage': 'R3_COMPLETE', 'latest_gate': results['gates'], 'next_operation': 'Freeze R4: stop treating empirical counts as known population probabilities; exact simultaneous sampling intervals and unselected population UNKNOWN.'})
    print(json.dumps({'tables': len(tables), 'gates': results['gates'], 'cases': [{k: r[k] for k in ['C', 'P', 'posterior_interval', 'width']} for r in cases], 'cost': results['cost']}, indent=2))
if __name__ == '__main__':
    main()
