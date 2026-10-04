from common import *
from fractions import Fraction as Q
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
import collections

def rational_bound(A, b, prep, film):
    candidates = []
    for j in range(A.shape[0]):
        row = A.getrow(j)
        co = [Q(float(x)) for x in row.data]
        if not co or any((x < 0 for x in co)):
            continue
        s = sum(co)
        rhs = Q(float(b[j])) - sum((x * (Q(float(prep[k])) + Q(float(film))) for (x, k) in zip(co, row.indices)))
        candidates.append((rhs / s, j))
    return min(candidates) if candidates else (None, None)

def verify_witness(coeff, lower, b):
    return bool(len(coeff) == len(lower) and all((Q(x) >= 0 for x in coeff)) and (sum((Q(a) * Q(x) for (a, x) in zip(coeff, lower))) > Q(b)))

def run():
    from item_analysis import load_rows
    rows = [r for r in load_rows() if r['split'] == 'test']
    names = read(ROOT / 'PREREG_R1B_BOUNDARY.json')['frontier_participants']
    by = collections.defaultdict(dict)
    for r in rows:
        by[r['task_id']][r['participant']] = r
    families = []
    for fam in sorted({r['family'] for r in rows}):
        xx = [v for v in by.values() if v[names[0]]['family'] == fam]
        m = np.array([[t[n]['L1'] == 'PASS' for n in names] for t in xx])
        families.append(dict(family=fam, items=len(xx), discriminating_items=int(np.sum(np.ptp(m.astype(int), axis=1) > 0)), participant_rates=dict(zip(names, m.mean(0))), participant_rate_range=float(np.ptp(m.mean(0))), resolution='POPULATION'))
    out = []
    excluded = []
    public = V4 / 'payload/public'
    cases = [r for r in read(V4 / 'payload/private/COHORT.json') if r['split'] == 'test'][:8]
    for case in cases:
        t = next((t for t in read(public / 'tasks' / f"{case['case_key']}.json") if t['family'] == 'molar_crown' and t['level'] == 'easy'))
        if t['status'] != 'READY':
            excluded.append(dict(task_id=t['task_id'], reason=t['site_error']))
            continue
        s = npz(public / t['geometry_file'])
        A = csr_matrix((s['a_data'], s['a_indices'], s['a_indptr']), shape=(int(s['a_rows']), len(s['xy'])))
        b = s['b'] - t['requirements']['clearance_mm']
        prep = np.full(len(s['xy']), t['preparation_height_mm'])
        film = t['requirements']['film_min_mm']
        (wall, j) = rational_bound(A, b, prep, film)
        if wall is None or wall <= 0:
            excluded.append(dict(task_id=t['task_id'], reason='no positive admissible critical wall or no nonnegative obstacle rows'))
            continue
        for sign in [-1, 1]:
            w = wall + sign * Q('0.000001')
            lo = [Q(float(v)) + Q(float(film)) + w for v in prep]
            row = A.getrow(j)
            co = [str(Q(float(x))) for x in row.data]
            lower = [str(lo[k]) for k in row.indices]
            rhs = str(Q(float(b[j])))
            witness = verify_witness(co, lower, rhs)
            changed_b = str(sum((Q(a) * Q(x) for (a, x) in zip(co, lower))) + 1)
            fault_rejected = not verify_witness(co, lower, changed_b)
            lp = linprog(np.zeros(len(prep)), A_ub=A, b_ub=b, bounds=[(float(x), None) for x in lo], method='highs', options={'primal_feasibility_tolerance': 1e-09, 'dual_feasibility_tolerance': 1e-09, 'threads': 1})
            out.append(dict(task_id=t['task_id'], offset_sign=sign, critical_wall_mm=float(wall), critical_wall_exact=str(wall), wall_mm=float(w), exact_core_feasible=not witness, LP_success=bool(lp.success), LP_status=int(lp.status), witness_row=j, witness_coefficients=co, witness_lower=lower, witness_rhs=rhs, wrong_rhs_rejected=fault_rejected if witness else None, agreement=bool(lp.success) == (not witness), resolution='PER_POINT', scope='Exact rational claim for stored PL wall/film/obstacle core, not physical or complete restoration feasibility'))
    result = dict(claim_type='capability', frontier_families=families, rows=out, dropout={'candidate_cases': len(cases), 'excluded': excluded, 'fraction': len(excluded) / len(cases)}, gate=all((r['agreement'] for r in out)) and all((r['wrong_rhs_rejected'] for r in out if not r['exact_core_feasible'])))
    dump(ROOT / 'raw/R1B_BOUNDARY.json', result)
    return result
if __name__ == '__main__':
    run()
