from common import *
from itertools import combinations
import time, copy
import numpy as np
from scipy.optimize import linprog
from field_engine.experimental._farkas_check import check_dual, check_primal

def nonempty_all(rows):
    assert rows and all((r['controls'] for r in rows)), 'EMPTY_CONTROL_SELECTION'
    return all((c['pass_gate'] for r in rows for c in r['controls']))

def linear(A, b):
    n = len(b)
    a = [list(r) + [v] for (r, v) in zip(A, b)]
    for j in range(n):
        i = next((i for i in range(j, n) if a[i][j]), None)
        if i is None:
            raise ValueError('singular')
        (a[i], a[j]) = (a[j], a[i])
        v = a[j][j]
        a[j] = [x / v for x in a[j]]
        for i in range(n):
            if i != j:
                v = a[i][j]
                a[i] = [x - v * y for (x, y) in zip(a[i], a[j])]
    return [r[-1] for r in a]

def exact_lp(A, b, c):
    feasible = []
    for ids in combinations(range(len(c)), len(b)):
        B = [[r[j] for j in ids] for r in A]
        try:
            v = linear(B, b)
        except ValueError:
            continue
        if min(v) < 0:
            continue
        x = [Q(0)] * len(c)
        for (i, w) in zip(ids, v):
            x[i] = w
        z = sum((v * w for (v, w) in zip(c, x)))
        y = linear(list(zip(*B)), [c[j] for j in ids])
        feasible.append((z, x, y, ids))
    assert feasible, 'NO_EXACT_FEASIBLE_BASIS'
    z = min((r[0] for r in feasible))
    (z, x, y, ids) = next((r for r in feasible if r[0] == z and all((sum((A[i][j] * r[2][i] for i in range(len(b)))) <= c[j] for j in range(len(c))))))
    return (z, x, y)

def constraints(U, L):
    A = [[Q(1)] * 3 + [Q(0)] * 3, [Q(0)] * 3 + [Q(1)] * 3, [p[0] for p in U] + [-p[0] for p in L], [p[1] for p in U] + [-p[1] for p in L]]
    return (A, [Q(1), Q(1), Q(0), Q(0)], [p[2] for p in U] + [-p[2] for p in L])

def cert(A, b, c, z, x, y, eps=Q(1, 1000)):
    n = len(c)
    AA = A + [[-v for v in r] for r in A] + [[-Q(i == j) for j in range(n)] for i in range(n)] + [c]
    bb = b + [-v for v in b] + [Q(0)] * n + [z - eps]
    slack = [c[j] - sum((A[i][j] * y[i] for i in range(len(b)))) for j in range(n)]
    weights = [max(-v, 0) for v in y] + [max(v, 0) for v in y] + slack + [Q(1)]
    assert check_primal('orthant', A, b, x).status == 'JA'
    assert check_dual('inequality', AA, bb, weights).status == 'NEJ'
    return dict(A=AA, b=bb, weights=weights, separation=-eps, primal=x, dual_equality=y, optimum_mm=z)

def main():
    t = time.perf_counter()
    state('X54_RUNNING', 'Frozen exact witness gate', 'Certify all stored winning pairs')
    faults = {}
    for (name, rows) in [('no_rows', []), ('no_controls', [dict(controls=[])])]:
        try:
            nonempty_all(rows)
            faults[name] = False
        except AssertionError:
            faults[name] = True
    faults['legacy_all_empty_returns_true'] = all([])
    rows = []
    full = []
    for case in [1, 2, 3]:
        p = BASE / 'LANE_X54_DEFORMABLE_CONTACT/raw' / f'geometry_{case:03d}_continuous.json'
        g = read(p, exact=True)
        controls = []
        for patch in g['patches']:
            U = [[Q(v) for v in p] for p in patch['continuous_upper_subtriangle']]
            L = [[Q(v) for v in p] for p in patch['continuous_lower_triangle']]
            (A, b, c) = constraints(U, L)
            (z, x, y) = exact_lp(A, b, c)
            w = cert(A, b, c, z, x, y)
            lp = linprog(np.array(c, float), A_eq=np.array(A, float), b_eq=np.array(b, float), bounds=(0, None), method='highs', options={'threads': 4})
            err = abs(float(z) - lp.fun)
            old_err = abs(float(z) - float(patch['vertical_gap_mm']))
            bad = w['weights'].copy()
            bad[-1] = 0
            reject = check_dual('inequality', w['A'], w['b'], bad).status != 'NEJ'
            bb = w['b'].copy()
            bb[-1] = z + Q(1, 1000)
            reject2 = check_dual('inequality', w['A'], bb, w['weights']).status != 'NEJ'
            row = dict(case=case, upper_fdi=patch['upper_fdi'], lower_fdi=patch['lower_fdi'], region=patch['region'], exact_optimum_mm=z, float_LP_error_mm=err, old_winner_error_mm=old_err, pass_gate=bool(lp.success and err <= 1e-08 and (old_err <= 1e-08) and reject and reject2), bad_dual_rejected=reject, false_lower_bound_rejected=reject2, exact_primal_dual_gap_mm=0)
            controls.append(row)
            full.append(dict(binding=dict(source_path=p, source_sha256=sha(p), upper_face=patch['upper_source_face'], lower_face=patch['lower_source_face'], region=patch['region']), triangles=dict(upper=U, lower=L), certificate=w))
        rows.append(dict(case=case, controls=controls))
    ok = nonempty_all(rows)
    assert ok and all(faults.values())
    write(ROOT / 'raw/X54_CERTIFICATES.json', full)
    out = dict(claim_type='capability', rows=rows, all_pass=ok, control_count=sum((len(r['controls']) for r in rows)), injections=faults, changed_conclusion=False, changed_certification=True, external_referent=dict(kind='published_code', locator='Field e590aa3:src/field_engine/experimental/_farkas_check.py', compared_quantity='exact strict rational separation and primal feasibility of stored winning triangle pair', refutes_us=True), resolution='PER_SURFACE_REGION', rigorous_enclosure='EXACT_STORED_WINNING_PAIR_ONLY; global winning-pair coverage, clipped-triangle upstream rounding, sampled area, normal projection, NCP and parameter boxes remain MISSING', round_scope={'R1': 'sampled minima remain sampling', 'R2': 'sampled minima remain sampling', 'R3': 'exact stored-pair minima; count emitted', 'R4': 'reuses R3 geometry; same pair certificates apply, no new force certification'}, cost_wall_s=time.perf_counter() - t)
    write(ROOT / 'raw/X54.json', out)
    state('X54_DONE', str(out['control_count']) + ' exact pair witnesses PASS', 'Correct X2 family size')
    print('X54', out['control_count'], ok)
if __name__ == '__main__':
    main()
