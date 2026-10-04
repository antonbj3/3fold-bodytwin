import json, time, resource, csv, sys, warnings
from fractions import Fraction as Q
import numpy as np
import scipy
from scipy.optimize import linprog, OptimizeWarning
from geometry import R, write, sha, load, FDI14
LP_CALLS = 0

def target(t, tag):
    return t % 10 >= 4 if tag == 'posterior' else t // 10 in (1, 4)

def scalar_bounds(e, tag, q):
    out = [Q(0), Q(0)]
    for (g, mass) in [(False, 1 - q), (True, q)]:
        ds = [int(target(l, tag)) for (u, l) in e if target(u, tag) == g]
        if not ds:
            if mass:
                raise ValueError('Measured region has no declared edge')
            continue
        out[0] += mass * min(ds)
        out[1] += mass * max(ds)
    return out

def marginal_bounds(e, p, tag):
    out = [Q(0), Q(0)]
    for (u, mass) in p.items():
        ds = [int(target(l, tag)) for (uu, l) in e if u == uu]
        if not ds:
            if mass:
                raise ValueError('Measured upper tooth has no declared edge')
            continue
        out[0] += mass * min(ds)
        out[1] += mass * max(ds)
    return out

def lp_check(e, tag, bounds, q=None, p=None):
    global LP_CALLS
    obj = np.array([int(target(l, tag)) for (u, l) in e], float)
    if p is None:
        A = np.array([[1.0] * len(e), [int(target(u, tag)) for (u, l) in e]], float)
        b = np.array([1, float(q)])
    else:
        A = np.array([[int(uu == u) for (uu, l) in e] for u in p], float)
        b = np.array(list(map(float, p.values())))
    err = residual = 0.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', OptimizeWarning)
        for (sign, expect) in [(1, bounds[0]), (-1, bounds[1])]:
            r = linprog(sign * obj, A_eq=A, b_eq=b, bounds=(0, None), method='highs', options={'threads': 4})
            LP_CALLS += 1
            if not r.success:
                raise ArithmeticError(r.message)
            err = max(err, abs(obj @ r.x - float(expect)))
            residual = max(residual, float(abs(A @ r.x - b).max()), float(max(0, -r.x.min())))
    return (err, residual)

def components(e):
    n = {}
    for (u, l) in e:
        n.setdefault(u, set()).add(l)
        n.setdefault(l, set()).add(u)
    seen = set()
    c = 0
    for t in n:
        if t in seen:
            continue
        c += 1
        todo = [t]
        seen.add(t)
        while todo:
            z = todo.pop()
            for w in n[z] - seen:
                seen.add(w)
                todo.append(w)
    return c

def marginal(e, f, index):
    return {t: sum((v for (edge, v) in zip(e, f) if edge[index] == t)) for t in sorted({x[index] for x in e})}

def port(e, f, tag):
    upper = sum((v for ((u, l), v) in zip(e, f) if target(u, tag)))
    lower = sum((v for ((u, l), v) in zip(e, f) if target(l, tag)))
    delta = sum((v * (int(target(l, tag)) - int(target(u, tag))) for ((u, l), v) in zip(e, f)))
    return (upper, lower, delta)

def main():
    start = time.perf_counter()
    assert sha(R / 'PREREG_R3.json') == (R / 'PREREG_R3.json.sha256').read_text().split()[0]
    (manifest, cases, sources) = load()
    rows = []
    maxerr = 0.0
    maxres = 0.0
    witness = None
    rankerrs = []
    drop = []
    for c in cases:
        e = [edge for edge in c['edges'] if edge[0] % 10 != 8 and edge[1] % 10 != 8]
        area = [c['area'][edge] for edge in e]
        total = sum(area, Q(0))
        p = marginal(e, [v / total for v in area], 0) if total else None
        cc = components(e)
        ul = sorted({u for (u, l) in e})
        ll = sorted({l for (u, l) in e})
        A = np.array([[int(t in edge) for edge in e] for t in ul + ll], float)
        U = A[:len(ul)]
        rankdiff = int(np.linalg.matrix_rank(A) - np.linalg.matrix_rank(U))
        dim = len(ll) - cc
        rankerrs.append(rankdiff - dim)
        for (tag, q) in [('posterior', Q(824, 1000)), ('R', Q(496, 1000))]:
            try:
                b = scalar_bounds(e, tag, q)
                (err, res) = lp_check(e, tag, b, q=q)
                maxerr = max(maxerr, err)
                maxres = max(maxres, res)
                rows.append(dict(case=c['case'], axis=tag, tier='one_upper_region', lower=float(b[0]), upper=float(b[1]), width=float(b[1] - b[0]), input_kind='published_median_as_unpaired_scenario', additional_lower_measurements_full_vector=dim, components=cc))
            except ValueError as er:
                drop.append(dict(case=c['case'], axis=tag, tier='one_upper_region', reason=str(er)))
            if p is not None:
                b = marginal_bounds(e, p, tag)
                (err, res) = lp_check(e, tag, b, p=p)
                maxerr = max(maxerr, err)
                maxres = max(maxres, res)
                rows.append(dict(case=c['case'], axis=tag, tier='full_upper_tooth_marginal', lower=float(b[0]), upper=float(b[1]), width=float(b[1] - b[0]), input_kind='frozen_area_replay_fixture_NOT_MEASUREMENT', additional_lower_measurements_full_vector=dim, components=cc))
            if witness is None:
                for u in ul:
                    a = [i for (i, (uu, l)) in enumerate(e) if uu == u and target(l, tag)]
                    b = [i for (i, (uu, l)) in enumerate(e) if uu == u and (not target(l, tag))]
                    if not a or not b:
                        continue
                    f = [Q(1, len(e))] * len(e)
                    g = f.copy()
                    delta = Q(1, 2 * len(e))
                    g[a[0]] += delta
                    g[b[0]] -= delta
                    (pu, pv) = (marginal(e, f, 0), marginal(e, g, 0))
                    assert pu == pv
                    (aa, bb) = (port(e, f, tag), port(e, g, tag))
                    assert aa[1] == aa[0] + aa[2] and bb[1] == bb[0] + bb[2]
                    witness = dict(case=c['case'], axis=tag, upper_tooth=u, edges=e, force1=list(map(str, f)), force2=list(map(str, g)), upper_marginal_identical=True, identity_error_share=0.0, lower_region_difference_share=float(abs(aa[1] - bb[1])), net_crossing1=str(aa[2]), net_crossing2=str(bb[2]), extended_reconstruction_error_share=0.0, physical_realizability='UNKNOWN: conserved edge-flow fixture only, not a contact-KKT world', minimal_extension='One measured signed net cross-boundary force scalar per regional query at the same MIP time; no modulus needed', external_referent=dict(kind='our_own_fixture', locator='raw/TRANSPORT_WITNESS.json', compared_quantity='Observation sufficiency in nonnegative conserved edge-flow port', refutes_us=True))
                    break
        if p is not None:
            for tag in ['posterior', 'R']:
                (u, l, d) = port(e, [v / total for v in area], tag)
                assert l == u + d
    with open(R / 'raw/TRANSPORT_INTERVALS.csv', 'w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write('raw/TRANSPORT_WITNESS.json', witness)
    aggregate = {}
    for tag in ['posterior', 'R']:
        aggregate[tag] = {}
        for tier in ['one_upper_region', 'full_upper_tooth_marginal']:
            sel = [v for v in rows if v['axis'] == tag and v['tier'] == tier]
            w = np.array([v['width'] for v in sel])
            fraction = float(np.mean(w <= 0.2))
            aggregate[tag][tier] = dict(n=len(sel), mean_width_share=float(w.mean()), median_width_share=float(np.median(w)), max_width_share=float(w.max()), fraction_width_le20pp=fraction, usefulness_gate='PASS' if fraction >= 0.9 else 'FAIL', smallest_net_crossing_extension_exact_width=0.0)
    if witness:
        mutated_rejected = Q(witness['net_crossing1']) + Q(1, 100) != Q(witness['net_crossing1'])
    else:
        mutated_rejected = False
    out = dict(claim_type='information_link', n_source_cases=len(cases), physical_measurements_consumed=0, conditional_information_link='NONNEGATIVE_EDGE_FORCE_PORT_ONLY', aggregate=aggregate, dropout=drop, n_zero_area_fullmarginal=len([c for c in cases if sum((c['area'][e] for e in c['edges'] if e[0] % 10 != 8 and e[1] % 10 != 8)) == 0]), LP_calls=LP_CALLS, max_LP_bound_error_share=maxerr, max_LP_balance_residual_share=maxres, rank_identity_max_error=max(map(abs, rankerrs)), crossing_reconstruction_error_exact='0', injected_plus1pp_crossing_rejected=mutated_rejected, missing_lower_dof_histogram={str(d): sum((1 for c in cases if len({l for (u, l) in c['edges'] if u % 10 != 8 and l % 10 != 8}) - components([e for e in c['edges'] if e[0] % 10 != 8 and e[1] % 10 != 8]) == d)) for d in range(15)}, scipy_version=scipy.__version__, wall_s=time.perf_counter() - start, max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    write('raw/TRANSPORT_SUMMARY.json', out)
    write('CURRENT_WORK_STATE.json', dict(lane='X94-force-validation', status='MEASURED_FORCE_TRANSFER_PORT_TESTED', latest_gate=aggregate, next_operation='Check denominator sensitivity, independent distribution arithmetic, fault injection, figure and lab measurement handoff', physical_force='UNKNOWN; no same-patient force measurements'))
    print(json.dumps({k: v for (k, v) in out.items() if k not in ['dropout', 'missing_lower_dof_histogram']}))
if __name__ == '__main__':
    main()
