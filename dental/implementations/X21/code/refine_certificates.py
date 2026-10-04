"""Exact rational contender refinement on preregistered panel only."""
import json, time, zipfile
from pathlib import Path
import numpy as np
from numba import njit
from contact import P, D, sha, write, state, load_case, boxes, buckets
from certificates import rational_min

@njit(cache=True)
def contenders(U, L, ul, ll, origin, h, shape, ptr, ids, initial, epsilon):
    skipped = np.full((29, 49), np.inf)
    seen = np.full(len(L), -1, np.int32)
    pairs = []
    checked = 0
    for i in range(len(U)):
        a = ul[i]
        x0 = max(0, int(np.floor((U[i, 0] - origin[0]) / h)))
        x1 = min(shape[0] - 1, int(np.floor((U[i, 3] - origin[0]) / h)))
        y0 = max(0, int(np.floor((U[i, 1] - origin[1]) / h)))
        y1 = min(shape[1] - 1, int(np.floor((U[i, 4] - origin[1]) / h)))
        for x in range(x0, x1 + 1):
            for y in range(y0, y1 + 1):
                bkt = x * shape[1] + y
                for k in range(ptr[bkt], ptr[bkt + 1]):
                    j = ids[k]
                    if seen[j] == i:
                        continue
                    seen[j] = i
                    b = ll[j]
                    if U[i, 0] > L[j, 3] or U[i, 3] < L[j, 0] or U[i, 1] > L[j, 4] or (U[i, 4] < L[j, 1]):
                        continue
                    if not np.isfinite(initial[a, b]):
                        continue
                    checked += 1
                    lower = np.nextafter(U[i, 2] - L[j, 5], -np.inf)
                    threshold = np.nextafter(initial[a, b] - epsilon, np.inf)
                    if lower < threshold:
                        pairs.append((i, j))
                    elif lower < skipped[a, b]:
                        skipped[a, b] = lower
    return (pairs, skipped, checked)

def brute_rational(U, L, ub, lb):
    best = None
    n = 0
    for (i, u) in enumerate(ub):
        jj = np.flatnonzero((u[0] <= lb[:, 3]) & (u[3] >= lb[:, 0]) & (u[1] <= lb[:, 4]) & (u[4] >= lb[:, 1]))
        for j in jj:
            (v, xy) = rational_min(U[i], L[j])
            n += 1
            if v is not None and (best is None or v < best):
                best = v
    return (best, n)

def run():
    pre = json.load(open(P / 'PREREG_R4.json'))
    m = json.load(open(P / 'raw/INPUT_MANIFEST.json'))
    epsilon = pre['target_width_mm']
    start = time.perf_counter()
    with zipfile.ZipFile(m['zip']) as z:
        for c in pre['cases']:
            path = P / 'raw/refined_certificates' / (c + '.json')
            if path.exists():
                continue
            state('R4_EXACT_REFINEMENT', 'R4 prereg frozen; R3 broad bounds retained', 'Exact possible-minimum triangle intersections ' + c)
            st = time.perf_counter()
            (meta, R, center, arches) = load_case(z, c)
            (U, ul, uc, uf, us) = arches[0]
            (L, ll, lc, lf, ls) = arches[1]
            ub = boxes(U)
            lb = boxes(L)
            source = json.load(open(P / 'raw/cases' / (c + '.json')))
            r3 = json.load(open(P / 'raw/certificates' / (c + '.json')))
            best = {}
            initial = np.full((29, 49), np.inf)
            witness = {}
            from fractions import Fraction
            for p in r3['pairs']:
                (a, b) = (p['upper_fdi'], p['lower_fdi'])
                q = p['exact_upper_fraction']
                best[a, b] = Fraction(int(q['numerator']), int(q['denominator']))
                initial[a, b] = p['minimum_gap_enclosure_mm'][1]
                witness[a, b] = dict(upper_source_face=p['upper_source_face'], lower_source_face=p['lower_source_face'], xy_fraction=p['exact_witness_xy_fraction'])
            origin = np.minimum(ub[:, :2].min(0), lb[:, :2].min(0)) - 0.8
            hi = np.maximum(ub[:, 3:5].max(0), lb[:, 3:5].max(0))
            shape = np.ceil((hi - origin) / 0.8).astype(np.int64) + 2
            (ptr, ids) = buckets(lb, origin, 0.8, shape)
            (cand, skipped, checks) = contenders(ub, lb, ul, ll, origin, 0.8, shape, ptr, ids, initial, epsilon)
            n_overlap = 0
            s = time.perf_counter()
            for (i, j) in cand:
                (a, b) = (int(ul[i]), int(ll[j]))
                (v, xy) = rational_min(U[i], L[j])
                if v is None:
                    continue
                n_overlap += 1
                if v < best[a, b]:
                    best[a, b] = v
                    witness[a, b] = dict(upper_source_face=int(uf[i]), lower_source_face=int(lf[j]), xy_fraction=[str(q) for q in xy])
            exact_s = time.perf_counter() - s
            out = []
            for p in source['pairs']:
                (a, b) = (p['upper_fdi'], p['lower_fdi'])
                q = best[a, b]
                exact_lo = float(np.nextafter(float(q), -np.inf))
                upper = float(np.nextafter(float(q), np.inf))
                lower = min(exact_lo, float(skipped[a, b]))
                nominal = p['minimum_projected_gap_mm']
                width = upper - lower
                out.append(dict(upper_fdi=a, lower_fdi=b, minimum_gap_enclosure_mm=[lower, upper], width_mm=width, nominal_minimum_mm=nominal, nominal_inside_tolerance=lower - 1e-07 <= nominal <= upper + 1e-07, gate=width <= 0.005 + 1e-09 and lower <= upper, exact_upper_fraction=dict(numerator=str(q.numerator), denominator=str(q.denominator)), witness=witness[a, b], injected_lower_plus1_rejected=lower + 1 > upper))
            (a, b) = source['numerical_control']['pair']
            s = time.perf_counter()
            (oracle, n) = brute_rational(U[ul == a], L[ll == b], ub[ul == a], lb[ll == b])
            control_s = time.perf_counter() - s
            pp = next((p for p in out if p['upper_fdi'] == a and p['lower_fdi'] == b))
            val = float(oracle)
            valid = pp['minimum_gap_enclosure_mm'][0] - 1e-09 <= val <= pp['minimum_gap_enclosure_mm'][1] + 1e-09
            write(path, dict(case_id=c, source_case_sha256=sha(P / 'raw/cases' / (c + '.json')), source_R3_sha256=sha(P / 'raw/certificates' / (c + '.json')), domain=pre['domain'], pairs=out, aabb_pairs_checked=checks, exact_contender_pairs=len(cand), exact_overlapping_contenders=n_overlap, exact_query_s=exact_s, full_rational_control=dict(pair=[a, b], minimum_mm=val, exact_pairs_tested=n, gate=valid, query_s=control_s), wall_s=time.perf_counter() - st))
            print(c, 'exact contenders', len(cand), 'query_s', round(exact_s, 2), 'total', round(time.perf_counter() - st, 2), 'elapsed', round(time.perf_counter() - start, 1), flush=True)
    summarize()

def summarize():
    rr = [json.load(open(p)) for p in sorted((P / 'raw/refined_certificates').glob('*.json'))]
    pairs = [p for r in rr for p in r['pairs']]
    width = [p['width_mm'] for p in pairs]
    gates = dict(cases=len(rr) == 12, enclosures=all((p['gate'] for p in pairs)), nominal_inside=all((p['nominal_inside_tolerance'] for p in pairs)), full_rational_controls=all((r['full_rational_control']['gate'] for r in rr)), injected_error_rejections=all((p['injected_lower_plus1_rejected'] for p in pairs)))
    result = dict(round='R4', claim_type='capability', external_referent=json.load(open(P / 'PREREG_R4.json'))['external_referent'], n_cases=len(rr), n_pairs=len(pairs), maximum_width_mm=max(width), median_width_mm=float(np.median(width)), gates=gates, primary_gate='PASS' if all(gates.values()) else 'FAIL', exact_contender_pairs=sum((r['exact_contender_pairs'] for r in rr)), all_aabb_pairs_checked=sum((r['aabb_pairs_checked'] for r in rr)), all_wrong_bounds_rejected=all((p['injected_lower_plus1_rejected'] for p in pairs)), cost=dict(wall_s=sum((r['wall_s'] for r in rr)), exact_query_s=sum((r['exact_query_s'] for r in rr)), independent_control_s=sum((r['full_rational_control']['query_s'] for r in rr))), limitations='12 fixed panel cases only;981 full-cohort cases retain broad R3 certificate. No physical pose/FDI/point/contact validation.')
    write(P / 'raw/RESULTS_R4.json', result)
    state('R4_DECIDED', result['primary_gate'], 'Package verified conditional geometric capability and independent measurement requirements')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    run()
