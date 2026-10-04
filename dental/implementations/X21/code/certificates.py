"""Conservative projected-minimum enclosure on the rounded mesh, not patient truth."""
import sys, json, time, zipfile
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from numba import njit
from contact import P, D, sha, write, state, load_case, boxes, buckets

@njit(cache=True)
def lower_bounds(U, L, ul, ll, origin, h, shape, ptr, ids):
    out = np.full((29, 49), np.inf)
    seen = np.full(len(L), -1, np.int32)
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
                    checked += 1
                    v = U[i, 2] - L[j, 5]
                    if v < out[a, b]:
                        out[a, b] = v
    return (out, checked)

def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]

def sub(a, b):
    return [a[0] - b[0], a[1] - b[1]]

def rational_min(U, L):
    u = [[F(float(x)) for x in v] for v in U]
    l = [[F(float(x)) for x in v] for v in L]
    poly = [v[:2] for v in u]
    sgn = 1 if cross(sub(l[1], l[0]), sub(l[2], l[0])) > 0 else -1
    for e in range(3):
        edge = sub(l[(e + 1) % 3], l[e])
        new = []
        if not poly:
            return (None, None)
        for (i, q) in enumerate(poly):
            p = poly[i - 1]
            dp = sgn * cross(edge, sub(p, l[e]))
            dq = sgn * cross(edge, sub(q, l[e]))
            pi = dp >= 0
            qi = dq >= 0
            if pi != qi:
                t = dp / (dp - dq)
                new.append([p[k] + t * (q[k] - p[k]) for k in range(2)])
            if qi:
                new.append(q)
        poly = new
    if not poly:
        return (None, None)

    def plane(T):
        a = [T[1][i] - T[0][i] for i in range(3)]
        b = [T[2][i] - T[0][i] for i in range(3)]
        det = cross(a, b)
        aa = (a[2] * b[1] - b[2] * a[1]) / det
        bb = (a[0] * b[2] - b[0] * a[2]) / det
        return (aa, bb, T[0][2] - aa * T[0][0] - bb * T[0][1])
    cu = plane(u)
    cl = plane(l)
    co = [cu[i] - cl[i] for i in range(3)]
    values = [co[0] * p[0] + co[1] * p[1] + co[2] for p in poly]
    idx = min(range(len(values)), key=lambda i: values[i])
    return (values[idx], poly[idx])

def independent_pair_lb(ub, lb):
    best = np.inf
    for U in ub:
        ok = (U[0] <= lb[:, 3]) & (U[3] >= lb[:, 0]) & (U[1] <= lb[:, 4]) & (U[4] >= lb[:, 1])
        if ok.any():
            best = min(best, float(np.min(U[2] - lb[ok, 5])))
    return float(np.nextafter(best, -np.inf))

def run():
    m = json.load(open(P / 'raw/INPUT_MANIFEST.json'))
    panel = set(m['numerical_panel'])
    done = {p.stem for p in (P / 'raw/certificates').glob('*.json')} if (P / 'raw/certificates').exists() else set()
    st = time.perf_counter()
    with zipfile.ZipFile(m['zip']) as z:
        for c in m['cases']:
            if c in done:
                continue
            state('R3_CERTIFICATES_RUNNING', 'PREREG_R3 frozen; conditional mesh proof only', 'All-face lower bound and exact rational upper witness ' + c, completed=len(done), expected=993)
            s = time.perf_counter()
            (meta, R, center, arches) = load_case(z, c)
            (U, ul, uc, uf, us) = arches[0]
            (L, ll, lc, lf, ls) = arches[1]
            ub = boxes(U)
            lb = boxes(L)
            origin = np.minimum(ub[:, :2].min(0), lb[:, :2].min(0)) - 0.8
            hi = np.maximum(ub[:, 3:5].max(0), lb[:, 3:5].max(0))
            shape = np.ceil((hi - origin) / 0.8).astype(np.int64) + 2
            (ptr, ids) = buckets(lb, origin, 0.8, shape)
            bs = time.perf_counter()
            (low, checked) = lower_bounds(ub, lb, ul, ll, origin, 0.8, shape, ptr, ids)
            bound_s = time.perf_counter() - bs
            source = json.load(open(P / 'raw/cases' / (c + '.json')))
            out = []
            for pair in source['pairs']:
                a = pair['upper_fdi']
                b = pair['lower_fdi']
                w = pair['witness']
                uv = np.array(w['upper_triangle_frame_mm'])
                lv = np.array(w['lower_triangle_frame_mm'])
                ui = int(np.searchsorted(uf, w['upper_source_face']))
                li = int(np.searchsorted(lf, w['lower_source_face']))
                if not np.array_equal(U[ui], uv) or not np.array_equal(L[li], lv):
                    raise ValueError('ROUND_MODEL_WITNESS_COORDINATE_DRIFT ' + c)
                (upper, xy) = rational_min(uv, lv)
                lower = float(np.nextafter(low[a, b], -np.inf))
                nominal = pair['minimum_projected_gap_mm']
                hi0 = float(np.nextafter(float(upper), np.inf)) if upper is not None else None
                valid = hi0 is not None and lower <= hi0 and (lower - 1e-07 <= nominal <= hi0 + 1e-07)
                out.append(dict(upper_fdi=a, lower_fdi=b, minimum_gap_enclosure_mm=[lower, hi0], nominal_minimum_mm=nominal, interval_width_mm=hi0 - lower if hi0 is not None else None, valid=valid, exact_upper_fraction={'numerator': str(upper.numerator), 'denominator': str(upper.denominator)} if upper is not None else None, exact_witness_xy_fraction=[str(q) for q in xy] if xy is not None else None, upper_source_face=w['upper_source_face'], lower_source_face=w['lower_source_face'], assured_separation_beyond_0_1mm=bool(lower > 0.1), witnessed_penetration_beyond_0_03mm=bool(hi0 is not None and hi0 < -0.03), proximity_witness_within_0_1mm=bool(hi0 is not None and hi0 <= 0.1), gap_offset_scenario_mm=[lower - 0.05, hi0 + 0.05] if hi0 is not None else None, physical_status='UNKNOWN'))
            control = None
            if c in panel:
                (a, b) = source['numerical_control']['pair']
                s0 = time.perf_counter()
                truth = independent_pair_lb(ub[ul == a], lb[ll == b])
                me = float(np.nextafter(low[a, b], -np.inf))
                err = abs(truth - me)
                target = next((x for x in out if x['upper_fdi'] == a and x['lower_fdi'] == b))
                control = dict(pair=[a, b], independent_lower_bound_mm=truth, bucket_lower_bound_mm=me, parity_mm=err, gate=err <= 1e-12, lower_bound_injected_plus100_rejected=target['minimum_gap_enclosure_mm'][0] + 100 > target['minimum_gap_enclosure_mm'][1], query_s=time.perf_counter() - s0)
            write(P / 'raw/certificates' / (c + '.json'), dict(case_id=c, source_case_sha256=sha(P / 'raw/cases' / (c + '.json')), rounded_mesh_coordinates='IEEE754 float64 X11-frame outputs treated as exact rational input', domain='All labelled nondegenerate-projected triangles, no anatomical/pose/physical certificate', pairs=out, n_aabb_pairs_scanned=int(checked), bound_query_s=bound_s, control=control, wall_s=time.perf_counter() - s))
            done.add(c)
            if len(done) % 20 == 0:
                print('certified', len(done), 'elapsed', round(time.perf_counter() - st, 1), 's', flush=True)
    summarize()

def summarize():
    files = sorted((P / 'raw/certificates').glob('*.json'))
    rr = [json.load(open(p)) for p in files]
    pairs = [x for r in rr for x in r['pairs']]
    valid = [p for p in pairs if p['valid']]
    width = np.array([p['interval_width_mm'] for p in valid])
    ctrl = [r['control'] for r in rr if r['control']]
    med = float(np.median(width))
    frac = float(np.mean(width <= 0.2))
    gates = dict(completeness=len(rr) == 993 and len(valid) == len(pairs), nominal_enclosed=all((p['valid'] for p in pairs)), numerical_controls=len(ctrl) == 12 and all((c['gate'] and c['lower_bound_injected_plus100_rejected'] for c in ctrl)), useful_median_width=med <= 0.1, useful_width_fraction=frac >= 0.9)
    result = dict(round='R3', claim_type='capability', external_referent=json.load(open(P / 'PREREG_R3.json'))['external_referent'], n_cases=len(rr), n_pairs=len(pairs), valid_enclosures=len(valid), interval_width_quantiles_mm=np.quantile(width, [0, 0.25, 0.5, 0.75, 0.9, 0.95, 1]).tolist(), fraction_width_le_0_20mm=frac, assured_separation_pairs=sum((p['assured_separation_beyond_0_1mm'] for p in pairs)), witnessed_penetration_pairs=sum((p['witnessed_penetration_beyond_0_03mm'] for p in pairs)), cases_with_witnessed_penetration=sum((any((p['witnessed_penetration_beyond_0_03mm'] for p in r['pairs'])) for r in rr)), proximity_witness_pairs=sum((p['proximity_witness_within_0_1mm'] for p in pairs)), gates=gates, primary_gate='PASS' if all(gates.values()) else 'FAIL', mathematical_domain='Rounded projected mesh with predicted labels, excluded degenerate projections; arithmetic does not validate patient geometry', control_maximum_parity_mm=max((c['parity_mm'] for c in ctrl)), all_wrong_lower_bound_controls_rejected=all((c['lower_bound_injected_plus100_rejected'] for c in ctrl)), costs=dict(total_wall_s=sum((r['wall_s'] for r in rr)), lower_bound_query_s=sum((r['bound_query_s'] for r in rr)), independent_control_s=sum((c['query_s'] for c in ctrl))), invalid_examples=[dict(case_id=r['case_id'], pair=p) for r in rr for p in r['pairs'] if not p['valid']][:10])
    write(P / 'raw/RESULTS_R3.json', result)
    state('R3_DECIDED', result['primary_gate'], 'Package exact conditional maps, clinical disagreement and next physical validation contract')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    summarize() if len(sys.argv) > 1 and sys.argv[1] == 'summarize' else run()
