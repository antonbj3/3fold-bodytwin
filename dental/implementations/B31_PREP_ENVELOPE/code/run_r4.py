"""Exact retest of frozen positive original/shoulder triangle pairs only."""
from exact_union import Q, dot, cross, clip
from common import *
from controls import triangle_lp
from collections import Counter

def sub(a, b):
    return tuple((x - y for (x, y) in zip(a, b)))

def normal(A):
    return cross(sub(A[1], A[0]), sub(A[2], A[0]))

def bary(A, p):
    n = normal(A)
    drop = next((i for (i, x) in enumerate(n) if x))
    ij = [i for i in range(3) if i != drop]
    (a, b, c) = [tuple((v[i] for i in ij)) for v in A]
    q = tuple((p[i] for i in ij))
    area = lambda x, y: x[0] * y[1] - x[1] * y[0]
    u = sub(b, a)
    v = sub(c, a)
    w = sub(q, a)
    det = area(u, v)
    beta = area(w, v) / det
    gamma = area(u, w) / det
    return (1 - beta - gamma, beta, gamma)

def cut_plane(A, n, b):
    e = [dot(n, p) - b for p in A]
    pts = [p for (p, v) in zip(A, e) if v == 0]
    for (i, j) in [(0, 1), (1, 2), (2, 0)]:
        if e[i] * e[j] < 0:
            t = e[i] / (e[i] - e[j])
            pts.append(tuple((A[i][k] + t * (A[j][k] - A[i][k]) for k in range(3))))
    return list(dict.fromkeys(pts))

def primal_valid(A, B, p, a, b):
    return sum(a) == sum(b) == 1 and min(a + b) >= 0 and all((sum((a[i] * A[i][k] for i in range(3))) == p[k] == sum((b[i] * B[i][k] for i in range(3))) for k in range(3)))

def exact_pair(Af, Bf):
    A = [tuple((Q(float(x)) for x in p)) for p in Af]
    B = [tuple((Q(float(x)) for x in p)) for p in Bf]
    na = normal(A)
    nb = normal(B)
    if not any(na) or not any(nb):
        return {'status': 'UNKNOWN_DEGENERATE_TRIANGLE'}
    direction = cross(na, nb)
    coplanar = not any(direction)
    if coplanar:
        if dot(na, sub(B[0], A[0])) != 0:
            return {'status': 'EXACT_DISJOINT_PARALLEL_PLANES'}
        drop = next((i for (i, x) in enumerate(na) if x))
        ij = [i for i in range(3) if i != drop]
        project = lambda p: tuple((p[i] for i in ij))
        orient = lambda a, b, p: (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
        BB = [project(p) for p in B]
        sign = 1 if orient(*BB) > 0 else -1
        poly = [(Q(1), Q(0), Q(0)), (Q(0), Q(1), Q(0)), (Q(0), Q(0), Q(1))]
        for i in range(3):
            row = tuple((-sign * orient(BB[i], BB[(i + 1) % 3], project(p)) for p in A))
            poly = clip(poly, row, Q(0), True)
            if not poly:
                return {'status': 'EXACT_DISJOINT_COPLANAR_CLIP'}
        aa = tuple((sum((p[k] for p in poly)) / len(poly) for k in range(3)))
        p = tuple((sum((aa[i] * A[i][k] for i in range(3))) for k in range(3)))
        bb = bary(B, p)
    else:
        ap = cut_plane(A, nb, dot(nb, B[0]))
        bp = cut_plane(B, na, dot(na, A[0]))
        if not ap or not bp:
            return {'status': 'EXACT_DISJOINT_PLANE_SECTION'}
        axis = next((i for (i, x) in enumerate(direction) if x))
        (alo, ahi) = (min((p[axis] for p in ap)), max((p[axis] for p in ap)))
        (blo, bhi) = (min((p[axis] for p in bp)), max((p[axis] for p in bp)))
        lo = max(alo, blo)
        hi = min(ahi, bhi)
        if lo > hi:
            return {'status': 'EXACT_DISJOINT_LINE_INTERVAL'}
        value = (lo + hi) / 2
        if len(ap) == 1:
            p = ap[0]
        else:
            a = min(ap, key=lambda p: p[axis])
            b = max(ap, key=lambda p: p[axis])
            t = (value - a[axis]) / (b[axis] - a[axis])
            p = tuple((a[k] + t * (b[k] - a[k]) for k in range(3)))
        aa = bary(A, p)
        bb = bary(B, p)
    assert primal_valid(A, B, p, aa, bb)
    strict = min(aa + bb) > 0
    bad = (p[0] + Q(1, 5), p[1], p[2])
    reject = not primal_valid(A, B, bad, aa, bb)
    assert reject
    return {'status': ('EXACT_COPLANAR_INTERIOR_OVERLAP' if coplanar else 'EXACT_STRICT_INTERIOR_CROSSING') if strict else 'EXACT_BOUNDARY_TOUCH', 'barycentric_added': [str(x) for x in aa], 'barycentric_source': [str(x) for x in bb], 'point_exact_mm': [str(x) for x in p], 'point_mm': [float(x) for x in p], 'exact_primal_residual': '0', 'injected_0.2mm_witness_shift_rejected': reject, 'resolution': 'PER_POINT'}

def main():
    tic = time.perf_counter()
    frozen = {r['key']: r for r in read(ROOT / 'FROZEN_PREDICTIONS_R4.json')['inputs']}
    rows = []
    for rec in read(ROOT / 'raw/R2_RESULTS.json'):
        r = {'key': rec['key'], 'family': rec['family'], 'resolution': 'PER_SURFACE_REGION', 'whole_crown_status': 'UNKNOWN', 'scope': 'only frozen positive original/shoulder pairs; global embedding UNKNOWN'}
        if 'candidate_path' not in rec:
            r.update(status='UNKNOWN_NO_CANDIDATE', missing=rec.get('reason'))
        else:
            assert sha(rec['candidate_path']) == frozen[rec['key']]['candidate_sha256']
            d = load_np(rec['candidate_path'])
            added = d['vertices'][d['faces'][d['roles'] == 2]]
            original = d['original']
            cp = ROOT / 'raw/r2' / (rec['key'] + '_surface_control.json')
            assert sha(cp) == frozen[rec['key']]['control_sha256']
            control = read(cp)
            pairs = [(w['added_facet'], w['source_facet']) for w in control['LP_raw'] if w['status'] == 'NUMERIC_INTERSECTION_WITNESS']
            assert [list(p) for p in pairs] == frozen[rec['key']]['numeric_positive_pairs']
            out = []
            for (ai, bi) in pairs:
                (A, B) = (added[ai], original[bi])
                q = exact_pair(A, B)
                h = triangle_lp(A, B)
                hit = q['status'] in ['EXACT_STRICT_INTERIOR_CROSSING', 'EXACT_COPLANAR_INTERIOR_OVERLAP', 'EXACT_BOUNDARY_TOUCH']
                q.update(added_facet=ai, source_facet=bi, added_triangle_mm=A, source_triangle_mm=B, numeric_control=h, LP_agrees=hit == (h['status'] == 'NUMERIC_INTERSECTION_WITNESS'))
                out.append(q)
            counts = dict(Counter((q['status'] for q in out)))
            proper = sum((counts.get(k, 0) for k in ['EXACT_STRICT_INTERIOR_CROSSING', 'EXACT_COPLANAR_INTERIOR_OVERLAP']))
            raw = ROOT / 'raw/r4' / (rec['key'] + '_exact_surface_pairs.json')
            dump(raw, {'key': rec['key'], 'candidate_sha256': rec['candidate_sha256'], 'frozen_surface_control_sha256': frozen[rec['key']]['control_sha256'], 'pairs': out})
            r.update(status='EXACT_EXPORTED_SHELL_CROSSING_COUNTERWITNESS' if proper else 'NO_PROPER_CROSSING_IN_FROZEN_POSITIVE_PAIRS', tested_pairs=len(out), proper_crossing_pairs=proper, status_counts=counts, LP_disagreements=sum((not q['LP_agrees'] for q in out)), exact_witness_file=str(raw), exact_witness_sha256=sha(raw), whole_crown_status='REJECT_EXPORTED_SHELL_EMBEDDING' if proper else 'UNKNOWN')
        rows.append(r)
        dump(ROOT / 'raw/R4_RESULTS.json', rows)
        print(r['key'], r['status'], r.get('proper_crossing_pairs'), flush=True)
    dump(ROOT / 'raw/R4_COST.json', {'seconds': time.perf_counter() - tic, **usage()})
    state('R4_DECIDED', {'completed': len(rows), 'requested': 18, 'proper_crossing_cases': sum((r['status'] == 'EXACT_EXPORTED_SHELL_CROSSING_COUNTERWITNESS' for r in rows)), 'complete_crowns': 0}, 'Preserve exact surface refutations; next construction must change shoulder Boolean boundary')
if __name__ == '__main__':
    main()
