from common_chain import *
import numpy as np
from fractions import Fraction as Q
from decimal import Decimal, localcontext
import math

def upward_sqrt(s, den=10 ** 12):
    n = (s.numerator * den * den + s.denominator - 1) // s.denominator
    k = math.isqrt(n)
    if k * k < n:
        k += 1
    u = Q(k, den)
    assert u * u >= s
    return u

def rational_point(tri, bary):
    return [sum((Q(float(tri[j, k])) * bary[j] for j in range(3))) for k in range(3)]

def make_witness(v, f, roles, proposal):
    import trimesh
    (i, j) = (proposal['exterior_face'], proposal['intaglio_face'])
    assert roles[i] == 0 and roles[j] == 1
    (a, b) = (v[f[i]], v[f[j]])
    approximate = trimesh.triangles.points_to_barycentric(b[None], np.array(proposal['inner_point_mm'])[None])[0]
    qb = [Q(float(max(0, x))).limit_denominator(10 ** 9) for x in approximate]
    qb = [x / sum(qb) for x in qb]
    pa = [Q(1, 3)] * 3
    p = rational_point(a, pa)
    q = rational_point(b, qb)
    s = sum(((x - y) ** 2 for (x, y) in zip(p, q)))
    u = upward_sqrt(s)
    return dict(exterior_face=i, intaglio_face=j, outer_barycentric=list(map(str, pa)), inner_barycentric=list(map(str, qb)), squared_distance_mm2=str(s), upper_distance_mm=str(u), upper_distance_decimal_mm=float(u))

def verify(v, f, roles, w):
    try:
        (i, j) = (w['exterior_face'], w['intaglio_face'])
        a = list(map(Q, w['outer_barycentric']))
        b = list(map(Q, w['inner_barycentric']))
        if roles[i] != 0 or roles[j] != 1 or len(a) != 3 or (len(b) != 3) or (min(a + b) < 0) or (sum(a) != 1) or (sum(b) != 1):
            return False
        p = rational_point(v[f[i]], a)
        q = rational_point(v[f[j]], b)
        s = sum(((x - y) ** 2 for (x, y) in zip(p, q)))
        u = Q(w['upper_distance_mm'])
        return s == Q(w['squared_distance_mm2']) and u >= 0 and (u * u >= s)
    except (ValueError, IndexError, KeyError):
        return False

def run():
    import copy, time
    verify_freeze(ROOT / 'PREREG_R2.json')
    t = time.perf_counter()
    rows = []
    for r in read(ROOT / 'raw/R1_ALL.json')['rows']:
        assert sha(r['mesh_path']) == r['mesh_sha256']
        with np.load(r['mesh_path'], allow_pickle=False) as a:
            (v, f, roles) = (a['vertices'], a['faces'], a['roles'])
        w = make_witness(v, f, roles, r['wall_candidate'])
        assert verify(v, f, roles, w)
        u = Q(w['upper_distance_mm'])
        s = Q(w['squared_distance_mm2'])
        budget = Q(1, 20)
        with localcontext() as ctx:
            ctx.prec = 80
            decimal_dist = (Decimal(s.numerator) / Decimal(s.denominator)).sqrt()
            decimal_upper = Decimal(u.numerator) / Decimal(u.denominator)
            assert decimal_dist <= decimal_upper
        faults = []
        for (field, badval) in [('upper_distance_mm', '0'), ('squared_distance_mm2', '1'), ('inner_barycentric', ['-1', '1', '1']), ('intaglio_face', w['exterior_face'])]:
            bad = copy.deepcopy(w)
            bad[field] = badval
            faults.append(dict(field=field, rejected=not verify(v, f, roles, bad)))
        assert all((x['rejected'] for x in faults))
        rows.append(dict(key=r['key'], mesh_sha256=r['mesh_sha256'], witness=w, verified=True, minimum_total_surface_motion_necessary_mm=float(max(Q(0), Q(1, 2) - u)), fixed_inner_50um_edit_status='FAIL_REPAIR' if u + budget < Q(1, 2) else 'UNKNOWN', both_surfaces_50um_edit_status='FAIL_REPAIR' if u + 2 * budget < Q(1, 2) else 'UNKNOWN', exact_fixed_inner_remaining_deficit_mm=str(Q(1, 2) - u - budget), faults=faults, decimal_independent_distance_mm=str(decimal_dist), numerical_candidate_difference_mm=abs(float(u) - r['wall_candidate']['distance_mm'])))
    a = [Q(1, 4), Q(3, 4)]
    b = [Q(1, 2), Q(1, 2)]
    assert sum(a) == sum(b)
    suff = dict(kind='our_own_fixture', summary='area-weighted mean thickness on two equal-area planar patches', resolution='PER_SURFACE_REGION', states_mm=[list(map(str, a)), list(map(str, b))], identical_summary_mm=str(sum(a) / 2), identity_error=0, downstream_minimum_mm=[str(min(a)), str(min(b))], downstream_difference_mm=str(min(b) - min(a)), wall_gate=[min(a) >= Q(1, 2), min(b) >= Q(1, 2)], minimum_extension='For this fixed wall-threshold question: certified global minimum; for edits: tracked local separation or full field. Mean is insufficient.')
    result = dict(claim_type='capability', round='R2', rows=rows, sufficiency=suff, seconds=time.perf_counter() - t, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, external_referent=read(ROOT / 'PREREG_R2.json')['external_referent'])
    dump(ROOT / 'raw/R2_EXACT_REPAIR.json', result)
    return result
if __name__ == '__main__':
    run()
