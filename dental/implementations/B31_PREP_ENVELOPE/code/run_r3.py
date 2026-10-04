"""Whole exported crown vs the finite preparation's entire axial withdrawal sweep.

Clip every triangle by swept-carrier halfspaces in exact rational barycentric coordinates.
SciPy HiGHS is the distinct bounded comparator on all potentially intersecting pairs.
"""
from geometry import *
from exact_union import *
from scipy.optimize import linprog
import warnings

def sweep_planes(c, base, length=Q(20)):
    (x, y, a) = map(lambda x: Q(float(x)), c)
    b = Q(float(base))
    r0 = H / 2 - S * H / 2 + S * Q(1, 20) - Q(51, 1000)
    a -= Q(1, 20)
    rad = r0 + S * (a - b)
    return [((Q(1), Q(0), S), x + r0 + S * a), ((Q(-1), Q(0), S), -x + r0 + S * a), ((Q(0), Q(1), S), y + r0 + S * a), ((Q(0), Q(-1), S), -y + r0 + S * a), ((Q(1), Q(0), Q(0)), x + rad), ((Q(-1), Q(0), Q(0)), -x + rad), ((Q(0), Q(1), Q(0)), y + rad), ((Q(0), Q(-1), Q(0)), -y + rad), ((Q(0), Q(0), Q(1)), a), ((Q(0), Q(0), Q(-1)), length - b)]

def exact_triangle_sweep(A, planes):
    aq = [tuple((Q(float(v)) for v in p)) for p in A]
    E = [tuple((dot(n, p) - b for p in aq)) for (n, b) in planes]
    vals = np.array([[float(x) for x in row] for row in E])
    cand = np.flatnonzero(vals.min(1) > 1e-09)
    if len(cand):
        j = int(cand[0])
        if min(E[j]) > 0:
            return ({'status': 'EXACT_DISJOINT_SINGLE_HALFSPACE', 'separator_plane': j, 'minimum_exact_residual': str(min(E[j]))}, vals)
    poly = [(Q(1), Q(0), Q(0)), (Q(0), Q(1), Q(0)), (Q(0), Q(0), Q(1))]
    for row in E:
        poly = clip(poly, row, Q(0), True)
        if not poly:
            return ({'status': 'EXACT_DISJOINT_CLIPPED_SIMPLEX'}, vals)
    alpha = poly[0]
    point = tuple((sum((alpha[i] * aq[i][j] for i in range(3))) for j in range(3)))
    assert sum(alpha) == 1 and min(alpha) >= 0 and all((dot(row, alpha) <= 0 for row in E))
    return ({'status': 'EXACT_COLLISION_WITNESS', 'barycentric': [str(x) for x in alpha], 'point_exact_mm': [str(x) for x in point], 'point_mm': [float(x) for x in point], 'primal_residual_exact': '0'}, vals)

def highs_control(vals):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = linprog(np.zeros(3), A_ub=vals, b_ub=np.zeros(len(vals)), A_eq=np.ones((1, 3)), b_eq=[1.0], bounds=[(0, 1)] * 3, method='highs', options={'threads': 1, 'primal_feasibility_tolerance': 1e-09})
    valid = bool(res.success and min(res.x) >= -1e-09 and (abs(sum(res.x) - 1) < 1e-08) and (max(vals @ res.x) <= 1e-08))
    return {'solver_status': res.status, 'valid_numeric_primal': valid, 'weights': res.x.tolist() if res.success else None, 'infeasible_status_is_not_rigorous': True}

def one(rec):
    tic = time.perf_counter()
    d = load_np(rec['candidate_path'])
    tri = d['vertices'][d['faces']]
    planes = [sweep_planes(c, float(d['base'])) for c in d['carrier_coords']]
    flags = np.zeros((len(tri), len(planes)), np.uint8)
    witnesses = []
    controls = []
    clipcount = 0
    for (j, P) in enumerate(planes):
        for (i, A) in enumerate(tri):
            (out, val) = exact_triangle_sweep(A, P)
            if out['status'] == 'EXACT_DISJOINT_SINGLE_HALFSPACE':
                flags[i, j] = 1
                continue
            clipcount += 1
            h = highs_control(val)
            exacthit = out['status'] == 'EXACT_COLLISION_WITNESS'
            flags[i, j] = 3 if exacthit else 2
            h.update(shell_facet=i, carrier=j, exact_hit=exacthit, comparator_agrees=exacthit == h['valid_numeric_primal'])
            controls.append(h)
            if exacthit:
                out.update(shell_facet=i, carrier=j, role=int(d['roles'][i]), source_triangle_mm=A, axial_displacement_mm=max(0.0, float(d['base']) - out['point_mm'][2]))
                witnesses.append(out)
    path = DATA / 'r3' / (rec['key'] + '_whole_path_flags.npz')
    path.parent.mkdir(exist_ok=True)
    np.savez_compressed(path, flags=flags, roles=d['roles'])
    raw = ROOT / 'raw/r3' / (rec['key'] + '_exact_path.json')
    dump(raw, {'key': rec['key'], 'candidate_sha256': rec['candidate_sha256'], 'flags_sha256': sha(path), 'witnesses': witnesses, 'controls': controls, 'claim_scope': 'exact digital triangle and rational swept frusta, axial20mm path; actual root/preparation/source uncertainty excluded'})
    return {'key': rec['key'], 'family': rec['family'], 'resolution': 'PER_SURFACE_REGION', 'candidate_sha256': rec['candidate_sha256'], 'whole_path_status': 'EXACT_COLLISION_COUNTERWITNESS' if witnesses else 'EXACT_CLEAR_FOR_FINITE_PREP_AXIS_MODEL', 'path_mm': 20, 'direction': [0, 0, 1], 'triangles': len(tri), 'swept_carriers': len(planes), 'whole_pairs_checked': len(tri) * len(planes), 'single_halfspace_certificates': int((flags == 1).sum()), 'exact_simplex_clips': clipcount, 'collision_pairs': len(witnesses), 'collision_roles': {str(k): sum((w['role'] == k for w in witnesses)) for k in [0, 1, 2]}, 'LP_controls': len(controls), 'LP_control_disagreements': sum((not c['comparator_agrees'] for c in controls)), 'flag_path': path, 'flag_sha256': sha(path), 'exact_witness_file': raw, 'exact_witness_sha256': sha(raw), 'seconds': time.perf_counter() - tic, 'physical_status': 'UNKNOWN_REAL_ROOT_FINISH_LINE_PRODUCT', 'whole_crown_status': 'UNKNOWN_SOURCE_SELF_INTERSECTION_AND_WALL_ARITHMETIC' if not witnesses else 'REJECT_CROWN_INSERTION_FOR_EXPORTED_GEOMETRY'}

def main():
    tic = time.perf_counter()
    rows = []
    for rec in read(ROOT / 'raw/R2_RESULTS.json'):
        if 'candidate_path' not in rec:
            r = {'key': rec['key'], 'family': rec['family'], 'whole_path_status': 'UNKNOWN_NO_CANDIDATE', 'missing': rec.get('reason'), 'resolution': 'PER_SURFACE_REGION'}
        else:
            assert sha(rec['candidate_path']) == rec['candidate_sha256']
            r = one(rec)
        rows.append(r)
        dump(ROOT / 'raw/R3_RESULTS.json', rows)
        state('R3_RUNNING', {'completed': len(rows), 'requested': 18}, 'Exact full-triangle sweep checks, then preserve remaining source/margin uncertainty')
        print(r['key'], r['whole_path_status'], r.get('collision_pairs'), r.get('LP_controls'), flush=True)
    dump(ROOT / 'raw/R3_COST.json', {'seconds': time.perf_counter() - tic, **usage()})
    state('R3_DECIDED', {'whole_path_clear': sum((r['whole_path_status'] == 'EXACT_CLEAR_FOR_FINITE_PREP_AXIS_MODEL' for r in rows)), 'requested': 18}, 'Package exact digital results and missing source/shoulder observations; fresh copy replay')
if __name__ == '__main__':
    main()
