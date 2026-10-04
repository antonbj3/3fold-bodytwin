from geometry import *
from carriers import explicit_check
from fractions import Fraction as Q
from scipy.optimize import linprog
import warnings

def triangle_lp(A, B, direction=None):
    A = np.asarray(A)
    B = np.asarray(B)
    nvar = 7 if direction is not None else 6
    E = np.zeros((5, nvar))
    E[:3, :3] = A.T
    E[:3, 3:6] = -B.T
    E[3, :3] = 1
    E[4, 3:6] = 1
    if direction is not None:
        E[:3, 6] = -np.asarray(direction)
    rhs = np.array([0, 0, 0, 1, 1.0])
    bounds = [(0, 1)] * 6 + ([(0, 20)] if direction is not None else [])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = linprog(np.zeros(nvar), A_eq=E, b_eq=rhs, bounds=bounds, method='highs', options={'threads': 1, 'primal_feasibility_tolerance': 1e-09})
    if not res.success:
        return {'status': 'UNKNOWN_NO_PRIMAL', 'solver_status': res.status, 'not_a_rigorous_separation': True}
    x = res.x
    err = float(np.max(abs(E @ x - rhs)))
    return {'status': 'NUMERIC_INTERSECTION_WITNESS' if err < 1e-08 and x.min() >= -1e-10 else 'UNKNOWN_RESIDUAL', 'weights': x, 'residual': err, 'point_A': x[:3] @ A, 'point_B': x[3:6] @ B, 'translation': float(x[6]) if direction is not None else 0.0, 'rigorous_reconstruction': 'not yet performed'}

def exact_volume_under_roof(tris, zlower=-1):
    total = Q(0)
    for t in tris:
        (a, b, c) = [[Q(float(v)) for v in p] for p in t]
        twice = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
        total += twice / 2 * ((a[2] + b[2] + c[2]) / 3 - Q(zlower))
    return total

def roof_fixture():
    x = np.arange(-2.0, 3.0)
    (xx, yy) = np.meshgrid(x, x, indexing='ij')
    v = np.c_[xx.ravel(), yy.ravel(), np.full(xx.size, -0.6)]
    idx = np.arange(len(v)).reshape(xx.shape)
    f = np.r_[np.stack([idx[:-1, :-1], idx[1:, :-1], idx[:-1, 1:]], -1).reshape(-1, 3), np.stack([idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]], -1).reshape(-1, 3)]
    outer = v.copy()
    outer[:, 2] = 0
    changed = v.copy()
    changed[np.flatnonzero((v[:, 0] == -1) & (v[:, 1] == 0))[0], 2] += 0.2
    changed[np.flatnonzero((v[:, 0] == 1) & (v[:, 1] == 0))[0], 2] -= 0.2
    return (outer[f], v[f], changed[f])

def run_controls():
    tic = time.perf_counter()
    (o, a, b) = roof_fixture()
    bulge_a = wall_bound(o, a)
    bulge_b = wall_bound(o, b)
    sa = b.copy()
    zlo = float(b[np.any(b[:, :, 0] == 1, axis=1), :, 2].min())
    zhi = float(b[:, :, 2].max())
    for tr in sa:
        for point in tr:
            if point[1] == 0 and abs(point[0]) == 1:
                point[2] = zlo if point[0] == -1 else zhi
    so = o.copy()
    so[:, :, 2] = np.where(so[:, :, 0] <= 0, 0, 0.4 * np.minimum(so[:, :, 0], 1))
    va = exact_volume_under_roof(sa)
    vb = exact_volume_under_roof(b)
    sample = np.array([[-2, -2, 0], [-2, 2, 0], [2, -2, 0.4], [2, 2, 0.4]], float)
    da = distances(source_mesh(sa), sample)
    db = distances(source_mesh(b), sample)
    wa = wall_bound(so, sa)
    wb = wall_bound(so, b)

    def gate(w):
        return 'FAIL' if w['sampled_min_mm'] < 0.5 else 'MODEL_PASS' if w['geometric_lower_mm'] >= 0.5 else 'UNKNOWN'
    ga = gate(wa)
    gb = gate(wb)
    suff = {'summary': 'exact volume under roof plus sampled thickness multiset', 'volume_A_mm3': float(va), 'volume_B_mm3': float(vb), 'exact_volume_difference': str(va - vb), 'machine_volume_identity_error_mm3': abs(float(va) - float(vb)), 'sampled_histogram_A_mm': np.sort(da), 'sampled_histogram_B_mm': np.sort(db), 'sampled_histogram_identity_error_mm': float(np.max(abs(np.sort(da) - np.sort(db)))), 'downstream_whole_wall_A': wa, 'downstream_whole_wall_B': wb, 'minimum_sampled_difference_mm': wa['sampled_min_mm'] - wb['sampled_min_mm'], 'decision_A': ga, 'decision_B': gb, 'source_resolution': 'PER_POINT', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'code/controls.py:roof_fixture; exact piecewise planar roof integration', 'compared_quantity': 'summary insufficiency only, not external dental truth', 'refutes_us': True}, 'minimum_extension': 'retain certified continuous wall minimum for this gate; spatial wall field for future regional constraints', 'brief_logical_limit': 'identical full actual thickness multisets cannot have different minima; identical sampled histograms can.'}
    checks = [{'name': 'continuous_wall_0.2mm_bulge', 'valid': gate(bulge_a) == 'MODEL_PASS', 'injected_result': gate(bulge_b), 'injected_rejected': gate(bulge_b) == 'FAIL', 'amplitude_mm': 0.2, 'baseline_wall': bulge_a, 'injected_wall': bulge_b}, {'name': 'sufficiency', 'valid': va == vb and np.array_equal(np.sort(da), np.sort(db)) and (ga != gb), 'injected_rejected': not np.array_equal(np.sort(da), np.sort(db) + 0.01), 'injected_value': 'histogram changed by 0.01mm'}]
    T = np.ones((7, 7, 4), bool)
    F = np.zeros_like(T)
    F[3, 3, 2] = 1
    C = np.array([[3, 3, 2]])
    good = explicit_check(T, F, C)
    badT = T.copy()
    badT[3, 3, 1] = False
    bad = explicit_check(badT, F, C)
    checks.append({'name': 'source_carrier_box_containment', 'valid': not good[0]['outside_boxes'], 'injected_rejected': bool(bad[0]['outside_boxes']), 'injected_value': 'one source box removed', 'actual_control': bad})
    A = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], float)
    B = np.array([[0.25, 0.25, -1], [0.25, 0.25, 1], [0.75, 0.25, 0]], float)
    hit = triangle_lp(A, B)
    sweep = triangle_lp(A, B + np.array([0, 0, 4]), [0, 0, -1])
    invalid = {'weights': np.zeros(6), 'reported': 'intersection'}
    validator = lambda w: bool(len(w) == 6 and min(w) >= 0 and (abs(sum(w[:3]) - 1) < 1e-08) and (abs(sum(w[3:]) - 1) < 1e-08) and (np.max(abs(w[:3] @ A - w[3:] @ B)) < 1e-08))
    checks.append({'name': 'independent_triangle_LP', 'valid': hit['status'] == 'NUMERIC_INTERSECTION_WITNESS' and validator(np.array(hit['weights'])), 'injected_rejected': not validator(invalid['weights']), 'injected_value': 'all barycentric weights zero', 'hit': hit, 'sweep': sweep})
    from exact_union import primitive, mesh_exact, volume2
    P = [primitive((Q(0), Q(0)), Q(2), Q(0), Q(1, 5)), primitive((Q(1, 5), Q(0)), Q(2), Q(0), Q(1, 5))]
    (_, certificate, _, _) = mesh_exact(P)
    expected = volume2((Q(0), Q(0)), Q(2), (Q(1, 5), Q(0)), Q(2), Q(0), Q(1, 5))
    checks.append({'name': 'exact_union_volume', 'valid': Q(certificate['exact_volume_mm3']) == expected and certificate['watertight'] and (certificate['exact_cone_failures'] == 0), 'injected_value': 'reported volume plus0.2mm3', 'injected_rejected': Q(certificate['exact_volume_mm3']) + Q(1, 5) != expected, 'actual_volume_mm3': certificate['volume_mm3']})
    from run_r3 import sweep_planes, exact_triangle_sweep, highs_control
    SP = sweep_planes([0, 0, 2], 0)
    clear = np.array([[10, 0, 1], [10.01, 0, 1], [10, 0.01, 1]])
    (before, _) = exact_triangle_sweep(clear, SP)
    (after, mat) = exact_triangle_sweep(clear - np.array([10, 0, 0]), SP)
    checks.append({'name': 'exact_whole_triangle_sweep', 'valid': before['status'].startswith('EXACT_DISJOINT') and after['status'] == 'EXACT_COLLISION_WITNESS' and highs_control(mat)['valid_numeric_primal'], 'injected_value': 'shell triangle translated into actual swept preparation', 'injected_rejected': after['status'] == 'EXACT_COLLISION_WITNESS', 'actual_witness': after})
    np.savez_compressed(DATA / 'SUFFICIENCY_AND_BULGE.npz', outer=o, original_inner=a, bulged_inner=b, suff_outer=so, suff_A=sa, suff_B=b, sample_points=sample)
    dump(ROOT / 'raw/SUFFICIENCY.json', suff)
    dump(ROOT / 'raw/CONTROLS.json', {'rows': checks, 'all_pass': all((c['valid'] and c['injected_rejected'] for c in checks)), 'seconds': time.perf_counter() - tic, 'resolution': 'PER_POINT', 'sufficiency_path': DATA / 'SUFFICIENCY_AND_BULGE.npz', 'sufficiency_sha256': sha(DATA / 'SUFFICIENCY_AND_BULGE.npz'), 'claim_scope': 'geometry falsifiability, not dental validation', 'LP_runtime_threads': 1})
    print(json.dumps(clean({'controls': checks, 'volume_identical': va == vb, 'histogram_identical': np.array_equal(da, db), 'wall_A': ga, 'wall_B': gb})), flush=True)
if __name__ == '__main__':
    run_controls()
