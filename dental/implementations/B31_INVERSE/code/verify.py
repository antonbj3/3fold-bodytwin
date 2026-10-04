from common import *
from inverse import build
from elastic_primitives import stress_tensor, geometry
from parameter_bounds import Response, ratio
from dual_bounds import dual_cell
from fractions import Fraction as F
from scipy.optimize import linprog
import time, warnings

def exact_cover(cover):
    rectangles = [tuple((F(x) for x in row[2:])) for row in cover]
    xs = sorted(set([x for (a, b, l, h) in rectangles for x in [a, b]]))
    if xs[0] != F(0) or xs[-1] != F(30):
        return False
    for (x, y) in zip(xs, xs[1:]):
        bands = sorted(((l, h) for (a, b, l, h) in rectangles if a <= x and b >= y))
        if not bands or bands[0][0] != F(2.0) or bands[-1][1] != F(15.6):
            return False
        if any((h != l1 for ((l, h), (l1, h1)) in zip(bands, bands[1:]))):
            return False
    return True

def bound_record_ok(stored, recomputed, actual):
    return bool(np.isfinite(stored) and abs(stored - recomputed) <= 1e-12 and (stored <= actual + 1e-10))

def run():
    start = time.perf_counter()
    checks = {}
    faults = {}
    locks = read(ROOT / 'INPUT_LOCK.json')
    checks['locked_inputs'] = all((sha(x['path']) == x['sha256'] for x in locks))
    checks['preregs_and_predictions'] = all((sha(p) == p.with_suffix('.sha256').read_text().split()[0] for p in list(ROOT.glob('PREREG_*.json')) + list(ROOT.glob('FROZEN_PREDICTIONS*.json'))))
    z = np.load(DATA / 'B_FIELDS.npz')
    fit = read(ROOT / 'raw/B_FIT.json')
    x = np.array(fit['coefficients_mm'])
    (models, V, T, modes, slaves, faces) = build()
    errs = []
    supererrs = []
    eq = []
    for model in models:
        for (label, xx) in [('reference', np.zeros(4)), ('optimized', x)]:
            a = model.solve(xx)
            saved = z[label + '__' + model.name + '__stress_basis']
            err = float(np.max(abs(a['stress'][:model.nc] - saved)) / max(np.max(abs(saved)), 1e-20))
            errs.append(err)
            eq.append(a['residual'])
            c = np.array([np.cos(np.deg2rad(17)), np.sin(np.deg2rad(17))])
            fr = model.fr @ c
            ur = np.zeros(a['Kr'].shape[0])
            ur[model.free] = a['fac'].solve(fr[model.free])
            u = np.asarray(model.Q @ ur)
            strain = np.einsum('eij,ej->ei', a['B'], u[model.edof])
            s = np.einsum('eij,ej->ei', model.C, strain)[:model.nc]
            predicted = saved @ c
            scale = max(np.max(abs(s)), 1e-20)
            supererrs.append(float(np.max(abs(s - predicted)) / scale))
            faults['wrong_lateral_tensor_rejected'] = bool(np.max(abs(s - (saved[:, :, 0] * c[0] - saved[:, :, 1] * c[1]))) / scale > 1e-09)
    checks['regenerated_field_max_rel_error'] = max(errs)
    checks['direct_load_max_rel_error'] = max(supererrs)
    checks['max_equilibrium_residual'] = max(eq)
    checks['fields_pass'] = max(errs) <= 1e-09 and max(supererrs) <= 1e-09 and (max(eq) <= 1e-08)
    del models, model, a
    cover_checks = {}
    shared_predicate_checks = {}
    for support in ['8600', '18000', 'rigid']:
        ref = Response(z['reference__' + support + '__stress_basis'], z['reference__' + support + '__vol'])
        opt = Response(z['optimized__' + support + '__stress_basis'], z['optimized__' + support + '__vol'])
        cover = read(ROOT / f'raw/D_COVER_{support}.json')
        errors = []
        violations = []
        record_checks = []
        for (lo, _, a, b, l, h) in cover:
            recomputed = dual_cell(ref, opt, a, b, l, h)
            errors.append(abs(lo - recomputed))
            actual = ratio(ref, opt, (a + b) / 2, (l + h) / 2)
            violations.append(max(0, lo - actual))
            record_checks.append(bound_record_ok(lo, recomputed, actual))
        good = exact_cover(cover)
        shared_predicate_checks[support] = all(record_checks)
        cover_checks[support] = dict(exact_rational_cover=good, stored_bound_error=max(errors), midpoint_violation=max(violations), cells=len(cover))
        faults['missing_cover_cell_rejected_' + support] = not exact_cover(cover[:-1])
        (lo, _, a, b, l, h) = cover[0]
        actual = ratio(ref, opt, (a + b) / 2, (l + h) / 2)
        faults['inflated_bound_rejected_' + support] = not bound_record_ok(actual + 0.01, dual_cell(ref, opt, a, b, l, h), actual)
    checks['continuous_covers'] = cover_checks
    checks['covers_shared_predicate'] = shared_predicate_checks
    checks['covers_pass'] = all(shared_predicate_checks.values()) and all((v['exact_rational_cover'] and v['stored_bound_error'] <= 1e-12 and (v['midpoint_violation'] <= 1e-10) for v in cover_checks.values()))
    C = read(ROOT / 'raw/C2.json')
    lp_errors = []
    for row in C['rows']:
        support = row['support']
        theta = np.deg2rad(row['angle_deg'])
        c = np.array([np.cos(theta), np.sin(theta)])
        m = row['m']
        p = []
        for label in ['reference', 'optimized']:
            s = z[label + '__' + support + '__stress_basis'] @ c
            p.append(np.maximum(np.linalg.eigvalsh(stress_tensor(s))[:, -1], 0))
        scale = max(max(p[0]), max(p[1]))
        a = z['reference__' + support + '__vol'] * (p[0] / scale) ** m
        b = z['optimized__' + support + '__vol'] * (p[1] / scale) ** m
        d = (b - a) / sum(a)
        k = row['kappa']
        w = np.where(d > 0, k, 1)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sol = linprog(-1000000.0 * d, bounds=[(1, k)] * len(d), options=dict(threads=1, dual_feasibility_tolerance=1e-09, primal_feasibility_tolerance=1e-09), method='highs')
        lp_errors.append(abs(-sol.fun / 1000000.0 - w @ d))
        faults['uniform_witness_rejected'] = bool(d.sum() < 0 and np.where(d > 0, 1 + (k - 1) * 1.01, 1) @ d > 0)
    checks['LP_max_error'] = max(lp_errors)
    checks['LP_pass'] = max(lp_errors) <= 1e-09
    tet = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    (inv, B, vol, sgn) = geometry(tet, np.array([[0, 1, 2, 3]]))
    u = tet @ np.diag([1.0, 2.0, 3.0])
    eps = B[0] @ u.ravel()
    checks['affine_patch_error'] = float(np.max(abs(eps - np.array([1, 2, 3, 0, 0, 0]))))
    faults['wrong_strain_factor_rejected'] = not np.array_equal(2 * eps, np.array([1, 2, 3, 0, 0, 0]))

    def admissible(certificates):
        return all((v == 'PASS' for v in certificates.values()))
    faults['unknown_hard_condition_rejected'] = not admissible(dict(wall='PASS', milling='UNKNOWN'))
    faults['failed_hard_condition_rejected'] = not admissible(dict(wall='FAIL', milling='PASS'))
    source = (BASE / 'LANE_X1C_CROWN_LITERATURE/inputs/PMC8558575.txt').read_text()
    numbers = read(ROOT / 'sources/SOURCE_NUMBERS.json')
    checks['published_numeric_source'] = all((format(v, 'g') in source for v in numbers['m']))

    def source_value_ok(v):
        return v == numbers['characteristic_force_N'][6] and format(v, 'g') in source
    faults['source_force_x10_rejected'] = not source_value_ok(10 * numbers['characteristic_force_N'][6])
    checks['all_faults_reject'] = all(faults.values())
    passed = all([checks['locked_inputs'], checks['preregs_and_predictions'], checks['fields_pass'], checks['covers_pass'], checks['LP_pass'], checks['affine_patch_error'] == 0, checks['published_numeric_source'], checks['all_faults_reject']])
    result = dict(pass_gate=passed, checks=checks, faults=faults, seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, scope='Implementation replay; analytic/discrete model only, not physical validation')
    save(ROOT / 'raw/VERIFICATION.json', result)
    print(json.dumps(clean(result)), flush=True)
    if not passed:
        raise RuntimeError('Verification failed; retain results, do not qualify design')
if __name__ == '__main__':
    run()
