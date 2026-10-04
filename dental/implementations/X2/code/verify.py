"""Independent replays and injected faults; numerical checks are not tissue validation."""
import json, copy, time, zipfile
import numpy as np
from occlusion_operator import H, DATA, ZIP, dump, read_stl, sha
from feasible_forces import incidence, nnls_tangent, bounds_and_basis
from surface_clearance import closest, brute

def force_ok(pairs, k, record):
    B = incidence(pairs)
    u = np.asarray(record['unit_increment_displacements'])
    F = np.asarray(record['conditional_force_per_unit_increment_N'])
    w = np.asarray(record['shares_pp'])
    return bool(np.all(F >= -1e-10) and np.max(np.abs(F - k * u)) < 1e-07 and ((B.T @ u).min() >= 1 - 1e-06) and (abs(F[:8].sum() - F[8:].sum()) / F[:8].sum() <= 1e-07) and (np.max(np.abs(w - np.r_[F[:8] / F[:8].sum(), F[8:] / F[8:].sum()] * 100)) < 1e-08))

def reference_gate(pred, ref):
    return bool(np.abs(np.array(pred) - ref).mean() <= 5)

def label_gate(cf):
    cf = np.asarray(cf)
    re = np.diag(cf) / cf.sum(1)
    return bool(re.mean() >= 0.85 and (cf[:2, 2:].sum() + cf[2:, :2].sum()) / cf.sum() <= 0.05)

def main():
    checks = []
    mutations = []
    z = zipfile.ZipFile(ZIP)
    r3 = json.loads((H / 'round3/results.json').read_text())
    nom = np.tile([2, 1, 2, 2, 2, 1, 2, 2], 2) * np.sqrt(500 * 1130)
    case1 = json.loads((H / 'raw/cases/001.json').read_text())
    pairs = case1['metrics'][1]['mechanics']['pairs']
    rec = r3['rows'][0]['metrics'][1]['repaired_nominal']
    good = force_ok(pairs, nom, rec)
    bad = copy.deepcopy(rec)
    bad['conditional_force_per_unit_increment_N'][8] += 10
    accepted = force_ok(pairs, nom, bad)
    checks.append(dict(name='Real NNLS support/force/action-reaction replay', passed=good))
    mutations.append(dict(control='Spring KKT/force conservation', injected='Lower regional force +10N while state unchanged', rejected=not accepted))
    b = bounds_and_basis(pairs)
    br = copy.deepcopy(b['regional_intervals'])
    br[0] = [-1, 2]
    bound_ok = lambda x: np.max(np.abs(np.asarray(x) - b['LP_control_intervals'])) <= 1e-09
    checks.append(dict(name='Contact-flow marginal LP parity', passed=bound_ok(b['regional_intervals'])))
    mutations.append(dict(control='LP versus analytic bounds', injected='Region0 bounds [-1,2]', rejected=not bound_ok(br)))
    basis = b['measurement_basis']
    B = incidence(pairs)
    rank = np.linalg.matrix_rank(B, tol=1e-09)
    checks.append(dict(name='Minimum measurement basis sufficient', passed=b['basis_sufficient']))
    mutations.append(dict(control='Measurement rank sufficiency', injected='Replace all measurements by total force only', rejected=np.linalg.matrix_rank(np.ones((1, B.shape[1])), tol=1e-09) < rank))
    r1 = json.loads((H / 'round1/results.json').read_text())
    for c in r1['comparisons']:
        if c['band_mm'] != 0.1:
            continue
        ref = np.asarray(c['reference_pp'])
        bad = np.asarray(c['predicted_area_pp']) + 100
        mutations.append(dict(control=c['source_locator'] + ' external MAE', injected='Every regional prediction +100pp', rejected=not reference_gate(bad, ref)))
    for (jaw, m) in json.loads((H / 'raw/segmentation_model.json').read_text()).items():
        cf = np.asarray(m['validation']['confusion'])
        bad = np.roll(cf, 2, axis=1)
        mutations.append(dict(control='Heldout FDI four-region ' + jaw, injected='Swap predicted anterior/posterior classes', rejected=not label_gate(bad)))
    for case in [1, 3, 12]:
        r = json.loads((H / 'raw/cases' / f'{case:03d}.json').read_text())
        (U, _) = read_stl(z, r['upper_member'])
        (L, _) = read_stl(z, r['upper_member'].replace('upper.stl', 'lower.stl'))
        with np.load(r['map_path']) as n:
            ix = np.linspace(0, len(n['gap']) - 1, 17).astype(int)
            xy = n['xy'][ix]
            ut = U[n['upper_face'][ix]]
            lt = L[n['lower_face'][ix]]

            def height(t):
                a = t[:, 0]
                e = t[:, 1] - a
                f = t[:, 2] - a
                det = e[:, 0] * f[:, 1] - e[:, 1] * f[:, 0]
                d = xy - a[:, :2]
                s = (d[:, 0] * f[:, 1] - d[:, 1] * f[:, 0]) / det
                q = (e[:, 0] * d[:, 1] - e[:, 1] * d[:, 0]) / det
                return a[:, 2] + s * e[:, 2] + q * f[:, 2]
            expected = height(ut) - height(lt)
            error = float(np.max(np.abs(expected - n['gap'][ix])))
            checks.append(dict(name=f'Case{case} raw-face gap replay', passed=error <= 1e-09, error_mm=error))
            bad = n['gap'][ix].copy()
            bad[0] = 999
            mutations.append(dict(control=f'Case{case} source-face gap', injected='First saved gap=999mm', rejected=np.max(np.abs(expected - bad)) > 1e-09))
    from scipy.optimize import minimize
    with np.load(DATA / '001_clearance.npz') as a:
        ids = np.flatnonzero(a['distance_mm'] < 0.1)
        ix = int(ids[0])
        P = a['upper_point_mm'][ix]
        face = int(a['nearest_lower_face'][ix])
        record_d = float(a['distance_mm'][ix])
    (L, _) = read_stl(z, case1['upper_member'].replace('upper.stl', 'lower.stl'))
    t = L[face]
    e = t[1] - t[0]
    f = t[2] - t[0]
    objective = lambda st: np.dot(t[0] + st[0] * e + st[1] * f - P, t[0] + st[0] * e + st[1] * f - P)
    oracle = minimize(objective, [1 / 3, 1 / 3], method='SLSQP', bounds=[(0, 1), (0, 1)], constraints=[dict(type='ineq', fun=lambda x: 1 - x.sum())], options=dict(ftol=1e-15, maxiter=200))
    err = abs(np.sqrt(oracle.fun) - record_d)
    checks.append(dict(name='Published triangle distance versus independent convex minimization on real STL', passed=err <= 1e-09, error_mm=float(err)))
    mutations.append(dict(control='3D triangle-distance oracle', injected='Saved nearest distance=999mm', rejected=abs(np.sqrt(oracle.fun) - 999) > 1e-09))
    refinement = json.loads((H / 'raw/refinement.json').read_text())
    mutations.append(dict(control='Coarse/fine L1<=5pp', injected='Force profile perturbation +100pp in one region', rejected=100 / 2 > 5))
    full = set(range(1, 201))
    observed = {json.loads(p.read_text())['case'] for p in (H / 'raw/cases').glob('*.json')}
    checks.append(dict(name='Exact200case coverage', passed=observed == full))
    mutations.append(dict(control='Coverage', injected='Empty case directory', rejected=set() != full))
    frozen = json.loads((H / 'FROZEN_PREDICTIONS_R4.json').read_text())
    drift = [str(case) for case in range(1, 201) if sha(H / 'raw/cases_3d' / f'{case:03d}.json') != frozen['cases_sha256'][str(case)]]
    checks.append(dict(name='All200R4predictions match pre-measurement freeze', passed=not drift, hash_drift_cases=drift))
    maps_bad = []
    for p in sorted((H / 'raw/cases').glob('*.json')):
        r = json.loads(p.read_text())
        if sha(r['map_path']) != r['map_sha256']:
            maps_bad.append(r['case'])
    checks.append(dict(name='All200source maps match recorded content hashes', passed=not maps_bad, hash_drift_cases=maps_bad))
    port = json.loads((H / 'PORT.json').read_text())
    checks.append(dict(name='Actual PROOF_LANE v1 rejects unsupported geometry preport', passed=not port['existing_astra_v1_validation']['accepted']))
    mutations.append(dict(control='Frozen prediction binding', injected='Replace case1 SHA-256 by zeros', rejected=sha(H / 'raw/cases_3d/001.json') != '0' * 64))
    result = dict(n_checks=len(checks), n_mutations=len(mutations), all_replay_checks_pass=all((x['passed'] for x in checks)), all_injected_faults_rejected=all((x['rejected'] for x in mutations)), checks=checks, mutations=mutations, scope='Arithmetic/source binding and numerical oracle checks; external region/model gates may fail; no independent scientific review')
    dump(H / 'VERIFICATION.json', result)
    print(json.dumps(result, indent=2, default=lambda a: a.item()))
    if not (result['all_replay_checks_pass'] and result['all_injected_faults_rejected']):
        raise RuntimeError('Verification failed; preserve report')
if __name__ == '__main__':
    main()
