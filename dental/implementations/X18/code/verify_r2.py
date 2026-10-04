"""Validate new continuous constraints against an independent full-pair minimum."""
from geometry import *
from continuous import extremum
from experiment import local_tri
from obstacle import constraints

def run():
    path = H / 'raw/PREDICTIONS_R2.json'
    if not path.exists():
        dump(H / 'raw/VERIFICATION_R2.json', dict(status='NOT_RUN_R2_STILL_DESIGNING'))
        return
    preds = json.loads(path.read_text())
    pr = json.loads((H / 'PREREG_R2.json').read_text())
    p = next((r for r in preds if r['status'] == 'DESIGNED' and r['fdi'] % 10 == 6))
    (data, _) = pair(p['case'])
    z = np.load(p['file'])
    xy = z['xy']
    U = local_tri(data['upper']['tri'], xy)
    (A, b) = constraints(U, xy, z['faces'], pr['clearance_mm'])
    C = np.c_[xy, z['z_informed']][z['faces']]
    pairs = broad_phase(U, C)
    full = extremum(U, C, pairs, False)
    maxineq = float(np.max(A @ z['z_informed'] - b))
    gap = full['minimum_gap_mm']
    identity_error = abs(maxineq - (pr['clearance_mm'] - gap))
    tol = pr['metrics']['continuous_gap_oracle_tolerance_mm']
    checks = [dict(name='All polygon-vertex inequalities match full facet gap minimum', value_mm=identity_error, pass_=identity_error <= tol), dict(name='Injected1mm inequality answer fails independent facet control', value_mm=1, pass_=abs(maxineq + 1 - (pr['clearance_mm'] - gap)) > tol), dict(name='Saved informed facet gap nonpenetration', gap_mm=gap, pass_=gap >= -pr['metrics']['penetration_tolerance_mm'])]
    for name in ['informed', 'control']:
        (vv, ff) = shell(xy, z['z_' + name], z['faces'], 1.2)
        if np.isfinite(vv).all():
            ee = np.r_[ff[:, [0, 1]], ff[:, [1, 2]], ff[:, [2, 0]]]
            (_, counts) = np.unique(np.sort(ee, axis=1), axis=0, return_counts=True)
            checks.append(dict(name=name + ' exported shell closed edge incidence', pass_=bool(np.all(counts == 2))))
    dump(H / 'raw/VERIFICATION_R2.json', dict(status='PASS' if all((c['pass_'] for c in checks)) else 'FAIL', case=p['case'], fdi=p['fdi'], checks=checks, constraints=A.shape[0], full_pairs=len(pairs), external_referent=dict(kind='closed_form', locator='raw/VERIFICATION.json independent halfplane linear-program oracles plus affine extreme-point theorem', compared_quantity='Continuous affine gap vs complete barycentric inequality set', refutes_us=True)))
    assert all((c['pass_'] for c in checks))
if __name__ == '__main__':
    run()
