import json, time
import numpy as np
from contact_model import P, BASE, ENGINE, sha, write, ncp, solve, load_module, local_path
compare = load_module('x54_measurement_port', P / 'code/measurement_port.py').compare

def main():
    start = time.perf_counter()
    checks = []

    def check(name, passed, **evidence):
        checks.append(dict(name=name, pass_gate=bool(passed), **evidence))
    for path in sorted(P.glob('PREREG*.json')):
        check('prereg_' + path.name, sha(path) == path.with_name(path.name + '.sha256').read_text().strip())
    src = json.loads((P / 'SOURCE_MANIFEST.json').read_text())
    for s in src['sources']:
        check('source_' + s['path'], sha(s['path']) == s['sha256'])
    for name in ['FROZEN_PREDICTIONS.json', 'FROZEN_PREDICTIONS_R2.json', 'FROZEN_PREDICTIONS_R3.json', 'FROZEN_PREDICTIONS_R4.json']:
        fp = P / name
        check('freeze_' + name, sha(fp) == (P / (name + '.sha256')).read_text().strip())
        f = json.loads(fp.read_text())
        for item in f['manifest']:
            check('prediction_hash_' + item['path'], sha(item['path']) == item['sha256'])
            arr = item.get('arrays')
            if arr:
                check('array_hash_' + arr['path'], sha(arr['path']) == arr['sha256'])
    maximum = 0.0
    count = 0
    rejected = 0
    for cn in [0.002, 0.04]:
        for ct in [0.005, 0.07]:
            for bn in [-0.1, 0.1]:
                for bt in [0.0, 0.001, 0.1]:
                    for mu in [0.0, 0.2, 0.4]:
                        W = np.diag([cn, ct, ct])
                        b = np.array([bn, bt, 0])
                        m = np.array([mu])
                        ln = max(-bn / cn, 0.0)
                        lt = -min(bt / ct, mu * ln) if ln > 0 else 0.0
                        facit = np.array([ln, lt, 0.0])
                        (lam, _, _) = ncp.solve_ncp_pgs(W, b, m, iters=100, tol=1e-10)
                        error = float(np.max(np.abs(lam - facit)))
                        maximum = max(maximum, error)
                        count += 1
                        injected = lam.copy()
                        injected[0] += 1.0
                        rejected += ncp.natural_residual(injected, W, b, m) > 1e-06
    check('closed_form_contact', maximum <= 1e-06, count=count, max_error_N=maximum, external_referent=dict(kind='closed_form', locator='Signorini/Coulomb diagonal contact, arXiv2304.06372 eq14', compared_quantity='reaction N', refutes_us=True))
    check('one_contact_injected_wrong_force', rejected == count, rejected=int(rejected), attempted=count)
    for case in range(1, 13):
        s = json.loads((P / 'rounds/R1' / f'case{case:03d}_d0.05.json').read_text())
        a = np.load(local_path(s['arrays']['path']))
        res = ncp.natural_residual(a['lambda_N'], a['W'], a['b'], a['mu'])
        check('NCP_live_array_' + str(case), res <= 1e-06, residual_N=res)
    g = json.loads((P / 'raw/geometry_001_h02.json').read_text())
    (s, a) = solve(g, 0.05, controls=True)
    old = json.loads((P / 'rounds/R1/case001_d0.05.json').read_text())
    arr = np.load(local_path(old['arrays']['path']))
    err = float(np.max(np.abs(a['lambda_N'] - arr['lambda_N'])))
    check('cold_contact_replay', err <= 0.0001, max_force_error_N=err)
    for (k, v) in s['controls'].items():
        check('control_' + k, v.get('pass_gate', v.get('rejected', False)))
    pred = json.loads((P / 'rounds/R2/case001_force100.json').read_text())
    n = len(pred['rows'])
    S = np.eye(n)
    value = np.array([r['vector_on_upper_N'][2] for r in pred['rows']])
    rec = dict(kind='independent_measurement', locator='our_own_fixture:interface_only', applicability_evidence_locator='our_own_fixture:interface_only', case=1, geometry_sha256=sha(P / 'raw/geometry_001_h02.json'), force_target_N=100.0, quantity='regional_vertical_force_N', observation_matrix=S, values_N=value)
    check('measurement_interface_zero_error', compare(rec, pred, rec['geometry_sha256'])['gate'], scope='our_own_fixture , no physical reference')
    rec['values_N'] = value.copy()
    rec['values_N'][0] += 20.0
    check('measurement_interface_injected_20N', not compare(rec, pred, rec['geometry_sha256'])['gate'], scope='our_own_fixture')
    rec['quantity'] = 'normal_pressure_MPa'
    try:
        compare(rec, pred, rec['geometry_sha256'])
        wrong = False
    except ValueError:
        wrong = True
    check('wrong_observation_quantity_rejected', wrong)
    r3b = json.loads((P / 'rounds/R3B/results.json').read_text())
    check('R3B_continuous_control_repair', r3b['all_pass'])
    result = dict(checks=checks, all_pass=all((x['pass_gate'] for x in checks)), count=len(checks), wall_s=time.perf_counter() - start, scientific_scope='Numerical/interface checks only; original physical/refinement/profile failures retained. No independent physical measurement performed.')
    write(P / 'VERIFICATION.json', result)
    print('verification', result['count'], 'checks', result['all_pass'], 'closedform max', maximum, flush=True)
    if not result['all_pass']:
        raise SystemExit(1)
if __name__ == '__main__':
    main()
