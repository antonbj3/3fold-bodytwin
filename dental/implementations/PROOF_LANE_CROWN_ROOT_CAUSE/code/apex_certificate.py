from diagnose import *
import sympy as sp
import warnings, datetime
from scipy.optimize import linprog

def certify(m):
    t = m['vertices'][m['faces'][m['roles'] == 1]]
    raw = -np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    norm = raw / np.linalg.norm(raw, axis=1)[:, None]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        lp = linprog(np.zeros(len(norm)), A_eq=np.r_[norm.T, np.ones((1, len(norm)))], b_eq=[0, 0, 0, 1], bounds=(0, None), method='highs', options={'threads': 1})
        cheby = linprog([0, 0, 0, 1], A_ub=np.r_[np.c_[norm, -np.ones(len(norm))], np.c_[-norm, -np.ones(len(norm))]], b_ub=np.r_[np.full(len(norm), 0.05), np.full(len(norm), -0.05)], bounds=[(None, None)] * 3 + [(0, None)], method='highs', options={'threads': 1})
    if not lp.success:
        return dict(status='UNKNOWN_NO_POSITIVE_DEPENDENCE', lp=lp.message, minimax=cheby.fun if cheby.success else None)
    ids = np.flatnonzero(lp.x > 1e-08)
    normals = []
    for i in ids:
        pts = [[sp.Rational(float(c)) for c in p] for p in t[i]]
        (a, b, c) = map(sp.Matrix, pts)
        nn = -(b - a).cross(c - a)
        normals.append(nn)
    M = sp.Matrix.hstack(*normals)
    cert = None
    for z in M.nullspace():
        if all((x <= 0 for x in z)):
            z = -z
        if all((x >= 0 for x in z)) and any((x > 0 for x in z)):
            assert M * z == sp.zeros(3, 1)
            cert = dict(triangle_ids=ids.tolist(), triangles=t[ids].tolist(), raw_normals=[[str(c) for c in n] for n in normals], weights=[str(c) for c in z], exact_zero_vector=True)
            break
    if cert is None:
        return dict(status='UNKNOWN_EXACT_RECONSTRUCTION', support=ids.tolist())
    bad = [F(x) for x in cert['weights']]
    bad[0] += 1
    ns = [[F(x) for x in n] for n in cert['raw_normals']]
    corruption = any((sum((y * n[j] for (y, n) in zip(bad, ns))) != 0 for j in range(3)))
    return dict(status='PROVED_NO_SINGLE_APEX_OFFSET_WITHIN10UM', certificate=cert, exact_minimax_residual_mm='1/20', strong_control_minimax_mm=float(cheby.fun), minimax_error_mm=float(abs(cheby.fun - 0.05)), fault_injected_bad_weight_rejected=corruption, proof='lambda>=0 and sum(lambda_i*n_i)=0 with nonzero normals imply not all n_i dot delta>0. Every0.05+/-0.01mm unit-normal displacement would be strictly positive. Delta0 attains max residual0.05, so exact minimax is0.05.', resolution='PER_POINT')

def run():
    p = R / 'PREREG_A2.json'
    if not p.exists():
        dump(p, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), claim_type='capability', capability='Certify whether changing least-squares solver could save one-apex offsets; extend fixed-exterior wall witness to18 teeth', obstacle='Least-squares failure alone is not proof that another objective fails', changed_operation='Positive dependence/Farkas certificate of exact facet cross products; independent Chebyshev LP', consumer='Reject solver tuning when topology must change', metrics=dict(max_allowed_residual_mm=0.01, uniform_apex_gap_mm=0.05, LP_parity_mm=1e-08, wall_min_mm=0.5), strongest_equally_informed_control='Explicit same-information Chebyshev LP', falsifiers=['LP residual<=.01', 'No exact positive dependence'], external_referent=dict(kind='closed_form', locator=str(B / 'LANE_NEXT_D_INSERTION_PROOF/code/cone.py'), compared_quantity='Farkas positive dependence, extended here to exact vertex cross products', refutes_us=True), full_cost=dict(preparation='Reuse18 meshes', fit='None', discovery='Exact certificates and witnesses timed', validation='LP and injected multiplier error', queries=18, fallback='UNKNOWN on failed reconstruction', historical='UNKNOWN')))
        p.with_suffix('.sha256').write_text(hashlib.sha256(p.read_bytes()).hexdigest() + '  ' + p.name + '\n')
        dump(R / 'DECOMPOSITION_A2.json', dict(idea='A positive dependence among facet normals prevents a common outward apex displacement', equation='sum(lambda_i*n_i)=0, lambda>=0 -> not(all n_i dot d>0)', operation='LP support discovery; rational nullspace; exact checker', representation='Stored binary64 vertices lifted to rationals', leaves=[dict(status='DERIVED_UNDER_ASSUMPTIONS', statement='Positive-dependence and minimax proof', stopping_argument='Exact finite surface model; clinical geometry is outside scope'), dict(status='UNKNOWN', statement='Scan uncertainty', stopping_argument='No scanner calibration')]))
    start = time.perf_counter()
    rs = read(B / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_PREDICTIONS_C.json')['rows']
    inp = {r['key']: r for r in read(B / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json')['records']}
    rows = []
    for rec in rs:
        if rec['method'] != 'exact_margin':
            continue
        m = load(rec['mesh_path'])
        (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
        native = load(inp[rec['key']]['private_path'])['target']
        prep = v[f[roles == 1]]
        out = dict(key=rec['key'], family=rec['family'], apex=certify(m), native_wall_witness=pair_witness(native, prep), generated_wall_witness=pair_witness(v[f[roles == 0]], prep))
        rows.append(out)
        dump(R / 'raw/A2_CERTIFICATES.json', rows)
        print(rec['key'], out['apex']['status'], out['native_wall_witness']['exact']['distance_squared_mm2'] == '0', flush=True)
    dump(R / 'RESULTS_A2.json', dict(claim_type='capability', rows=rows, seconds=time.perf_counter() - start, external_referent=read(p)['external_referent']))
if __name__ == '__main__':
    run()
