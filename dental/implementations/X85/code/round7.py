from common import *
from scipy.optimize import linprog
from regional import predict
from round2 import load_item, simulated_calibration
from round01 import full_control
from geometry import basis, roof_checks, ball_access, shell, export_stl
import time

def inverse_signed(model, B, weights, allow, extra_A=None, extra_b=None):
    J = np.array(model['J_N_per_mm'])
    f0 = np.array(model['f0'])
    eta = model['eta_N']
    beta = np.array(model['beta_N_per_mm'])
    D = np.c_[J, -J]
    br = np.r_[beta, beta]
    G = []
    b = []
    for i in range(2):
        G.extend([D[i] + br, -D[i] + br])
        b.extend([9 - f0[i] - eta, f0[i] - eta - 7])
    G.extend(-D + br)
    b.extend(f0 - eta - 1e-06)
    G.extend(np.c_[-B, B])
    b.extend(allow)
    if extra_A is not None:
        G.extend(np.c_[extra_A, -extra_A])
        b.extend(extra_b)
    sol = linprog(np.r_[weights @ B, weights @ B], A_ub=np.array(G), b_ub=np.array(b), bounds=[(0, 0.05)] * 4, method='highs', options={'threads': 1})
    if not sol.success:
        return dict(status='INFEASIBLE_OR_UNKNOWN', reason=sol.message)
    a = sol.x[:2] - sol.x[2:]
    p = predict(model, a)
    iv = np.array(p['force_interval_N'])
    good = p['positive_contact_path'] and np.all(iv[:2, 0] >= 7 - 1e-07) and np.all(iv[:2, 1] <= 9 + 1e-07)
    return dict(status='CONDITIONAL_SIGNED_TARGET' if good else 'UNKNOWN', coefficients_mm=a, at_selected_height=p, positive_height_requires_new_manufacture=bool(np.any(a > 0)), mean_absolute_design_change_mm=float(sol.fun / sum(weights)))

def run():
    start = time.perf_counter()
    worlds = read(ROOT / 'raw/R6_WORLDS.json')
    r6p = read(ROOT / 'FROZEN_PREDICTIONS_R6.json')['predictions']
    (t, _, _) = load_item(read(ROOT / 'PREREG_R4.json')['cohort'][0])
    B = basis(t['xy'])
    inner = t['preparation_z'] + t['requirements']['film_min_mm']
    z = inner + t['requirements']['wall_mm'] + 0.05
    (_, _, A, C, P, f0) = simulated_calibration(t, z, inner, B)
    answers = []
    for w in worlds:
        u = np.array(w['calibration']['final_coefficient_error_mm'])
        ans = inverse_signed(w['calibration'], B, t['weights'], z - inner - t['requirements']['wall_mm'] - B @ u)
        answers.append(ans)
    if not (ROOT / 'FROZEN_PREDICTIONS_R7.json').exists():
        freeze(ROOT / 'FROZEN_PREDICTIONS_R7.json', dict(answers=answers, independent_QP_queries_before_freeze=0, physical_measurement='NOT_RUN'))
    from review_integrity import assert_frozen_payload
    assert_frozen_payload(ROOT, 'FROZEN_PREDICTIONS_R7.json', 'answers', clean(answers))
    checks = []
    q = v6_contact()
    for (i, (w, ans)) in enumerate(zip(worlds, answers)):
        if ans['status'] != 'CONDITIONAL_SIGNED_TARGET':
            checks.append(dict(world=i, enclosed=False, injected_plus3_N_rejected=False, status=ans['status']))
            continue
        a = np.array(ans['coefficients_mm'])
        truea = a + np.array(w['final_error_mm'])
        truth = full_control(C, f0, P @ truea, A, A.T @ f0)
        iv = np.array(ans['at_selected_height']['force_interval_N'])
        checks.append(dict(world=i, enclosed=bool(np.all(truth >= iv[:, 0] - 1e-07) and np.all(truth <= iv[:, 1] + 1e-07)), injected_plus3_N_rejected=bool(truth[0] + 3 > iv[0, 1]), actual_force_N=truth, positive_height_requires_new_manufacture=ans['positive_height_requires_new_manufacture']))
    first = answers[0]
    representative = None
    if first['status'] == 'CONDITIONAL_SIGNED_TARGET':
        final = z + B @ np.array(first['coefficients_mm'])
        geom = roof_checks(t, final, inner)
        tool = ball_access(t['xy'], t['faces'], final)
        with np.load(V4 / 'payload/private/references' / (t['case_key'] + '.npz'), allow_pickle=False) as ref:
            rg = t['ceiling'] - ref[t['family']]
        representative = dict(inverse=first, geometry=geom, milling=tool, source_contact=q.compare(t['xy'], t['faces'], t['ceiling'] - final, rg), note='One declared spring fixture on a source geometry; no measured loaded-gap/force binding')
        from plugin import generate
        cal = worlds[0]['calibration']
        task = dict(t, height_operation='signed_generation', regional_calibration=cal, actuation_basis=B, unadjusted_z=z, inner_z=inner, geometry_sha256=cal['geometry_sha256'], basis_sha256=cal['basis_sha256'])
        design = generate(task)
        representative['plugin_status'] = design['status']
        representative['plugin_diagnostics'] = design.get('diagnostics')
        if design['status'] == 'DESIGN':
            dest = ROOT / 'exports/R7_signed_generation'
            dest.mkdir(parents=True, exist_ok=True)
            (v, f) = shell(t['xy'], t['faces'], np.array(design['outer_vertices'])[:, 2], inner)
            export_stl(dest / 'signed_roof.stl', v, f)
            dump(dest / 'GenCAD_design.json', design)
            representative['export'] = dict(path=str(dest / 'signed_roof.stl'), sha256=sha(dest / 'signed_roof.stl'), clinical_eligibility='UNKNOWN')
    summary = dict(worlds=len(checks), targets_obtained=sum((a['status'] == 'CONDITIONAL_SIGNED_TARGET' for a in answers)), enclosed=sum((c['enclosed'] for c in checks)), faults_rejected=sum((c['injected_plus3_N_rejected'] for c in checks)), positive_height_worlds=sum((c.get('positive_height_requires_new_manufacture', False) for c in checks)), seconds=time.perf_counter() - start, status='CONDITIONAL_SIGNED_DESIGN_PASS' if all((c['enclosed'] and c['injected_plus3_N_rejected'] for c in checks)) else 'FAIL')
    dump(ROOT / 'raw/R7_CONTROLS.json', checks)
    dump(ROOT / 'raw/R7_SUMMARY.json', summary)
    dump(ROOT / 'raw/R7_REPRESENTATIVE.json', representative)
    state('R7_DECIDED', str(summary), 'Package negative removal-only result and signed-generation recovery; acquire actual specimen calibration next')
    print(json.dumps(summary))
if __name__ == '__main__':
    run()
