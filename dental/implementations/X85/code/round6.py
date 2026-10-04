from common import *
from height_uncertainty import fit_bounded_heights
from regional import inverse
from round2 import load_item, simulated_calibration
from geometry import basis
from round01 import full_control
import time, copy

def run():
    start = time.perf_counter()
    pr = read(ROOT / 'PREREG_R6.json')
    p = pr['parameters']
    (t, _, _) = load_item(read(ROOT / 'PREREG_R4.json')['cohort'][0])
    B = basis(t['xy'])
    inner = t['preparation_z'] + t['requirements']['film_min_mm']
    z = inner + t['requirements']['wall_mm'] + 0.05
    (model, raw, A, C, P, f0) = simulated_calibration(t, z, inner, B)
    J = np.array(model['J_N_per_mm'])
    rng = np.random.default_rng(p['seed'])
    pred = []
    worlds = []
    u = np.full(2, p['final_coefficient_error_mm'])
    eps = p['force_channel_error_N']
    delta = p['probe_height_error_mm']
    s = p['probe_halfstep_mm']
    for i in range(p['worlds']):
        probes = clean(copy.deepcopy(raw))
        dp = rng.uniform(-delta, delta, 2)
        dm = rng.uniform(-delta, delta, 2)
        e0 = rng.uniform(-eps, eps, 6)
        ep = rng.uniform(-eps, eps, (2, 6))
        em = rng.uniform(-eps, eps, (2, 6))
        probes.update(height_error_mm=delta, baseline_force_N=(f0 + e0).tolist(), plus_force_N=np.array([f0 + J[:, j] * (s + dp[j]) + ep[j] for j in range(2)]).tolist(), minus_force_N=np.array([f0 + J[:, j] * (-s + dm[j]) + em[j] for j in range(2)]).tolist())
        cal = fit_bounded_heights(probes, u)
        ans = inverse(cal, B, t['weights'], z - inner - t['requirements']['wall_mm'] - B @ u)
        a = np.array(ans['coefficients_mm']) if ans['status'] == 'CONDITIONAL_TARGET_DESIGN' else np.zeros(2)
        du = rng.uniform(-u, u)
        worlds.append(dict(probes=probes, final_error_mm=du, calibration=cal))
        pred.append(dict(world=i, inverse=ans, actual_coefficients_mm=a + du))
    if not (ROOT / 'FROZEN_PREDICTIONS_R6.json').exists():
        freeze(ROOT / 'FROZEN_PREDICTIONS_R6.json', dict(predictions=pred, heldout_full_force_QP_queries_before_freeze=0, physical_measurement='NOT_RUN'))
    from review_integrity import assert_frozen_payload
    assert_frozen_payload(ROOT, 'FROZEN_PREDICTIONS_R6.json', 'predictions', clean(pred))
    checks = []
    for (i, (world, prediction)) in enumerate(zip(worlds, pred)):
        ans = prediction['inverse']
        truth = full_control(C, f0, P @ np.array(prediction['actual_coefficients_mm']), A, A.T @ f0)
        if ans['status'] == 'CONDITIONAL_TARGET_DESIGN':
            iv = np.array(ans['at_selected_height']['force_interval_N'])
            checks.append(dict(world=i, enclosed=bool(np.all(truth >= iv[:, 0] - 1e-07) and np.all(truth <= iv[:, 1] + 1e-07)), injected_plus3_N_rejected=bool(truth[0] + 3 > iv[0, 1]), maximum_actual_force_deviation_N=float(np.max(abs(truth - ans['at_selected_height']['force_N']))), joint_radius_N=ans['at_selected_height']['joint_l2_error_radius_N']))
        else:
            checks.append(dict(world=i, enclosed=False, injected_plus3_N_rejected=False, inverse_status=ans['status']))
    bad = clean(copy.deepcopy(raw))
    bad['height_error_mm'] = s
    try:
        fit_bounded_heights(bad, u)
        rejected = False
    except ValueError:
        rejected = True
    summary = dict(worlds=len(checks), enclosed=sum((r['enclosed'] for r in checks)), faults_rejected=sum((r['injected_plus3_N_rejected'] for r in checks)), zero_denominator_rejected=rejected, maximum_joint_radius_N=max((r.get('joint_radius_N', 0) for r in checks)), seconds=time.perf_counter() - start, status='CONDITIONAL_BOUNDED_HEIGHTS_PASS' if all((r['enclosed'] and r['injected_plus3_N_rejected'] for r in checks)) and rejected else 'FAIL')
    dump(ROOT / 'raw/R6_WORLDS.json', worlds)
    dump(ROOT / 'raw/R6_CONTROLS.json', checks)
    dump(ROOT / 'raw/R6_SUMMARY.json', summary)
    state('R6_DECIDED', str(summary), 'Freeze actual lab predictions only after receiving calibrated same-specimen proberuns; no physical data currently')
    print(json.dumps(summary))
if __name__ == '__main__':
    run()
