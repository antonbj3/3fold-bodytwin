from common import *
import time
from scipy.optimize import minimize

def full_control(C, f0, h, A, w):
    sol = minimize(lambda f: 0.5 * (f - f0) @ C @ (f - f0) - h @ f, f0, jac=lambda f: C @ (f - f0) - h, bounds=[(0, None)] * len(f0), constraints=[dict(type='eq', fun=lambda f: A.T @ f - w, jac=lambda f: A.T)], method='SLSQP', options=dict(ftol=1e-13, maxiter=500))
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.x

def run():
    t = time.perf_counter()
    u = np.array([1.0, 1.0, -2.0])
    v = np.array([1.0, -1.0, 0.0])
    f0 = np.array([8.0, 8.0, 16.0])
    h = np.array([1 / 64, 0.0, 0.0])
    hh = np.array([1 / 64, 1 / 64, 0.0])
    ans = []
    summ = []
    errors = []
    for b in [64.0, 128.0]:
        K = 64 * np.outer(u, u) + b * np.outer(v, v)
        ff = f0 + K @ h
        C = np.linalg.pinv(K) + np.ones((3, 3)) / 3 / 256
        control = full_control(C, f0, h, np.ones((3, 1)), np.array([32.0]))
        ans.append(ff)
        summ.append([float(f0[:2].sum()), *(K @ hh).tolist()])
        errors.append(float(np.max(abs(control - ff))))
    r0 = dict(summary_identity_error=float(np.max(abs(np.array(summ[0]) - summ[1]))), bitidentical_summary=np.array_equal(summ[0], summ[1]), summaries=summ, regional_edit_mm=h, forces_N=ans, downstream_max_region_force_difference_N=float(np.max(abs(ans[0] - ans[1]))), control_error_N=errors, status='TOOTH_SUMMARY_INSUFFICIENT', smallest_extension='Response for the single differential height mode of the two regions, and region-resolved force observation', resolution='PER_SURFACE_REGION', external_force_validation=False)
    dump(ROOT / 'raw/R0_SUFFICIENCY.json', r0)
    hf = x82()
    rows = []
    for p in read(ROOT / 'PREREG_R1.json')['cohort']:
        s = read(p)
        fdi = s['predicted_fdi']
        i = fdi.index(16)
        f0 = np.array(s['baseline_force_N'])
        target = float((f0[fdi.index(15)] + f0[fdi.index(17)]) / 2)
        out = hf.inverse(s, 16, target, 1.0)
        r = dict(case=s['case'], fdi=16, target_N=target, target_band_N=[target - 1, target + 1], baseline_force_N=float(f0[i]), inverse=out, region_transfer='UNKNOWN_NO_REGION_CALIBRATION', resolution='PER_TOOTH')
        if out.get('at_selected_height'):
            (A, N, H, f0, w, hb, rho, eta) = hf.validate(s)
            G = np.linalg.inv(H)
            C = N @ G @ N.T + (np.eye(len(f0)) - N @ N.T) / 750
            h = np.array(out['at_selected_height']['height_change_mm'])
            ff = full_control(C, f0, h, A, w)
            r['control_error_N'] = float(np.max(abs(ff - out['at_selected_height']['force_N'])))
            (lo, hi) = out['at_selected_height']['force_interval_N'][i]
            r['injected_plus3_N_rejected'] = not lo <= ff[i] + 3 <= hi
        rows.append(r)
    ref = read(X82 / 'inputs/F5367_R4_state.json')
    example = hf.predict(ref, hf.make_edit(ref, ref['predicted_fdi'].index(16), -0.02))
    result = dict(rows=rows, conditional_target_passes=sum((r['inverse']['status'] == 'CONDITIONAL_TARGET_BAND' for r in rows)), requested=len(rows), example_16_minus20um=example, physical_measurement='NOT_RUN', geometry_transfer_passes=0, seconds=time.perf_counter() - t)
    dump(ROOT / 'raw/R1_UNIFORM_INVERSE.json', result)
    state('R0_R1_DECIDED', 'Tooth summary fails exact regional sufficiency; uniform inverse rerun', 'Read external measurement facit, then freeze regional actuation and geometry certificates')
    (ROOT / 'history/HANDOFF_R0_R1.md').write_text('R0: identical tooth response and baseline, exact identity error0; regional edit changes max region force1N. Add the differential regional response column and spatial force observation. R1 reuses X82 inverse; see raw/R1_UNIFORM_INVERSE.json. No matched crown-force binding is available. Next: region-labelled occlusal roof basis with explicit calibration/error contract and continuous geometry gates.\n')
    dump(ROOT / 'ATTEMPTS.json', [dict(round='R0', outcome=r0['status'], evidence='raw/R0_SUFFICIENCY.json', next_operation=r0['smallest_extension']), dict(round='R1', outcome='UNIFORM_INVERSE_REUSED_REGIONAL_TRANSFER_UNKNOWN', evidence='raw/R1_UNIFORM_INVERSE.json', next_operation='Calibrate regional actuation basis on actual crown/fixture')])
    print(json.dumps(dict(R0_difference_N=r0['downstream_max_region_force_difference_N'], R1_passes=result['conditional_target_passes'], requested=result['requested'])))
if __name__ == '__main__':
    run()
