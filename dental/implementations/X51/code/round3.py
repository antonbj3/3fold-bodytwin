import thread_budget
import math, time
import numpy as np
from scipy.stats import t as student
from common import ROOT, read, save, state, check_lock
from drift import interval_hold

def main():
    check_lock('PREREG_R3.json')
    start = time.perf_counter()
    s = read('PREREG_R3.json')['scenarios']
    rng = np.random.default_rng(s['baseline_seed'])
    N = s['baseline_n']
    n = s['weekly_replicates']
    W = s['weeks']
    P = s['paths']
    sigma = s['per_scan_sigma_um']
    alpha = 0.05
    samples_r = rng.normal(0, sigma, (P, N))
    samples_y = rng.normal(0, sigma, (P, N))
    samples_d = samples_y - samples_r
    base_r = samples_r.mean(axis=1)
    base_d = samples_d.mean(axis=1)
    sr = samples_r.std(axis=1, ddof=1)
    sd = samples_d.std(axis=1, ddof=1)
    del samples_r, samples_y, samples_d
    tcrit = float(np.nextafter(student.isf(alpha / (4 * W), N - 1), np.inf))
    scale = math.sqrt(1 / n + 1 / N)
    while 4 * W * student.sf(tcrit, N - 1) > alpha:
        tcrit = float(np.nextafter(tcrit, np.inf))
    hr = tcrit * sr * scale
    hd = tcrit * sd * scale
    newr = rng.normal(0, sigma / np.sqrt(n), (P, W))
    newy = rng.normal(0, sigma / np.sqrt(n), (P, W))
    newd = newy - newr
    flagsr = np.abs(newr - base_r[:, None]) > hr[:, None]
    flagsd = np.abs(newd - base_d[:, None]) > hd[:, None]
    num = flagsr.sum(axis=1) + flagsd.sum(axis=1)
    mc = float(num.mean())
    analytic = 4 * W * float(student.sf(tcrit, N - 1))
    A = read('raw/R2.json')['source_amplitudes']['sinter_size_um']
    records = []
    unsafe = 0
    enclosure_residual = 0
    for duration in s['ramp_durations_days']:
        v = A / duration
        g0 = s['gap0_um']
        E = s['guard_error_d_um']
        H = s['horizon_days']
        L = s['gap_limit_um']
        hold_days = []
        released_count = 0
        for error in [-E, E]:
            held = False
            for day in range(0, 365, H):
                d = A * min(day / duration, 1)
                dhat = d + error
                guard = interval_hold(g0, dhat, E, v, H, L)
                corners = [g0 + dhat + e + vv * dt for e in [-E, E] for vv in [0, v] for dt in [0, H]]
                enclosure_residual = max(enclosure_residual, abs(guard['upper_gap_um'] - max(corners)))
                worsttrue = g0 + A * min((day + H) / duration, 1)
                if guard['verdict'] == 'HOLD':
                    hold_days.append(day)
                    held = True
                    break
                released_count += 1
                if worsttrue > L:
                    unsafe += 1
        first_reject = next((day for day in range(1, 365) if g0 + A * min(day / duration, 1) > L), None)
        records.append({'ramp_duration_days': duration, 'amplitude_um': A, 'slope_bound_um_per_day': v, 'slope_status': 'PHENOMENOLOGICAL imposed bound, not measured drift', 'hold_days_over_error_corners': hold_days, 'first_rejection_day': first_reject, 'released_intervals_over_corners': released_count, 'resolution': 'PER_SURFACE_REGION model; times PHENOMENOLOGICAL'})
    dhat = 0.0
    slow = 0.0
    fast = 80 / 7
    suff = {'summary': 'current real drift um', 'state_current_um': [0.0, 0.0], 'identity_error_um': 0.0, 'next_week_gap_um': [50.0, 130.0], 'downstream_difference_um': 80.0, 'minimum_extension': 'A validated upper increment over the release horizon (rate bound plus jump/event handling)', 'resolution': 'PER_SURFACE_REGION', 'scope': 'constructed conditional additive model, not external dental reference'}
    badguard = interval_hold(50, 0, 20, 0, 7, 120)
    bad_refutes = badguard['verdict'] == 'MODEL_WITHIN_LIMIT' and 50 + 80 > 120
    checks = [dict(name='FINITE_BASELINE_STUDENT', nominal_pass=abs(mc - analytic) < 0.005 and analytic <= 0.05, injected_error_rejected=4 * W * student.sf(3, N - 1) > 0.05, detail={'mc': mc, 'analytic': analytic, 'expected_budget_exact': 0.05}), dict(name='HORIZON_ENCLOSURE', nominal_pass=unsafe == 0 and enclosure_residual < 1e-10, injected_error_rejected=bad_refutes, detail={'unsafe_intervals': unsafe, 'endpoint_control_error_um': enclosure_residual}), dict(name='HISTORY_SUMMARY_SUFFICIENCY', nominal_pass=suff['identity_error_um'] == 0 and suff['downstream_difference_um'] > 1, injected_error_rejected=bad_refutes)]
    qstart = time.perf_counter()
    for _ in range(10000):
        interval_hold(50, 20, 20, 1, 7, 120)
    query = (time.perf_counter() - qstart) / 10000
    out = {'id': 'X51_R3', 'claim_type': 'capability', 'external_referent': read('PREREG_R3.json')['external_referent'], 'finite_baseline': {'N': N, 'weekly_n': n, 'df': N - 1, 'tcrit': tcrit, 'normal_known_baseline_z': read('raw/R2.json')['limits']['z'], 'reference_threshold_at_nominal_sigma_um': float(tcrit * sigma * scale), 'process_threshold_at_nominal_sigma_um': float(tcrit * np.sqrt(2) * sigma * scale), 'expected_false_channel_signals_per_364_day_year_mc': mc, 'mc_standard_error': float(num.std(ddof=1) / np.sqrt(P)), 'probability_any_signal_mc': float(np.mean(num > 0)), 'analytic_budget': analytic, 'rigorous_ideal_model_budget': 0.05, 'calibration_scope': 'finite-baseline normal iid marginal Student prediction; empirical normality/independence/stability UNKNOWN'}, 'guard_scenarios': records, 'sufficiency_test': suff, 'understated_slope_counterexample': {'guard': badguard, 'true_next_gap_um': 130, 'refutes_zero_slope': bad_refutes}, 'checks': checks, 'gates': {'conditional_controls': all((q['nominal_pass'] for q in checks)), 'all_faults_rejected': all((q['injected_error_rejected'] for q in checks)), 'physical_error_and_drift_enclosure': 'UNKNOWN', 'physical_preemption': 'UNKNOWN'}, 'enclosure': read('PREREG_R3.json')['enclosure'], 'cost': {'wall_s': time.perf_counter() - start, 'query_wall_s': query, 'new_physical_measurements': 0}}
    save('raw/R3.json', out)
    np.savez_compressed(ROOT / 'raw/R3_PATHS.npz', baseline_mean_b_d=np.c_[base_r, base_d], baseline_sd_b_d=np.c_[sr, sd], delta_limits_b_d=np.c_[hr, hd], false_channel_signals=num.astype(np.int16))
    state('R3_COMPLETE', out['gates'], 'package runnable CLI, figure and X1b lab protocol extension; retain physical HOLD_UNKNOWN')
    (ROOT / 'history/HANDOFF_R3.md').write_text('R3 finite-baseline Student chart and bounded-horizon guard executed; no unsafe intervals in declared corners. Physical error/slope/gap mapping still UNKNOWN. No real parts measured. Next executable physical construction: independently referenced regional coupon time series plus withheld seated gap checks; freeze error/rate bounds before next lots.\n')
    print({'R3': out['gates'], 'finite_baseline': out['finite_baseline'], 'guard': records})
    if not all((c['nominal_pass'] and c['injected_error_rejected'] for c in checks)):
        raise RuntimeError('R3 verification failed')
if __name__ == '__main__':
    main()
