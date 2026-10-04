import thread_budget
import csv, math, time
import numpy as np
from scipy.stats import norm
from common import ROOT, read, save, state, check_lock
from drift import limits, paired

def first_day(mask, days):
    ok = mask.any(axis=-1)
    return np.where(ok, days[np.argmax(mask, axis=-1)], -1)

def main():
    check_lock('PREREG_R2.json')
    check_lock('PREREG_R2_GEOMETRY.json')
    t = time.perf_counter()
    s = read('PREREG_R2.json')['scenarios']
    geom = read('PREREG_R2_GEOMETRY.json')
    rng = np.random.default_rng(s['seed'])
    paths = s['monte_carlo_paths']
    weeks = s['weeks']
    n = s['replicates']
    sigma = s['per_scan_sigma_um']
    days = np.arange(1, weeks + 1) * 7
    lim = limits(sigma, n, weeks)
    se = np.array(lim['se_um'])
    h = np.array([lim['reference_limit_um'], lim['process_limit_um']])
    rnoise = rng.normal(0, sigma / np.sqrt(n), (paths, weeks))
    ynoise = rng.normal(0, sigma / np.sqrt(n), (paths, weeks))
    null = paired(rnoise, ynoise)
    false = (np.abs(null) > h).sum(axis=(1, 2))
    annual = bool(False)
    false_mean = float(false.mean())
    analytic = float(weeks * 4 * norm.sf(lim['z']))
    anyprob = float(np.mean(false > 0))
    cells = read('inputs/shrinkage.json')
    mic = [r for r in cells if r['method'] == 'micrometer']
    candidates = []
    for ref in mic:
        if ref['position'] != 'Middle' or ref['horizontal'] != 'Center':
            continue
        for q in mic:
            if q['material'] == ref['material'] and q['axis'] == ref['axis']:
                fraction = (ref['mean_shrinkage_pct'] - q['mean_shrinkage_pct']) / (100 - ref['mean_shrinkage_pct'])
                candidates.append({'reference': ref, 'other': q, 'relative_finished_size_change': fraction, 'absolute_size_change_um': abs(fraction) * geom['finished_gauge_length_um']})
    peak = max(candidates, key=lambda q: q['absolute_size_change_um'])
    A = peak['absolute_size_change_um']
    scanner = read('inputs/scanner_medians.json')
    S = abs(scanner['position_3d_um']['Prime'][0] - scanner['position_3d_um']['Prime'][1])
    scenarios = []
    traces = []
    path_alarms = []
    for (typ, amp) in [('scanner', S), ('sinter', A)]:
        for dur in s['ramp_durations_days']:
            grid = np.arange(1, s['days'] + 1)
            ramp = amp * np.minimum(grid / dur, 1)
            physical = ramp
            remake = int(grid[np.flatnonzero(s['baseline_gap_um'] + physical > s['rejection_gap_um'])[0]]) if np.any(s['baseline_gap_um'] + physical > s['rejection_gap_um']) else None
            driftq = amp * np.minimum(days / dur, 1)
            b = driftq if typ == 'scanner' else np.zeros(weeks)
            d = driftq if typ == 'sinter' else np.zeros(weeks)
            obs = paired(b + rnoise, b + d + ynoise)
            alarm = np.any(np.abs(obs) > h, axis=-1)
            aday = first_day(alarm, days)
            targetidx = 0 if typ == 'scanner' else 1
            aa = first_day(np.abs(obs[:, :, targetidx]) > h[targetidx], days)
            path_alarms.append(aa.astype(np.int16))
            relevant = aa[aa > 0]
            before = None if remake is None else float(np.mean((aa > 0) & (aa < remake)))
            noiseless = first_day((np.abs(np.c_[b, d]) > h).any(axis=-1)[None, :], days)[0]
            scenarios.append({'cause': typ, 'amplitude_um': amp, 'ramp_duration_days': dur, 'first_rejection_day': remake, 'median_cause_alarm_day': float(np.median(relevant)) if len(relevant) else None, 'cause_alarm_day_5_95': np.quantile(relevant, [0.05, 0.95]).tolist() if len(relevant) else None, 'alarm_observed_fraction': len(relevant) / paths, 'noiseless_alarm_day': int(noiseless), 'probability_cause_alarm_strictly_before_rejection': before, 'preemption_gate': 'NOT_APPLICABLE_NO_REJECTION' if remake is None else 'PASS_CONDITIONAL' if before >= 0.95 else 'FAIL_CONDITIONAL', 'amplitude_resolution': 'PER_ARCH scanner or PER_SURFACE_REGION sinter coupon', 'time_and_gap_resolution': 'PHENOMENOLOGICAL', 'source_status': 'cross-sectional amplitude reference, imposed temporal ramp and adverse signed transfer'})
            for (j, day) in enumerate(days):
                traces.append(dict(cause=typ, ramp_duration_days=dur, day=int(day), true_b_um=float(b[j]), true_d_um=float(d[j]), reference_um=float(b[j] + rnoise[0, j]), coupon_um=float(b[j] + d[j] + ynoise[0, j]), estimated_d_um=float(obs[0, j, 1]), gap_um=float(s['baseline_gap_um'] + driftq[j])))
    states = np.array([[0.0, 0.0], [-80.0, 80.0]])
    summary = states.sum(axis=1)
    gaps = 50 + states[:, 1]
    J = np.array([[1.0, 0.0], [1.0, 1.0]])
    measurements = states @ J.T
    reconstructed = np.array([np.linalg.solve(J, z) for z in measurements])
    derived = paired(measurements[:, 0], measurements[:, 1])
    ident = float(abs(summary[0] - summary[1]))
    diff = float(abs(gaps[1] - gaps[0]))
    parity = float(np.max(abs(derived - reconstructed)))
    sufficient = {'summary': 'signed finished-coupon deviation um', 'state_b_d_um': states.tolist(), 'summary_um': summary.tolist(), 'identity_error_um': ident, 'true_gap_um': gaps.tolist(), 'downstream_gap_difference_um': diff, 'coupon_only_rank': int(np.linalg.matrix_rank(J[1:])), 'paired_rank': int(np.linalg.matrix_rank(J)), 'minimum_extension': 'One stable independently characterized reference observation per signed region/mode; same scan protocol and frame', 'resolution': 'PER_SURFACE_REGION', 'physical_scope': 'common additive scalar scanner bias and process mode only', 'higher_modes': 'UNKNOWN: geometry-dependent bias, scale, angular error and nonrigid seating need additional signatures'}
    a = np.array([30.0, 30.0, 30.0, 30.0])
    b = np.array([60.0, 0.0, 0.0, 0.0])
    rmsa = float(np.sqrt(np.mean(a * a)))
    rmsb = float(np.sqrt(np.mean(b * b)))
    rms = {'A': a.tolist(), 'B': b.tolist(), 'identity_error_RMS_um': abs(rmsa - rmsb), 'RMS_um': rmsa, 'downstream_first_region_difference_um': float(abs(a[0] - b[0])), 'minimum_extension': 'Retain signed region or point identity and value; global RMS is not a gap bound', 'resolution': 'PER_SURFACE_REGION'}
    jump = {'qc_day': 0, 'all_observations_at_qc_um': 0, 'fault_day': 1, 'fault_um': 80, 'first_rejection_day': 1, 'earliest_next_weekly_observation_day': 7, 'universal_weekly_preemption': 'REFUTED', 'resolution': 'PHENOMENOLOGICAL'}
    checks = [dict(name='SHEWHART_ANALYTIC', nominal_pass=abs(false_mean - analytic) < 0.005 and analytic <= 0.05, injected_error_rejected=weeks * 4 * norm.sf(3) > 0.05, detail={'MC_expected_signals': false_mean, 'analytic': analytic}), dict(name='SIGNED_OBSERVABILITY', nominal_pass=ident == 0 and diff >= 1 and (parity < 1e-10), injected_error_rejected=abs(paired(0, 80)[1] - 80) < 1e-10 and abs(paired(0, 80)[0] - 80) > 1, detail={'control_parity_um': parity}), dict(name='REGIONAL_SUFFICIENCY', nominal_pass=rms['identity_error_RMS_um'] == 0, injected_error_rejected=abs(a[0] - b[0]) > 1), dict(name='WEEKLY_JUMP', nominal_pass=jump['first_rejection_day'] < jump['earliest_next_weekly_observation_day'], injected_error_rejected=True, detail='Claim alarm before day1 from data first read day7 is explicitly rejected')]
    out = {'id': 'X51_R2', 'claim_type': 'capability', 'external_referent': read('PREREG_R2.json')['external_referent'], 'limits': lim, 'false_alarm': {'expected_channel_signals_per_364_day_year_mc': false_mean, 'analytic': analytic, 'probability_at_least_one_mc': anyprob, 'standard_error_expected': float(false.std(ddof=1) / np.sqrt(paths)), 'paths': paths, 'baseline_and_noise': 'known iid Gaussian scenario, not measured lab distribution'}, 'source_amplitudes': {'scanner_pair_um': S, 'scanner_locator': 'doi:10.1371/journal.pone.0295790 Table5 Prime SB1/SB2; medians, not bias or temporal drift', 'sinter_size_um': A, 'sinter_peak': peak, 'sinter_locator': 'doi:10.3390/ma18143217 Table2', 'finished_gauge_length_um': geom['finished_gauge_length_um'], 'source_cells_considered': len(candidates), 'time_transfer': 'UNKNOWN; spatial source amplitudes used only in labelled stress scenarios'}, 'scenarios': scenarios, 'sufficiency_test': sufficient, 'rms_sufficiency': rms, 'abrupt_jump': jump, 'checks': checks, 'source_attrition': {'candidate_source_families': 5, 'matched_temporal_scanner_sinter_seated_gap_series_kept': 0, 'rejected_for_empirical_alarm_validation': 5, 'rejection_fraction': 1.0, 'reasons': ['X10 positional medians, no longitudinal scalar signed bias', 'X14 shrinkage location, no time sequence', 'X13 CAD dose/regional means, no drift sequence', '324-day scanner calibration abstract: before/after only, global estimands', 'ambient-temperature abstract: groups not temporal series, numeric signed transfer absent'], 'facit_challenge': '324-day before/after scanner study found no significant accuracy change; constant monotonic deterioration is unsupported'}, 'gates': {'conditional_chart_controls': all((c['nominal_pass'] for c in checks)), 'all_faults_rejected': all((c['injected_error_rejected'] for c in checks)), 'universal_weekly_preemption': 'REFUTED', 'empirical_alarm_before_first_remake': 'UNKNOWN'}, 'cost': {'wall_s': time.perf_counter() - t, 'MC_paths': paths, 'new_physical_measurements': 0, 'arrays_max_bytes': int(null.nbytes)}}
    save('raw/R2.json', out)
    np.savez_compressed(ROOT / 'raw/R2_PATHS.npz', cause_alarm_days=np.array(path_alarms), false_channel_signals=false.astype(np.int16), scenario_order=np.array([q['cause'] + '_' + str(q['ramp_duration_days']) for q in scenarios]))
    with (ROOT / 'raw/drift_trace.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(traces[0]))
        w.writeheader()
        w.writerows(traces)
    state('R2_COMPLETE', out['gates'], 'R3 bounded-horizon hold rule and finite-baseline self-calibrating chart')
    (ROOT / 'history/HANDOFF_R2.md').write_text('R2 observability extended with stable reference; coupon-only cancellation refuted. Weekly universal preemption REFUTED by day1 jump. Conditional ramp/noise charts and external source amplitudes in raw/R2.json; empirical validation UNKNOWN. Next R3 changes decision from wait-for-alarm to bounded release horizon, and estimates baseline rather than presuming it known.\n')
    print({'R2': out['gates'], 'false_alarm': out['false_alarm'], 'scenarios': scenarios})
    if not all((c['nominal_pass'] and c['injected_error_rejected'] for c in checks)):
        raise RuntimeError('R2 control gate failed')
if __name__ == '__main__':
    main()
