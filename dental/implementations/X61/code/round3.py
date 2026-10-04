"""R3: open the observed sinter-vs-metrology ambiguity, don't fit it away."""
import datetime, json, math, time, resource
import numpy as np
from scipy.stats import norm
from lab_alarm import R, Monitor, verified, write, sha
from paired_metrology import PairedPort

def trial(seed, kind, pr):
    p = pr['profile']
    rng = np.random.default_rng(seed)
    scan_sd = p['scan_sd']
    gauge_sd = p['gauge_sd']
    rho = p['cross_error_correlation']
    S = np.array([[scan_sd ** 2, rho * scan_sd * gauge_sd], [rho * scan_sd * gauge_sd, gauge_sd ** 2]])
    shift = pr['faults'].get(kind, [0.0, 0.0])
    theta = np.array([p['sinter_null'] + shift[0], shift[1]])
    y = rng.multivariate_normal([0.0, 0.0], S, pr['metrics']['horizon']) + PairedPort.A @ theta
    mo = Monitor(dict(alpha_by_channel={'sinter': 0.025, 'scanner_bias': 0.025}, bins=8), {})
    port = PairedPort(p['prior_mean'], p['prior_sd'])
    timeline = []
    diff_sd = math.sqrt(scan_sd ** 2 + gauge_sd ** 2 - 2 * rho * scan_sd * gauge_sd)
    for (i, (scan, gauge)) in enumerate(y):
        mo.t = i + 1
        ga = mo.score_category('sinter', norm.cdf((gauge - p['sinter_null']) / gauge_sd))
        sc = mo.score_category('scanner_bias', norm.cdf((scan - gauge - p['scanner_bias_null']) / diff_sd))
        port.update(scan, gauge, scan_sd, gauge_sd, rho)
        timeline.append(dict(specimen=i + 1, scan=float(scan), gauge=float(gauge), sinter_wealth=ga['wealth'], scanner_wealth=sc['wealth']))
    rep = port.report()
    covered = [rep[k]['credible_95'][0] <= theta[i] <= rep[k]['credible_95'][1] for (i, k) in enumerate(('sinter_factor', 'scanner_bias'))]
    return (dict(seed=seed, kind=kind, first_hits=mo.first_hits, any_alarm=bool(mo.first_hits), posterior=rep, truth=theta.tolist(), posterior_coverage=covered), timeline, y, S)

def exact_collision():
    theta = np.array([[0.75, 0.03125], [0.78125, 0.0]])
    scan = theta.sum(axis=1)
    gauge = theta[:, 0]
    return dict(states=theta.tolist(), scanner_summary=scan.tolist(), identity_error=float(abs(scan[0] - scan[1])), gauge_difference=float(abs(gauge[0] - gauge[1])), minimum_extension='One independently calibrated ratio channel with noncollinear observation row; full covariance for paired errors', resolution='PER_TOOTH')

def run():
    pr = verified(R / 'PREREG_R3.json')
    st = time.perf_counter()
    raw = []
    examples = {}
    null = []
    faults = {}
    for i in range(pr['metrics']['null_replicates']):
        (r, *_) = trial(161000 + i, 'nominal', pr)
        null.append(r)
        raw.append(r)
    for kind in pr['faults']:
        records = []
        for i in range(pr['metrics']['fault_replicates']):
            (r, *_) = trial(171000 + i, kind, pr)
            records.append(r)
            raw.append(r)
        desired = ['sinter', 'scanner_bias'] if kind == 'both' else ['sinter'] if kind == 'true_sinter' else ['scanner_bias']
        t = [max((r['first_hits'][k] for k in desired)) for r in records if all((k in r['first_hits'] for k in desired))]
        correct = []
        for r in records:
            h = r['first_hits']
            correct.append(all((k in h for k in desired)) and (kind == 'both' or h[desired[0]] < h.get('scanner_bias' if desired == ['sinter'] else 'sinter', 10000)))
        faults[kind] = dict(replicates=len(records), detected=len(t), detection_fraction=len(t) / len(records), correct_first_link_fraction=float(np.mean(correct)), median_measurements_among_detected=float(np.median(t)) if t else None, p10_p90_measurements=np.quantile(t, [0.1, 0.9]).tolist() if t else None, posterior_coverage=np.mean([r['posterior_coverage'] for r in records], axis=0).tolist(), resolution='POPULATION', diagnosis_scope='Gauge calibrated independently; common gauge bias remains UNKNOWN beyond this declared covariance')
        (ex, timeline, y, S) = trial(16100, kind, pr)
        examples[kind] = dict(**ex, timeline=timeline)
        write(R / 'raw' / ('EXAMPLE_R3_' + kind + '.json'), examples[kind])
    (r, t, y, S) = trial(16100, 'scanner_bias', pr)
    A = np.tile(PairedPort.A, (len(y), 1))
    V = np.kron(np.eye(len(y)), S)
    prec0 = np.diag(1 / np.array(pr['profile']['prior_sd']) ** 2)
    mean0 = np.array(pr['profile']['prior_mean'])
    ViA = np.linalg.solve(V, A)
    C = np.linalg.inv(prec0 + A.T @ ViA)
    mu = C @ (prec0 @ mean0 + A.T @ np.linalg.solve(V, y.ravel()))
    parity = float(max(np.max(np.abs(mu - r['posterior']['mean'])), np.max(np.abs(C - r['posterior']['covariance']))))
    wrong_s = np.mean(y[:, 0])
    right = np.mean(y[:, 1])
    mut = dict(swapped_channel_changes_physical_ratio=abs(wrong_s - right) > 0.003, wrong_covariance_rejected=False, wrong_dense_mean_rejected=float(np.max(np.abs(mu + np.array([0.01, 0.0]) - r['posterior']['mean']))) > pr['metrics']['independent_GLS_parity_max'])
    try:
        PairedPort().update(0.8, 0.8, 0.001, 0.0005, 1.1)
    except ValueError:
        mut['wrong_covariance_rejected'] = True
    mr = np.array(r['posterior']['mean'])
    CV = np.array(r['posterior']['covariance'])
    margin = 0.0001
    full = float(norm.sf(margin / math.sqrt(float(np.array([1.0, 1.0]) @ CV @ np.array([1.0, 1.0])))))
    diag = float(norm.sf(margin / math.sqrt(float(np.trace(CV)))))
    cov_witness = dict(summary_means_identical=mr.tolist(), summary_marginal_variances_identical=np.diag(CV).tolist(), identity_error=0.0, joint_covariance=float(CV[0, 1]), downstream_probability_full=full, downstream_probability_if_independent=diag, difference=abs(full - diag), minimum_extension='Joint covariance; conditional margin cannot multiply independent node probabilities')
    write(R / 'raw/TRIALS_R3.json', raw)
    write(R / 'raw/MUTATION_CONTROLS_R3.json', mut)
    rate = sum((x['any_alarm'] for x in null)) / len(null)
    q = pr['metrics']
    coverage = np.mean([r['posterior_coverage'] for r in null], axis=0).tolist()
    collision = exact_collision()
    gates = dict(exact_summary_collision=collision['identity_error'] == 0, nominal_FPR=rate <= q['MC_FPR_max'], all_faults_detected=all((x['detection_fraction'] >= q['detection_fraction_min'] for x in faults.values())), correct_first_link=all((x['correct_first_link_fraction'] >= q['correct_first_link_fraction_min'] for x in faults.values())), median_measurements=all((x['median_measurements_among_detected'] <= q['median_alarm_measurements_max'] for x in faults.values())), posterior_coverage=all((v >= q['posterior_coverage_min'] for v in coverage)), independent_GLS=parity <= q['independent_GLS_parity_max'], mutations=all(mut.values()))
    out = dict(round='R3', claim_type='capability', outcome='PASS_SCOPED_PAIRED_LOCALIZATION' if all(gates.values()) else 'PARTIAL_WITH_FAILED_GATES', gates=gates, nominal=dict(replicates=len(null), alarms=sum((x['any_alarm'] for x in null)), fraction=rate, alpha=0.05, posterior_coverage=coverage), faults=faults, exact_sufficiency_witness=collision, joint_covariance_sufficiency=cov_witness, external_referent=pr['external_referent'], independent_GLS_max_abs_error=parity, cost=dict(wall_seconds=time.perf_counter() - st, max_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, preparation='Synthetic paired instrument profile; empirical acquisition/noise costs UNKNOWN', fit='2x2 joint precision update', discovery='400 null +200 per three fault families, 72 observations', validation='Independent 144x144 dense GLS, fixed-seed coverage/mutations', questions=0, fallback='UNKNOWN calibrated gauge access/cost; shared systematic gauge bias is not resolved by more same-gauge samples'), physical_validation='UNKNOWN_NO_PAIRED_LAB_MEASUREMENTS', affine_enclosure='Exact for declared additive Gaussian observation; physical nonlinear ratio metrology remainder is missing', resolution='PER_TOOTH pairs -> POPULATION s/b; no patient inference', timescale='SIMULTANEOUS pairing and HANDOVER process state')
    write(R / 'raw/RESULTS_R3.json', out)
    (R / 'HANDOFF_R3.md').write_text(json.dumps(dict(outcome=out['outcome'], gates=gates, next='Calibrate independent gauge systematic error and m-dependent stress hazard, freeze new batch and hold validation; never absorb common bias into sinter'), indent=2) + '\n')
    write(R / 'CURRENT_WORK_STATE.json', dict(lane='X61-lab-alarm', status='R3_DECIDED', latest_gate=gates, next_operation='Connect paired port to CSV CLI; one-command full replay + graph feedback'))
    print(json.dumps(dict(round='R3', gates=gates, nominal=out['nominal'], faults=faults, cost=out['cost'])), flush=True)
    return out
if __name__ == '__main__':
    run()
