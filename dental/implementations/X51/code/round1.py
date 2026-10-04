import thread_budget
import math, time
import numpy as np
from scipy.stats import binom, weibull_min
from common import ROOT, read, save, state, check_lock
from acceptance import upper, n_zero, plan, endpoint_rows, known_shape_lower, known_shape_plan, weibull_fit, profile_lower
import x25_fatigue as aft

def explicit_binomial(n, c, p):
    return sum((math.comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(c + 1)))

def main():
    check_lock('PREREG_R1.json')
    check_lock('FROZEN_PREDICTIONS.json')
    spec = read('PREREG_R1.json')
    t = time.perf_counter()
    checks = []

    def check(name, nominal, error_rejected, detail=None):
        checks.append(dict(name=name, nominal_pass=bool(nominal), injected_error_rejected=bool(error_rejected), detail=detail))
    examples = [{'n': n, 'upper95': upper(n, 0), 'hand': 1 - 0.05 ** (1 / n)} for n in [3, 12, 24, 25, 29, 58, 59]]
    e = max((abs(r['upper95'] - r['hand']) for r in examples))
    check('CP_HAND', e < 1e-10, abs(upper(3, 0) - upper(3, 1)) > 0.01, e)
    zs = [n_zero(0.1), n_zero(0.1, groups=4), n_zero(0.05)]
    check('X25_MINIMUM_COUNTS', [r['n'] for r in zs] == [29, 42, 59], upper(28, 0) > 0.1)
    p = plan(**{k: spec['scenarios'][k] for k in ['p_bad', 'p_good']}, power=0.8)
    ref = explicit_binomial(p['n'], p['failures_allowed'], 0.1)
    check('EXACT_CONSUMER_AND_PRODUCER_RISK', abs(ref - p['false_acceptance_at_p_bad']) < 1e-10 and ref <= 0.05 and (p['power_at_p_good'] >= 0.8), explicit_binomial(p['n'], p['failures_allowed'] + 1, 0.1) > 0.05, p)
    rows = read('inputs/specimens.json')
    groups = sorted({r['configuration'] for r in rows})
    fits = read('inputs/x25_fits.json')
    replay = []
    maxfit = maxload = 0.0
    for (family, label) in [('lognormal', 'lognormal'), ('Weibull', 'Weibull')]:
        f = aft.fit(rows, groups, family)
        old = fits[label]
        maxfit = max(maxfit, abs(f['nll'] - old['nll']))
        for g in groups:
            for q in [0.1, 0.5, 0.9]:
                val = aft.force_at_risk(f, q, g)
                expected = aft.force_at_risk(old, q, g)
                maxload = max(maxload, abs(val - expected))
                replay.append(dict(family=family, configuration=g, risk=q, force_N=val, source_force_N=expected, resolution='PER_TOOTH', model_assumptions='same censored AFT law and source regime'))
    check('X25_CENSORED_LOAD_REPLAY', maxload < 1e-06 and maxfit < 1e-06, abs(aft.force_at_risk(fits['naive'], 0.1, groups[0]) - replay[0]['force_N']) > 1e-06, {'max_load_error_N': maxload, 'max_nll_error': maxfit})
    cell = [r for r in rows if r['configuration'] == 'MUSLA 04020' and r['max_load_N'] == 100]
    ep = endpoint_rows([dict(value=r['cycles'], failure=r['failure']) for r in cell], aft.NSTAR)
    check('X25_THREE_RUNOUTS', ep['good'] == 3 and abs(ep['upper95_conservative_all_unknown_bad'] - (1 - 0.05 ** (1 / 3))) < 1e-10, endpoint_rows([dict(value=1, failure=0), dict(value=1, failure=0), dict(value=1, failure=0)], aft.NSTAR)['upper95_conservative_all_unknown_bad'] == 1, ep)
    strengths = read('inputs/strength_rows.json')
    table = []
    max_mle = 0.0
    for g in sorted({r['group'] for r in strengths}):
        x = np.array([r['strength_MPa'] for r in strengths if r['group'] == g])
        f = weibull_fit(x)
        sc = weibull_min.fit(x, floc=0)
        max_mle = max(max_mle, abs(f['m'] - sc[0]) / sc[0], abs(f['eta'] - sc[2]) / sc[2])
        table.append({'group': g, 'n': len(x), 'sample_mean_MPa': float(x.mean()), 'sample_sd_MPa': float(x.std(ddof=1)), 'm_fitted': f['m'], 'eta_MPa': f['eta'], 'mean_profile': profile_lower(x, quantity='mean'), 'q05_profile': profile_lower(x, quantity='q05'), 'threshold300_exact_risk': endpoint_rows([dict(value=v, failure=1) for v in x], 300), 'ISO_conformity': 'UNKNOWN; raw groups are different surface/fixture conditions; full standard protocol not verified', 'resolution': 'POPULATION'})
    check('WEIBULL_MLE_STRONG_CONTROL', max_mle < 1e-05, abs(weibull_fit([100, 100, 100, 500])['m'] - weibull_fit([100, 200, 300, 400])['m']) > 0.1, max_mle)
    x = np.sort(np.array([r['strength_MPa'] for r in strengths if r['group'] == table[0]['group']]))
    r = len(x) // 2
    type2 = np.r_[x[:r], np.repeat(x[r - 1], len(x) - r)]
    known = known_shape_lower(type2, 5, censoring='TYPE_II', r=r, shape_independently_known=True)
    censfit = weibull_fit(type2, np.r_[np.ones(r, int), np.zeros(len(x) - r, int)])
    complete = known_shape_lower([2, 2, 2], 1, shape_independently_known=True)
    from scipy.stats import chi2
    hand = 12 / chi2.isf(0.05, 6)
    faults = []
    for kw in [dict(censoring='TYPE_I', shape_independently_known=True), dict(censoring='COMPLETE', shape_independently_known=False)]:
        try:
            known_shape_lower(x, 5, **kw)
            faults.append(False)
        except ValueError:
            faults.append(True)
    check('KNOWN_SHAPE_AND_CENSORING_CONTRACT', abs(complete['eta_lower'] - hand) < 1e-10 and censfit['status'] == 'FITTED', all(faults), {'hand_error': abs(complete['eta_lower'] - hand)})
    a = np.array([100.0, 100.0, 100.0, 100.0])
    b = np.array([40.0, 120.0, 120.0, 120.0])
    err = abs(a.mean() - b.mean())
    difference = float(np.mean(b < 75) - np.mean(a < 75))
    sufficient = {'summary': 'sample mean MPa', 'state_A': a.tolist(), 'state_B': b.tolist(), 'identity_error_MPa': float(err), 'downstream': 'empirical fraction below frozen75MPa', 'difference': difference, 'resolution': 'POPULATION', 'external_referent': {'kind': 'closed_form', 'locator': 'Bernoulli event count and arithmetic mean definitions; ISO6872 Table1 is a mean criterion', 'compared_quantity': 'same mean versus different threshold-failure fractions', 'refutes_us': True}, 'minimal_extension': 'For one fixed proof stress retain count below that stress plus unresolved censor count. General quantile questions need distribution or additional thresholds.'}
    check('SUMMARY_SUFFICIENCY', err == 0 and difference > 0.0, abs(a.mean() - (b + 1).mean()) > 0, sufficient)
    kp = [known_shape_plan(m, ratio) for m in [3, 5, 10] for ratio in [1.05, 1.1, 1.2]]
    out = {'id': 'X51_R1', 'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'external_referent': spec['external_referent'], 'zero_failure_examples': examples, 'minimum_plans': zs, 'producer_power_plan': p, 'x25_replay': replay, 'x25_endpoint_cell': ep, 'strength_groups': table, 'known_shape_planning': kp, 'typeII_reanalysis': {'group': table[0]['group'], 'n': len(x), 'fixed_r': r, 'records': type2.tolist(), 'assumed_known_m5_result': known, 'unknown_shape_censored_fit': censfit, 'scope': 'retrospective censor reconstruction; m5 scenario not externally established for these bars'}, 'sufficiency_test': sufficient, 'checks': checks, 'gates': {'exact_replay_and_controls': all((q['nominal_pass'] for q in checks)), 'all_faults_rejected': all((q['injected_error_rejected'] for q in checks)), 'empirical_crown_lot_conformity': 'UNKNOWN', 'unknown_shape_finite_sample_coverage': 'UNKNOWN'}, 'attrition': {'strength_rows_total': len(strengths), 'kept_for_model_reanalysis': len(strengths), 'rejected_for_crown_ISO_conformity': len(strengths), 'rejection_fraction': 1.0, 'reason': 'Different specimen/fixture/population and incomplete protocol; bar data do not certify crowns', 'implant_rows_total': len(rows), 'kept_for_likelihood': len(rows), 'unknown_endpoints': sum((not aft.endpoint_known(r) for r in rows)), 'endpoint_missing_fraction': sum((not aft.endpoint_known(r) for r in rows)) / len(rows)}, 'cost': {'wall_s': time.perf_counter() - t, 'new_physical_measurements': 0}}
    save('raw/R1.json', out)
    state('R1_COMPLETE', out['gates'], 'R2 paired-reference drift observability and weekly chart')
    (ROOT / 'history/HANDOFF_R1.md').write_text('R1 exact test inversion and X25 replay complete. See raw/R1.json; all original gates retained. Unknown-shape profile is approximate, not a finite-sample certificate. Next: R2 paired-reference weekly drift chart; same mean does not fix tail risk.\n')
    print(json_summary(out))
    if not all((q['nominal_pass'] and q['injected_error_rejected'] for q in checks)):
        raise RuntimeError('R1 verification gate failed')

def json_summary(out):
    import json
    return json.dumps({'R1': out['gates'], 'plan': out['producer_power_plan'], 'source_groups': len(out['strength_groups'])}, indent=2)
if __name__ == '__main__':
    main()
