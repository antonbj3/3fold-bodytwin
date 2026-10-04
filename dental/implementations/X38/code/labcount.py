"""Fixed-plan laboratory specimen count, horizon bounds and two-design power."""
import argparse, json, math, sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.stats import beta, binom, chi2, hypergeom, nct, norm, t
NIST = 'https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/exacbici.htm'
FISHER = 'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html'
TTEST = 'https://www.statsmodels.org/stable/generated/statsmodels.stats.power.TTestIndPower.html'

def probability(p, name):
    if not math.isfinite(p) or not 0 < p < 1:
        raise ValueError(name + ' must be strictly between0 and1')

def integer(n, name, minimum=0):
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < minimum:
        raise ValueError(name + ' must be an integer >= ' + str(minimum))

def upper(k, n, alpha):
    integer(k, 'failures')
    integer(n, 'n')
    probability(alpha, 'tail alpha')
    if k > n:
        raise ValueError('failures exceed n')
    if n == 0 or k == n:
        return 1.0
    if k == 0:
        return float(-np.expm1(math.log(alpha) / n))
    return float(beta.ppf(1 - alpha, k + 1, n - k))

def plan(risk=0.1, confidence=0.95, max_failures=0, family_tests=1, true_risk=None):
    probability(risk, 'target risk')
    probability(confidence, 'confidence')
    integer(max_failures, 'maximum accepted failures')
    integer(family_tests, 'family tests', 1)
    alpha = (1 - confidence) / family_tests
    if max_failures == 0:
        n = max(1, math.floor(math.log(alpha) / math.log1p(-risk)) + 1)
    else:
        n = max_failures + 1
        while n <= 100000 and upper(max_failures, n, alpha) >= risk:
            n += 1
        if n > 100000:
            raise ValueError('sample search budget exceeded100000')
    while upper(max_failures, n, alpha) >= risk:
        n += 1
    while n > max_failures + 1 and upper(max_failures, n - 1, alpha) < risk:
        n -= 1
    x = dict(claim_type='capability', n_per_condition=n, max_accepted_failures=max_failures, target_risk_strictly_below=risk, family_confidence=confidence, family_tests=family_tests, one_sided_tail_alpha=alpha, upper_if_acceptance_boundary=upper(max_failures, n, alpha), previous_n_upper=upper(max_failures, n - 1, alpha), fixed_plan=True, resolution='POPULATION', unit='probability', meaning='Conditional minimum if at most the declared number fail. This is not statistical power or probability of success.', assumptions=['iid specimens from one condition/batch sampling regime', 'all outcomes known at prespecified cycle horizon', 'freeze sample size and acceptance boundary before testing', 'do not repeatedly peek or select load after outcomes'], external_referent=dict(kind='closed_form', locator=NIST, compared_quantity='exact one-sided binomial upper confidence limit', refutes_us=True))
    if true_risk is not None:
        if not math.isfinite(true_risk) or not 0 <= true_risk <= 1:
            raise ValueError('true risk must be in[0,1]')
        x['assumed_true_risk'] = true_risk
        x['probability_of_acceptance_at_assumed_risk'] = float(binom.cdf(max_failures, n, true_risk))
    return x

def endpoint_summary(rows, horizon, configuration, force, confidence=0.95, family_tests=1):
    if not math.isfinite(horizon) or horizon <= 0 or (not math.isfinite(force)) or (force <= 0):
        raise ValueError('positive finite cycle horizon and force required')
    probability(confidence, 'confidence')
    integer(family_tests, 'family tests', 1)
    if not isinstance(rows, list) or not rows:
        raise ValueError('rows must be a nonempty JSON list')
    ids = []
    chosen = []
    excluded = []
    for row in rows:
        i = row['specimen_id']
        if i in ids:
            raise ValueError('duplicate specimen ID ' + str(i))
        ids.append(i)
        if row['configuration'] != configuration or row['max_load_N'] != force:
            excluded.append(dict(specimen_id=i, reason='different prespecified configuration/load'))
            continue
        if row.get('cycles_unit') != 'cycles' or row.get('force_unit') != 'N' or (not row.get('locator')):
            raise ValueError('selected record missing unit or locator')
        cycles = row['cycles']
        failure = row['failure']
        censored = row['censored']
        if not isinstance(cycles, (float, int)) or isinstance(cycles, bool) or (not math.isfinite(cycles)) or (cycles <= 0):
            raise ValueError('invalid cycle observation')
        if failure not in (0, 1) or censored not in (0, 1) or failure + censored != 1:
            raise ValueError('inconsistent event/censor flags')
        event = 'failure' if failure and cycles <= horizon else 'survivor' if cycles >= horizon else 'unknown_early_censor'
        chosen.append(dict(specimen_id=i, endpoint=event, observed_cycles=cycles, locator=row['locator'], resolution='PER_TOOTH'))
    if not chosen:
        raise ValueError('no specimens at chosen condition/load')
    n = len(chosen)
    k = sum((x['endpoint'] == 'failure' for x in chosen))
    u = sum((x['endpoint'] == 'unknown_early_censor' for x in chosen))
    alpha = (1 - confidence) / family_tests
    return dict(claim_type='capability', configuration=configuration, force_N=force, horizon_cycles=horizon, n=n, known_failures=k, complete_survivors=n - k - u, unknown_early_censors=u, sample_identification_interval=[k / n, (k + u) / n], population_upper_one_sided=upper(k + u, n, alpha), best_completion_population_upper=upper(k, n, alpha), family_confidence=confidence, family_tests=family_tests, resolution='POPULATION', observation_resolution='PER_TOOTH', time_scale='SIMULTANEOUS', early_censors_discarded=0, dropout_fraction_selected=u / n, screening=dict(total=len(rows), selected=n, excluded=len(excluded), exclusion_fraction=len(excluded) / len(rows), reasons={'different prespecified configuration/load': len(excluded)}), endpoints=chosen, excluded_records=excluded, assumptions=['common prespecified load/protocol', 'iid sampling', 'fixed horizon; unresolved endpoints are kept as unknown', 'no parametric lifetime/censor-independence closure used'], external_referent=dict(kind='published_dataset', locator=sorted({x['locator'].split('#')[0] for x in chosen}), compared_quantity='per-specimen time/event and known horizon endpoint', refutes_us=True))

def fisher_pvalues(n):
    integer(n, 'n per group', 2)
    if n > 1000:
        raise ValueError('exact binary power limited to1000/group')
    x = np.arange(n + 1)[:, None]
    y = np.arange(n + 1)[None, :]
    return np.minimum(1.0, 2 * hypergeom.cdf(np.minimum(x, y), 2 * n, x + y, n))

def binary_power(n, p_a, p_b, alpha=0.05):
    probability(p_a, 'pA')
    probability(p_b, 'pB')
    probability(alpha, 'alpha')
    reject = fisher_pvalues(n) <= alpha
    mass_a = binom.pmf(np.arange(n + 1), n, p_a)
    mass_b = binom.pmf(np.arange(n + 1), n, p_b)
    return float(mass_a @ (reject @ mass_b))

def continuous_power(n, effect, alpha=0.05):
    integer(n, 'n per group', 2)
    probability(alpha, 'alpha')
    if not math.isfinite(effect) or effect <= 0:
        raise ValueError('effect must be positive abs(meanA-meanB)/common_SD')
    df = 2 * n - 2
    cut = t.isf(alpha / 2, df)
    delta = effect * math.sqrt(n / 2)
    value = float(nct.sf(cut, df, delta) + nct.cdf(-cut, df, delta))
    if math.isfinite(value) and 0 <= value <= 1:
        return value
    if not math.isfinite(delta):
        raise ValueError('effect/noncentrality exceeds numerical range')
    (value, error) = quad(lambda y: norm.pdf(y) * chi2.cdf(df * ((y + delta) / cut) ** 2, df), -10, 10, epsabs=1e-10, epsrel=1e-10, limit=200)
    if not math.isfinite(value) or error > 1e-09 or (not 0 <= value <= 1):
        raise ValueError('noncentral-t fallback precision unresolved')
    return float(value)

def compare(endpoint, power=0.8, alpha=0.05, p_a=None, p_b=None, effect=None, max_n=500):
    probability(power, 'desired power')
    probability(alpha, 'alpha')
    integer(max_n, 'max n', 2)
    if endpoint == 'binary':
        probability(p_a, 'pA')
        probability(p_b, 'pB')
        if p_a == p_b:
            raise ValueError('pA and pB must differ')
        if max_n > 1000:
            raise ValueError('binary power search capped1000/group')
        query = lambda n: binary_power(n, p_a, p_b, alpha)
        assumptions = ['independent unpaired iid groups', 'fixed common horizon; every binary outcome known', 'prespecified two-sided Fisher test', 'pA and pB are planning assumptions, not measured risks']
        source = FISHER
        quantity = 'exact Fisher test rejection probability under two independent binomial groups'
    elif endpoint == 'continuous':
        query = lambda n: continuous_power(n, effect, alpha)
        assumptions = ['independent unpaired groups', 'normal population common SD', 'prespecified two-sided pooled t-test', 'effect=absolute mean difference / common SD supplied by lab, not calibrated here']
        source = TTEST
        quantity = 'noncentral-t two-sample test power'
    else:
        raise ValueError('unknown endpoint')
    curve = []
    selected = None
    for n in range(2, max_n + 1):
        value = query(n)
        curve.append(dict(n_per_group=n, power=value))
        if value >= power:
            selected = n
            break
    x = dict(claim_type='capability', endpoint=endpoint, status='PLANNED' if selected else 'SEARCH_LIMIT', n_per_group=selected, total_n=2 * selected if selected else None, desired_power=power, alpha=alpha, test_sidedness='two-sided', achieved_power=curve[-1]['power'], search_curve=curve, maximum_previous_power=max([q['power'] for q in curve[:-1]], default=0), search_max_n=max_n, pA=p_a, pB=p_b, standardized_effect=effect, assumptions=assumptions, resolution='POPULATION', unit='probability', uncertainty='Exact conditional mathematical power; population assumptions have no empirical confidence band', external_referent=dict(kind='published_code', locator=source, compared_quantity=quantity, refutes_us=True))
    return x

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='command', required=True)
    pl = sub.add_parser('plan')
    pl.add_argument('--risk', required=True, type=float)
    pl.add_argument('--confidence', default=0.95, type=float)
    pl.add_argument('--max-failures', default=0, type=int)
    pl.add_argument('--family-tests', default=1, type=int)
    pl.add_argument('--true-risk', type=float)
    ep = sub.add_parser('endpoints')
    ep.add_argument('--data', required=True)
    ep.add_argument('--horizon-cycles', required=True, type=float)
    ep.add_argument('--configuration', required=True)
    ep.add_argument('--force-N', required=True, type=float)
    ep.add_argument('--confidence', default=0.95, type=float)
    ep.add_argument('--family-tests', default=1, type=int)
    co = sub.add_parser('compare')
    co.add_argument('--endpoint', required=True, choices=['binary', 'continuous'])
    co.add_argument('--power', default=0.8, type=float)
    co.add_argument('--alpha', default=0.05, type=float)
    co.add_argument('--p-a', type=float)
    co.add_argument('--p-b', type=float)
    co.add_argument('--effect', type=float)
    co.add_argument('--max-n', default=500, type=int)
    for p in [pl, ep, co]:
        p.add_argument('--output')
    a = ap.parse_args()
    try:
        if a.command == 'plan':
            x = plan(a.risk, a.confidence, a.max_failures, a.family_tests, a.true_risk)
        elif a.command == 'endpoints':
            x = endpoint_summary(json.loads(Path(a.data).read_text()), a.horizon_cycles, a.configuration, a.force_N, a.confidence, a.family_tests)
        else:
            if a.endpoint == 'binary' and (a.p_a is None or a.p_b is None):
                raise ValueError('binary effect requires both --p-a and --p-b')
            if a.endpoint == 'continuous' and a.effect is None:
                raise ValueError('continuous comparison requires --effect')
            x = compare(a.endpoint, a.power, a.alpha, a.p_a, a.p_b, a.effect, a.max_n)
        if a.output:
            Path(a.output).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')
        print(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False))
        return 0 if x.get('status') != 'SEARCH_LIMIT' else 3
    except (ValueError, TypeError, KeyError, OSError) as e:
        x = dict(status='REFUSED', reason=str(e))
        if a.output:
            Path(a.output).write_text(json.dumps(x, indent=2) + '\n')
        print(json.dumps(x))
        return 2
if __name__ == '__main__':
    sys.exit(main())
