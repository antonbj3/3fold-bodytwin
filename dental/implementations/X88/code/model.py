"""Study-bound pathway statistics; no treatment-selection rule."""
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import product
from math import comb
import numpy as np
from scipy.optimize import brentq
from scipy.stats import beta, binomtest, fisher_exact
CATEGORIES = ('normal', 'widening', 'lesion')

def cp(k, n, alpha=0.05):
    if not (isinstance(k, int) and isinstance(n, int) and (0 <= k <= n) and (n > 0)):
        raise ValueError('Integer 0 <= k <= n, n > 0 required')
    return [0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1)), 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))]

def literal_binomial_cdf(k, n, p):
    return sum((comb(n, j) * p ** j * (1 - p) ** (n - j) for j in range(k + 1)))

def cp_control(k, n, alpha=0.05):
    lo = 0.0 if k == 0 else brentq(lambda p: literal_binomial_cdf(k - 1, n, p) - (1 - alpha / 2), 0.0, 1.0, xtol=1e-14)
    hi = 1.0 if k == n else brentq(lambda p: literal_binomial_cdf(k, n, p) - alpha / 2, 0.0, 1.0, xtol=1e-14)
    return [lo, hi]

def counts(source):
    rows = []
    for g in CATEGORIES:
        k = source['CBCT_completed_counts']['fallback_RCT'][g]
        n = k + source['CBCT_completed_counts']['P'][g]
        alternatives = [(k, n), (k + 1, n + 1)]
        rows.append(dict(category=g, k=k, n=n, p=k / n, CI95=cp(k, n), simultaneous_CI95=cp(k, n, 0.05 / 3), missing_point_envelope=[min((a / b for (a, b) in alternatives)), max((a / b for (a, b) in alternatives))], missing_simultaneous_CI95_envelope=[min((cp(a, b, 0.05 / 3)[0] for (a, b) in alternatives)), max((cp(a, b, 0.05 / 3)[1] for (a, b) in alternatives))], resolution='POPULATION', subject_prediction='UNKNOWN', source_locator='doi:10.1111/iej.14144 Table 4'))
    return rows

def r1_result(source, external):
    rows = counts(source)
    checks = []
    for r in rows:
        for a in [0.05, 0.05 / 3]:
            checks.append(max((abs(x - y) for (x, y) in zip(cp(r['k'], r['n'], a), cp_control(r['k'], r['n'], a)))))
    pooled = 25 / 86
    ext = []
    for n in external['confidence_denominators']:
        ext.append(dict(k=1, n=n, p=1 / n, CI95=cp(1, n), frozen_pooled_probability=pooled, exact_binomial_p=float(binomtest(1, n, pooled).pvalue), two_cohort_Fisher_p=float(fisher_exact([[25, 61], [1, n - 1]]).pvalue), binomial_screen='REJECT_FIXED_PRIMARY_PROBABILITY' if binomtest(1, n, pooled).pvalue < 0.05 else 'NOT_REFUTED', inference='Descriptive transport challenge only: endpoint/protocol/selection not fully matched.'))
    witness = {'states': ['widening', 'lesion'], 'summary': [1, 1], 'summary_name': 'CBCT abnormal flag', 'identity_error': 0, 'machine_identical': True, 'downstream_risk': [rows[1]['p'], rows[2]['p']], 'downstream_difference': rows[2]['p'] - rows[1]['p'], 'downstream_difference_exact': str(Fraction(9, 19) - Fraction(10, 27)), 'resolution': 'POPULATION', 'external_referent': 'doi:10.1111/iej.14144 Table 4', 'minimum_extension': 'Retain widening versus lesion (three-category CBCT label); not a causal patient-state proof.'}
    return dict(round='R1', claim_type='information_link', rows=rows, pooled=dict(k=25, n=86, p=pooled, CI95=cp(25, 86), resolution='POPULATION'), missing=dict(n=1, fraction=1 / 86, known_conversion=True, handling='Retain all 3 possible category assignments; no imputation.'), sufficiency=witness, interval_control_max_error=max(checks), external_comparisons=ext, external_conditional_CBCT_calibration='UNKNOWN', empirical_category_times_min=None, clinical_treatment_choice=None, gates={'source_counts': 'PASS', 'CP95_control': 'PASS' if max(checks) <= 1e-10 else 'FAIL', 'binary_CBCT_summary_sufficient': 'FAIL', 'external_transport': 'UNKNOWN_UNMATCHED_SELECTION_AND_ENDPOINT', 'source_CBCT_timing_identified': 'UNKNOWN'}, inherited_negative_gate={'locator': 'LANE_XBREAK_HUNT_4/raw/R6_RESULTS.json, independently checked in LANE_XREVIEW_BATCH27', 'robust_LOO_relative_Brier_improvement_ge_5percent': 'FAIL (4.221-7.023% across missing allocations); not reclassified.'})

def exact_scenario(p, base, extra):
    """Enclose affine arithmetic exactly over declared boxes using rationals.

    Probability endpoints are exact encodings of the computed floats, not a
    proof of biological coverage or transcendental beta-quantile rounding.
    """
    ps = [Fraction.from_float(float(x)) for x in p]
    bs = [Fraction(str(x)) for x in base]
    es = [Fraction(str(x)) for x in extra]
    corners = [b + q * e for (b, q, e) in product(bs, ps, es)]
    (lo, hi) = (min(corners), max(corners))
    return dict(lower=float(np.nextafter(float(lo), -np.inf)), upper=float(np.nextafter(float(hi), np.inf)), exact_lower=str(lo), exact_upper=str(hi), corner_values=[str(x) for x in corners], guarantee='Exact rational corner enclosure on supplied box; no rigorous beta-quantile rounding enclosure supplied.', resolution='PHENOMENOLOGICAL')

def poisson_binomial_exact(probabilities):
    pmf = [Fraction(1)]
    for p in probabilities:
        q = Fraction(p)
        new = [Fraction(0)] * (len(pmf) + 1)
        for (j, v) in enumerate(pmf):
            new[j] += v * (1 - q)
            new[j + 1] += v * q
        pmf = new
    assert sum(pmf) == 1
    return pmf

def enumeration_control(probabilities):
    pmf = [Fraction(0)] * (len(probabilities) + 1)
    for bits in product([0, 1], repeat=len(probabilities)):
        value = Fraction(1)
        for (b, p) in zip(bits, probabilities):
            value *= p if b else 1 - p
        pmf[sum(bits)] += value
    return pmf

def tail(pmf, reserve):
    return sum(pmf[reserve + 1:], Fraction(0))

def reserve_quantile(pmf, tail_limit=Fraction(1, 20)):
    return next((j for j in range(len(pmf)) if tail(pmf, j) <= tail_limit))

def timing_sufficiency():
    worlds = [dict(normal=[40, 60], lesion=[100, 120]), dict(normal=[100, 120], lesion=[40, 60])]
    summaries = []
    differences = []
    for w in worlds:
        values = np.array(sorted(w['normal'] + w['lesion']), dtype=np.float64)
        summaries.append([float(values.mean()), float(values.std(ddof=1))])
        differences.append(float(np.mean(w['normal']) - np.mean(w['lesion'])))
    err = max((abs(x - y) for (x, y) in zip(*summaries)))
    assert summaries[0] == summaries[1] and err == 0.0
    return dict(worlds=worlds, summaries=summaries, summary='Pooled mean, sample SD, counts and exact time multiset', identity_error_min=err, machine_identical=True, downstream_normal_minus_lesion_min=differences, downstream_difference_min=abs(differences[0] - differences[1]), minimum_extension='Category-conditioned time sum/first moment for this mean query; individual paired times for prediction intervals.', external_referent=dict(kind='our_own_fixture', locator='code/model.py timing_sufficiency', compared_quantity='summary sufficiency, not clinical timing', refutes_us=False), resolution='PHENOMENOLOGICAL', source_moments_reproduced=False)

def r2_result(r1, base_min, extra_min, base_cost, extra_cost):
    scenarios = []
    for r in r1['rows']:
        box = r['missing_simultaneous_CI95_envelope']
        scenarios.append(dict(category=r['category'], source_probability_box=box, time=exact_scenario(box, base_min, extra_min), cost=exact_scenario(box, base_cost, extra_cost), point_time=base_min[0] + r['p'] * extra_min[0], point_cost=base_cost[0] + r['p'] * extra_cost[0], empirical_category_time=None, target_site_risk=None))
    panel = ['normal'] * 4 + ['widening'] * 2 + ['lesion'] * 2
    lookup = {r['category']: r for r in r1['rows']}
    ps = [Fraction(lookup[g]['k'], lookup[g]['n']) for g in panel]
    lo = [Fraction.from_float(lookup[g]['missing_simultaneous_CI95_envelope'][0]) for g in panel]
    hi = [Fraction.from_float(lookup[g]['missing_simultaneous_CI95_envelope'][1]) for g in panel]
    mid = poisson_binomial_exact(ps)
    low = poisson_binomial_exact(lo)
    high = poisson_binomial_exact(hi)
    control = enumeration_control(ps)
    pooled = poisson_binomial_exact([Fraction(25, 86)] * len(panel))
    q = reserve_quantile(high)
    pointq = reserve_quantile(mid)
    panel_out = dict(categories=panel, n=len(panel), mean_conversions_exact=str(sum(ps)), mean_conversions=float(sum(ps)), exact_PMF=[str(x) for x in mid], P_at_least_one=float(1 - mid[0]), point_reserve95=pointq, uncertainty_box_reserve95=q, pooled_point_reserve95=reserve_quantile(pooled), risk_exceeds_point_reserve_box=[float(tail(low, pointq)), float(tail(high, pointq))], risk_exceeds_robust_reserve_box=[float(tail(low, q)), float(tail(high, q))], upper_tail_exact_at_robust_reserve=str(tail(high, q)), enumeration_control_error_exact=str(max((abs(a - b) for (a, b) in zip(mid, control)))), resolution='PHENOMENOLOGICAL', source_probability_resolution='POPULATION', assumptions='Independent future conversions conditional on fixed category probabilities inside input box. No patient/site transport guarantee. 95% tail conditional on box; 95% confidence box is separate, not joint 95% coverage.', rigorous_probability_box_enclosure='Monotone coupling Y_i=1[U_i<p_i] encloses any exceedance probability by endpoint PMFs; coefficients/tails are exact rationals. Computed CP box itself has no rigorous special-function rounding enclosure.')
    same_mean = [Fraction(1, 2)] * 8
    different = [Fraction(1)] * 4 + [Fraction(0)] * 4
    mean1 = sum(same_mean)
    mean2 = sum(different)
    pmf1 = poisson_binomial_exact(same_mean)
    pmf2 = poisson_binomial_exact(different)
    suff = dict(states=[['1/2'] * 8, ['1'] * 4 + ['0'] * 4], summary_expected_conversions=[str(mean1), str(mean2)], identity_error=0, machine_identical=True, downstream_reserve95=[reserve_quantile(pmf1), reserve_quantile(pmf2)], downstream_tail_above4=[float(tail(pmf1, 4)), float(tail(pmf2, 4))], minimum_extension='For this pair, Bernoulli variance distinguishes the states; for all reserve quantiles retain the PMF or individual probabilities, mean+variance is not claimed sufficient.', external_referent=dict(kind='closed_form', locator='Product of Bernoulli probability-generating functions, explicit exact enumeration in code/model.py', compared_quantity='Poisson-binomial exceedance distribution', refutes_us=False), clinical_external_facit=False, resolution='PHENOMENOLOGICAL')
    return dict(round='R2', claim_type='information_link', scenarios=scenarios, panel=panel_out, timing_sufficiency=timing_sufficiency(), reserve_sufficiency=suff, gates={'scenario_arithmetic_enclosure': 'PASS', 'exact_reserve_control': 'PASS' if mid == control else 'FAIL', 'clinical_time_cost_prediction': 'UNKNOWN', 'target_site_transport': 'UNKNOWN', 'expected_conversions_sufficient_for_reserve': 'FAIL'}, planning_change='Executable conditional time/cost and extra-procedure-count scenarios, with study/protocol/endpoint retained; no inferred treatment choice or booking instruction.', minimum_new_measurement='Same-tooth local preoperative CBCT category, assignment, protocol, bleeding outcome, final treatment, first/all visit durations and cost; then held-out-site calibration.')
