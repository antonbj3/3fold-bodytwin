"""Controls use real inputs and each has a deliberately rejected mutant."""
from copy import deepcopy
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import product
import re
from bs4 import BeautifulSoup
import numpy as np
from scipy.stats import poisson_binom
from model import cp_control, exact_scenario, poisson_binomial_exact, enumeration_control, tail

def table_check(source, html):
    soup = BeautifulSoup(html, 'html.parser')
    names = {'Normal': 'normal', 'Abnormal (slight enlargement)': 'widening', 'Abnormal (Lesion)': 'lesion'}
    actual = {}
    mismatches = []
    for tr in soup.select('tr'):
        cols = [x.get_text(' ', strip=True) for x in tr.find_all(['th', 'td'], recursive=False)]
        if cols and cols[0] in names:
            nums = [int(re.match('\\d+', s).group()) for s in cols[1:5]]
            actual[names[cols[0]]] = dict(total_reported=nums[0], P=nums[1], direct_RCT=nums[2], fallback_RCT=nums[3])
            if nums[0] != sum(nums[1:]):
                mismatches.append(dict(category=names[cols[0]], reported_total=nums[0], sum_group_cells=sum(nums[1:])))
    ok = len(actual) == 3
    for (g, r) in actual.items():
        for b in ['P', 'fallback_RCT', 'direct_RCT']:
            ok = ok and source['CBCT_completed_counts'][b][g] == r[b]
    ok = ok and source['N_completed']['P'] + source['N_completed']['fallback_RCT'] == 86
    ok = ok and sum(source['CBCT_completed_counts']['fallback_RCT'].values()) == 25
    return (ok, dict(parsed=actual, source_total_discrepancies=mismatches, handling='Use assigned 86 and chosen within-treatment cells, not Table4 total header171 or normal total72. Preserve source discrepancies.'))

def independent_check(external, xml):
    soup = BeautifulSoup(xml, 'xml')
    pars = [p.get_text(' ', strip=True) for p in soup.find_all('p')]
    group = next((p for p in pars if p.startswith('In Group II')))
    ok = external['n_randomized_full'] == 33 and external['n_full_received'] == 31 and (external['n_uncontrolled_bleeding'] == 1) and (external['n_partial_pulp_necrosis'] == 1) and ('except one' in group) and ('partially necrosed' in group) and ('10 min' in group) and (external['confirmed_RCT_conversion_counts'] is None)
    return ok

def atlas_safe(atlas):
    return all((r['CBCT_periapical_disease_category'] is None and r['conversion_probability'] is None and (not r['same_subject_as_RCT']) and (r['join_kind'] == 'ANATOMY_CONTEXT_ONLY') for r in atlas['rows']))

def decimal_corners(p, base, extra):
    with localcontext() as ctx:
        ctx.prec = 90
        return [Decimal(str(b)) + Decimal.from_float(float(q)) * Decimal(str(e)) for (b, q, e) in product(base, p, extra)]

def scenario_check(result, p, base, extra):
    vals = decimal_corners(p, base, extra)
    lo = Decimal.from_float(result['lower'])
    hi = Decimal.from_float(result['upper'])
    return all((lo <= v <= hi for v in vals))

def run(source, external, table_html, xml, atlas, r1, r2, base_min, extra_min, base_cost, extra_cost):
    rows = []
    (good, detail) = table_check(source, table_html)
    mutant = deepcopy(source)
    mutant['CBCT_completed_counts']['fallback_RCT']['normal'] += 1
    rows.append(dict(name='source_table_counts', valid_pass=good, mutant_rejected=not table_check(mutant, table_html)[0], detail=detail))
    good = independent_check(external, xml)
    mutant = deepcopy(external)
    mutant['n_uncontrolled_bleeding'] = 2
    rows.append(dict(name='independent_endpoint', valid_pass=good, mutant_rejected=not independent_check(mutant, xml)))
    interval = r1['rows'][0]['CI95']
    expected = cp_control(5, 39)
    good = max((abs(a - b) for (a, b) in zip(interval, expected))) <= 1e-10
    mutant = [interval[0] + 0.01, interval[1]]
    rows.append(dict(name='CP_interval', valid_pass=good, mutant_rejected=max((abs(a - b) for (a, b) in zip(mutant, expected))) > 1e-10))
    p = r1['rows'][0]['missing_simultaneous_CI95_envelope']
    checks = []
    for (b, e) in [(base_min, extra_min), (base_cost, extra_cost)]:
        env = exact_scenario(p, b, e)
        good = scenario_check(env, p, b, e)
        mutant = deepcopy(env)
        mutant['upper'] = env['lower'] - 1.0
        checks.append(dict(valid_pass=good, mutant_rejected=not scenario_check(mutant, p, b, e)))
    rows.append(dict(name='scenario_corner_enclosure', valid_pass=all((x['valid_pass'] for x in checks)), mutant_rejected=all((x['mutant_rejected'] for x in checks))))
    ps = [Fraction(5, 39)] * 4 + [Fraction(10, 27)] * 2 + [Fraction(9, 19)] * 2
    pmf = poisson_binomial_exact(ps)
    oracle = enumeration_control(ps)
    mutant = pmf.copy()
    mutant[0] += Fraction(1, 1000)
    rows.append(dict(name='contingency_PMF', valid_pass=pmf == oracle, mutant_rejected=mutant != oracle))
    published = poisson_binom.pmf(np.arange(9), np.array([float(p) for p in ps]))
    error = max((abs(float(x) - y) for (x, y) in zip(pmf, published)))
    mutant_error = max((abs(float(x) - y) for (x, y) in zip(mutant, published)))
    rows.append(dict(name='published_scipy_PMF', valid_pass=error <= 1e-12, mutant_rejected=mutant_error > 1e-12, max_error=error, external_referent=dict(kind='published_code', locator='https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.stats.poisson_binom.html', compared_quantity='source panel PMF', refutes_us=False)))
    lookup = {r['category']: r for r in r1['rows']}
    probs = [Fraction.from_float(lookup[g]['missing_simultaneous_CI95_envelope'][1]) for g in r2['panel']['categories']]
    high = poisson_binomial_exact(probs)
    q = r2['panel']['uncertainty_box_reserve95']
    good = tail(high, q) <= Fraction(1, 20) and (q == 0 or tail(high, q - 1) > Fraction(1, 20))
    rows.append(dict(name='reserve_tail', valid_pass=good, mutant_rejected=tail(high, 0) > Fraction(1, 20)))
    good = atlas_safe(atlas)
    mutant = deepcopy(atlas)
    mutant['rows'][0]['conversion_probability'] = 0.29
    rows.append(dict(name='atlas_risk_join', valid_pass=good, mutant_rejected=not atlas_safe(mutant)))
    w = r2['timing_sufficiency']
    good = w['summaries'][0] == w['summaries'][1]
    mutant = deepcopy(w)
    mutant['summaries'][1][0] += 0.0001
    rows.append(dict(name='machine_identity', valid_pass=good, mutant_rejected=mutant['summaries'][0] != mutant['summaries'][1]))
    assert all((r['valid_pass'] and r['mutant_rejected'] for r in rows)), rows
    return dict(checks=rows, n_checks=len(rows), n_valid_pass=int(sum((r['valid_pass'] for r in rows))), n_mutants_rejected=int(sum((r['mutant_rejected'] for r in rows))), independent_scientific_review=False)
