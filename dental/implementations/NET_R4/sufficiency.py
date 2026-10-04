"""Exact adversarial states, explicitly not measurements."""
import copy, math
from fractions import Fraction as Q
from common import *

def summary(x):
    q = list(map(Q, x))
    mean = sum(q) / len(q)
    variance = sum(((a - mean) ** 2 for a in q)) / len(q)
    return [Q(len(q)), min(q), max(q), mean, variance]

def paired_witness(a, b):
    sa = summary(a['first']) + summary(a['tenth'])
    sb = summary(b['first']) + summary(b['tenth'])
    da = [Q(x) - Q(y) for (x, y) in zip(a['first'], a['tenth'])]
    db = [Q(x) - Q(y) for (x, y) in zip(b['first'], b['tenth'])]
    error = max((abs(x - y) for (x, y) in zip(sa, sb)))
    diff = abs(min(da) - min(db))
    return {'pass': error == 0 and diff > 0, 'state_A': a, 'state_B': b, 'summary_A_exact': [str(v) for v in sa], 'summary_B_exact': [str(v) for v in sb], 'identity_error': float(error), 'bitwise_identical': all((float(x).hex() == float(y).hex() for (x, y) in zip(sa, sb))), 'minimum_individual_loss_A_N': float(min(da)), 'minimum_individual_loss_B_N': float(min(db)), 'downstream_difference_N': float(diff), 'marginal_multisets_identical': all((sorted(a[k]) == sorted(b[k]) for k in a)), 'minimal_extension': 'Specimen ID joining first/tenth measurements and state/protocol; covariance alone only resolves variance, not tails.', 'resolution_level': 'PHENOMENOLOGICAL', 'uses_source_means_as_fixture_centres_only': True}

def distribution_witness(a, b):
    sa = summary(a)
    sb = summary(b)
    err = max((abs(x - y) for (x, y) in zip(sa, sb)))
    da = sum((x > 2 for x in a))
    db = sum((x > 2 for x in b))
    return {'pass': err == 0 and da != db, 'state_A': a, 'state_B': b, 'summary_names': ['n', 'min', 'max', 'mean', 'population_variance'], 'summary_A_exact': [str(x) for x in sa], 'summary_B_exact': [str(x) for x in sb], 'identity_error': float(err), 'bitwise_identical': all((float(x).hex() == float(y).hex() for (x, y) in zip(sa, sb))), 'downstream_count_above_2_A': da, 'downstream_count_above_2_B': db, 'downstream_difference': abs(da - db), 'minimal_extension': 'For this query retain exceedance count above 2; for arbitrary thresholds retain empirical CDF or raw values.', 'resolution_level': 'PHENOMENOLOGICAL', 'unit': 'fixture_units'}

def run():
    a = {'first': ['319.9', '339.9'], 'tenth': ['243.7', '263.7']}
    b = {'first': a['first'].copy(), 'tenth': list(reversed(a['tenth']))}
    p = paired_witness(a, b)
    d = distribution_witness([0, 0, 0, 3, 3, 4], [0, 0, 1, 1, 4, 4])
    changed = copy.deepcopy(b)
    changed['first'][0] = '319.8'
    tests = [{'fault': 'same pairing removes downstream difference', 'rejected': not paired_witness(a, a)['pass']}, {'fault': 'changed summary invalidates exact witness', 'rejected': not paired_witness(a, changed)['pass']}, {'fault': 'same distribution removes downstream difference', 'rejected': not distribution_witness(d['state_A'], d['state_A'])['pass']}, {'fault': 'different maximum invalidates summary identity', 'rejected': not distribution_witness(d['state_A'], [0, 0, 1, 1, 4, 5])['pass']}]
    r = read('SOURCE_OBSERVATIONS.json')['records'][0]['reported']
    s0 = Q(r['first_SD_N'])
    s1 = Q(r['tenth_SD_N'])
    out = {'claim_type': 'capability', 'pairing': p, 'range_and_moments': d, 'conditional_SD_difference_enclosure': {'lower_N': float(abs(s0 - s1)), 'upper_N': float(s0 + s1), 'proof': 'For centered vectors u,v, reverse and forward triangle inequalities give |sd(u)-sd(v)| <= sd(u-v) <= sd(u)+sd(v). Same paired n and sample normalization required.', 'assumptions': 'Printed SDs used exactly; source rounding and measurement error not enclosed. This does not bound individual losses or infer paired covariance.', 'resolution_level': 'POPULATION', 'source': read('SOURCE_OBSERVATIONS.json')['records'][0]['reference']}, 'external_referent': {'kind': 'our_own_fixture', 'locator': 'sufficiency.py:paired_witness and distribution_witness', 'compared_quantity': 'Exact equality of marginal summaries versus individual losses and tails', 'refutes_us': True}, 'controls': tests, 'all_pass': p['pass'] and d['pass'] and all((t['rejected'] for t in tests)), 'rigorous_enclosure': 'Exact rational arithmetic on fixture states; conditional norm enclosure for printed SDs only. Physical transfer enclosure MISSING. No affine sensitivity claimed.'}
    write('SUFFICIENCY.json', out)
    return out
if __name__ == '__main__':
    print(run()['all_pass'])
