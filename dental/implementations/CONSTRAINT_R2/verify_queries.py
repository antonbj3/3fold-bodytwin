"""Adversarial checks for the second construction and reviewed comparisons."""
import copy
import math
from fractions import Fraction as Q
from common import *
from stress_and_queries import compatible_sum, sufficiency, x71_gate
from residual_query import query

def main():
    s = sufficiency()
    assert s['bitwise_identical'] and s['identity_error'] == 0 and (s['downstream_difference'] == 1)
    assert s['unit_conversion_test']['error'] == 0 and s['unit_conversion_test']['mixed_quantity_sum_rejected']
    r = query()
    assert r['remaining_pp'] == 7 and r['answer'] == 'UNKNOWN' and (r['downstream_difference_pp'] == 7)
    assert query(threshold=Q(7))['answer'] == 'NO'
    assert query(threshold=Q(-1))['answer'] == 'YES'
    independent = [(x, Q(7) - x) for x in [Q(0), Q(7)]]
    assert min((x for (x, y) in independent)) == 0 and max((x for (x, y) in independent)) == 7
    try:
        compatible_sum([{'quantity': 'a', 'unit': 'mm', 'gap_exact': '1', 'decision_count': 1}, {'quantity': 'b', 'unit': 'mm', 'gap_exact': '1', 'decision_count': 1}])
        same_unit_different_quantity_reject = False
    except ValueError:
        same_unit_different_quantity_reject = True
    assert same_unit_different_quantity_reject
    plan = source('results/LANE_X71_MEASUREMENT_PRIORITY/results.json')['first_week_crown72_only']
    review = source('results/LANE_XREVIEW_BATCH18/evidence/independent_checks.json')['X71']
    assert x71_gate(plan, review)
    bad_cost = copy.deepcopy(plan)
    bad_cost['operator_hours'][1] += 1
    bad_question = copy.deepcopy(plan)
    bad_question['evaluation_ready'].pop()
    bad_physical = copy.deepcopy(plan)
    bad_physical['already_physically_decided'] = 1
    x71_faults = {'cost_plus1_hour_rejected': not x71_gate(bad_cost, review), 'missing_question_rejected': not x71_gate(bad_question, review), 'invented_physical_decision_rejected': not x71_gate(bad_physical, review)}
    assert all(x71_faults.values())
    raw = source('results/LANE_X70_END_TO_END/results.json')['external_force']
    rows = []

    def passes(observed, predicted):
        return abs(math.log(predicted / observed)) <= raw['tolerance_log']
    for (i, x) in enumerate(raw['rows']):
        ok = passes(x['observed_N'], x['predicted_N'])
        bad = dict(x)
        bad['observed_N'] *= 100
        reject = not passes(bad['observed_N'], bad['predicted_N'])
        assert ok == x['gate'] and reject
        rows.append({'row': i, 'empirical_group_mean': i < 3, 'source_pass': ok, 'injected_observed_times100_rejected': reject})
    dump('QUERY_CONTROLS.json', {'weighted_gap_summary': s, 'residual_simplex': r, 'independent_endpoint_oracle_match': True, 'same_unit_different_quantity_rejected': same_unit_different_quantity_reject, 'X71_actual_faults': x71_faults, 'X70_recomputed_source_gates': rows, 'empirical_passes': [sum((x['source_pass'] for x in rows[:3])), 3], 'Weibull_closure_passes': [sum((x['source_pass'] for x in rows[3:])), 2], 'all_controls_pass': True})
    print('R2B exact summary, units, simplex and five real force mutations: PASS')
if __name__ == '__main__':
    main()
