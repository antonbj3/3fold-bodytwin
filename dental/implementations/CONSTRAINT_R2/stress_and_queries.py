"""Quantity-safe stress views and explicit failure of scalar sufficiency."""
from collections import defaultdict
from fractions import Fraction as Q
import hashlib
import json
from common import *

def normalized(e):
    g = e.get('gap', {})
    if g.get('value') is None:
        return None
    if e.get('comparison_error'):
        return (Q(str(e['comparison_error']['value'])) / Q(str(e['comparison_tolerance']['value'])), 'source_log_error_over_source_tolerance')
    if g.get('normalization_scale'):
        return (Q(g['exact_decimal']) / Q(str(g['normalization_scale']['value'])), 'absolute_residual_over_source_tolerance')
    return None

def compatible_sum(terms):
    units = {(t['quantity'], t['unit']) for t in terms}
    if len(units) > 1:
        raise ValueError('INCOMMENSURATE_QUANTITIES_OR_UNITS')
    return sum((Q(t['gap_exact']) * t['decision_count'] for t in terms), Q(0))

def x71_gate(plan, review):
    return sorted(plan['evaluation_ready']) == sorted(review['ready_question_ids']) and plan['operator_hours'] == review['cost_operator_h'] and (plan['already_physically_decided'] == review['physically_decided'] == 0)

def sufficiency():
    a = {'edge1': Q(2), 'edge2': Q(1)}
    b = {'edge1': Q(1), 'edge2': Q(2)}
    sa = sum(a.values())
    sb = sum(b.values())
    ra = a['edge2']
    rb = b['edge2']
    scale_a = Q(1, 2) / Q(1, 10)
    scale_b = Q(500) / Q(100)
    try:
        compatible_sum([{'quantity': 'length', 'unit': 'mm', 'gap_exact': '1', 'decision_count': 1}, {'quantity': 'force', 'unit': 'N', 'gap_exact': '1', 'decision_count': 1}])
        rejects = False
    except ValueError:
        rejects = True
    return {'summary': 'sum of gap times one decision per edge; same two edges and equal weights', 'state_A': {k: str(v) for (k, v) in a.items()}, 'state_B': {k: str(v) for (k, v) in b.items()}, 'summary_A': float(sa), 'summary_B': float(sb), 'identity_error': float(abs(sa - sb)), 'bitwise_identical': float(sa).hex() == float(sb).hex(), 'intervention': 'resolve edge1 completely', 'remaining_A': float(ra), 'remaining_B': float(rb), 'downstream_difference': float(abs(ra - rb)), 'decision_A': 'PASS', 'decision_B': 'FAIL', 'unit': 'source_tolerance_multiple', 'resolution_level': 'PHENOMENOLOGICAL', 'outcome': 'WEIGHTED_SUM_INSUFFICIENT', 'minimal_extension': 'For this fixed consumer retain gap by edge identity. For the network also retain consumer AND prerequisites, sign, specimen, region, protocol and time.', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'stress_and_queries.py:sufficiency', 'compared_quantity': 'weighted-score invariance versus remaining decision gap', 'refutes_us': True}, 'unit_conversion_test': {'before_exact': str(scale_a), 'after_exact': str(scale_b), 'error': float(abs(scale_a - scale_b)), 'mm_to_um_factor': 1000, 'mixed_quantity_sum_rejected': rejects}}

def build(net):
    unresolved = [e for e in net['edges'] if e['status'] in ['OPEN', 'UNKNOWN']]
    rows = []
    all_terms = []
    for (name, v) in net['variables'].items():
        incident = [e for e in unresolved if name in e['between']]
        bands = defaultdict(list)
        contributions = []
        missing = []
        not_normalized = []
        total = Q(0)
        for e in incident:
            g = e['gap']
            if g['value'] is None:
                missing.append(e['id'])
                continue
            term = {'edge': e['id'], 'quantity': g['quantity'], 'unit': g['unit'], 'gap_exact': g['exact_decimal'], 'gap': g['value'], 'gap_resolution': g['resolution_level'], 'decision_count': e['decision_count'], 'decision_ids': e['decision_ids'], 'scope': g['scope']}
            bands[g['quantity'], g['unit']].append(term)
            norm = normalized(e)
            if norm is None:
                not_normalized.append(e['id'])
            else:
                (value, metric) = norm
                contribution = value * e['decision_count']
                total += contribution
                contributions.append({'edge': e['id'], 'metric': metric, 'normalized_gap': float(value), 'decision_count': e['decision_count'], 'weighted_exact': str(contribution)})
        dimensional = []
        for ((qty, unit), ts) in bands.items():
            value = compatible_sum(ts)
            dimensional.append({'quantity': qty, 'unit': unit + ' * decision', 'weighted_sum': float(value), 'weighted_exact': str(value), 'terms': ts})
        control = Q(0)
        for e in unresolved:
            if name in e['between']:
                n = normalized(e)
                if n:
                    control += n[0] * len(set(e['decision_ids']))
        assert control == total
        rows.append({'variable': name, 'desc': v['desc'], 'resolution_level': v['resolution_level'], 'incident_open_unknown_edges': [e['id'] for e in incident], 'decision_ids': sorted({k for e in incident for k in e['decision_ids']}), 'dimensional_stress': dimensional, 'known_normalized_stress': float(total) if contributions else None, 'known_normalized_stress_exact': str(total) if contributions else None, 'contributions': contributions, 'unidentified_gap_edges': missing, 'no_source_scale_edges': not_normalized, 'normalized_edge_coverage': [len(contributions), len(incident)], 'global_rank_identified': False, 'rank_reason': 'Incomplete/unscaled gaps and distinct physical quantities prevent a physical total ranking. Score is a partial source-criterion workload view.', 'no_incident_open_edges': not incident})
    rows.sort(key=lambda r: (r['known_normalized_stress'] is None, -(r['known_normalized_stress'] or 0), r['variable']))
    for (i, r) in enumerate(rows):
        r['display_order'] = i + 1
    by_quantity = defaultdict(list)
    for r in rows:
        for q in r['dimensional_stress']:
            by_quantity[q['quantity'], q['unit']].append({'variable': r['variable'], 'weighted_sum': q['weighted_sum'], 'term_edges': [t['edge'] for t in q['terms']]})
    dim_ranks = []
    for ((qty, unit), rs) in sorted(by_quantity.items()):
        dim_ranks.append({'quantity': qty, 'unit': unit, 'ranking': sorted(rs, key=lambda x: (-x['weighted_sum'], x['variable'])), 'scope': 'Known contribution only; missing full-chain gaps do not contribute zero.'})
    old = source('results/PROOF_LANE_CONSTRAINT_NET_DENTAL/CONSTRAINT_STRESS_MAP_DENTAL.json')
    ports = {r['measurement']: r['quantity_variables'] for r in old['x71_comparison']['mapped_measurements']}
    x = source('results/LANE_X71_MEASUREMENT_PRIORITY/results.json')
    rev = source('results/LANE_XREVIEW_BATCH18/evidence/independent_checks.json')['X71']
    plan = x['first_week_crown72_only']
    ids = plan['ids']
    assert x71_gate(plan, rev)
    mapping = []
    for (idx, candidate) in enumerate(x['candidate_measurements']):
        m = candidate['measurement']
        vs = ports.get(m, [])
        es = [e for e in unresolved if set(vs) & set(e['between'])]
        mapping.append({'X71_rank': idx + 1, 'measurement': m, 'name': candidate['name'], 'variable_ids': vs, 'incident_edges': [e['id'] for e in es], 'missing_numeric_edges': [e['id'] for e in es if e['gap']['value'] is None], 'interpretation': 'Observational correspondence only. Touching a gap does not show this measurement resolves its unknown law or mismatched population.'})
    mapped_top = []
    for r in [r for r in rows if r['known_normalized_stress'] is not None][:10]:
        mapped_top.append({'variable': r['variable'], 'measurement_ids': [m for (m, vs) in ports.items() if r['variable'] in vs], 'known_partial_score': r['known_normalized_stress'], 'coverage': r['normalized_edge_coverage']})
    return {'schema': 'dental-constraint-stress/r2', 'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'resolution_level': 'PHENOMENOLOGICAL', 'metric': 'Within each named quantity/unit: sum(abs(scoped residual)*number of direct K decisions). Separate partial normalized view uses frozen source tolerances only.', 'ranking_status': 'PARTIAL_ONLY_FULL_RANK_NOT_IDENTIFIED', 'knowledge_debt': 'Measure matched per-edge residuals and decision-specific acquisition effect/cost; source tolerance multiples are not empirical utility.', 'all_variables_sorted': rows, 'rankings_by_quantity': dim_ranks, 'deduplication_warning': 'Inherited K edges may share an observation; this requested edge sum measures incidence, not independent evidence. No variance fusion.', 'x71_comparison': {'status': 'FINAL_REVIEWED_PLAN_COMPARED; NOT_A_SHARED_OPTIMUM', 'selected_measurements': ids, 'evaluable_questions': plan['evaluation_ready'], 'operator_hours': plan['operator_hours'], 'cost_resolution_level': 'PHENOMENOLOGICAL', 'physically_decided': 0, 'full_chains_closed': 0, 'independent_review_match': True, 'source': evidence('results/LANE_X71_MEASUREMENT_PRIORITY/results.json', '/first_week_crown72_only', 'operator_hour;decision', 'PHENOMENOLOGICAL'), 'measurement_mapping': mapping, 'top_partial_score_mapping': mapped_top, 'conclusion': 'The gap sum does not establish a competing measurement plan: many residuals are unidentified, and source populations are not X71 specimens. X71 prioritizes matched scan+film+fracture with prerequisites and cost.'}, 'sufficiency': sufficiency()}

def main():
    s = build(read(HERE / 'CONSTRAINT_NET_DENTAL_R2.json'))
    dump('CONSTRAINT_STRESS_MAP_DENTAL_R2.json', s)
    dump('SUFFICIENCY_R2B.json', s['sufficiency'])
    print(json.dumps({'quantity_rankings': len(s['rankings_by_quantity']), 'top_partial': [(r['variable'], r['known_normalized_stress']) for r in s['all_variables_sorted'][:5]], 'X71_review_match': s['x71_comparison']['independent_review_match'], 'summary_identity_error': s['sufficiency']['identity_error'], 'downstream_difference': s['sufficiency']['downstream_difference']}))
if __name__ == '__main__':
    main()
