from fractions import Fraction as Q
from common import *

def build(net):
    rows = []
    for (id, v) in net['variables'].items():
        es = [e for e in net['edges'] if e['status'] in ['OPEN', 'UNKNOWN'] and id in e['between']]
        terms = []
        for e in es:
            g = e['gap']
            value = None
            if g['value'] is None:
                continue
            if e.get('comparison_error'):
                value = Q(str(e['comparison_error']['value'])) / Q(str(e['comparison_tolerance']['value']))
            elif g.get('normalization_scale'):
                value = Q(g['exact_decimal']) / Q(str(g['normalization_scale']['value']))
            if value is not None:
                terms.append({'edge': e['id'], 'normalized_exact': str(value), 'decision_count': e['decision_count'], 'weighted_exact': str(value * e['decision_count'])})
        rows.append({'variable': id, 'resolution_level': v['resolution_level'], 'open_unknown_edge_ids': [e['id'] for e in es], 'numeric_edges': [e['id'] for e in es if e['gap']['value'] is not None], 'physical_blocker_edges': [e['id'] for e in es if e['gap']['kind'] == 'BLOCKED_ON_PHYSICAL_MEASUREMENT'], 'nonphysical_exception_edges': [e['id'] for e in es if e['gap']['value'] is None and e['gap']['kind'] != 'BLOCKED_ON_PHYSICAL_MEASUREMENT'], 'measurement_ids': sorted({m for e in es for m in e.get('measurement_obligation', {}).get('measurement_ids', [])}), 'inherited_normalized_terms': terms, 'partial_normalized_sum_exact': str(sum((Q(t['weighted_exact']) for t in terms), Q())) if terms else None, 'global_rank_identified': False})
    out = {'schema': 'constraint-stress-map-r3', 'claim_type': 'capability', 'rows': rows, 'ranking': 'NO_PHYSICAL_TOTAL_ORDER', 'explanation': 'Partial normalized source-criterion sums preserved only as diagnostic context. Missing physical observations are grouped into typed AND obligations, never assigned zero stress.', 'lab_groups': 'LAB_MEASUREMENT_COVERAGE.json', 'sufficiency': 'SUFFICIENCY.json', 'rigorous_enclosure': 'Exact rational arithmetic for the stated normalized scores and finite set counts; physical model enclosure absent. No affine sensitivity reported.'}
    dump('CONSTRAINT_STRESS_MAP_DENTAL_R3.json', out)
    return out
if __name__ == '__main__':
    build(read(HERE / 'CONSTRAINT_NET_DENTAL_R3.json'))
