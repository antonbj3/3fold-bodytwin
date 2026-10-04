"""Exact, explicitly conditional completion of the reviewed partial force vector."""
import argparse
from fractions import Fraction as Q
from common import *

def query(tooth=18, threshold=Q(7, 2)):
    net = read(HERE / 'CONSTRAINT_NET_DENTAL_R2.json')
    e = next((e for e in net['edges'] if e['id'] == 'D-E-TANDLAST-PARTIAL-REFERENCE'))
    known = sum((Q(str(ev_value(r))) for r in e['gap']['operands']), Q(0))
    remainder = Q(100) - known
    if Q(str(e['gap']['value'])) != remainder:
        raise ValueError('WRONG_REMAINDER')
    if tooth not in e['missing_fdi']:
        raise ValueError('Query must name an unobserved coordinate, 18 or28')
    other = next((t for t in e['missing_fdi'] if t != tooth))
    a = {tooth: Q(0), other: remainder}
    b = {tooth: remainder, other: Q(0)}
    assert sum(a.values()) == sum(b.values()) == remainder
    decision = 'YES' if threshold < 0 else 'NO' if threshold >= remainder else 'UNKNOWN'
    return {'claim_type': 'capability', 'known_sum_pp': float(known), 'remaining_pp': float(remainder), 'query': f'Is unobserved coordinate FDI{tooth} greater than {float(threshold)} percentage points?', 'answer': decision, 'coordinate_interval_pp': [0, float(remainder)], 'endpoint_witnesses': [{str(k): float(v) for (k, v) in s.items()} for s in [a, b]], 'summary_identity_error_pp': 0, 'bitwise_identical': float(sum(a.values())).hex() == float(sum(b.values())).hex(), 'downstream_difference_pp': float(abs(a[tooth] - b[tooth])), 'resolution_level': 'POPULATION', 'anatomical_support_level': 'PER_TOOTH', 'minimal_extension': 'One independent named force-share observation (or a justified allocation constraint); conservation then fixes the other coordinate.', 'rigorous_enclosure': 'Exact rational simplex endpoints under nonnegativity and normalization100pp. Does not prove these completions mechanically realizable or measured.', 'source_referent': {'kind': 'external_review', 'locator': str(DENT / 'results/LANE_XREVIEW_BATCH18/REVIEW_SOL_FALT_TANDLAST_20261003.json'), 'compared_quantity': '14 observed FDI sum93pp on16-FDI support; unobserved mass7pp', 'refutes_us': True}, 'endpoint_witness_referent': {'kind': 'our_own_fixture', 'locator': 'residual_query.py:query', 'compared_quantity': 'Algebraically feasible residual allocations only', 'refutes_us': True}, 'physical_allocation': 'UNKNOWN', 'same_subject_transfer': 'UNKNOWN', 'measurement_performed': False}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tooth', type=int, default=18)
    ap.add_argument('--threshold', default='3.5')
    a = ap.parse_args()
    r = query(a.tooth, Q(a.threshold))
    dump('RESIDUAL_QUERY_R2B.json', r)
    print(json.dumps(r, ensure_ascii=False, indent=2))
if __name__ == '__main__':
    main()
