"""Exact original-identity binding, including the batch21 source-substitution correction."""
import copy
from fractions import Fraction as Q
from common import *
from extract_sources import extract

def check_record(record, expected):
    errors = []
    for field in expected:
        if not exact_equal(record.get(field), expected[field]):
            errors.append('FROZEN_SOURCE_BINDING:' + field)
    if set(record) != set(expected):
        errors.append('RECORD_FIELD_CENSUS')
    return errors

def check_proposal(p, original):
    errors = []
    for k in ['source_file', 'source_sha256', 'key_path', 'edge_id', 'unit', 'measured', 'value']:
        if type(p.get(k)) != type(original.get(k)) or p.get(k) != original.get(k):
            errors.append('FROZEN_PROPOSAL_BINDING:' + k)
    return errors

def run():
    obs = read('SOURCE_OBSERVATIONS.json')
    expected = extract(False)
    errors = []
    tests = []
    if len(obs['records']) != 3:
        errors.append('SOURCE_CENSUS')
    for (r, e) in zip(obs['records'], expected['records']):
        errors += check_record(r, e)
    (a, b, c) = obs['records']
    arithmetic = {'Sagheb': Q(a['reported']['first_mean_N']) - Q(a['reported']['tenth_mean_N']) == Q('76.2'), 'Cunali': Q(b['reported']['replica_mean_um']) - Q(b['reported']['ct_mean_um']) == Q('60.82'), 'Choi_TX': sum((x['group'] == 'TX' and x['max_N'] == 285 and (x['events'] == [False] * 3) for x in c['reported']['fatigue_rows'])) == 1, 'Choi_EV': sum((x['group'] == 'EV' and x['max_N'] == 316 and (x['events'] == [False] * 3) for x in c['reported']['fatigue_rows'])) == 1}
    if not all(arithmetic.values()):
        errors.append('DIRECT_TABLE_ARITHMETIC')
    for (i, e) in enumerate(expected['records']):
        for (label, mut) in [('source_path', lambda r: r['reference'].__setitem__('source_file', 'another_same_value_copy.json')), ('source_hash', lambda r: r['reference'].__setitem__('source_sha256', '0' * 64)), ('table_locator', lambda r: r['reference'].__setitem__('locator', 'wrong table')), ('population', lambda r: r.__setitem__('population', {'n': 999})), ('protocol', lambda r: r['protocol'].__setitem__('id', 'different')), ('wrong_quantity', lambda r: r.__setitem__('quantity', 'marginal_gap')), ('unit', lambda r: r.__setitem__('unit', 'mm')), ('scientific_promotion', lambda r: r.__setitem__('full_edge_validated', True)), ('value', lambda r: r.__setitem__('reported', {'value': 801})), ('resolution', lambda r: r.__setitem__('resolution_level', 'PER_POINT'))]:
            p = copy.deepcopy(e)
            mut(p)
            errs = check_record(p, e)
            tests.append({'record': e['id'], 'injected_fault': label, 'rejected': bool(errs), 'errors': errs})
    original = source('original_proposals')['proposals']
    delivered = source('r3_proposals')['rows']
    proposal_ok = []
    for (row, p) in zip(delivered, original):
        proposal_ok.append(not check_proposal(row['proposal'], p))
    if len(original) != 25 or len(delivered) != 25 or (not all(proposal_ok)):
        errors.append('R3_ORIGINAL_PROPOSAL_BINDING')
    for i in [9, 10]:
        p = copy.deepcopy(delivered[i]['proposal'])
        p['source_file'] = 'frozen_packaged_X58_alias_with_identical_value_and_hash.json'
        err = check_proposal(p, original[i])
        tests.append({'record': delivered[i]['id'], 'injected_fault': 'batch21_same_value_source_substitution', 'rejected': bool(err), 'errors': err})
    out = {'gate': 'R4A', 'pass': not errors and all(arithmetic.values()) and all((t['rejected'] for t in tests)), 'errors': errors, 'arithmetic_control': arithmetic, 'new_source_bindings': 3, 'full_edge_promotions': 0, 'r3_original_proposal_identity_checks': sum(proposal_ok), 'r3_scope': 'Original proposal identity pinned to frozen R3 audit; underlying R3 numerical findings inherited from independent batch21 review. No claim to rerun all prior physics.', 'fault_controls': tests, 'fault_count': len(tests), 'external_referent': {'kind': 'external_review', 'locator': entry('r3_review')['original_path'], 'compared_quantity': 'Exact original source identity guard for P10/P11', 'refutes_us': True}}
    write('SOURCE_GATE.json', out)
    state('R4A_COMPLETE' if out['pass'] else 'R4A_FAILED', 'SOURCE_GATE.json', 'Construct context-dependent ranges and test summary insufficiency')
    return out
if __name__ == '__main__':
    print({k: v for (k, v) in run().items() if k not in ['fault_controls']})
