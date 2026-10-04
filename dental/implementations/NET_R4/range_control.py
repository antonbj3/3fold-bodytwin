"""Independent numerical oracle and boundary tests, not a competing algorithm."""
from decimal import Decimal
import copy
from common import *
from validate import check_observation

def run(n):
    tests = []
    for (k, v) in n['variables'].items():
        for r in v['plausible_range']['cases']:
            for (label, val, expected) in [('lower', r['lower'], True), ('upper', r['upper'], True), ('below', str(Decimal(r['lower']) - 1), False), ('above', str(Decimal(r['upper']) + 1), False)]:
                o = {'variable': k, 'case_id': r['id'], 'value': float(val), 'unit': r['unit'], 'reference': r['source'], **{f: r[f] for f in ['quantity', 'protocol_id', 'resolution_level', 'value_role', 'context']}}
                oracle = Decimal(r['lower']) <= Decimal(str(o['value'])) <= Decimal(r['upper'])
                got = not check_observation(n, o)
                tests.append({'case': r['id'], 'value': float(val), 'position': label, 'candidate_accept': got, 'independent_decimal_accept': oracle, 'expected_accept': expected, 'pass': got == oracle == expected})
    for (k, v) in n['variables'].items():
        r = v['plausible_range']
        if r['status'] == 'MATHEMATICAL_DOMAIN':
            for (val, expected) in [(float(r['lower']), True), (float(r['upper']), True), (float(r['upper']) + 1, False)]:
                o = {'variable': k, 'value': val, 'unit': v['unit'], 'quantity': v['desc'], 'resolution_level': v['resolution_level'], 'value_role': 'mathematical_domain_only'}
                got = not check_observation(n, o)
                oracle = Decimal(r['lower']) <= Decimal(str(val)) <= Decimal(r['upper'])
                tests.append({'case': k, 'value': val, 'candidate_accept': got, 'independent_decimal_accept': oracle, 'expected_accept': expected, 'pass': got == oracle == expected})
    ablation = []
    for (i, x) in [(0, 80100), (3, 801000)]:
        o = copy.deepcopy(n['r4_observations'][i])
        o['value'] = x
        old = o['unit'] == n['variables'][o['variable']]['unit']
        errs = check_observation(n, o)
        ablation.append({'observation': o['id'], 'injected_value': x, 'unit': o['unit'], 'unit_only_accepts': old, 'range_rejects': 'OUT_OF_PLAUSIBLE_RANGE' in errs, 'scope': 'Ablation of numeric intake, not a claim that complete R3 source-hash checks accept altered frozen source values.'})
    out = {'all_pass': all((t['pass'] for t in tests)) and all((x['range_rejects'] for x in ablation)), 'boundary_and_domain_tests': tests, 'count': len(tests), 'ablation': ablation, 'interpretation': 'Independent membership oracle matches; no algorithmic speed or accuracy gain claimed.'}
    write('RANGE_CONTROL.json', out)
    return out
