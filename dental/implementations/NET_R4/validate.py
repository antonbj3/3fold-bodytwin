"""Fail closed scalar ingestion and exact inherited/scoped source contracts."""
import argparse, math, copy
from fractions import Fraction as Q
from common import *
from extract_sources import extract
from gate_sources import check_record, check_proposal
from build_net import observations

def finite(x):
    return type(x) in [int, float] and math.isfinite(x)

def check_observation(net, o):
    errors = []
    k = o.get('variable')
    v = net['variables'].get(k)
    if v is None:
        return ['UNKNOWN_VARIABLE']
    r = v.get('plausible_range', {})
    if not finite(o.get('value')):
        return ['NONFINITE_OR_NONNUMERIC_VALUE']
    if o.get('unit') != v['unit']:
        errors.append('UNIT_MISMATCH')
    if r.get('status') in ['UNKNOWN', 'PARTIAL_DOMAIN', 'NOT_SCALAR']:
        if r.get('lower') is not None and Q(str(o['value'])) < Q(r['lower']):
            errors.append('OUT_OF_DOMAIN')
        return errors + ['PLAUSIBLE_RANGE_UNKNOWN_OR_NOT_SCALAR']
    if r.get('status') == 'MATHEMATICAL_DOMAIN':
        if o.get('value_role') != 'mathematical_domain_only':
            errors.append('EMPIRICAL_PLAUSIBILITY_UNKNOWN')
        if o.get('quantity') != v['desc'] or o.get('resolution_level') != v['resolution_level']:
            errors.append('SCOPE_MISMATCH')
        chosen = r
    else:
        chosen = next((x for x in r.get('cases', []) if x['id'] == o.get('case_id')), None)
        if chosen is None:
            return errors + ['RANGE_CASE_UNKNOWN']
        for f in ['quantity', 'protocol_id', 'resolution_level', 'value_role', 'context']:
            if type(o.get(f)) != type(chosen.get(f)) or o.get(f) != chosen.get(f):
                errors.append('SCOPE_MISMATCH:' + f)
        if o.get('reference') != chosen['source']:
            errors.append('RANGE_SOURCE_MISMATCH')
    if chosen.get('lower') is None or chosen.get('upper') is None:
        return errors + ['FINITE_PLAUSIBLE_RANGE_MISSING']
    x = Q(str(o['value']))
    lo = Q(chosen['lower'])
    hi = Q(chosen['upper'])
    if x < lo or x > hi:
        errors.append('OUT_OF_PLAUSIBLE_RANGE')
    return errors

def validate(n):
    verify()
    p = source('r3_net')
    errs = []
    rc = read('RANGE_CONTRACTS.json')['variables']
    e0 = {e['id']: e for e in p['edges']}
    es = {e['id']: e for e in n['edges']}
    if set(es) != set(e0) or len(es) != len(n['edges']):
        errs.append('EDGE_CENSUS')
    for (id, e) in es.items():
        if id not in e0:
            continue
        if not exact_equal({k: v for (k, v) in e.items() if k != 'r4_reference_binding'}, e0[id]):
            errs.append('INHERITED_EDGE_CHANGED:' + id)
    for (k, v) in p.items():
        if k not in ['variables', 'edges', 'schema_version', 'seed', 'claim_type', 'review_state'] and (not exact_equal(n.get(k), v)):
            errs.append('INHERITED_TOP_LEVEL_CHANGED:' + k)
    if set(n['variables']) != set(rc):
        errs.append('VARIABLE_CENSUS')
    for (k, v) in n['variables'].items():
        if v.get('plausible_range') != rc.get(k):
            errs.append('FROZEN_RANGE_CHANGED:' + k)
        if k in p['variables'] and {f: x for (f, x) in v.items() if f != 'plausible_range'} != p['variables'][k]:
            errs.append('INHERITED_VARIABLE_CHANGED:' + k)
    expected = extract(False)['records']
    if len([e for e in es.values() if 'r4_reference_binding' in e]) != 3:
        errs.append('REFERENCE_CENSUS')
    for r in expected:
        b = es.get(r['edge_id'], {}).get('r4_reference_binding', {})
        errs += [r['edge_id'] + ':' + x for x in check_record(b.get('record', {}), r)]
        if b.get('supports_entire_edge') is not False or b.get('target_transfer_status') != 'UNKNOWN':
            errs.append('FALSE_EDGE_PROMOTION')
        if b.get('binding_resolution') != r['resolution_level'] or b.get('timescale') != r['timescale'] or b.get('reference_port') != r['variable']:
            errs.append('REFERENCE_BINDING_SCOPE')
        if b.get('review_state') != 'PENDING_INDEPENDENT_REVIEW' or b.get('evidence_status') != 'VERIFIED_PUBLISHED_OBSERVATION_SCOPED' or b.get('inherited_direction_is_not_validated_by_reference') is not True:
            errs.append('REFERENCE_STATUS_PROMOTION')
    expected_obs = observations(expected)
    if n.get('r4_observations') != expected_obs:
        errs.append('OBSERVATION_SOURCE_CONTRACT')
    checks = [{'id': o['id'], 'errors': check_observation(n, o)} for o in n.get('r4_observations', [])]
    for x in checks:
        errs += [x['id'] + ':' + e for e in x['errors']]
    physical = [e['id'] for e in es.values() if e.get('gap', {}).get('kind') == 'BLOCKED_ON_PHYSICAL_MEASUREMENT']
    nonphys = [e['id'] for e in es.values() if e['status'] in ['OPEN', 'UNKNOWN'] and e.get('gap', {}).get('value') is None and (e['id'] not in physical)]
    missing_ranges = [k for (k, v) in rc.items() if v['status'] in ['UNKNOWN', 'PARTIAL_DOMAIN']]
    return {'integrity_pass': not errs, 'errors': errs, 'range_checks': checks, 'range_observations_accepted': sum((not x['errors'] for x in checks)), 'range_observations_rejected': sum((bool(x['errors']) for x in checks)), 'physical_blocker_edges': physical, 'nonphysical_strict_failures': nonphys, 'numeric_plausibility_unknown_variables': missing_ranges, 'source_scoped_variables_not_universal': [k for (k, v) in rc.items() if v['cases']], 'full_physical_chain_validated': False, 'strict_ready': not errs and (not physical) and (not nonphys) and (not missing_ranges), 'review_state': 'PENDING_INDEPENDENT_REVIEW'}

def mutations(n):
    tests = []

    def obs_test(label, o, code):
        err = check_observation(n, o)
        tests.append({'fault': label, 'expected': code, 'errors': err, 'rejected': code in err})
    for o in n['r4_observations']:
        case = next((x for x in n['variables'][o['variable']]['plausible_range']['cases'] if x['id'] == o['case_id']))
        for (side, x) in [('above', float(case['upper']) + 1), ('below', float(case['lower']) - 1)]:
            q = copy.deepcopy(o)
            q['value'] = x
            obs_test(o['id'] + ' ' + side, q, 'OUT_OF_PLAUSIBLE_RANGE')
    for (field, value) in [('unit', 'wrong'), ('quantity', 'marginal_gap'), ('protocol_id', 'target_lab'), ('resolution_level', 'PER_POINT'), ('context', {}), ('reference', {})]:
        q = copy.deepcopy(n['r4_observations'][3])
        q[field] = value
        code = 'UNIT_MISMATCH' if field == 'unit' else 'RANGE_SOURCE_MISMATCH' if field == 'reference' else 'SCOPE_MISMATCH:' + field
        obs_test('wrong ' + field, q, code)
    for x in [True, None, float('nan'), float('inf'), '801']:
        q = copy.deepcopy(n['r4_observations'][0])
        q['value'] = x
        obs_test('invalid numeric ' + str(x), q, 'NONFINITE_OR_NONNUMERIC_VALUE')
    for (k, v) in n['variables'].items():
        if v['plausible_range']['status'] in ['UNKNOWN', 'PARTIAL_DOMAIN', 'NOT_SCALAR']:
            obs_test('unsupported scalar ' + k, {'variable': k, 'value': 801, 'unit': v['unit']}, 'PLAUSIBLE_RANGE_UNKNOWN_OR_NOT_SCALAR')
    for (i, value) in [(0, 80100), (3, 801000)]:
        q = copy.deepcopy(n['r4_observations'][i])
        q['value'] = value
        obs_test('unit-correct absurd magnitude ' + q['id'], q, 'OUT_OF_PLAUSIBLE_RANGE')
    for (label, mut, code) in [('range widened', lambda x: x['variables']['preload']['plausible_range']['cases'][0].__setitem__('upper', '100000'), 'FROZEN_RANGE_CHANGED:preload'), ('edge promoted', lambda x: x['edges'][0].__setitem__('status', 'VALIDATED'), 'INHERITED_EDGE_CHANGED:'), ('reference state', lambda x: next((e for e in x['edges'] if e['id'] == 'D-E-PRELOAD-HISTORY'))['r4_reference_binding']['record']['protocol'].__setitem__('fixture_separation_mm', '0'), 'FROZEN_SOURCE_BINDING:protocol')]:
        q = copy.deepcopy(n)
        mut(q)
        err = validate(q)['errors']
        tests.append({'fault': label, 'expected': code, 'errors': err, 'rejected': any((code in e for e in err))})
    return {'tests': tests, 'count': len(tests), 'all_rejected': all((x['rejected'] for x in tests))}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--net', default='CONSTRAINT_NET_DENTAL_R4.json')
    ap.add_argument('--strict', action='store_true')
    ap.add_argument('--mutations', action='store_true')
    ap.add_argument('--observation')
    a = ap.parse_args()
    n = read(a.net)
    if a.observation:
        base = validate(n)
        r = {'errors': base['errors'] + check_observation(n, read(a.observation))}
        print(json.dumps(r))
        raise SystemExit(bool(r['errors']))
    r = validate(n)
    if a.mutations:
        r['fault_controls'] = mutations(n)
    write('VALIDATION.json' if not a.strict else 'STRICT_VALIDATION.json', r)
    print(json.dumps({k: v for (k, v) in r.items() if not isinstance(v, (list, dict))}))
    raise SystemExit(0 if (r['strict_ready'] if a.strict else r['integrity_pass']) and (not a.mutations or r['fault_controls']['all_rejected']) else 1)
if __name__ == '__main__':
    main()
