"""Validate syntax, exact provenance, review-scoped TIGHT and graph binding.

No empirical validity is inferred from these integrity checks.
"""
import argparse, copy, json, math
from collections import Counter
from common import *
LEVELS = ['PER_POINT', 'PER_SURFACE_REGION', 'PER_TOOTH', 'PER_ARCH', 'POPULATION', 'PHENOMENOLOGICAL']
SCHEMA = {'$schema': 'https://json-schema.org/draft/2020-12/schema', 'type': 'object', 'required': ['seed', 'objective', 'schema_note', 'source_files', 'variables', 'edges', 'chains', 'claim_type'], 'properties': {'claim_type': {'const': 'capability'}, 'variables': {'type': 'object', 'minProperties': 1, 'additionalProperties': {'type': 'object', 'required': ['desc', 'unit', 'resolution_level', 'working_bindings', 'chains', 'knowledge_debt'], 'properties': {'unit': {'type': 'string', 'minLength': 1}, 'resolution_level': {'enum': LEVELS}, 'working_bindings': {'type': 'array', 'minItems': 1}}}}, 'edges': {'type': 'array', 'minItems': 1, 'items': {'type': 'object', 'required': ['id', 'between', 'constraint', 'status', 'evidence', 'consumer_chains', 'resolution_level', 'timescale', 'working_bindings', 'knowledge_debt', 'rigorous_enclosure', 'scope', 'external_referent'], 'properties': {'id': {'type': 'string'}, 'between': {'type': 'array', 'minItems': 2, 'uniqueItems': True, 'items': {'type': 'string'}}, 'status': {'enum': ['TIGHT', 'OPEN', 'UNKNOWN']}, 'resolution_level': {'enum': LEVELS}, 'timescale': {'enum': ['SIMULTANEOUS', 'HANDOVER']}, 'evidence': {'type': 'array', 'minItems': 1, 'items': {'type': 'object', 'required': ['source_file', 'key', 'value', 'source_sha256', 'resolution_level', 'unit'], 'properties': {'key': {'type': 'string', 'pattern': '^/'}, 'resolution_level': {'enum': LEVELS}}}}, 'working_bindings': {'type': 'array', 'minItems': 1}}}}, 'chains': {'type': 'array', 'minItems': 52, 'maxItems': 52, 'items': {'type': 'object', 'required': ['id', 'variables', 'edge_ids', 'working_ids']}}}}

def exact(a, b):
    return json.dumps(a, sort_keys=True, ensure_ascii=False, allow_nan=False) == json.dumps(b, sort_keys=True, ensure_ascii=False, allow_nan=False)

def validate(net, live=False, curation=True, live_values=False):
    import jsonschema
    errors = []
    checked = 0
    cache = {}
    evidence_seen = set()
    hash_drift = {}

    def fail(code, where, detail):
        errors.append({'code': code, 'where': where, 'detail': str(detail)})
    for er in jsonschema.Draft202012Validator(SCHEMA).iter_errors(net):
        fail('SCHEMA', '/'.join(map(str, er.path)), er.message)
    if errors:
        return {'pass': False, 'errors': errors, 'evidence_rows_checked': 0}

    def ev(row, where):
        nonlocal checked
        checked += 1
        try:
            origin = row['source_file']
            record = SOURCES[origin]
            p = Path(origin) if live or live_values else HERE / record['snapshot']
            if origin not in cache:
                cache[origin] = (sha(p), read(p))
            (h, d) = cache[origin]
            if h != row['source_sha256'] or h != record['sha256']:
                hash_drift[origin] = {'frozen_sha256': record['sha256'], 'current_sha256': h}
                if not live_values:
                    fail('EVIDENCE_HASH', where, origin)
            actual = pointer(d, row['key'])
            if not exact(actual, row['value']):
                fail('EVIDENCE_VALUE', where, row['key'])
            evidence_seen.add((origin, row['key']))
        except (OSError, KeyError, IndexError, ValueError, TypeError) as e:
            fail('EVIDENCE_LOOKUP', where, e)
    var = net['variables']
    edge = {e['id']: e for e in net['edges']}
    chain = {c['id']: c for c in net['chains']}
    if len(edge) != len(net['edges']):
        fail('DUPLICATE_EDGE', 'edges', 'IDs')
    if set(chain) != {f'K{i:02}' for i in range(1, 53)}:
        fail('CHAIN_COVERAGE', 'chains', sorted(chain))
    view = src('results/GRAPH_WORKING_VIEW_20260923/WORKING_VIEW.json')
    known = {**view['nodes'], **view['planning_nodes']}

    def binding(b, where):
        if b['id'] not in known:
            fail('WORKING_NODE', where, b['id'])
        ev(b['source_claim'], where)
        if b['source_claim']['key'].split('/')[-2] != b['id']:
            fail('WORKING_NODE_ID', where, b['id'])
    for (v, o) in var.items():
        for b in o['working_bindings']:
            binding(b, 'variable:' + v)
        if not o['chains']:
            fail('UNBOUND_VARIABLE', v, 'no chain')
        for c in o['chains']:
            if c not in chain or v not in chain[c]['variables']:
                fail('VARIABLE_CHAIN', v, c)
        if not o['knowledge_debt'].get('replacement_measurement'):
            fail('DEBT', v, 'missing replacement')
    for e in net['edges']:
        id_ = e['id']
        for v in e['between']:
            if v not in var:
                fail('DANGLING_QUANTITY', id_, v)
        for b in e['working_bindings']:
            binding(b, id_)
        for row in e['evidence']:
            ev(row, id_)
        for c in e['consumer_chains']:
            if c not in chain or id_ not in chain[c]['edge_ids']:
                fail('EDGE_CONSUMER', id_, c)
        if e['status'] == 'TIGHT':
            try:
                r = e['review']
                d = src(r['record'])
                if not str(pointer(d, r['decision_key'])).startswith('ACCEPT'):
                    fail('TIGHT_REVIEW', id_, d.get('decision'))
                if SOURCES[r['record']]['sha256'] != r['record_sha256']:
                    fail('REVIEW_HASH', id_, r['record'])
                actual = SOURCES[r['reviewed_result']]['sha256']
                if live or live_values:
                    if sha(r['reviewed_result']) != actual:
                        fail('TIGHT_LIVE_RESULT_HASH', id_, r['reviewed_result'])
                if actual != d['result_sha256'] or actual != r['reviewed_result_sha256']:
                    fail('TIGHT_RESULT_HASH', id_, actual)
                er = e['external_referent']
                if er['kind'] in ['our_own_fixture', 'external_review'] or not er['locator'] or (not er['compared_quantity']):
                    fail('TIGHT_EXTERNAL_FACIT', id_, er)
                ec = e['external_comparison']
                if ec['executed'] is not True:
                    fail('TIGHT_NOT_EXECUTED', id_, ec)
                if ec['result_file'] != r['reviewed_result'] or ec['source_result_sha256'] != actual:
                    fail('TIGHT_COMPARISON_IDENTITY', id_, ec)
                for key in ec['quantity_paths']:
                    if not any((x['source_file'] == ec['result_file'] and x['key'] == key for x in e['evidence'])):
                        fail('TIGHT_QUANTITY_NOT_EVIDENCED', id_, key)
            except (KeyError, TypeError) as err:
                fail('TIGHT_PROVENANCE_MISSING', id_, err)
    for c in chain.values():
        for eid in c['edge_ids']:
            if eid not in edge or c['id'] not in edge[eid]['consumer_chains']:
                fail('CHAIN_EDGE', c['id'], eid)
        for v in c['variables']:
            if v not in var or c['id'] not in var[v]['chains']:
                fail('CHAIN_VARIABLE', c['id'], v)
        for nid in c['working_ids']:
            if nid not in known:
                fail('CHAIN_WORKING_NODE', c['id'], nid)
    if curation:
        contract = read(HERE / 'CURATION_CONTRACT.json')
        for expected in contract['edges']:
            e = edge.get(expected['id'], {})
            for (key, v) in expected.items():
                if not exact(e.get(key), v):
                    fail('CURATED_SCOPE', expected['id'], key)
    return {'pass': not errors, 'errors': errors, 'evidence_rows_checked': checked, 'unique_evidence_rows': len(evidence_seen), 'source_files_checked': len(cache), 'variables': len(var), 'edges': len(edge), 'chains': len(chain), 'status_counts': dict(Counter((e['status'] for e in net['edges']))), 'check_mode': 'LIVE_VALUES_WITH_REVIEW_HASHES' if live_values else 'LIVE_SOURCES' if live else 'FROZEN_SOURCES', 'whole_file_hash_drift': hash_drift, 'scientific_validity_inferred': False}

def mutations(net):
    probes = []

    def run(name, change, expected):
        n = copy.deepcopy(net)
        change(n)
        r = validate(n, curation=False)
        probes.append({'name': name, 'rejected': not r['pass'], 'expected_detector': expected, 'detectors': sorted({e['code'] for e in r['errors']})})
    tight = next((i for (i, e) in enumerate(net['edges']) if e['status'] == 'TIGHT'))
    run('false_numeric_value', lambda n: n['edges'][tight]['evidence'][0].update(value=999999.0), 'EVIDENCE_VALUE')
    run('unknown_promoted_to_TIGHT', lambda n: n['edges'][0].update(status='TIGHT'), 'TIGHT_EXTERNAL_FACIT')
    run('missing_review', lambda n: n['edges'][tight].update(review=None), 'TIGHT_PROVENANCE_MISSING')
    run('wrong_reviewed_hash', lambda n: n['edges'][tight]['review'].update(reviewed_result_sha256='0' * 64), 'TIGHT_RESULT_HASH')
    run('fixture_disguised_as_facit', lambda n: n['edges'][tight]['external_referent'].update(kind='our_own_fixture'), 'TIGHT_EXTERNAL_FACIT')
    run('dangling_quantity', lambda n: n['edges'][0]['between'].append('DOES_NOT_EXIST'), 'DANGLING_QUANTITY')
    run('missing_evidence_locator', lambda n: n['edges'][tight]['evidence'][0].update(key='/absent'), 'EVIDENCE_LOOKUP')
    run('missing_unit', lambda n: n['variables']['cement_gap'].pop('unit'), 'SCHEMA')
    run('wrong_working_node', lambda n: n['variables']['cement_gap']['working_bindings'][0].update(id='DENT-NONEXISTENT'), 'WORKING_NODE')
    run('missing_chain', lambda n: n['chains'].pop(), 'SCHEMA')
    run('wrong_time_scale', lambda n: n['edges'][tight].update(timescale='INSTANT_HEALING'), 'SCHEMA')
    for (name, change) in [('wrong_quantity', lambda n: n['edges'][tight]['external_referent'].update(compared_quantity='patient bone modulus')), ('wrong_assay_handover', lambda n: next((e for e in n['edges'] if e['id'] == 'D-E-DRILL-CULTURE')).update(timescale='SIMULTANEOUS'))]:
        n = copy.deepcopy(net)
        change(n)
        r = validate(n)
        probes.append({'name': name, 'rejected': not r['pass'], 'expected_detector': 'CURATED_SCOPE', 'detectors': sorted({e['code'] for e in r['errors']})})
    return {'pass': all((p['rejected'] and p['expected_detector'] in p['detectors'] for p in probes)), 'probes': probes, 'fraction_rejected': sum((p['rejected'] for p in probes)) / len(probes)}

def main():
    a = argparse.ArgumentParser()
    a.add_argument('--net', default='CONSTRAINT_NET_DENTAL.json')
    a.add_argument('--live', action='store_true')
    a.add_argument('--live-values', action='store_true')
    a.add_argument('--mutations', action='store_true')
    args = a.parse_args()
    dump('CONSTRAINT_NET_DENTAL.schema.json', SCHEMA)
    net = read(HERE / args.net)
    r = validate(net, args.live, live_values=args.live_values)
    if args.mutations:
        r['mutation_tests'] = mutations(net)
    dump('VALIDATION_LIVE_VALUES.json' if args.live_values else 'VALIDATION_LIVE.json' if args.live else 'VALIDATION.json', r)
    print(json.dumps({k: v for (k, v) in r.items() if k != 'mutation_tests'}))
    raise SystemExit(0 if r['pass'] and r.get('mutation_tests', {'pass': True})['pass'] else 1)
if __name__ == '__main__':
    main()
