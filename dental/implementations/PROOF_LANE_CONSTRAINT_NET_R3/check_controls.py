import copy
from common import *
from audit_proposals import evaluate
from lab_plan import query, direct_control
from validate import validate

def run():
    audit = read(HERE / 'PROPOSAL_AUDIT.json')
    rows = []
    for r in audit['rows']:
        p = copy.deepcopy(r['proposal'])
        p['value'] += 1
        out = evaluate(p, r['semantic_decision'])
        rows.append({'id': r['id'], 'fault': 'source value +1', 'rejected': not out['hash_value_unit_pass'], 'errors': out['errors']})
    for (field, value) in [('source_sha256', '0' * 64), ('unit', 'kg'), ('key_path', '/missing')]:
        r = audit['rows'][9]
        p = copy.deepcopy(r['proposal'])
        p[field] = value
        out = evaluate(p, r['semantic_decision'])
        rows.append({'id': field, 'fault': 'corrupt proposal ' + field, 'rejected': not out['hash_value_unit_pass'], 'errors': out['errors']})
    for r in audit['rows']:
        if not r['proposal']['measured']:
            p = copy.deepcopy(r['proposal'])
            p['measured'] = True
            out = evaluate(p, r['semantic_decision'])
            rows.append({'id': r['id'], 'fault': 'mislabel derived value measured:true', 'rejected': not out['accepted_numeric_release_statistic']})
    c = read(HERE / 'MEASUREMENT_CONTRACTS.json')
    plan_checks = []
    q = query(['M02'], c)
    plan_checks.append({'fault': 'omit M01 prerequisite', 'rejected': q['missing_protocol_prerequisites'] == {'M02': ['M01']} and (not q['all_listed_acquisition_types_present'])})
    q = query(['M01', 'M02'], c)
    plan_checks.append({'fault': 'claim full pilot without fracture', 'rejected': 'D-E-X71-PLAN' not in q['all_listed_acquisition_types_present']})
    q = query(['M01', 'M02', 'M04'], c)
    claimed = q['all_listed_acquisition_types_present'] + ['D-E-K43']
    plan_checks.append({'fault': 'claim sinter/process link from scan+film only', 'rejected': sorted(claimed) != direct_control(['M01', 'M02', 'M04'], c)})
    try:
        query(['MADE_UP'], c)
        rejected = False
    except ValueError:
        rejected = True
    plan_checks.append({'fault': 'unknown measurement', 'rejected': rejected})
    s = read(HERE / 'SUFFICIENCY.json')
    net = read(HERE / 'VALIDATION.json')['fault_controls']
    r = {'proposal_faults': rows, 'planning_faults': plan_checks, 'network_faults': net, 'join_faults': s['actual_join_faults'], 'all_rejected': all((x['rejected'] for x in rows + plan_checks + s['actual_join_faults'])) and net['all_rejected'], 'count': len(rows) + len(plan_checks) + len(s['actual_join_faults']) + net['count']}
    dump('CONTROLS.json', r)
    return r
if __name__ == '__main__':
    print({k: v for (k, v) in run().items() if k in ['all_rejected', 'count']})
