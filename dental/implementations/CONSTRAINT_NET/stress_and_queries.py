"""Exact consumer-blocker sets; no probabilities or physical margins invented."""
import argparse, itertools, json
from collections import Counter
from common import *
X71D = 'results/LANE_X71_MEASUREMENT_PRIORITY/inputs/DECISIONS.json'
X71M = 'results/LANE_X71_MEASUREMENT_PRIORITY/inputs/MEASUREMENTS.json'
MEASUREMENT_PORTS = {'M01': ['scan_surface', 'pose', 'signed_error', 'asbuilt_error'], 'M02': ['cement_gap', 'marginal_gap'], 'M03': ['cement_gap', 'heat_flux'], 'M04': ['fracture_mean', 'measured_response'], 'M05': ['fracture_mean', 'fracture_quantile', 'weibull_m', 'weibull_sigma0'], 'M06': ['failure_risk', 'material_history'], 'M07': ['patch_force', 'tooth_force', 'tooth_displacement'], 'M08': ['preload', 'assembly_gap', 'contact_traction'], 'M09': ['micromotion', 'tooth_force', 'contact_traction'], 'M10': ['hu_ref', 'bone_density', 'bone_E', 'cortex_thickness'], 'M11': ['anatomy_surface', 'wall_sigma', 'pose', 'nerve_distance', 'dentin_thickness', 'enamel_thickness'], 'M12': ['pdl_E', 'tooth_force', 'tooth_displacement'], 'M13': [], 'M14': ['cement_setting', 'cement_gap', 'marginal_gap'], 'M15': ['cement_gap', 'marginal_gap']}
MEASUREMENT_PORTS['M03'] = ['cement_gap']

def stress(net):
    edges = net['edges']
    out = []
    for (v, o) in net['variables'].items():
        es = [e for e in edges if v in e['between'] and e['status'] != 'TIGHT']
        consumers = sorted({k for e in es for k in e['consumer_chains']})
        out.append({'variable': v, 'desc': o['desc'], 'unit': o['unit'], 'resolution_level': o['resolution_level'], 'unresolved_edges': [e['id'] for e in es], 'consumer_decisions': consumers, 'unresolved_edge_count': len(es), 'decision_count': len(consumers), 'measurement_ids': [m for (m, vs) in MEASUREMENT_PORTS.items() if v in vs], 'count_scope': 'Distinct direct K consumer questions in this curated net, not global physical or clinical importance'})
    out.sort(key=lambda x: (-x['decision_count'], -x['unresolved_edge_count'], x['variable']))
    maxerr = 0
    for row in out:
        seen = []
        for c in net['chains']:
            found = False
            for e in edges:
                if e['status'] in ('OPEN', 'UNKNOWN') and row['variable'] in e['between'] and (c['id'] in e['consumer_chains']):
                    found = True
            if found:
                seen.append(c['id'])
        maxerr = max(maxerr, abs(len(seen) - row['decision_count']))
        if set(seen) != set(row['consumer_decisions']):
            raise AssertionError('Stress consumer mismatch')
    return {'cell': 'DENT-CONSTRAINT-STRESS-MAP', 'generated_by': 'stress_and_queries.py', 'n_variables_total': len(net['variables']), 'n_edges_total': len(edges), 'status_counts': dict(Counter((e['status'] for e in edges))), 'ranking_metric': 'Count distinct direct K consumers of incident OPEN/UNKNOWN edges; ties retain all IDs', 'resolution_level': 'PHENOMENOLOGICAL', 'interpretation': 'Structural priority heuristic; no physical effect sizes or probability weights', 'replacement_measurement': 'Validate question-specific dependencies and acquire matched physical observations before using rank as experiment utility', 'all_variables_sorted': out, 'top_stress_points': out[:12], 'selfcheck_after_update': {'independent_count_max_error': maxerr, 'all_checks_pass': maxerr == 0}, 'sufficiency': 'Scalar rank refuted as measurement selector; retain consumer incidence and full AND blocker sets', 'x71_comparison': compare_x71(out)}

def compare_x71(rows):
    measurements = src(X71M)
    decisions = src(X71D)
    map_rows = []
    for (i, m) in enumerate(measurements):
        vs = MEASUREMENT_PORTS[m['id']]
        map_rows.append({'measurement': m['id'], 'name': m['name'], 'quantity_variables': vs, 'binding_status': 'PROPOSED_OBSERVABLE_MATCH_NOT_CALIBRATION', 'source': evidence(X71M, '/' + str(i), 'PHENOMENOLOGICAL', 'operator_hour;instrument_hour;count', role='prospective_measurement_contract'), 'scope_limit': 'No validated numerical observation map. ' + m['notes']})
    quantitative = [r for r in rows if r['unit'] != '1' and (not r['unit'].startswith('SE3'))][:10]
    overlap = [r for r in quantitative if r['measurement_ids']]
    return {'status': 'PARTIAL_OBSERVABLE_AGREEMENT; FINAL_PRIORITY_ORDER_UNAVAILABLE_AT_FREEZE', 'source_state': src('results/LANE_X71_MEASUREMENT_PRIORITY/CURRENT_WORK_STATE.json'), 'input_lock': src('results/LANE_X71_MEASUREMENT_PRIORITY/INPUT_LOCK.json'), 'mapped_measurements': map_rows, 'top10_quantitative_variables': quantitative, 'overlap_count': len(overlap), 'denominator': len(quantitative), 'contradiction': 'No contradiction of a final X71 priority order can be tested: no such order was available in the frozen input. Structural degree is not an acquisition optimum.', 'coverage_missing': ['Optical ray/translucency/cure variables required by X71 M13 are outside K01–K52 core quantities here.', 'A full rheology parameter field and chemical cement cure history are not present. M03 association only targets the film consumer.', 'X71 feasibility, costs and sample compatibility are prospective assumptions, not measurements.']}

def unlocked(blockers, resolved):
    return sorted((k for (k, rs) in blockers.items() if set(rs) <= set(resolved)))

def sufficiency():
    a = {'D1': ['e1'], 'D2': ['e2']}
    b = {'D1': ['e1', 'e2'], 'D2': ['e1', 'e2']}
    sa = (float(len({e for es in a.values() for e in es})), float(len(a)))
    sb = (float(len({e for es in b.values() for e in es})), float(len(b)))
    ua = unlocked(a, {'e1'})
    ub = unlocked(b, {'e1'})
    return {'claim_type': 'capability', 'summary': '[incident unresolved edge count, affected decision count]', 'state_A': a, 'state_B': b, 'summary_A': sa, 'summary_B': sb, 'bitwise_identical': all((x.hex() == y.hex() for (x, y) in zip(sa, sb))), 'identity_error': max((abs(x - y) for (x, y) in zip(sa, sb))), 'intervention': 'resolve e1', 'unlocked_A': ua, 'unlocked_B': ub, 'downstream_difference': len(ua) - len(ub), 'unit': 'answerable_consumer_count', 'resolution_level': 'PHENOMENOLOGICAL', 'outcome': 'SCALAR_SUMMARY_INSUFFICIENT', 'minimal_sufficient_extension': 'For this deterministic model: consumer identity plus remaining AND blocker set. For physical decisions also specimen/protocol/uncertainty/criterion.', 'proof': 'Answerable(d,S) iff blockers(d) subset S. Incidence plus AND semantics determines every intervention; counts do not.', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'stress_and_queries.py:sufficiency', 'compared_quantity': 'Set logic counterexample, not an external physical measurement', 'refutes_us': True}}

def x71_portfolio(ids):
    ms = {m['id']: m for m in src(X71M)}
    ds = src(X71D)
    ids = set(ids)
    if not ids <= set(ms):
        raise ValueError('Unknown measurement ID')
    missing = {p for i in ids for p in ms[i]['prerequisites'] if p not in ids}
    supplied = {s for i in ids for s in ms[i]['supplies']}
    unlocked_ids = [d['id'] for d in ds if set(d['requires']) <= supplied and (not d['other_blockers'])]
    other = {d['id']: {'missing_measurements': sorted(set(d['requires']) - supplied), 'other_blockers': d['other_blockers']} for d in ds if d['id'] not in unlocked_ids}
    hours = [sum((ms[i]['operator_hours_interval'][j] for i in ids)) for j in (0, 1)]
    return {'measurements': sorted(ids), 'prerequisites_satisfied': not missing, 'missing_prerequisites': sorted(missing), 'supplied_observables': sorted(supplied), 'conditionally_answerable': unlocked_ids if not missing else [], 'count': len(unlocked_ids) if not missing else 0, 'remaining': other, 'operator_hours_interval': hours, 'cost_resolution': 'PHENOMENOLOGICAL', 'cost_replacement': 'Lab quote and timed dry run; sums assume sequential nonshared operator work', 'physical_validation': False, 'scope': 'Prospective dependency feasibility only; no actual measurements taken and no protocol compatibility optimization'}

def build_bindings(net):
    rows = []
    for (kind, objects) in [('variable', net['variables'].items()), ('edge', ((e['id'], e) for e in net['edges']))]:
        for (id_, obj) in objects:
            nodes = [b['id'] for b in obj['working_bindings']]
            rows.append({'id': 'DENT-CONSTRAINT-' + kind.upper() + '-' + id_.upper().replace('_', '-'), 'type': 'CLAIM', 'status': 'OPEN', 'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'claim': obj.get('desc', obj.get('constraint')), 'constraint_kind': kind, 'constraint_id': id_, 'annotation_for': nodes[0], 'working_ids': nodes, 'constraint_status': obj.get('status', 'VARIABLE_DEFINITION'), 'source_contract': obj, 'depends_on': [], 'proposed_semantic_edges': [{'to': n, 'relation': 'quantity_constraint_context', 'status': 'PENDING_INDEPENDENT_REVIEW', 'numeric_fusion': False} for n in nodes], 'evidence': ['results/PROOF_LANE_CONSTRAINT_NET_DENTAL/CONSTRAINT_NET_DENTAL.json'], 'source_status_preserved': True, 'no_scientific_admission': True})
    for c in net['chains']:
        rows.append({'id': 'DENT-CONSTRAINT-CHAIN-' + c['id'], 'type': 'CLAIM', 'status': 'OPEN', 'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'claim': c['desc'], 'annotation_for': c['working_ids'][0], 'working_ids': c['working_ids'], 'northstar_id': c['id'], 'variables': c['variables'], 'constraint_edges': c['edge_ids'], 'source_contract': c, 'depends_on': [], 'no_scientific_admission': True})
    payload = ''.join((json.dumps(r, ensure_ascii=False, separators=(',', ':')) + '\n' for r in rows))
    (HERE / 'constraint_net_dental_20261003.jsonl').write_text(payload)
    dest = DENT / 'notes/expansion/constraint_net_dental_20261003.jsonl'
    if dest.exists() and dest.read_text() != payload:
        raise RuntimeError('Expansion already exists with different content; do not overwrite another session')
    if not dest.exists():
        dest.write_text(payload)
    return {'rows': len(rows), 'variables': len(net['variables']), 'edges': len(net['edges']), 'chains': len(net['chains']), 'path': str(dest), 'sha256': sha(dest), 'status': 'PROPOSAL_NOT_NATIVE_ADMISSION'}

def main():
    a = argparse.ArgumentParser()
    a.add_argument('--chain')
    a.add_argument('--measurements', nargs='*')
    a.add_argument('--build', action='store_true')
    args = a.parse_args()
    net = read(HERE / 'CONSTRAINT_NET_DENTAL.json')
    if args.chain:
        c = next((c for c in net['chains'] if c['id'] == args.chain))
        print(json.dumps({'chain': c, 'relations': [e for e in net['edges'] if e['id'] in c['edge_ids']]}, ensure_ascii=False, indent=2))
        return
    if args.measurements is not None:
        print(json.dumps(x71_portfolio(args.measurements), ensure_ascii=False, indent=2))
        return
    out = stress(net)
    dump('CONSTRAINT_STRESS_MAP_DENTAL.json', out)
    dump('SUFFICIENCY_R1.json', sufficiency())
    dump('BINDINGS_RECEIPT.json', build_bindings(net))
    print(json.dumps({'top_variables': [(x['variable'], x['decision_count']) for x in out['top_stress_points']], 'x71_quantitative_overlap': [out['x71_comparison']['overlap_count'], out['x71_comparison']['denominator']], 'summary_identity_error': sufficiency()['identity_error'], 'downstream_difference': sufficiency()['downstream_difference']}))
if __name__ == '__main__':
    main()
