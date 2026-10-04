"""Queryable AND-obligations. Counts are conditional observation coverage, never physical closure."""
import argparse, copy, math
from collections import Counter
from common import *

def query(selected, contracts):
    selected = set(selected)
    groups = contracts['groups']
    unknown = selected - set(groups)
    if unknown:
        raise ValueError('UNKNOWN_MEASUREMENT:' + ','.join(sorted(unknown)))
    missing = {i: sorted(set(groups[i]['prerequisites']) - selected) for i in selected if set(groups[i]['prerequisites']) - selected}
    touched = []
    ready = []
    remaining = {}
    for (id, c) in contracts['edges'].items():
        if c['status'] != 'BLOCKED_ON_PHYSICAL_MEASUREMENT':
            continue
        req = set(c['measurement_ids'])
        if req & selected:
            touched.append(id)
        if req <= selected and (not missing):
            ready.append(id)
        else:
            remaining[id] = sorted(req - selected)
    return {'selected': sorted(selected), 'missing_protocol_prerequisites': missing, 'touched_edges': sorted(touched), 'all_listed_acquisition_types_present': sorted(ready), 'remaining': remaining, 'actual_edges_closed': 0, 'specimen_and_state_compatibility': 'UNKNOWN_UNTIL_OBSERVATION_RECORDS', 'known_inventory_conflicts': ['M05 full72 and M06 destructive retention24 cannot use unchanged same72'] if {'M05', 'M06'} <= selected else [], 'scope': 'ACQUISITION_TYPE_COVERAGE_ONLY', 'code_model_holdout_prerequisites': 'RETAINED_PER_EDGE', 'counts_resolution': 'PHENOMENOLOGICAL', 'interpretation': 'Potential evaluability only; all requisite measurements on compatible specimen/state, sufficient uncertainty and remaining code/model validations still required.'}

def direct_control(selected, c):
    out = []
    sel = list(selected)
    if any((p not in sel for i in sel for p in c['groups'][i]['prerequisites'])):
        return []
    for (id, e) in c['edges'].items():
        if e['status'] == 'BLOCKED_ON_PHYSICAL_MEASUREMENT' and all((any((i == j for j in sel)) for i in e['measurement_ids'])):
            out.append(id)
    return sorted(out)

def matched_pilot(rows):
    required = ['M01', 'M02', 'M04']
    by = {}
    for r in rows:
        if r['measurement'] not in required or r['unit'] != 'record' or r['state_chain'] != 'scan>seated>fractured' or (r['reference_valid'] is not True):
            continue
        by.setdefault((r['specimen'], r['protocol']), set()).add(r['measurement'])
    return sorted((k for (k, v) in by.items() if set(required) <= v))

def sufficiency():

    def rows(ids):
        return [{'measurement': m, 'specimen': s, 'protocol': 'pilot12', 'unit': 'record', 'state_chain': 'scan>seated>fractured', 'reference_valid': True} for (m, ss) in ids.items() for s in ss]
    a = rows({m: list(range(1, 13)) for m in ['M01', 'M02', 'M04']})
    b = rows({'M01': list(range(1, 13)), 'M02': list(range(1, 13)), 'M04': list(range(13, 25))})
    summary = lambda r: [24] + [sum((x['measurement'] == m for x in r)) for m in ['M01', 'M02', 'M04']]
    sa = summary(a)
    sb = summary(b)
    da = len(matched_pilot(a))
    db = len(matched_pilot(b))
    tests = []
    for (field, value) in [('specimen', 999), ('protocol', 'unmatched'), ('unit', 'mm'), ('state_chain', 'fractured>scanned'), ('reference_valid', False)]:
        x = copy.deepcopy(a)
        next((r for r in x if r['measurement'] == 'M04'))[field] = value
        tests.append({'fault': field, 'injected_value': value, 'matched_count': len(matched_pilot(x)), 'rejected': len(matched_pilot(x)) == da - 1})
    sites = source('results/LANE_X58_CANAL_WALL_SPREAD/raw/PER_SITE_LOCATED.json')
    va = [s['old_minus_new_mid_mm'] for s in sites]
    vb = va.copy()
    imax = max(range(len(va)), key=lambda i: va[i])
    iz = va.index(0.0)
    (vb[imax], vb[iz]) = (vb[iz], vb[imax])
    display = lambda v: [max(v), math.fsum((abs(x) for x in v)) / len(v)]
    ua = display(va)
    ub = display(vb)
    return {'claim_type': 'capability', 'count_witness': {'summary_names': ['manufactured_slots', 'M01_count', 'M02_count', 'M04_count'], 'state_A_summary': sa, 'state_B_summary': sb, 'identity_error': max((abs(a - b) for (a, b) in zip(sa, sb))), 'bitwise_identical': all((float(a).hex() == float(b).hex() for (a, b) in zip(sa, sb))), 'downstream_A_same_specimen_triples': da, 'downstream_B_same_specimen_triples': db, 'downstream_difference': da - db, 'minimal_extension': 'Specimen/protocol ID and observation-state chain per record; count alone cannot determine matched experiments.', 'resolution_level': 'PHENOMENOLOGICAL', 'data': {'state_A': a, 'state_B': b}, 'external_referent': {'kind': 'our_own_fixture', 'locator': 'lab_plan.py:sufficiency', 'compared_quantity': 'exact count invariance versus matched observation triples', 'refutes_us': True}, 'protocol_anchor': logical('results/LANE_X71_MEASUREMENT_PRIORITY/inputs/MEASUREMENTS.json')}, 'located_summary_witness': {'state_A_summary': ua, 'state_B_summary': ub, 'identity_error': max((abs(a - b) for (a, b) in zip(ua, ub))), 'bitwise_identical': all((a.hex() == b.hex() for (a, b) in zip(ua, ub))), 'full_multiset_identical': sorted(va) == sorted(vb), 'queried_site_index': imax, 'permuted_site_index': iz, 'downstream_A_mm': va[imax], 'downstream_B_mm': vb[imax], 'downstream_difference_mm': abs(va[imax] - vb[imax]), 'summary_resolution': 'POPULATION', 'downstream_resolution': 'PER_TOOTH', 'minimal_extension': 'Preserve signed delta bound to case, FDI and pose; active PER_POINT wall witness supports this minimum query.', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'lab_plan.py:sufficiency uses frozen X58 site deltas but permutes their assignment', 'compared_quantity': 'aggregate-invariant site permutation vs located clearance', 'refutes_us': True}}, 'actual_join_faults': tests, 'all_join_faults_rejected': all((t['rejected'] for t in tests))}

def build():
    c = read(HERE / 'MEASUREMENT_CONTRACTS.json')
    group_rows = []
    for (id, g) in c['groups'].items():
        bundle = {id}
        for _ in range(len(c['groups'])):
            nxt = bundle | {p for x in bundle for p in c['groups'][x]['prerequisites']}
            if nxt == bundle:
                break
            bundle = nxt
        q = query(bundle, c)
        single = query([id], c)
        involved = sorted((e for (e, r) in c['edges'].items() if id in r.get('measurement_ids', [])))
        single_req = sorted((e for (e, r) in c['edges'].items() if r.get('measurement_ids') == [id]))
        group_rows.append({'id': id, 'name': g['name'], 'origin': g['origin'], 'x71_matches': g['x71_matches'], 'directly_touched_edges': involved, 'directly_touched_count': len(involved), 'sole_listed_acquisition_edges': single_req, 'sole_listed_acquisition_count': len(single_req), 'bundle_including_protocol_prerequisites': sorted(bundle), 'bundle_acquisition_types_present': q['all_listed_acquisition_types_present'], 'bundle_conditional_count': len(q['all_listed_acquisition_types_present']), 'actual_edges_closed': 0, 'resolution_level': 'PHENOMENOLOGICAL', 'specimen': g['specimen'], 'measurement': g['measurement'], 'operator_hours_interval': g['operator_hours_interval'], 'cost_status': g['cost_status']})
    plans = {k: query(v, c) for (k, v) in {'X71_crown_pilot': ['M01', 'M02', 'M04'], 'X71_all_declared_protocols': [f'M{i:02}' for i in range(1, 16)], 'all_R3_proposed_acquisitions': list(c['groups'])}.items()}
    controls = []
    for (id, g) in c['groups'].items():
        selected = next((r['bundle_including_protocol_prerequisites'] for r in group_rows if r['id'] == id))
        q = query(selected, c)
        controls.append({'bundle': selected, 'candidate': q['all_listed_acquisition_types_present'], 'direct_control': direct_control(selected, c), 'equal': q['all_listed_acquisition_types_present'] == direct_control(selected, c)})
    x71 = source('results/LANE_X71_MEASUREMENT_PRIORITY/results.json')['first_week_crown72_only']
    result = {'claim_type': 'capability', 'groups': group_rows, 'plans': plans, 'X71_original': x71, 'X71_comparison_note': 'X71 five D-questions and R3 edge counts have DIFFERENT denominators. No new costs or empirical yield inferred.', 'controls': controls, 'all_controls_equal': all((r['equal'] for r in controls)), 'physical_closure_count': 0, 'global_lab_optimum': 'UNKNOWN_COST_ACCESS_AND_JOINT_SPECIMEN_COMPATIBILITY', 'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    dump('LAB_MEASUREMENT_COVERAGE.json', result)
    dump('SUFFICIENCY.json', sufficiency())
    return result
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--measurements', nargs='+')
    a = ap.parse_args()
    if a.measurements:
        print(json.dumps(query(a.measurements, read(HERE / 'MEASUREMENT_CONTRACTS.json')), ensure_ascii=False, indent=2))
    else:
        print({k: len(v['all_listed_acquisition_types_present']) for (k, v) in build()['plans'].items()})
