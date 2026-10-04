import copy
from common import *

def build():
    verify()
    n = copy.deepcopy(source('r3_net'))
    obs = read('SOURCE_OBSERVATIONS.json')['records']
    ranges = read('RANGE_CONTRACTS.json')['variables']
    n['schema_version'] = 'dental-constraint-net-r4'
    n['seed'] = 'PROOF_LANE-constraint-net-r4'
    n['claim_type'] = 'capability'
    n['review_state'] = 'PENDING_INDEPENDENT_REVIEW'
    n['r4_parent'] = {'file': entry('r3_net'), 'independent_review': entry('r3_review'), 'unchanged_physical_status': True}
    n['variables']['precementation_internal_gap'] = {'desc': 'Internal dry/PVS replica gap before cementation; source observation port, distinct from assembled cement film and marginal opening', 'unit': 'um', 'resolution_level': 'PER_SURFACE_REGION', 'value': None, 'value_status': 'NO_GLOBAL_VALUE_USE_SCOPED_EDGE_OBSERVATIONS', 'instance_key': ['specimen_or_population', 'region_or_point', 'protocol', 'time', 'source_lineage'], 'knowledge_debt': {'status': 'UNKNOWN', 'replacement_measurement': 'Paired region/state observations with repeatability and controlled seating.'}}
    for (k, v) in n['variables'].items():
        v['plausible_range'] = ranges[k]
    by = {e['id']: e for e in n['edges']}
    for r in obs:
        e = by[r['edge_id']]
        e['r4_reference_binding'] = {'evidence_status': 'VERIFIED_PUBLISHED_OBSERVATION_SCOPED', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'record': r, 'reference_port': r['variable'], 'binding_resolution': r['resolution_level'], 'anatomical_support': r.get('anatomical_support_level', 'PER_TOOTH'), 'timescale': r['timescale'], 'supports_entire_edge': False, 'inherited_direction_is_not_validated_by_reference': True, 'numerical_reference_matches_gap': r['id'] != 'CHOI_FATIGUE', 'target_transfer_status': 'UNKNOWN'}
    n['r4_observations'] = observations(obs)
    n['r4_status_note'] = 'Only new reference-facet availability changes. All R3 edge status, gap values, measurement obligations and directions remain exact; reference does not validate these directions.'
    n['r4_range_note'] = read('RANGE_CONTRACTS.json')['policy']
    write('CONSTRAINT_NET_DENTAL_R4.json', n)
    return n

def observations(records):
    out = []
    (s, c, ch) = records

    def add(id, var, value, case, ref):
        r = next((x for x in read('RANGE_CONTRACTS.json')['variables'][var]['cases'] if x['id'] == case))
        out.append({'id': id, 'variable': var, 'value': float(value), 'unit': r['unit'], 'case_id': case, **{k: r[k] for k in ['quantity', 'protocol_id', 'resolution_level', 'value_role', 'context']}, 'reference': ref})
    for seq in ['first', 'tenth']:
        add('preload_' + seq, 'preload', s['reported'][seq + '_mean_N'], 'SAGHEB_' + seq.upper(), s['reference'])
    add('tightening', 'tightening_torque', 25, 'SAGHEB_TORQUE', s['reference'])
    for method in ['replica', 'ct']:
        add('gap_' + method, 'precementation_internal_gap', c['reported'][method + '_mean_um'], 'CUNALI_AMANN_MO', c['reference'])
    for group in ['TX', 'EV']:
        add('static_' + group, 'fracture_mean', ch['reported']['static_mean_N'][group], 'CHOI_STATIC_' + group, ch['reference'])
    for row in ch['reported']['fatigue_rows']:
        add('fatigue_' + row['group'] + '_' + str(row['loading_level_percent']), 'fatigue_load', row['max_N'], 'CHOI_' + row['group'], ch['reference'])
    return out
if __name__ == '__main__':
    build()
