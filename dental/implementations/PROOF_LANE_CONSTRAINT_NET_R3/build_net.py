import copy
from fractions import Fraction as Q
from common import *
from audit_proposals import run as audit
CORRECTED = ['D-E-CONTRAST-' + str(i).zfill(2) for i in range(1, 9)] + ['D-E-CEMENT-TRANSFER-FAIL', 'D-E-BRIDGE-SYSTEM']

def residual(ev, quantity, scope):
    v = Q(str(ev['value']))
    return {'value': float(abs(v)), 'signed_value': float(v), 'exact_decimal': str(abs(v)), 'unit': ev['unit'], 'quantity': quantity, 'resolution_level': ev['resolution_level'], 'kind': 'OBSERVED_RELEASE_DIFFERENCE', 'operation': 'reported_residual', 'operands': [ev], 'scope': scope, 'measured_physical_gap': False, 'not_complete': True}

def build():
    verify_frozen()
    a = audit()
    n = copy.deepcopy(source('results/PROOF_LANE_CONSTRAINT_NET_R2/CONSTRAINT_NET_DENTAL_R2.json'))
    contracts = read(HERE / 'MEASUREMENT_CONTRACTS.json')
    n.update({'schema_version': 'dental-constraint-net-r3', 'seed': 'PROOF_LANE-constraint-net-r3', 'r3_parent_sha256': lock()[logical('results/PROOF_LANE_CONSTRAINT_NET_R2/CONSTRAINT_NET_DENTAL_R2.json')]['sha256'], 'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'gap_policy': 'Finite correctly sourced residual OR explicit physical measurement obligation. Nonphysical exceptions remain strict failures. Source edge status preserved; gap.kind/blocker_status are the new availability status.'})
    by = {e['id']: e for e in n['edges']}
    for (id, c) in contracts['edges'].items():
        e = by[id]
        e['r2_gap'] = copy.deepcopy(e['gap'])
        e['gap']['kind'] = c['status']
        e['gap']['status'] = c['status']
        e['blocker_status'] = c['status']
        e['measurement_obligation'] = c
        e['gap']['reason'] = 'R3 independent classification; missing source quantity remains null.'
        if c['status'] == 'BLOCKED_ON_PHYSICAL_MEASUREMENT':
            e['gap']['required_measurement'] = c['measurement']
            e['gap']['specimen'] = c['specimen']
            e['gap']['would_decide'] = c['would_decide']
            e['gap']['measurement_ids'] = c['measurement_ids']
        else:
            e['gap']['required_operation'] = c['required_operation']
    for row in a['rows']:
        p = row['proposal']
        rule = row['semantic_decision']
        e = by[p['edge_id']]
        e.setdefault('r3_verified_context', []).append({'proposal_id': row['id'], 'evidence': evidence(p['source_file'], p['key_path'], p['unit'], rule['resolution_level'], 'context_not_full_gap'), 'quantity_status': rule['kind'], 'reason': rule['reason'], 'supports_full_edge_physical_closure': False})
    e = by['D-E-CANAL-RELEASE']
    e['r2_gap'] = copy.deepcopy(e['gap'])
    sitepath = 'results/LANE_X58_CANAL_WALL_SPREAD/raw/PER_SITE_LOCATED.json'
    sites = source(sitepath)
    e['gap'] = residual(evidence('results/LANE_X58_CANAL_WALL_SPREAD/results.json', '/R3/max_located_old_minus_new_mm', 'mm', 'POPULATION'), 'Maximum old minus new release clearance on located sites', 'Display aggregate only; no anatomical-truth or individual safety claim.')
    e['blocker_status'] = 'NUMERIC_RELEASE_DIFFERENCE_WITH_PHYSICAL_TRUTH_UNKNOWN'
    e['gap_components'].append(residual(evidence('results/LANE_X58_CANAL_WALL_SPREAD/results.json', '/R3/mean_located_abs_change_mm', 'mm', 'POPULATION'), 'Mean absolute release clearance change', 'Population aggregate, not per-point error.'))
    e['located_release_field'] = []
    for (i, s) in enumerate(sites):
        ev = evidence(sitepath, f'/{i}/old_minus_new_mid_mm', 'mm', 'PER_TOOTH', 'signed_located_release_difference')
        e['located_release_field'].append({'site_index': i, 'evidence': ev, 'instance': {'case': s['case'], 'fdi': s['fdi'], 'pose': s['pose']}, 'spatial_witness': {'old': s['old'].get('witness'), 'new': s['new'].get('witness')}, 'point_artifact': s['point_artifact'], 'binding_resolution': 'PER_POINT', 'observable_resolution': 'PER_TOOTH', 'timescale': 'SIMULTANEOUS', 'truth_error': 'UNKNOWN'})
    e['physical_truth_obligation'] = {'status': 'BLOCKED_ON_PHYSICAL_MEASUREMENT', 'measurement_ids': ['M11', 'A19'], 'measurement': contracts['groups']['A19']['measurement'], 'specimen': contracts['groups']['A19']['specimen'], 'would_decide': 'Independent anatomical boundary error; release differences alone do not decide this.'}
    n['r3_proposal_audit'] = 'PROPOSAL_AUDIT.json'
    n['measurement_groups'] = contracts['groups']
    n['r3_no_new_physical_observations'] = True
    dump('CONSTRAINT_NET_DENTAL_R3.json', n)
    dump('CURRENT_WORK_STATE.json', {'milestone': 'R3A_COMPLETE_R3B_BUILT', 'latest_gate': '25 exact source checks, 2 release statistics accepted; 8 declared values excluded', 'current': 'validate numeric and physical-obligation branches, recount lab groups', 'next_operation': 'actual mutations, specimen matching sufficiency, stress map', 'review_state': 'PENDING_INDEPENDENT_REVIEW'})
    return n
if __name__ == '__main__':
    print('Built', len(build()['edges']), 'edges')
