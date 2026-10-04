"""Exact source verification with frozen, manually scoped semantic decisions."""
from common import *
from collections import Counter
PROP = str(DENT.parent / 'bodytwin/results/ASSEMBLY_DENTAL_GAP_PROPOSALS/PROPOSALS_V1.json')

def evaluate(p, rule):
    errors = []
    try:
        row = lock()[p['source_file']]
        if row['sha256'] != p['source_sha256']:
            errors.append('SOURCE_HASH')
        actual = pointer(source(p['source_file']), p['key_path'])
        if type(actual) != type(p['value']) or actual != p['value']:
            errors.append('SOURCE_VALUE')
    except Exception as exc:
        errors.append('SOURCE_READ:' + str(exc))
        actual = None
    if p['unit'] != rule['expected_unit']:
        errors.append('UNIT')
    if p['key_path'] != rule['key_path'] or p['edge_id'] != rule['edge_id']:
        errors.append('ESTIMAND_BINDING')
    eligible = not errors and rule['disposition'] == 'ACCEPT_RELEASE_DELTA' and (p['measured'] is True)
    return {'id': rule['id'], 'proposal': p, 'actual_source_value': actual, 'hash_value_unit_pass': not errors, 'errors': errors, 'semantic_decision': rule, 'accepted_numeric_release_statistic': eligible, 'measured_physical_gap': False}

def run():
    props = source(PROP)
    rules = read(HERE / 'SEMANTIC_RULES.json')['rows']
    rows = [evaluate(p, r) for (p, r) in zip(props['proposals'], rules)]
    assert len(rows) == len(props['proposals']) == 25
    direct = []
    for r in rows:
        p = r['proposal']
        f = HERE / lock()[p['source_file']]['snapshot']
        v = json.loads(f.read_text())
        for part in p['key_path'][1:].split('/'):
            v = v[int(part)] if isinstance(v, list) else v[part]
        direct.append(v == p['value'] and sha(f) == p['source_sha256'])
    refs = {'D-E-K04': ('results/LANE_NIGHT_I_DECISION_CERT/CALIBRATION_ACQUISITION_SPEC.json', '/existing_component_parameters/canal/status'), 'D-E-K12': ('results/LANE_X24_CBCT_HU_CALIBRATION/results.json', ''), 'D-E-K26': ('results/LANE_X68_INSERTION_TORQUE/results.json', '/coverage'), 'D-E-K27': ('results/LANE_X56_MICROMOTION/results.json', '/external_referent'), 'D-E-K38': ('results/LANE_NIGHT_C_TOLERANCE_STACK/results.json', '/attempts/3/physical_validity'), 'D-E-K43': ('results/PROOF_LANE_GENCAD_V5B/results.json', '/physical/measured_seated_film'), 'D-E-K45': ('results/LANE_X64_MANUFACTURING/results.json', '/selection/all_candidates/5/reason'), 'D-E-K48': ('results/LANE_X55_METROLOGY/results.json', '/r1/known_deformations'), 'D-E-SINTER-FIT': ('results/DEMO48_PACKAGE/demos/X14/results.json', '/geometric_reconstruction'), 'D-E-MICROMOTION-OBSERVATION': ('results/LANE_X56_MICROMOTION/results.json', '/external_referent'), 'D-E-X71-PLAN': ('results/LANE_X71_MEASUREMENT_PRIORITY/results.json', '/first_week_crown72_only')}
    near = []
    for p in props['_meta']['unfilled_edges']:
        if not p['near_miss_rejected']:
            continue
        (path, key) = refs[p['edge_id']]
        ev = evidence(path, key, p['unit'], 'PHENOMENOLOGICAL', 'reason_to_reject_not_numeric_evidence')
        if p['edge_id'] == 'D-E-K12':
            ev['value_excerpt'] = 'HU_ref calibration; no same-specimen density/modulus measured. Full record at source pointer.'
        near.append({'edge_id': p['edge_id'], 'decision': 'REJECT_AS_GAP', 'independent_source_check': ev, 'reason': p['reason']})
    result = {'claim_type': 'capability', 'rows': rows, 'counts': {'proposals': len(rows), 'source_exact_pass': sum((r['hash_value_unit_pass'] for r in rows)), 'accepted_numeric_statistics': sum((r['accepted_numeric_release_statistic'] for r in rows)), 'not_accepted_as_numeric_gap': sum((not r['accepted_numeric_release_statistic'] for r in rows)), 'measured_false': sum((not r['proposal']['measured'] for r in rows)), 'measured_false_admitted': sum((r['accepted_numeric_release_statistic'] and (not r['proposal']['measured']) for r in rows)), 'near_misses_rejected': len(near)}, 'independent_pointer_control': {'all_match': all(direct)}, 'near_misses': near, 'dropout_reasons': dict(Counter((r['semantic_decision']['reason'] for r in rows if not r['accepted_numeric_release_statistic']))), 'source_description_corrections': ['X65 three proposed rows are strict_eligible=false; strict_n=3 belongs to other rows.', 'X67 entry 5 is 1.04 degC but not minimum of the 16-value vector (minimum0.3).', 'X58 median patient p95 is not sigma.', 'X31 radial magnitudes are not signed individual deviations.']}
    dump('PROPOSAL_AUDIT.json', result)
    return result
if __name__ == '__main__':
    print(run()['counts'])
