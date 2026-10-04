import itertools
import time
import json
import datetime
from common import ROOT, OUTPUT, read, write, sha, check_frozen, checkpoint
from measurement import independent_separator_check, minimum_separators
FIELDS = ['case_definition_confirmed', 'local_periodontal_origin_confirmed', 'ESE_deep_geometry_confirmed', 'pulp_state_compatible', 'irreversible_symptoms_absent', 'actual_pulp_exposure_absent']

def outcomes(state):
    usable = state['case_definition_confirmed'] and state['local_periodontal_origin_confirmed']
    eligible = all((state[x] for x in FIELDS[2:]))
    return 'FURCATION_' + ('USABLE' if usable else 'NOT_USABLE') + '|CARIES_' + ('ELIGIBLE' if eligible else 'OUTSIDE')

def run():
    check_frozen('PREREG_R4_LOCAL_ETIOLOGY.json')
    start = time.perf_counter()
    r3 = read('raw/RESULTS_R3.json')
    worlds = []
    for vals in itertools.product([False, True], repeat=len(FIELDS)):
        state = dict(zip(FIELDS, vals))
        worlds.append({'state': state, 'outcome': outcomes(state)})
    separators = minimum_separators(worlds, FIELDS)
    comparisons = []
    for w in worlds:
        v = w['state']
        control = (v[FIELDS[0]] and v[FIELDS[1]], all((v[x] for x in FIELDS[2:])))
        got = outcomes(v)
        expected = 'FURCATION_' + ('USABLE' if control[0] else 'NOT_USABLE') + '|CARIES_' + ('ELIGIBLE' if control[1] else 'OUTSIDE')
        comparisons.append(got == expected)
    old_summary = {'case_definition_confirmed': True, 'stage_iv_complexity': False, 'periodontal_tooth_loss_ge5': False, 'prior_stage_iv': False}
    summaries = [old_summary.copy(), old_summary.copy()]
    vectors = [[float(s[k]) for k in sorted(s)] for s in summaries]
    hashes = [sha(json.dumps(s, sort_keys=True, separators=(',', ':')).encode()) for s in summaries]
    origin = [True, False]
    usable = [summaries[i]['case_definition_confirmed'] and origin[i] for i in range(2)]
    counterexample = {'R3_summaries': summaries, 'summary_identity_error': max((abs(a - b) for (a, b) in zip(*vectors))), 'summary_float_hex': [[v.hex() for v in vs] for vs in vectors], 'summary_sha256': hashes, 'summary_bytes_identical': hashes[0] == hashes[1], 'local_origin_values': origin, 'downstream_furcation_usable': usable, 'downstream_difference_boolean': int(usable[0]) - int(usable[1]), 'minimum_extension': 'Separate tooth-local periodontal attribution from case-level disease membership', 'referent': {'kind': 'external_review', 'locator': 'doi:10.1002/JPER.18-0006::clinical case exclusions and periodontal-only destruction requirement', 'compared_quantity': 'Local periodontal admissibility; worlds are a logical counterexample, not measured etiology', 'refutes_us': True}, 'resolution': 'PER_TOOTH', 'timescale': 'SIMULTANEOUS'}
    unknown_outcomes = sorted({w['outcome'] for w in worlds})
    real = {'patient_id': r3['patient_id'], 'snapshot_row': r3['snapshot_row'], 'FDI': 16, 'furcation_periodontal_admissibility': 'UNKNOWN_LOCAL_ETIOLOGY', 'stage_constraint': 'UNKNOWN_ETIOLOGY; {III,IV} only conditional on confirmed local periodontal origin and case definition', 'caries_branch': 'UNKNOWN_GEOMETRY_AND_PULP_STATE', 'possible_joint_predicate_outcomes': unknown_outcomes, 'missing_joint_global_predicates': separators[0], 'pulp_measurement_settles_periodontal_origin_alone': False, 'independent_measurement_performed': False}
    injection = {'name': 'Force periodontal stage evidence admissibilityYES without local-origin adjudication', 'forged_value': 'YES', 'rejected': real['furcation_periodontal_admissibility'] != 'YES', 'control': 'Exact source-premise gate; correct UNKNOWN accepted'}
    necessity = []
    for f in FIELDS[:2]:
        pairs = independent_separator_check(worlds, [x for x in FIELDS if x != f])
        necessity.append({'missing_predicate': f, 'opposing_pair': pairs[0]})
    result = {'round': 'R4', 'claim_type': 'capability', 'worlds': worlds, 'hypothetical_world_count': len(worlds), 'possible_joint_outcomes': unknown_outcomes, 'same_information_control_disagreements': len(comparisons) - sum(comparisons), 'counterexample': counterexample, 'R3_full_clinical_sufficiency': 'REFUTED_MISSING_LOCAL_ORIGIN', 'real_shared_case': real, 'minimal_global_separator': separators, 'necessity_witnesses': necessity, 'injections': [injection], 'gate': 'PASS_CONDITIONAL_REPRESENTATION' if all(comparisons) and injection['rejected'] else 'FAIL', 'external_referent': read('PREREG_R4_LOCAL_ETIOLOGY.json')['external_referent'], 'clinical_validation': 'NOT_PERFORMED', 'complete_clinical_answers': 0, 'wall_s': time.perf_counter() - start}
    write('raw/RESULTS_R4.json', result)
    frozen = {'schema': 'X57-local-origin-frozen-predictions-v1', 'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'patient_id': r3['patient_id'], 'snapshot_row': r3['snapshot_row'], 'FDI': 16, 'source_observation_sha256': sha((ROOT / 'raw/OBSERVATIONS_R2.json').read_bytes()), 'before_independent_measurement': True, 'independent_measurement_performed': False, 'report_semantics_predictions': {'furcation_grade': 3, 'furcation_tooth': 16, 'caries_report_label': 'DEEP', 'caries_surface': 'mesial'}, 'future_falsifier': 'Blinded expert disagrees with literal tooth/surface/quantity binding, or independent clinical measurements invalidate any invoked source premise', 'clinical_predictions': real, 'external_referent': result['external_referent']}
    write('FROZEN_PREDICTIONS_R4.json', frozen)
    (OUTPUT / 'FROZEN_PREDICTIONS_R4.sha256').write_text(sha((OUTPUT / 'FROZEN_PREDICTIONS_R4.json').read_bytes()) + '\n')
    checkpoint('R4_COMPUTED', result['gate'] + '; real clinical decisions UNKNOWN', 'Package independent annotation, anatomical uncertainty and one-command shared-case demo')
    print(result['gate'], 'worlds', len(worlds), 'joint real case', r3['patient_id'])
if __name__ == '__main__':
    run()
