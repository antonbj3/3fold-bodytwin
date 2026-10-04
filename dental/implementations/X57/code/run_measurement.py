import datetime
import resource
import time
from common import ROOT, OUTPUT, read, write, sha, check_frozen, checkpoint
from measurement import make_worlds, minimum_separators, minimum_answer_certificates, independent_separator_check, answer

def run():
    check_frozen('PREREG_R3_MINIMUM_MEASUREMENT.json')
    start = time.perf_counter()
    cpu = time.process_time()
    r2 = read('raw/RESULTS_R2.json')
    row = r2['joint_localized_snapshot_rows'][0]
    case = next((x for x in read('raw/PATIENT_MODELS_R2.json') if x['snapshot_row'] == row))
    prediction = {'schema': 'X57-case-measurement-prediction-v1', 'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'patient_id': case['patient_id'], 'snapshot_row': row, 'FDI': 16, 'source': 'raw/PATIENT_MODELS_R2.json', 'source_sha256': sha((ROOT / 'raw/PATIENT_MODELS_R2.json').read_bytes()), 'external_referents': read('PREREG_R3_MINIMUM_MEASUREMENT.json')['external_referent'], 'predictions_before_independent_measurement': {'BD10': 'Conditional stage III/IV if periodontitis case and journal furcation III confirmed; other etiology remains possible', 'BD08': 'Generic deep-caries label does not establish ESE barrier/inner-quarter criteria or compatible pulpal state; branch eligibility UNKNOWN'}, 'independent_measurement_performed': False, 'frozen_expected_model_certificate_predicate_counts': {'BD10_global': 4, 'BD10_answerIII': 4, 'BD10_answerIV_minimum': 2, 'BD10_caseexcluded': 1, 'BD08_global': 4, 'BD08_eligible': 4, 'BD08_outside_minimum': 1}, 'code_hashes': {p.name: sha(p.read_bytes()) for p in (ROOT / 'code').glob('*.py')}}
    write('FROZEN_PREDICTIONS.json', prediction)
    (OUTPUT / 'FROZEN_PREDICTIONS.sha256').write_text(sha((OUTPUT / 'FROZEN_PREDICTIONS.json').read_bytes()) + '\n')
    chains = {}
    failure_count = 0
    for chain in ['BD10', 'BD08']:
        (fields, worlds) = make_worlds(chain)
        separators = minimum_separators(worlds, fields)
        certs = []
        for w in worlds:
            solutions = minimum_answer_certificates(worlds, fields, w['state'])
            certs.append({'hypothetical_state': w['state'], 'outcome': w['outcome'], 'minimum_predicate_count': len(solutions[0]), 'certificates': solutions, 'source': 'finite logical completions, not observed patient data'})
        necessity = []
        for separator in separators:
            for f in separator:
                omitted = [x for x in separator if x != f]
                conflicts = independent_separator_check(worlds, omitted)
                necessity.append({'deleted_predicate': f, 'opposite_outcome_pair': conflicts[0] if conflicts else None})
                failure_count += not bool(conflicts)
            failure_count += bool(independent_separator_check(worlds, separator))
        chains[chain] = {'fields': fields, 'worlds': worlds, 'world_count': len(worlds), 'current_answer': answer(chain, {}), 'global_minimum_separators': separators, 'necessity_witnesses': necessity, 'per_answer_certificates': certs, 'resolution': 'PER_ARCH' if chain == 'BD10' else 'PER_TOOTH', 'timescale': 'SIMULTANEOUS'}
    result = {'round': 'R3', 'claim_type': 'capability', 'patient_id': case['patient_id'], 'snapshot_row': row, 'FDI': 16, 'chains': chains, 'certificate_failures': int(failure_count), 'certificate_gate': 'PASS' if not failure_count else 'FAIL', 'external_referent': prediction['external_referents'], 'independent_measurement_performed': False, 'clinical_answer': 'UNKNOWN', 'logical_enclosure': 'Exhaustive on declared clinical predicates; not a bound on actual clinical/anatomical error', 'cost': {'wall_s': time.perf_counter() - start, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'fit_s': 0, 'clinical_measurement': 'UNKNOWN_COST_NOT_PERFORMED', 'questions': 0}}
    write('raw/RESULTS_R3.json', result)
    checkpoint('R3_COMPUTED', result['certificate_gate'] + ' finite measurement-sufficiency; clinical UNKNOWN', 'Run external vectors, exact summary counterexamples, source audit and falsifying injections; package demo')
    print('R3', case['patient_id'], 'row', row, 'FDI16', 'certificate_gate', result['certificate_gate'])
if __name__ == '__main__':
    run()
