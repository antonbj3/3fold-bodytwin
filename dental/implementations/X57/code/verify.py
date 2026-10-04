import copy
import itertools
import math
import time
from common import read, write, sha, ROOT, OUTPUT
from clinical_rules import severity_rbl, severity_cal, severity_pd, lesion_class, interval_image, validate_observation, periodontal
from measurement import make_worlds, minimum_separators, independent_separator_check, minimum_answer_certificates, answer
from extract import source_rows, extract_row, case_id

def raised(fn):
    try:
        fn()
    except (ValueError, AssertionError):
        return True
    return False

def verify():
    start = time.perf_counter()
    vectors = read('sources/REFERENCE_VECTORS.json')
    table = []
    functions = {'RBL': severity_rbl, 'CAL': severity_cal, 'PD': severity_pd}
    for v in vectors['periodontal']:
        got = functions[v['quantity']](v['value'])
        table.append({'id': v['id'], 'source': v['source'], 'expected': v['expected_severity_band'], 'got': got, 'ok': got == v['expected_severity_band']})
    for v in vectors['caries']:
        got = lesion_class(v['depth_fraction'], v['barrier_present'])
        table.append({'id': v['id'], 'source': v['source'], 'expected': v['expected_class'], 'got': got, 'ok': got == v['expected_class']})
    rbl_a = [0.125, 0.375]
    rbl_b = [0.25, 0.25]
    rdt_a = [0.0, 0.5]
    rdt_b = [0.25, 0.25]
    rbl_mean_a = math.fsum(rbl_a) / 2
    rbl_mean_b = math.fsum(rbl_b) / 2
    rdt_mean_a = math.fsum(rdt_a) / 2
    rdt_mean_b = math.fsum(rdt_b) / 2
    sufficient = {'referent': {'kind': 'our_own_fixture', 'locator': 'code/verify.py::exact_binary_counterexamples', 'compared_quantity': 'Exact aggregation sufficiency only; no measured patient anatomy', 'refutes_us': True}, 'RBL': {'states': [rbl_a, rbl_b], 'unit': 'root-length fraction', 'summary': 'regional mean', 'summary_resolution': 'PER_TOOTH', 'input_resolution': 'PER_SURFACE_REGION', 'timescale': 'SIMULTANEOUS', 'summary_values': [rbl_mean_a, rbl_mean_b], 'float_hex': [rbl_mean_a.hex(), rbl_mean_b.hex()], 'identity_error': abs(rbl_mean_a - rbl_mean_b), 'downstream_severity_bands': [severity_rbl(max(rbl_a)), severity_rbl(max(rbl_b))], 'downstream_difference_ordinal_bands': severity_rbl(max(rbl_a)) - severity_rbl(max(rbl_b)), 'minimum_extension_for_this_query': 'Worst-region RBL plus its source region; other stage/case modifiers still required'}, 'RDT': {'states': [rdt_a, rdt_b], 'unit': 'mm', 'summary': 'regional mean', 'summary_resolution': 'PER_TOOTH', 'input_resolution': 'PER_SURFACE_REGION', 'timescale': 'SIMULTANEOUS', 'summary_values': [rdt_mean_a, rdt_mean_b], 'float_hex': [rdt_mean_a.hex(), rdt_mean_b.hex()], 'identity_error': abs(rdt_mean_a - rdt_mean_b), 'downstream_zero_barrier_present': [min(rdt_a) == 0, min(rdt_b) == 0], 'downstream_difference_boolean': int(min(rdt_a) == 0) - int(min(rdt_b) == 0), 'minimum_extension_for_this_query': 'Minimum remaining barrier thickness and region. ESE class also needs local total dentin/depth and hard/firm character; physiology needs pulp predicates'}, 'gate': 'PASS' if rbl_mean_a == rbl_mean_b and rdt_mean_a == rdt_mean_b and (severity_rbl(max(rbl_a)) != severity_rbl(max(rbl_b))) and ((min(rdt_a) == 0) != (min(rdt_b) == 0)) else 'FAIL'}
    interval_checks = []
    for (quantity, breaks, fn) in [('RBL', [0.15, 0.33], severity_rbl), ('PD', [6.0], severity_pd), ('CAL', [3.0, 5.0], severity_cal)]:
        intervals = [(0, breaks[0]), (breaks[0], breaks[-1]), (breaks[-1], 1 if quantity == 'RBL' else 10)]
        for (lo, hi) in intervals:
            exact = interval_image([lo, hi], fn, breaks)
            dense = [lo + (hi - lo) * i / 100 for i in range(101)]
            for b in breaks:
                dense.extend([b, math.nextafter(b, -math.inf), math.nextafter(b, math.inf)])
            got = set((fn(v) for v in dense if lo <= v <= hi))
            interval_checks.append({'quantity': quantity, 'interval': [lo, hi], 'exact_set': exact, 'sample_outcomes_subset': got <= set(exact), 'rigorous_basis': 'All boundary points and open cells of piecewise-constant source predicate; sampling is secondary check'})
    obs = read('raw/OBSERVATIONS_R2.json')
    manual = {'O-3021adc41560c1e9': {'teeth': [22, 47], 'resolution': 'PER_SURFACE_REGION', 'reason': '22D,47D source tokens explicitly indicate distal surfaces; manual region label refined from source after residual audit'}, 'O-ace6bdd732052fa4': {'teeth': [25], 'resolution': 'PER_TOOTH', 'reason': '25CBCT source token'}, 'O-9588594c0214a012': {'teeth': [11, 21], 'resolution': 'PER_SURFACE_REGION', 'reason': 'Source explicitly says mesial pocket'}}
    sample = read('raw/REVIEW_SAMPLE_R1.json')
    r1bad = []
    audited = []
    current = {x['observation_id']: x for x in obs}
    for original in sample:
        oid = original['observation_id']
        expected = manual.get(oid, {'teeth': original['teeth'], 'resolution': original['resolution'], 'reason': 'Producer inspected literal observation and neighboring FDI scope'})
        old_ok = all((original[k] == expected[k] for k in ['teeth', 'resolution']))
        now = current[oid]
        new_ok = all((now[k] == expected[k] for k in ['teeth', 'resolution']))
        if not old_ok:
            r1bad.append(oid)
        audited.append({'observation_id': oid, 'expected': expected, 'R1_ok': old_ok, 'R2_ok': new_ok, 'reviewer': 'producer', 'independent_review': False})
    jointrows = read('raw/RESULTS_R2.json')['joint_localized_snapshot_rows']
    for x in obs:
        if x['snapshot_row'] in jointrows and x['observation_id'] not in {v['observation_id'] for v in audited}:
            expected_surface = 'furcation' if x['quantity'] == 'FURCATION' else 'mesial'
            ok = x['teeth'] == [16] and x['surface'] == expected_surface and (x['resolution'] == 'PER_SURFACE_REGION')
            audited.append({'observation_id': x['observation_id'], 'expected': {'teeth': [16], 'surface': expected_surface, 'resolution': 'PER_SURFACE_REGION'}, 'R2_ok': ok, 'reviewer': 'producer', 'independent_review': False})
    audit = {'reviewer': 'producer', 'independent_review': False, 'R1_n': len(sample), 'R1_wrong': r1bad, 'R1_fraction_correct': (len(sample) - len(r1bad)) / len(sample), 'R1_gate': 'PASS' if (len(sample) - len(r1bad)) / len(sample) >= 0.95 else 'FAIL', 'R2_n': len(audited), 'R2_fraction_correct': sum((x['R2_ok'] for x in audited)) / len(audited), 'R2_gate': 'PASS' if sum((x['R2_ok'] for x in audited)) / len(audited) >= 0.95 else 'FAIL', 'rows': audited, 'sample_scope': 'Original fixed 30 plus newly linked shared-tooth observations; does not estimate all-report sensitivity'}
    (rows, volumes, csv_hash) = source_rows()
    source_checks = []
    for o in obs:
        r = rows[o['snapshot_row'] - 2]
        text = r['Oral Check'].replace('\\n', ' ').replace('\n', ' ')
        loc = o['locator']
        ok = case_id(r['Filename']) == o['patient_id'] and sha(r['Oral Check'].encode()) == loc['original_column_sha256'] and (csv_hash == loc['csv_sha256']) and (text[loc['span_start']:loc['span_end']] == o['match_text'])
        source_checks.append({'observation_id': o['observation_id'], 'ok': ok})
    pd = next((x for x in obs if x['quantity'] == 'PD'))
    fg = next((x for x in obs if x['quantity'] == 'FURCATION'))
    badunit = copy.deepcopy(pd)
    badunit['unit'] = 'm'
    badgrade = copy.deepcopy(fg)
    badgrade['interval'] = [4, 4]
    (fields, worlds) = make_worlds('BD10')
    (cfields, cworlds) = make_worlds('BD08')
    complete_min = minimum_separators(worlds, fields)[0]
    wrong_case = copy.deepcopy(obs[0])
    wrong_case['patient_id'] = 'M-WRONG'
    forged = copy.deepcopy(read('FROZEN_PREDICTIONS.json'))
    forged['FDI'] = 17
    import json
    original_pred = read('FROZEN_PREDICTIONS.json')
    canonical = lambda p: (json.dumps(p, ensure_ascii=False, indent=2, allow_nan=False) + '\n').encode()
    good_frozen_hash = sha(canonical(original_pred)) == sha((OUTPUT / 'FROZEN_PREDICTIONS.json').read_bytes())
    from clinical_rules import caries
    near = next((x for x in obs if x['quantity'] == 'NEAR_PULP_REPORTED'))
    true_near = caries([near])
    forged_near = copy.deepcopy(true_near)
    forged_near['selective_excavation_eligibility'] = 'ELIGIBLE_RESEARCH_BRANCH'
    injections = [{'name': 'RBL10 measured value changed to0.50', 'rejected': severity_rbl(0.5) != vectors['periodontal'][0]['expected_severity_band'], 'control': 'external source vector'}, {'name': 'PD given in meters without unit conversion', 'rejected': raised(lambda : validate_observation(badunit)), 'control': 'dimension validation'}, {'name': 'Furcation grade4', 'rejected': raised(lambda : validate_observation(badgrade)), 'control': 'source category domain'}, {'name': 'Claim complete stage without case/IV predicates', 'rejected': len(answer('BD10', {})['possible_outcomes']) > 1, 'control': 'logical outcome image'}, {'name': 'Cold(-) used as a Boolean vascular pulp measurement', 'rejected': raised(lambda : answer('BD08', {'pulp_state_compatible': 'cold(-)'})), 'control': 'measurement type and semantics'}, {'name': 'Delete stage-IV history predicate from separator', 'rejected': bool(independent_separator_check(worlds, complete_min[:-1])), 'control': 'opposite-outcome pairs'}, {'name': 'Delete exposure predicate from caries separator', 'rejected': bool(independent_separator_check(cworlds, cfields[:-1])), 'control': 'opposite-outcome pairs'}, {'name': 'Wrong patient binding to real source observation', 'rejected': case_id(rows[wrong_case['snapshot_row'] - 2]['Filename']) != wrong_case['patient_id'], 'control': 'source case replay'}, {'name': 'Forged FDI in frozen future-measurement packet', 'rejected': good_frozen_hash and sha(canonical(forged)) != sha((OUTPUT / 'FROZEN_PREDICTIONS.json').read_bytes()), 'control': 'same canonical frozen bytes, accepts correct baseline'}, {'name': 'Near-pulp label asserted eligible without barrier', 'rejected': forged_near != caries([near]) and true_near == caries([near]), 'control': 'pipeline reproduction, accepts correct baseline'}]
    expected = read('FROZEN_PREDICTIONS.json')['frozen_expected_model_certificate_predicate_counts']
    cert_checks = []
    for chain in ['BD10', 'BD08']:
        (fields, worlds) = make_worlds(chain)
        count = len(minimum_separators(worlds, fields)[0])
        cert_checks.append({'quantity': chain + '_global', 'expected': expected[chain + '_global'], 'got': count, 'ok': count == expected[chain + '_global']})
        for w in worlds:
            n = len(minimum_answer_certificates(worlds, fields, w['state'])[0])
            if chain == 'BD10':
                key = 'BD10_answerIII' if w['outcome'] == 'III' else 'BD10_caseexcluded' if w['outcome'] == 'NOT_CONFIRMED_PERIODONTITIS_CASE' else 'BD10_answerIV_minimum'
            else:
                key = 'BD08_eligible' if w['outcome'] == 'ELIGIBLE_RESEARCH_BRANCH' else 'BD08_outside_minimum'
            cert_checks.append({'quantity': key, 'expected': expected[key], 'got': n, 'ok': n == expected[key]})
    passed = all((x['ok'] for x in table)) and sufficient['gate'] == 'PASS' and all((x['sample_outcomes_subset'] for x in interval_checks)) and all((x['ok'] for x in source_checks)) and all((x['rejected'] for x in injections)) and all((x['ok'] for x in cert_checks)) and (audit['R2_gate'] == 'PASS')
    result = {'status': 'PASS' if passed else 'FAIL', 'external_table': table, 'sufficiency': sufficient, 'rigorous_interval_image': interval_checks, 'source_scope_audit': audit, 'source_replay': source_checks, 'injections': injections, 'frozen_certificate_checks': cert_checks, 'wall_s': time.perf_counter() - start, 'clinical_validation': 'NOT_PERFORMED', 'linear_sensitivity': 'NOT_REPORTED'}
    write('raw/VERIFICATION.json', result)
    print('VERIFICATION', result['status'], 'R1 audit', audit['R1_gate'], 'R2 audit', audit['R2_gate'])
    if not passed:
        raise SystemExit(1)
if __name__ == '__main__':
    verify()
