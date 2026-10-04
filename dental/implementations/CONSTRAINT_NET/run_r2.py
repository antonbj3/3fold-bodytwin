import copy, json, time
from common import *
from stress_and_queries import X71D, X71M, x71_portfolio, unlocked, sufficiency

def enumerate_contracts(decisions, measurements):
    mids = [m['id'] for m in measurements]
    symbols = sorted({s for m in measurements for s in m['supplies']} | {s for d in decisions for s in d['requires']})
    sb = {s: 1 << i for (i, s) in enumerate(symbols)}
    mb = {s: 1 << i for (i, s) in enumerate(mids)}
    supply_masks = [sum((sb[s] for s in set(m['supplies']))) for m in measurements]
    required_masks = [sum((sb[s] for s in set(d['requires']))) for d in decisions]
    prereq_masks = [sum((mb[s] for s in set(m['prerequisites']))) for m in measurements]
    mismatch = 0
    valid = 0
    hist = {}
    examples = {}
    for bits in range(1 << len(mids)):
        chosen = [i for i in range(len(mids)) if bits & 1 << i]
        selected = {mids[i] for i in chosen}
        supplied = {s for i in chosen for s in measurements[i]['supplies']}
        candidate_valid = all((set(measurements[i]['prerequisites']) <= selected for i in chosen))
        candidate = [d['id'] for d in decisions if set(d['requires']) <= supplied and (not d['other_blockers'])] if candidate_valid else []
        sm = 0
        for i in chosen:
            sm |= supply_masks[i]
        oracle_valid = all((bits & prereq_masks[i] == prereq_masks[i] for i in chosen))
        oracle = []
        if oracle_valid:
            for (j, d) in enumerate(decisions):
                if required_masks[j] & sm == required_masks[j] and len(d['other_blockers']) == 0:
                    oracle.append(d['id'])
        mismatch += int(candidate != oracle or candidate_valid != oracle_valid)
        valid += candidate_valid
        if candidate_valid:
            hist[str(len(candidate))] = hist.get(str(len(candidate)), 0) + 1
    return {'portfolios': 1 << len(mids), 'prerequisite_valid': valid, 'set_vs_bitmask_mismatches': mismatch, 'answerable_count_histogram': hist, 'all_controls_pass': mismatch == 0, 'scope': 'Only exact logic of frozen X71 inputs; costs, specimen compatibility, model acceptance and physical success are not inferred'}

def arithmetic():
    d = src('results/LANE_X36_META_REGRESSION/results.json')
    rows = []
    tol = read(HERE / 'PREREG_R2.json')['criteria']['contrast_arithmetic_absolute_tolerance']
    for (i, c) in enumerate(d['within_study_contrasts']):
        scale = 10 if c['dataset'] == 'cement' else 1
        predicted = scale * (c['high_outcome'] - c['low_outcome']) / (c['high_x'] - c['low_x'])
        error = abs(predicted - c['estimate'])
        rows.append({'edge_id': f'D-E-CONTRAST-{i + 1:02}', 'recomputed': predicted, 'reported': c['estimate'], 'absolute_error': error, 'unit': c['unit'], 'resolution_level': c['resolution_level'], 'pass': error <= tol, 'injected_plus_one_rejected': abs(predicted - (c['estimate'] + 1)) > tol})
    return {'rows': rows, 'tolerance': tol, 'pass': all((r['pass'] and r['injected_plus_one_rejected'] for r in rows)), 'scope': 'Independent re-arithmetic of reviewed primary-table contrasts; does not make paired-arm covariance or transfer error known'}

def main():
    start = time.perf_counter()
    ds = src(X71D)
    ms = src(X71M)
    exhaustive = enumerate_contracts(ds, ms)
    cases = [[], ['M01'], ['M02'], ['M01', 'M02'], ['M01', 'M02', 'M03'], ['M01', 'M02', 'M04'], ['M07'], ['M10'], ['M11']]
    examples = [x71_portfolio(ids) for ids in cases]
    counter = sufficiency()
    injected_count_detected = len(unlocked(counter['state_A'], {'e1'})) + 1 != len(counter['unlocked_A'])
    forces = [[1, 3], [2, 2]]
    total = [sum(f) for f in forces]
    coeff = [1, 2]
    response = [sum((a * b for (a, b) in zip(f, coeff))) for f in forces]
    force_test = {'states_N': forces, 'summary_total_N': total, 'identity_error_N': abs(total[0] - total[1]), 'response_coefficients_MPa_per_N': coeff, 'responses_MPa': response, 'downstream_difference_MPa': abs(response[0] - response[1]), 'resolution_level': 'PER_SURFACE_REGION', 'leaf_status': 'DERIVED_UNDER_ASSUMPTIONS', 'minimal_extension': 'Registered per-patch force vector (one independent share for two patches with known total)', 'rigorous_enclosure': {'input_domain': 'F_i >= 0; F1+F2=4 N', 'exact_response_range_MPa': [4, 8], 'proof': 's=F1+2F2=4+F2, 0<=F2<=4; endpoint extrema attained'}, 'external_referent': {'kind': 'our_own_fixture', 'locator': 'run_r2.py', 'compared_quantity': 'Declared two-patch linear operator; not measured dental stress', 'refutes_us': True}}
    control = arithmetic()
    report = {'claim_type': 'capability', 'outcome': 'EXACT_CONSUMER_INCIDENCE_EXECUTABLE_PHYSICAL_CALIBRATION_UNKNOWN', 'exhaustive': exhaustive, 'examples': examples, 'source_arithmetic': control, 'sufficiency_extension': {'counterexample_retained': counter, 'wrong_downstream_count_rejected': injected_count_detected, 'force_summary_counterexample': force_test}, 'runtime_seconds': time.perf_counter() - start, 'physical_measurements': 0, 'physical_certificates': 0, 'remaining_obstruction': 'Actual observation values, conditional model closures and matched uncertainty remain missing', 'external_referent': {'kind': 'external_review', 'locator': str(path(X71D)), 'compared_quantity': 'Frozen X71 prospective dependency clauses checked without changing them', 'refutes_us': True}, 'checks_pass': exhaustive['all_controls_pass'] and control['pass'] and injected_count_detected}
    dump('RESULTS_R2.json', report)
    print(json.dumps({'checks_pass': report['checks_pass'], 'enumeration': exhaustive, 'examples': [(x['measurements'], x['count']) for x in examples], 'runtime_seconds': report['runtime_seconds']}))
    if not report['checks_pass']:
        raise SystemExit(1)
if __name__ == '__main__':
    main()
