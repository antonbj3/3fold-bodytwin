import json, csv, copy, math, inspect, time
from fractions import Fraction as Q
import numpy as np
import scipy.optimize
from geometry import R, sha, write
from validate import w1_bounds, KNOTS, REF
from force_port import project
checks = []

def check(name, val):
    checks.append(dict(name=name, pass_=bool(val)))
    if not val:
        raise AssertionError(name)

def exact_w1(x, values):
    x = sorted(map(Q, x))
    n = len(x)
    knots = list(map(Q, ['0', '.05', '.10', '.25', '.50', '.75', '.90', '.95', '1']))
    vals = [Q(str(v)) / 100 for v in values]
    rank = Q(1, 178)
    rd = Q(1, 2000)
    cuts = sorted(set([Q(k, n) for k in range(n + 1)] + [max(Q(0), min(Q(1), u + s * rank)) for u in knots for s in [-1, 1]]))
    lo = hi = Q(0)
    for (a, b) in zip(cuts[:-1], cuts[1:]):
        u = (a + b) / 2
        v = x[min(int(u * n), n - 1)]
        l = [i for (i, k) in enumerate(knots) if k <= u - rank]
        h = [i for (i, k) in enumerate(knots) if k >= u + rank]
        left = max(Q(0), vals[l[-1]] - rd) if l else Q(0)
        right = min(Q(1), vals[h[0]] + rd) if h else Q(1)
        lo += (b - a) * max(left - v, v - right, Q(0))
        hi += (b - a) * max(abs(v - left), abs(v - right))
    return [lo, hi]

def directed(q, lower):
    f = float(q)
    if Q(f) > q if lower else Q(f) < q:
        f = float(np.nextafter(f, -math.inf if lower else math.inf))
    return f

def main():
    start = time.perf_counter()
    for name in ['PREREG_R1.json', 'PREREG_R2.json', 'PREREG_R3.json', 'PREREG_R4.json', 'FROZEN_PREDICTIONS.json', 'FROZEN_PREDICTIONS_R4.json']:
        check('frozen_hash_' + name, sha(R / name) == (R / (name + '.sha256')).read_text().split()[0])
    frozen = json.load(open(R / 'FROZEN_PREDICTIONS.json'))
    for (p, h) in frozen['source_sha256'].items():
        check('source_hash_' + p, sha(p) == h)
    validate = json.load(open(R / 'raw/EXTERNAL_VALIDATION.json'))
    split = json.load(open(R / 'raw/PATIENT_PARTITION.json'))['partition']
    test = set(split['test'])
    check('train_cal_test_disjoint', all((not set(a) & set(b) for (i, a) in enumerate(split.values()) for b in list(split.values())[i + 1:])))
    rows = [r for r in csv.DictReader(open(R / 'raw/GEOMETRY_PROFILES.csv')) if r['case'] in test and r['jaw'] == 'upper' and (r['area_usable'] == 'True')]
    cases = {}
    for r in rows:
        cases.setdefault(r['case'], {})[int(r['fdi'])] = float(r['area_share'])
    exact = {}
    original_vals = {'R': [25.2, 33.6, 38.3, 45.6, 49.6, 55.5, 60.2, 65.8, 77.0], 'posterior': [39.8, 54.5, 62.4, 75.8, 82.4, 87.3, 90.4, 94.0, 99.6]}
    eps = validate['DKW']['external_epsilon'] + validate['DKW']['local_epsilon']
    for tag in ['R', 'posterior']:
        x = [sum((v for (t, v) in p.items() if (t // 10 in (1, 4) if tag == 'R' else t % 10 >= 4))) for p in cases.values()]
        q = exact_w1(x, original_vals[tag])
        f = w1_bounds(x, vals=REF[tag])
        check('independent_rational_W1_' + tag, max((abs(float(a) - b) for (a, b) in zip(q, f))) <= 1e-12)
        exact[tag] = dict(lower_exact=str(q[0]), upper_exact=str(q[1]), outward_float_enclosure=[directed(q[0], True), directed(q[1], False)], population_DKW_W1_enclosure=[max(0, float(q[0]) - eps), min(1, float(q[1]) + eps)], population_caveat='95% joint DKW under iid/same observable only; systematic sensor/cohort/FDI/pose errors UNKNOWN')
    tr = json.load(open(R / 'raw/TRANSPORT_SUMMARY.json'))
    check('real_7796_LP_queries', tr['LP_calls'] == 7796 and tr['max_LP_bound_error_share'] <= 1e-08)
    from transport import port, lp_check, marginal_bounds, marginal
    witness = json.load(open(R / 'raw/TRANSPORT_WITNESS.json'))
    e = witness['edges']
    forces = list(map(Q, witness['force1']))
    (upper, lower, cross) = port(e, forces, witness['axis'])
    p = marginal(e, forces, 0)
    bounds = marginal_bounds(e, p, witness['axis'])
    bounds[1] += Q(1, 100)
    (actual_fault_error, _) = lp_check(e, witness['axis'], bounds, p=p)
    check('injected_plus1pp_LP_bound_fails', actual_fault_error > 1e-08)
    observation = dict(context=dict(evidence_kind='SIMULATED', time_relation='SIMULTANEOUS', unit='force_share_same_total', same_pose_id=witness['case'], common_denominator_id='sum-edge-force-fixture', incidence_sha256=sha(R / 'raw/TRANSPORT_WITNESS.json')), upper_region_interval_share=[str(upper)] * 2, signed_net_crossing_interval_share=[str(cross)] * 2)
    a = project(observation)
    check('port_exact_conservation', a['status'] == 'CONDITIONAL_BOUND' and list(map(Q, a['lower_region_interval_share'])) == [lower, lower])
    bad = copy.deepcopy(observation)
    bad['signed_net_crossing_interval_share'] = [str(cross + Q(1, 100))] * 2
    b = project(bad)
    check('injected_crossing_measurement_fails_truth', b['status'] == 'CONDITIONAL_BOUND' and (not Q(b['lower_region_interval_share'][0]) <= lower <= Q(b['lower_region_interval_share'][1])))
    uncal = copy.deepcopy(observation)
    uncal['context']['evidence_kind'] = 'MEASURED'
    check('uncalibrated_sensor_refused', project(uncal)['status'] == 'UNKNOWN')
    wrong = copy.deepcopy(observation)
    wrong['context']['time_relation'] = 'HANDOVER'
    check('nonsimultaneous_force_refused', project(wrong)['status'] == 'UNKNOWN')
    err = copy.deepcopy(observation)
    err['upper_region_interval_share'] = ['1', '0']
    check('reversed_interval_rejected', project(err)['status'] == 'REJECTED')
    uncertain = copy.deepcopy(observation)
    uncertain['upper_region_interval_share'] = list(map(str, [upper - Q(1, 100), upper + Q(1, 100)]))
    uncertain['signed_net_crossing_interval_share'] = list(map(str, [cross - Q(1, 200), cross + Q(1, 200)]))
    ans = project(uncertain)
    check('analytic_uncertainty_inclusion', all((Q(ans['lower_region_interval_share'][0]) <= q + d <= Q(ans['lower_region_interval_share'][1]) for q in map(Q, uncertain['upper_region_interval_share']) for d in map(Q, uncertain['signed_net_crossing_interval_share']))))
    write('raw/MEASUREMENT_INPUT_SIMULATED.json', observation)
    write('raw/MEASUREMENT_OUTPUT_SIMULATED.json', a)
    ferr = validate['categorical_diagnostics'][0]
    check('Ferrato_raw_sum93_preserved', abs(ferr['raw_mass'] - 0.93) < 1e-12 and ferr['normalization_gate'] == 'FAIL')
    check('injected_Ferrato_sum100_fails_source_check', abs(1 - ferr['raw_mass']) > 0.0005 * 14)
    check('R1_force_mutation_rejected', json.load(open(R / 'raw/FORCE_WITNESS.json'))['injected_plus1N_rejected'])
    for (tag, val) in validate['injected_profiles_rejected'].items():
        check('quantile_false_profile_' + tag, val)
    check('R4_missing_label_mutation_rejected', json.load(open(R / 'raw/SUBGROUP_DENOMINATOR.json'))['mutated_missing_field_rejected'])
    source = inspect.getsourcefile(scipy.optimize.linprog)
    write('raw/RIGOROUS_W1.json', dict(source_quantile_unit='force share, externalTable1', exact_integrals=exact, rank_slack='1/178', rounding_share='1/2000', bound_type='Nonparametric rigorous deterministic enclosure under given quantile observation contract; monotonicity envelope, not synthetic external individual data'))
    write('raw/VERIFICATION.json', dict(n_checks=len(checks), additional_real_LP_fault_check_calls=2, LP_injected_bound_error_share=actual_fault_error, all_pass=all((v['pass_'] for v in checks)), checks=checks, linprog_source_path=source, linprog_source_sha256=sha(source), wall_s=time.perf_counter() - start, scope='Producer self-verification only; not independent scientific review'))
    print(json.dumps(dict(n_checks=len(checks), all_pass=True, exact_W1_float={k: v['outward_float_enclosure'] for (k, v) in exact.items()}, physical_measured_force='NONE')))
if __name__ == '__main__':
    main()
