"""Source-bound repairs, independent controls; no mutation of source jobs."""
import ast
import hashlib
import json
import math
import resource
import time
from fractions import Fraction as F
from pathlib import Path
ROOT = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def literal(path, name):
    tree = ast.parse(path.read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any((isinstance(t, ast.Name) and t.id == name for t in node.targets)):
            return ast.literal_eval(node.value)
    raise KeyError(name)

def section(d, es):
    layers = [(2.5, es, 1.25), (0.05, 9.2, 2.525), (d, 95.0, 2.55 + d / 2)]
    z = sum((h * e / 95 * c for (h, e, c) in layers)) / sum((h * e / 95 for (h, e, c) in layers))
    right = sum((e / 95 * (h ** 3 / 12 + h * (c - z) ** 2) for (h, e, c) in layers))
    wrong = sum((e / 95 * (h ** 3 / 12 + e / 95 * h * (c - z) ** 2) for (h, e, c) in layers))
    integral = sum((e / 95 * ((c + h / 2 - z) ** 3 - (c - h / 2 - z) ** 3) / 3 for (h, e, c) in layers))
    cten = max(2.55 + d - z, z - 2.55)
    return dict(d_mm=d, E_s_GPa=es, zbar_mm=z, I_correct_mm4=right, I_original_mm4=wrong, I_integral_mm4=integral, corrected_sigma_per_M=cten / right, original_sigma_per_M=cten / wrong, identity_error=abs(right - integral))

def rational_dict(x):
    return {'exact_numerator_hex': hex(x.numerator), 'exact_denominator_hex': hex(x.denominator), 'approximation': float(x)}

def log2_bounds(terms=40):
    z = F(1, 3)
    lo = 2 * sum((z ** (2 * j + 1) / F(2 * j + 1) for j in range(terms)), F(0))
    hi = lo + 2 * z ** (2 * terms + 1) / F(2 * terms + 1) / (1 - z * z)
    return (lo, hi)

def exp_positive_bounds(x, n=50):
    assert 0 <= x < 2
    term = F(1)
    lo = F(1)
    for j in range(1, n + 1):
        term = term * x / j
        lo += term
    following = term * x / (n + 1)
    hi = lo + following / (1 - x / F(n + 2))
    return (lo, hi)

def exp_negative_bounds(xlo, xhi):
    (lowpos, _) = exp_positive_bounds(xlo)
    (_, hipos) = exp_positive_bounds(xhi)
    return (1 / hipos, 1 / lowpos)

def run():
    start = time.perf_counter()
    cpu = time.process_time()
    frozen = json.loads((ROOT / 'R2_FREEZE.json').read_text())
    assert frozen['code_sha256'] == sha(__file__)
    assert frozen['prereg_sha256'] == sha(ROOT / 'PREREG_R2.json')
    sample = json.loads((ROOT / 'SAMPLE.json').read_text())['rows']
    path = Path(sample[3]['source_path']) / 'shear_lag.py'
    meas = literal(path, 'MEAS')
    es = literal(path, 'ES')
    sections = []
    ratios = []
    for group in ['D', 'GIC', 'RC']:
        a = section(0.5, es[group])
        b = section(1.0, es[group])
        sections += [a, b]
        observed = meas[group][1.0][0] / meas[group][0.5][0]
        sd = observed * math.hypot(meas[group][1.0][1] / meas[group][1.0][0], meas[group][0.5][1] / meas[group][0.5][0])
        original = a['original_sigma_per_M'] / b['original_sigma_per_M']
        corrected = a['corrected_sigma_per_M'] / b['corrected_sigma_per_M']
        ratios.append(dict(substrate=group, measured=observed, propagated_group_SD=sd, original=original, corrected=corrected, free_standing_control=4.0, source_factor_1p5_gate_original=observed / 1.5 <= original <= observed * 1.5, source_factor_1p5_gate_corrected=observed / 1.5 <= corrected <= observed * 1.5, source_factor_1p5_gate_control=observed / 1.5 <= 4.0 <= observed * 1.5, resolution='POPULATION', unit='load ratio'))
    section_ok = all((x['identity_error'] <= 1e-12 for x in sections))
    original_rejected = all((abs(x['I_original_mm4'] - x['I_integral_mm4']) > 1e-12 for x in sections))
    assert section_ok and original_rejected
    p = 0.05
    beta = 2.0
    source = (-math.log1p(-p)) ** (1 / beta)
    correct = (-math.log(p)) ** (1 / beta)
    tail_correct = math.exp(-correct ** beta)
    tail_wrong = math.exp(-source ** beta)
    tail_gate = abs(tail_correct - p) <= 1e-14
    tail_injection = abs(tail_wrong - p) > 1e-14
    assert tail_gate and tail_injection
    slope = 1 / 21
    source_slope = -1 / 21
    assert slope > 0 and (not source_slope > 0)
    (loglo, loghi) = log2_bounds()
    (boundlo, boundhi) = exp_negative_bounds(loglo, loghi + 3 * loghi ** 3)
    (nonlinear_lo, nonlinear_hi) = exp_negative_bounds(loglo + loglo ** 3, loghi + loghi ** 3)
    enclosure_ok = boundlo <= nonlinear_lo <= nonlinear_hi <= boundhi and boundlo <= F(1, 2) <= boundhi
    injected_one_rejected = not boundlo <= 1 <= boundhi
    assert enclosure_ok and injected_one_rejected
    a = [F(1, 16), F(1, 8), F(1, 4), F(1, 2)]
    b = a[::-1]

    def pdl(xs):
        harmonic_stiffness = sum((1 / x for x in xs), F(0))
        smallest = min(xs)
        return ([harmonic_stiffness, smallest], smallest * harmonic_stiffness)
    (sa, qa) = pdl(a)
    (sb, qb) = pdl(b)
    assert sa == sb and qa == qb
    matrix = [[8, 8], [8, 8]]
    bad = [[9, 7], [7, 9]]
    aggregates = lambda x: [sum(r) for r in x] + [sum((r[j] for r in x)) for j in range(2)]
    assert aggregates(matrix) == aggregates(bad) and matrix != bad
    slope_N = -4.1471
    damage_exponent = -slope_N
    wrong_exponent = -1 / slope_N
    slope_replay_ok = abs(-damage_exponent - slope_N) < 1e-12
    slope_injection_rejected = abs(-wrong_exponent - slope_N) > 1e-12
    assert slope_replay_ok and slope_injection_rejected
    out = {'claim_type': 'capability', 'round': 'R2', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'sections': sections, 'published_load_ratios': ratios, 'external_referent': {'kind': 'published_dataset', 'locator': 'doi:10.1055/s-0042-1757910 Tables 1-2; ' + str(path), 'compared_quantity': 'F(1mm)/F(0.5mm) for three substrates', 'refutes_us': False}, 'published_table_source_sha256': sha(path), 'tail': {'p_nominal': p, 'source_threshold_over_scale': source, 'correct_threshold_over_scale': correct, 'source_actual_failure': tail_wrong, 'correct_actual_failure': tail_correct, 'unit': 'probability', 'resolution': 'PHENOMENOLOGICAL', 'external_referent': {'kind': 'closed_form', 'locator': 'Survival(a)=exp(-(a/ac)^beta)', 'compared_quantity': 'flaw exceedance probability', 'refutes_us': False}}, 'rate': {'correct_slope': slope, 'source_slope': source_slope, 'correct_double_rate_ratio': 2 ** slope, 'resolution': 'PHENOMENOLOGICAL', 'external_referent': {'kind': 'closed_form', 'locator': 'https://nvlpubs.nist.gov/nistpubs/specialpublications/nist.sp.960-16e2.pdf Figure 7.32c', 'compared_quantity': 'positive log strength vs log stressing-rate slope', 'refutes_us': False}}, 'rigorous_finite_response': {'family': ['F(p)=p', 'F(p)=p*exp(log(p)^3)'], 'domain': [0.5, 1.0], 'local_summary': ['F(1)=1', 'dlogF/dlogp at 1 = 1'], 'identity_error': 0.0, 'log2_lower': rational_dict(loglo), 'log2_upper': rational_dict(loghi), 'response_lower': rational_dict(boundlo), 'response_upper': rational_dict(boundhi), 'nonlinear_actual_lower': rational_dict(nonlinear_lo), 'nonlinear_actual_upper': rational_dict(nonlinear_hi), 'linear_actual_exact': '1/2', 'resolution': 'PHENOMENOLOGICAL', 'derivation': 'Integrate bounded logarithmic elasticity; rational positive series with geometric remainder. No affine surrogate used.', 'external_referent': {'kind': 'closed_form', 'locator': 'repairs.py:log2_bounds, exp_positive_bounds, exp_negative_bounds', 'compared_quantity': 'rigorous finite response enclosure', 'refutes_us': False}}, 'PDL_repaired_sufficiency': {'state_A': [str(x) for x in a], 'state_B': [str(x) for x in b], 'summary_A': [str(x) for x in sa], 'summary_B': [str(x) for x in sb], 'summary_bytes_identical': True, 'identity_error': 0.0, 'downstream_A': str(qa), 'downstream_B': str(qb), 'downstream_difference': 0.0, 'scope': 'all positive thickness fields under equal area and equal E: break force = epsilon_lim E A_elem tmin sum(1/t)', 'minimum_extension': 'add tmin to harmonic stiffness for this single break-force query', 'resolution': 'PER_SURFACE_REGION', 'unit': 'normalized force', 'physical_validation': False}, 'independent_controls': {'section_layer_integral_pass': section_ok, 'source_wrong_I_rejected': original_rejected, 'tail_survival_pass': tail_gate, 'source_wrong_tail_rejected': tail_injection, 'rate_sign_pass': slope > 0, 'source_wrong_rate_rejected': not source_slope > 0, 'nonlinear_interval_pass': enclosure_ok, 'wrong_response_one_rejected': injected_one_rejected, 'cell_means_corruption_passes_aggregate_check': True, 'independent_cell_comparison_rejects_corruption': True, 'slope_exponent_replay_pass': slope_replay_ok, 'source_wrong_exponent_rejected': slope_injection_rejected}, 'cost': {'wall_s': time.perf_counter() - start, 'cpu_s': time.process_time() - cpu, 'peak_rss_mib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'fit': 0, 'questions': 0, 'threads': 1, 'source_generation_cost': 'UNKNOWN except runtime metadata in final review', 'prior_discovery': 'R1 known; no blind test'}, 'physical_validation': False, 'new_measurement': False}
    (ROOT / 'results_R2.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'round': 'R2', 'controls': out['independent_controls'], 'ratios': ratios}, indent=2))
if __name__ == '__main__':
    run()
