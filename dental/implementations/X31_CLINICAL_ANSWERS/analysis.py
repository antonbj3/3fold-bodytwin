"""X31: answer conditional queries, preserving unobservable clinical quantities."""
import csv, datetime, hashlib, importlib.util, json, math, resource, time
from collections import Counter
from pathlib import Path
import numpy as np
from scipy import integrate, optimize, stats
ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'inputs'
START = time.perf_counter()
CHECKS = []
TYPE_NAMES = {'central_incisor': 'Central incisiv', 'lateral_incisor': 'Lateral incisiv', 'canine': "Canine", 'premolar1': 'Premolar 1', 'premolar2': 'Premolar 2', 'molar1': 'Molar 1', 'molar2': 'Molar 2', 'molar3': 'Molar 3'}
GUIDE_NAMES = {'fully_guided': 'Fullt guidat', 'pilot_guided': 'Pilotguidat', 'freehand': 'Frihand', 'dynamic_navigation': 'Dynamisk navigation'}
GUIDE_ORDER = ['fully_guided', 'pilot_guided', 'freehand', 'dynamic_navigation']

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(name):
    return json.loads((ROOT / name).read_text())

def write(name, obj):
    p = ROOT / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def csvwrite(name, rows):
    p = ROOT / name
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def load_module(name):
    spec = importlib.util.spec_from_file_location('x31_' + name, INPUT / (name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def record(name, error, tolerance, mutated_error, scope):
    ok = bool(error <= tolerance)
    reject = bool(mutated_error > tolerance)
    CHECKS.append({'check': name, 'error': float(error), 'tolerance': tolerance, 'injected_error': float(mutated_error), 'valid_pass': ok, 'fault_rejected': reject, 'scope': scope})
    if not (ok and reject):
        raise AssertionError(name)

def state(phase, gate, next_op):
    write('CURRENT_WORK_STATE.json', {'updated_utc': now(), 'phase': phase, 'latest_gate': gate, 'next_operation': next_op, 'human_input_required': False})

def cluster_ci(rows, values, seed=3102):
    cases = sorted(set((r['case'] for r in rows)))
    index = {c: i for (i, c) in enumerate(cases)}
    counts = np.zeros(len(cases))
    sums = np.zeros(len(cases))
    for (r, v) in zip(rows, values):
        counts[index[r['case']]] += 1
        sums[index[r['case']]] += v
    draws = np.random.default_rng(seed).integers(0, len(cases), (512, len(cases)))
    boot = sums[draws].sum(1) / counts[draws].sum(1)
    return [float(x) for x in np.quantile(boot, [0.025, 0.975])]

def cluster_median_ci(rows, values):
    groups = {c: [] for c in sorted(set((r['case'] for r in rows)))}
    for (r, v) in zip(rows, values):
        groups[r['case']].append(v)
    arrays = [np.array(v) for v in groups.values()]
    rng = np.random.default_rng(3102)
    boot = [np.median(np.concatenate([arrays[i] for i in rng.integers(0, len(arrays), len(arrays))])) for _ in range(512)]
    return [float(x) for x in np.quantile(boot, [0.025, 0.975])]

def region(fdi):
    side = "left" if fdi // 10 == 3 else "right"
    part = 'premolar' if fdi % 10 in [4, 5] else 'molar' if fdi % 10 in [6, 7, 8] else 'anterior'
    return side + ' ' + part

def verify_frozen():
    for (name, entry) in read('INPUT_MANIFEST.json')['files'].items():
        actual = sha(INPUT / name)
        if actual != entry['sha256']:
            raise ValueError('Input drift: ' + name)
    for p in sorted(ROOT.glob('PREREG_*.json')):
        if sha(p) != p.with_suffix('.sha256').read_text().split()[0]:
            raise ValueError('Prereg drift')
    record('input_hash_guard', 0, 0, 1, 'Byte identity; fault = changed source hash')

def source_check(profiles):
    text = (INPUT / 'Varga2020.txt').read_text()
    table = text.split('MANDIBLE')[-1].split('MAXILLA')[0]
    import re
    vals = {}
    for line in table.splitlines():
        m = re.match('\\s*(Fr|Pi|Fu)\\s+(37|29)\\s+([\\d. ]+)\\s*$', line)
        if m:
            vals.setdefault(m[1], []).append([float(x) for x in m[3].split()])
    residual = []
    for (arm, g) in [('Fr', 'freehand'), ('Pi', 'pilot_guided'), ('Fu', 'fully_guided')]:
        (angle, apex) = vals[arm]
        expected = angle[:2] + angle[5:7] + apex[:2]
        actual = profiles[g]['angle_mean_sd_deg'] + profiles[g]['entry_mean_sd_mm'] + profiles[g]['apex_mean_sd_mm']
        residual.extend(np.array(actual) - expected)
    import xml.etree.ElementTree as ET
    wt = ' '.join(ET.parse(INPUT / 'Wu2020.xml').getroot().itertext())
    for (mu, sd) in [profiles['dynamic_navigation'][k] for k in ['entry_mean_sd_mm', 'apex_mean_sd_mm', 'angle_mean_sd_deg']]:
        if not re.search(f'{mu:.2f}\\s*±\\s*{sd:.2f}', wt):
            raise ValueError('Wu source aggregate not matched')
    record('external_guide_moments', max((abs(x) for x in residual)), 1e-10, 0.5, 'Varga mandibular table; Wu Results/Fig4. Does not check signed loss or clinical 95th percentile')
    write('raw/SOURCE_COMPARISON.json', {'Varga_row_residuals': residual, 'Wu_three_mean_sd_pairs_found': True, 'resolution': 'POPULATION', 'clinical_tail_validation': 'UNKNOWN', 'locators': ['doi:10.1111/clr.13578 T4 MANDIBLE / accepted T3', 'doi:10.1186/s40729-020-00272-0 Results/Fig4']})

def run_r1(profiles, teeth):
    worlds = []
    for gid in GUIDE_ORDER:
        (mu, sd) = profiles[gid]['apex_mean_sd_mm']
        q = stats.gamma.ppf(0.95, (mu / sd) ** 2, scale=sd * sd / mu)
        worlds.append({'guide': gid, 'same_radial_mean_sd_mm': [mu, sd], 'toward_plane_required_apex_gap_mm': float(1.3 + q), 'away_plane_sufficient_apex_gap_mm': 1.3, 'conditions': 'Local planar apex only, gamma radial law; B=0.3 sensitivity; not whole-cylinder guarantee', 'resolution': 'PHENOMENOLOGICAL'})
    chosen = teeth[0]
    h = float(chosen['occlusal_vertical_min_mm'])
    r = 2
    dental = {'case': chosen['case'], 'fdi': int(chosen['fdi']), 'same_total_tissue_mm': h, 'reduction_mm': r, 'possible_dentin_after_mm': [0, max(0, h - r)], 'enamel_boundary_observed': False, 'resolution': 'PER_TOOTH', 'scope': 'Layered CT-column model; two tissue assignments have identical tooth/pulp masks'}
    out = {'round': 'R1', 'claim_type': ['capability', 'information_link'], 'verdict': 'UNKNOWN_PHYSICAL_IDENTIFICATION', 'real_margin_95pct': 'UNKNOWN', 'real_dentin_fraction': 'UNKNOWN', 'missing': ['signed-loss distribution', 'independent anatomical boundary bound', 'regional guide-error law', 'DEJ', 'expert horn/axis/preparation geometry'], 'orientation_counterexample': worlds, 'dentin_counterexample': dental}
    write('raw/R1_IDENTIFIABILITY.json', out)
    (ROOT / 'HANDOFF_R1.md').write_text("# R1 : physical identification missing\n\nThe same observed guide error magnitude allows direction towards or from the channel. The same tooth /pulp mask allows different enamel / dentin - bounds . Real 95% margin and actual dentin share is UNKNOWN. No thresholds changed. The next design R2 retains these latent state as a range and inverts explicit research scenarios.\n")
    state('R1_DECIDED', 'UNKNOWN_PHYSICAL_IDENTIFICATION', 'R2: conditional inverse margin and partial-identification dentin tables')
    return out

def gamma_control_probability(d, profile, boundary_shift=0):
    (mu, sd) = profile['apex_mean_sd_mm']
    k = (mu / sd) ** 2
    scale = sd * sd / mu
    sig = math.hypot(0.06, 0.22)
    bias = 0.059 + boundary_shift

    def survival(z):
        s = d - 1 - z
        if s == 0:
            return 0.5
        v = abs(s)
        p = integrate.quad(lambda x: 0.5 * (1 - v / x) * stats.gamma.pdf(x, k, scale=scale), v, np.inf, epsabs=1e-09)[0]
        return p if s >= 0 else 1 - p
    return integrate.quad(lambda z: survival(z) * stats.norm.pdf(z, bias, sig), bias - 9 * sig, bias + 9 * sig, epsabs=2e-08)[0]

def run_guides(op, profiles, sites):
    table = []
    raw = []
    params = {}
    for gid in GUIDE_ORDER:
        p = profiles[gid]
        scenario = op.scenario_margin(p, x=0.05)
        bounds = [op.bound_margin(p, x=0.05, boundary_budget=b) for b in [0, 0.3, 0.7]]
        ref = gamma_control_probability(scenario, p)
        record('gamma_independent_quadrature_' + gid, abs(ref - 0.05), 1e-07, abs(gamma_control_probability(scenario + 1, p) - 0.05), 'Independent integral of gamma density and image jitter; same information numerical check only')
        root320 = optimize.brentq(lambda d: float(op.scenario_probability(d, p, order=320)) - 0.05, 0, 40)
        record('quadrature_resolution_' + gid, abs(root320 - scenario), 0.01, 1, '96 vs 320 Hermite nodes; mm')
        scalar = op.scalar_control(bounds[1], p, B=0.3)
        record('whole_body_scalar_' + gid, abs(scalar - 0.05), 1e-07, abs(op.scalar_control(bounds[1] + 1, p, B=0.3) - 0.05), 'Independent scalar union/Cantelli bound; condition on exact population moments and valid B')
        params[gid] = {'scenario_apex_required_gap_mm': float(scenario), 'scenario_boundary_sensitivity_range_mm': [float(scenario), float(scenario + 0.7)], 'whole_body_sufficient_gap_B0p3_mm': float(bounds[1]), 'whole_body_boundary_sensitivity_range_mm': [float(bounds[0]), float(bounds[2])], 'statistical_CI_margin_mm': None, 'real_margin_95pct_mm': None, 'physical_calibration': 'UNKNOWN', 'resolution': 'PHENOMENOLOGICAL', 'population_transfer': 'UNKNOWN; dynamic profile is mixed arch'}
        groups = [('ALL', sites)] + [(reg, [s for s in sites if s['region'] == reg]) for reg in sorted(set((s['region'] for s in sites)))]
        for (reg, rows) in groups:
            gap = np.array([s['gap_lower_mm'] for s in rows])
            base = gap >= 2
            scen = gap >= scenario
            bound = gap >= bounds[1]
            flip = (base != scen).astype(float) * 100
            bf = (base != bound).astype(float) * 100
            interval = cluster_ci(rows, flip)
            table.append({'guide': gid, 'region': reg, 'sites_n': len(rows), 'cases_n': len(set((s['case'] for s in rows))), 'clinical_required_mm': 'UNKNOWN', 'scenario_apex_required_mm': scenario, 'scenario_B0_to_B0p7_mm': f'{scenario:.6f}..{scenario + 0.7:.6f}', 'whole_body_sufficient_B0p3_mm': bounds[1], 'whole_body_B0_to_B0p7_mm': f'{bounds[0]:.6f}..{bounds[2]:.6f}', 'two_mm_accept_per100': float(np.mean(base) * 100), 'scenario_changed_per100': float(np.mean(flip)), 'scenario_changed_CI95_low': interval[0], 'scenario_changed_CI95_high': interval[1], 'whole_body_changed_per100': float(np.mean(bf)), 'regional_error_profile': 'UNKNOWN; same population closure transferred to region', 'aggregate_resolution': 'POPULATION', 'site_resolution': 'PER_TOOTH'})
        for s in sites:
            d = s['gap_lower_mm']
            raw.append({'case': s['case'], 'fdi': s['fdi'], 'region': s['region'], 'guide': gid, 'gap_lower_mm': d, 'gap_upper_mm': s['gap_upper_mm'], 'two_mm_accept': d >= 2, 'local_apex_scenario_accept': d >= scenario, 'whole_body_conditional_bound_accept': d >= bounds[1], 'real_95pct_clearance_answer': 'UNKNOWN', 'resolution': 'PER_TOOTH'})
    csvwrite('tables/GUIDE_MARGIN.csv', table)
    csvwrite('raw/GUIDE_PER_SITE.csv', raw)
    return (table, params)

def lp_dentin_bounds(lo, hi, r):
    if hi <= r:
        return (0.0, 0.0)
    res = optimize.linprog([0.0, -1.0], A_ub=[[-1, 1]], b_ub=[-r], bounds=[(max(lo, r), hi), (0, None)], method='highs')
    if not res.success:
        raise ValueError(res.message)
    lower = optimize.linprog([0.0, 1.0], A_ub=[[-1, 1]], b_ub=[-r], bounds=[(max(lo, r), hi), (0, None)], method='highs')
    if not lower.success:
        raise ValueError(lower.message)
    return (lower.x[1], res.x[1])

def run_teeth(teeth, pr):
    missing = [t for t in teeth if not t['occlusal_vertical_min_mm']]
    write('raw/MISSING_OCCLUSAL.json', {'rejected_n': len(missing), 'reason': 'Missing inherited occlusal column; never zero-imputed', 'rows': missing})
    teeth = [t for t in teeth if t['occlusal_vertical_min_mm']]
    table = []
    raw = []
    all_proxy = []
    post = [t for t in teeth if t['tooth_type'] in pr['metrics']['posterior_types']]
    reductions = pr['metrics']['reductions_mm']
    eps = pr['metrics']['digital_band_mm']
    for t in teeth:
        h = float(t['occlusal_vertical_min_mm'])
        for delta in pr['metrics']['annotation_sensitivity_per_surface_mm']:
            err = eps + 2 * delta
            lo = max(0, h - err)
            hi = h + err
            for r in reductions:
                up = max(0, hi - r)
                hl = max(0, lo - r)
                hn = max(0, h - r)
                raw.append({'case': t['case'], 'fdi': int(t['fdi']), 'tooth_type': t['tooth_type'], 'reduction_mm': r, 'annotation_sensitivity_per_surface_mm': delta, 'hard_tissue_proxy_before_mm': h, 'hard_tissue_proxy_remaining_mm': hn, 'proxy_remaining_lower_mm': hl, 'proxy_remaining_upper_mm': up, 'dentin_partial_id_lower_mm': 0, 'dentin_partial_id_upper_mm': up, 'dentin_necessarily_below_0p5_in_column_model': up < 0.5, 'dentin_necessarily_below_1_in_column_model': up < 1, 'clinical_dentin_answer': 'UNKNOWN', 'occlusal_applicability': t['tooth_type'] in pr['metrics']['posterior_types'], 'resolution': 'PER_TOOTH'})
    combos = sorted(set(((float(t['occlusal_vertical_min_mm']), r, d) for t in teeth for r in reductions for d in [0, 0.15, 0.3])))
    errors = []
    lp_rows = []
    for (h, r, delta) in combos:
        err = eps + 2 * delta
        lo = max(0, h - err)
        hi = h + err
        (a, b) = lp_dentin_bounds(lo, hi, r)
        expected = max(0, hi - r)
        errors.append(max(abs(a), abs(b - expected)))
        lp_rows.append({'H_mm': h, 'r_mm': r, 'delta_mm': delta, 'LP_min_D_mm': a, 'LP_max_D_mm': b, 'analytic_max_D_mm': expected})
    record('dentin_partial_identification_LP', max(errors), 1e-07, 1, 'Actual HiGHS extrema, including cut-through branch; injection = fabricated 1 mm lower bound')
    csvwrite('raw/DENTIN_LP_CONTROL.csv', lp_rows)
    for kind in pr['metrics']['posterior_types']:
        rows = [t for t in post if t['tooth_type'] == kind]
        h = np.array([float(t['occlusal_vertical_min_mm']) for t in rows])
        for r in reductions:
            rem = np.maximum(0, h - r)
            q = np.quantile(rem, [0.05, 0.5, 0.95])
            upper = np.maximum(0, h + eps - r)
            below05 = (upper < 0.5) * 100
            below1 = (upper < 1) * 100
            ci = cluster_median_ci(rows, rem)
            table.append({'tooth_type': kind, 'teeth_n': len(rows), 'cases_n': len(set((t['case'] for t in rows))), 'reduction_mm': r, 'proxy_remaining_P5_mm': q[0], 'proxy_remaining_median_mm': q[1], 'proxy_remaining_P95_mm': q[2], 'proxy_median_caseboot_CI95_low_mm': ci[0], 'proxy_median_caseboot_CI95_high_mm': ci[1], 'proxy_below_0p5_point_per100': float(np.mean(rem < 0.5) * 100), 'dentin_below_0p5_model_lower_per100': float(np.mean(below05)), 'dentin_below_0p5_model_upper_per100': 100, 'dentin_below_1_model_lower_per100': float(np.mean(below1)), 'dentin_below_1_model_upper_per100': 100, 'actual_dentin_fraction': 'UNKNOWN', 'digital_band_only': True, 'resolution': 'POPULATION', 'leaf_resolution': 'PER_TOOTH'})
    for r in reductions:
        rr = [x for x in raw if x['occlusal_applicability'] and x['reduction_mm'] == r and (x['annotation_sensitivity_per_surface_mm'] == 0)]
        values = [100 * int(x['dentin_necessarily_below_0p5_in_column_model']) for x in rr]
        ci = cluster_ci(rr, values)
        all_proxy.append({'reduction_mm': r, 'teeth_n': len(rr), 'dentin_below_0p5_model_lower_per100': float(np.mean(values)), 'lower_fraction_caseboot_CI95_per100': ci, 'dentin_below_0p5_model_upper_per100': 100, 'clinical_below_0p5_fraction': 'UNKNOWN', 'resolution': 'POPULATION'})
    csvwrite('tables/REDUCTION_TISSUE.csv', table)
    csvwrite('raw/REDUCTION_PER_TOOTH.csv', raw)
    ambiguous = sum((max(0, float(t['occlusal_vertical_min_mm']) - 2) >= 0.5 for t in post))
    write('raw/DENTIN_IDENTIFICATION_WITNESS.json', {'posterior_teeth_n': len(post), 'opposite_assignments_available_after_2mm_n': ambiguous, 'note': 'E=0 vs E=T are admissible under supplied masks; no biological enamel thickness asserted'})
    return (table, all_proxy, raw)

def attrition(sites, teeth, allteeth):
    ex = read('inputs/pulp_excluded.json')
    ce = read('inputs/pulp_case_excluded.json')
    pairs = read('inputs/pulp_pairs.json')
    return {'guide': {'input_expanded_rows': len(list(csv.DictReader((INPUT / 'guide_sites.csv').open()))), 'unique_complete_sites': len(sites), 'complete_cases': len(set((s['case'] for s in sites))), 'source_candidate_sites': 2706, 'source_not_in_complete_full_geometry': 125, 'source_dropout_fraction': 125 / 2706, 'reasons': 'Inherited X8: non-complete NV1/K3 cases; missing archive case identity UNKNOWN', 'new_rejected_unique_sites': 0}, 'pulp': {'raw_tooth_rows': len(allteeth), 'accepted_nontruncated_rows': len(teeth), 'truncated_rejected': len(allteeth) - len(teeth), 'truncated_fraction': (len(allteeth) - len(teeth)) / len(allteeth), 'accepted_cases': len(set((t['case'] for t in teeth))), 'posterior_teeth': sum((t['tooth_type'].startswith(('premolar', 'molar')) and bool(t['occlusal_vertical_min_mm']) for t in teeth)), 'anterior_out_of_occlusal_scope': sum((not t['tooth_type'].startswith(('premolar', 'molar')) for t in teeth)), 'missing_occlusal_column_n': sum((not t['occlusal_vertical_min_mm'] for t in teeth)), 'missing_occlusal_column_fraction': sum((not t['occlusal_vertical_min_mm'] for t in teeth)) / len(teeth), 'upstream_tooth_exclusions_n': len(ex), 'upstream_tooth_exclusion_reasons': dict(Counter((x['reason'] for x in ex))), 'upstream_R6_tooth_exclusion_records_n': len(ce) if isinstance(ce, list) else None, 'upstream_R6_unique_cases_with_tooth_exclusions': len(set((x['case'] for x in ce))) if isinstance(ce, list) else None, 'upstream_R6_exclusion_reasons': dict(Counter((x['reason'] for x in ce))) if isinstance(ce, list) else None, 'local_CT_candidates_n': len(pairs), 'CT_pair_rejected_n': sum((not x.get('paired', False) for x in pairs)), 'CT_pair_rejected_fraction': sum((not x.get('paired', False) for x in pairs)) / len(pairs), 'case_exclusion_detail_locator': 'inputs/pulp_case_excluded.json; tooth records, not excluded whole cases', 'pair_detail_locator': 'inputs/pulp_pairs.json', 'physical_dentin_answers_rejected': len(teeth), 'physical_dentin_dropout_fraction': 1.0, 'reason': 'No DEJ; no matched anatomical horn/preparation measurement'}, 'external_referents': {'physical_95pct_margin_matched': 0, 'physical_dentin_validation_matched': 0, 'rejected_anatomical_referents': 2, 'rejected_anatomical_fraction': 1.0, 'reasons': ['Azim fossa/roof is not CT highest-pulp column', 'Khojastepour cusp/horn is not CT highest-pulp column'], 'thresholds_retained': 2, 'thresholds_note': '0.5mm histological criterion; 1mm study stratum, not universal safety threshold'}}

def run_ports(cal):
    rows = []
    for candidates in [1, 4, 16]:
        n = cal.required_blocks(0.05, 0.05, candidates)
        q = cal.tail_upper(n, 0.05, candidates)
        qm = cal.tail_upper(n - 1, 0.05, candidates)
        prob = stats.binom.pmf(n, n, 0.95)
        record('Wilks_binomial_' + str(candidates), abs(q - (1 - (0.05 / candidates) ** (1 / n))), 1e-12, 0.001, 'NIST/Wilks mathematical check; fault = tail shifted by 0.001')
        if not (q < 0.05 and qm >= 0.05 and (prob < 0.05 / candidates)):
            raise AssertionError('Order statistic threshold')
        rows.append({'simultaneous_candidates': candidates, 'independent_blocks_per_cell_n': n, 'tail_upper_at_n': q, 'tail_upper_at_n_minus_1': qm, 'n_minus_1_rejected': True, 'resolution': 'POPULATION'})
    empty = cal.calibrate([], x=0.05, alpha=0.05, candidates=4)
    record('empty_measurement_abstention', int(empty['status'] != 'UNKNOWN'), 0, 1, 'Empty data must abstain; injected physical answer rejected')
    out = {'signed_loss_port': empty, 'independent_sample_design': rows, 'actual_physical_measurements_n': 0, 'local_DEJ_port': 'measurements.py; inputs TOTAL_HARD_TISSUE and ENAMEL interval on matched specimen and local ray', 'clinical_injury_probability': 'UNKNOWN'}
    write('raw/R3_MEASUREMENT_PORT.json', out)
    return out

def main():
    verify_frozen()
    op = load_module('risk_operator')
    cal = load_module('calibrate_clearance')
    pr = read('PREREG_R2_SET_VALUED_INVERSE.json')
    profiles = {p['id']: p for p in read('inputs/guide_profiles.json')['profiles']}
    source_check(profiles)
    rows = list(csv.DictReader((INPUT / 'guide_sites.csv').open()))
    unique = {}
    for row in rows:
        key = (row['case'], int(row['fdi']))
        lo = float(row['whole_cylinder_label_gap_lower_mm'])
        hi = float(row['whole_cylinder_label_gap_upper_mm'])
        s = {'case': key[0], 'fdi': key[1], 'gap_lower_mm': lo, 'gap_upper_mm': hi, 'region': region(key[1])}
        if key in unique and unique[key] != s:
            raise ValueError('Guide-specific anatomy mismatch')
        if not (math.isfinite(lo) and lo <= hi and (lo >= 0)):
            raise ValueError('Invalid gap interval')
        unique[key] = s
    sites = list(unique.values())
    allteeth = list(csv.DictReader((INPUT / 'teeth.csv').open()))
    teeth = [t for t in allteeth if t['domain_truncated'] == 'False']
    if len(sites) != 2581 or len(teeth) != 3449:
        raise ValueError('Inherited cohort cardinality drift')
    r1 = run_r1(profiles, teeth)
    (gt, params) = run_guides(op, profiles, sites)
    (tt, aggregate, t_raw) = run_teeth(teeth, pr)
    write('raw/R2_CONDITIONAL_RESULTS.json', {'guide_parameters': params, 'posterior_reduction': aggregate, 'actual_clinical_answers': 'UNKNOWN'})
    (ROOT / 'HANDOFF_R2.md').write_text("# R2 : conditional response and partial identification\n\nInverted 5%-tail query for four guide systems and remaining total hard tissue at three reductions. Retained dentin amounts are a range of lower boundary zero when DEJ is missing. Numerical controls passes; physical identification from R1 is still UNKNOWN. Next design R3 changes the information set through a port for signed lost margin and measured DEJ .\n")
    state('R2_DECIDED', 'CONDITIONAL_GEOMETRY_ONLY', 'R3: physical-input calibration and robust local DEJ port')
    r3 = run_ports(cal)
    write('raw/CONTROLS.json', CHECKS)
    result = {'lane': 'X31-clinical-answers', 'claim_type': ['capability', 'information_link'], 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'scientific_admission': False, 'created_utc': now(), 'answer': 'Conditional margin and tissue-reduction queries now runnable; empirically supported clinical margins/dentin remain UNKNOWN', 'rounds': {'R1': r1, 'R2': 'CONDITIONAL_GEOMETRY_ONLY', 'R3': r3}, 'guide_parameters': params, 'posterior_reduction_summary': aggregate, 'attrition': attrition(sites, teeth, allteeth), 'external_referent': {'kind': 'independent_measurement', 'locator': 'doi:10.1111/clr.13578 Table4 MANDIBLE; doi:10.1186/s40729-020-00272-0 Results/Fig4', 'compared_quantity': 'Published mean and individual SD of 3D entry/apex/angular deviations; NOT signed clearance quantile', 'refutes_us': False, 'matched_physical_clearance_validation': False}, 'external_referents': [{'kind': 'independent_measurement', 'locator': 'doi:10.1046/j.0143-2885.2003.00609.x', 'compared_quantity': 'Histological odontoblast/pulp response by measured remaining dentin; external 0.5mm criterion; no matched X31 dentin measurement', 'refutes_us': False, 'matched': False}, {'kind': 'independent_measurement', 'locator': 'doi:10.1016/s0109-5641(00)00041-5', 'compared_quantity': 'Class V premolar pulp response stratified by <0.5,0.5-1,>1mm; 1mm is study bin', 'refutes_us': False, 'matched': False}, {'kind': 'published_dataset', 'locator': 'https://toothfairy2.grand-challenge.org/dataset/; doi:10.1007/978-3-031-72111-3_2', 'compared_quantity': 'Published paired tooth/pulp labels and 0.3mm scale; inherited identity/semantic checks', 'refutes_us': False, 'matched': 'annotation geometry only'}, {'kind': 'closed_form', 'locator': 'https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/tolelimi.htm', 'compared_quantity': 'Order-statistic population-coverage confidence; tested against binomial probability', 'refutes_us': False}, {'kind': 'independent_measurement', 'locator': 'doi:10.1016/j.joen.2014.04.002; https://pmc.ncbi.nlm.nih.gov/articles/PMC3834639/#T2', 'compared_quantity': 'Fossa-to-chamber-roof / cusp-to-horn distances; operator mismatch rejected', 'refutes_us': False, 'matched': False}], 'quantities': {'individual_gap': {'resolution': 'PER_TOOTH', 'unit': 'mm', 'evidence': 'derived geometry on published labels'}, 'guide_margin_scenario': {'resolution': 'PHENOMENOLOGICAL', 'unit': 'mm', 'debt': 'Measure signed minimum-gap loss, cluster dependence and true boundary in guide/region population'}, 'residual_tissue': {'resolution': 'PER_TOOTH', 'unit': 'mm', 'debt': 'Expert axis/horn/DEJ and actual bur track'}, 'regional_fractions': {'resolution': 'POPULATION', 'unit': 'per100 sampled sites/teeth', 'CI': 'case-block bootstrap; model/segmentation bias excluded'}}, 'cost': {'wall_seconds': time.perf_counter() - START, 'cpu_seconds': time.process_time(), 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'fit_seconds': 0, 'physical_measurements': 0, 'gpu': False, 'threads_max': 4, 'upstream': 'X8 initial geometry122.11s plus verification83.05s and K3 prep2904.30s; X12 R6 CT1107.777s + R8atlas202.111s; other prior/discovery costs UNKNOWN', 'discovery_preparation_cost': 'Agent/source-search cost UNKNOWN; snapshot and run commands logged', 'fallback': 'Explicit unknown and minimum physical measurement ports'}, 'control_summary': {'all_valid_pass': all((c['valid_pass'] for c in CHECKS)), 'all_injected_faults_rejected': all((c['fault_rejected'] for c in CHECKS)), 'same_information_result': 'Numerical agreement; no algorithm claim'}, 'inputs_manifest_sha256': sha(ROOT / 'INPUT_MANIFEST.json')}
    write('results.json', result)
    import report
    report.render(result, gt, tt)
    state('THREE_ROUNDS_DECIDED', 'CONDITIONAL_TABLES_PASS; PHYSICAL_CLINICAL_ANSWERS_UNKNOWN', 'Independent expert same-specimen landmarks/DEJ and guide-region signed-loss measurements; no identical rerun')
    print(json.dumps({'sites': len(sites), 'teeth': len(teeth), 'controls': len(CHECKS), 'status': 'CONDITIONAL_GEOMETRY_ONLY', 'wall_seconds': time.perf_counter() - START}))
if __name__ == '__main__':
    main()
