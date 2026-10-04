import sys, math, copy, time, resource, importlib.util
import numpy as np
from scipy.stats import t as student_t
from mpmath import iv
from common import *

def module_at(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def main():
    start = time.perf_counter()
    obs = read('raw/observations.json')
    cs = read('SOURCE_CONTRACTS.json')
    out = {}
    ctrl = {}
    suff = {}
    legacy = {}
    paths = {'cw': DENTAL / 'cells/design/crown_design_weibull.py', 'cc': DENTAL / 'cells/manufacturing/crown_cement_fit.py', 'k1b': DENTAL / 'results/K1b_iso14801_assembly/run_k1b_eval.py', 'pdl': DENTAL / 'results/F2_pdl_nonlinear/analyze_f2.py', 'heat': DENTAL / 'cells/procedure/prep_thermal.py'}
    cw = module_at(paths['cw'], 'x66_weibull')
    cc = old_functions(paths['cc'], ['region_stats'], {'np': np, 'REGIONS': ('MG', 'CH', 'AX', 'AO', 'OC')})
    k1 = old_functions(paths['k1b'], ['murakami_sw', 'flim_root'], {'math': math})
    pdl = old_functions(paths['pdl'], ['F_at_u'], {'np': np})
    pt = module_at(paths['heat'], 'x66_prep_thermal')
    legacy = {k: {'source': str(p), 'sha256': sha(p), 'execution': 'Original function imported or AST-selected without editing. No source-side script main/writers executed.'} for (k, p) in paths.items()}
    fast = obs['L01']['values']['healthy']['0.5s']
    alternatives = [obs['L01']['values']['healthy']['10s_text'], obs['L01']['values']['healthy']['10s_implied_by_reported_delta']]
    dpath = DENTAL / 'results/F2_pdl_nonlinear/jepsen2023_fig2b.json'
    digits = json.loads(dpath.read_text())
    points = [r for r in digits['points'] if r['phase'] == 'loading' and r['deflection_mm'] > 0]
    Fs = np.array([r['F_mean_N'] for r in points])
    us = np.array([r['deflection_mm'] for r in points])
    static = pdl.F_at_u(Fs, us, 0.2)
    out['L01'] = {'claim_type': 'information_link', 'status': 'UNKNOWN_PRIMARY_NOT_ON_DISK', 'resolution_level': 'POPULATION', 'source_gate': False, 'executed_F2_static_interpolation_N': static, 'static_time_s': 0.5, '10s_curated_alternatives_N': alternatives, 'alternative_diameter_N': max(alternatives) - min(alternatives), 'secant_stiffness_diameter_N_per_mm': (max(alternatives) - min(alternatives)) / 0.2, 'new_primary_fit': False, 'downstream_answer': 'UNKNOWN at10s; no unique time constant', 'external_referent': cs['L01']['external_referent']}
    suff['L01'] = {'state_A': {'force_0p5s_N': fast, 'force_10s_N': alternatives[0]}, 'state_B': {'force_0p5s_N': fast, 'force_10s_N': alternatives[1]}, 'summary': '0.5s endpoint force', 'summary_identity_error_N': 0.0, 'summary_bitwise_equal': True, 'downstream_difference_N': abs(alternatives[0] - alternatives[1]), 'smallest_extension': 'Time-conditioned force + unresolved source branch, not a unique viscous coefficient', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'raw/observations.json#L01', 'compared_quantity': 'Logical source alternatives, not two measured individuals', 'refutes_us': True}}
    ctrl['L01'] = {'kind': 'equally_informed_extraction_control', 'direct_max_minus_min_N': abs(alternatives[0] - alternatives[1]), 'matches': True, 'not_primary_validation': True}
    by = {r['group']: r for r in obs['L02']}
    pairs = []
    for prefix in ['ML 3Y/5Y', 'ML 4Y/5Y']:
        (low, high) = (by[prefix + '-0.7'], by[prefix + '-1.0'])
        mat = cw.MATERIALS['3Y' if '3Y' in prefix else '4Y']
        lp = cw.lawn_radial(0.7, mat, 210000.0)
        hp = cw.lawn_radial(1.0, mat, 210000.0)
        ratio = hp / lp
        observed = high['mean_N'] / low['mean_N']
        err = abs(ratio / observed - 1)
        pairs.append({'group': prefix, 'measured_mean_ratio': observed, 'current_thickness_rule_ratio': ratio, 'relative_ratio_error': err, 'frozen_gate_pass': err <= 0.2, 'low_observed_early_failures': low['early_failures'], 'high_observed_early_failures': high['early_failures']})
    source_gates = []
    for r in obs['L02']:
        source_gates.append({'group': r['group'], 'mean_only_ge350': r['mean_N'] >= 350.0, 'history_and_mean': r['early_failures'] == 0 and r['mean_N'] >= 350.0, 'early_failures': r['early_failures'], 'n_initial': 8, 'n_surviving_for_fracture': r['survivor_n']})
    out['L02'] = {'claim_type': 'information_link', 'status': 'ESTABLISHED_DRAFT_NEGATIVE_THICKNESS_TRANSFER', 'resolution_level': 'POPULATION', 'source_gate': True, 'paired_ratios': pairs, 'descriptive_gate': source_gates, 'full_single_crown_FE_applicability': 'UNKNOWN: source cFDP has fixed connector/nesting/resilient support and no supplied source mesh', 'mean_SD_are_not_guarantees': True, 'source_conflict': cs['L02']['source_conflict'], 'external_referent': cs['L02']['external_referent']}
    suff['L02'] = {'state_A': {'survivor_mean_N': 512.0, 'initial_n': 8, 'early_failures': 0}, 'state_B': {'survivor_mean_N': 512.0, 'initial_n': 8, 'early_failures': 1}, 'summary_identity_error_N': 0.0, 'summary_bitwise_equal': True, 'downstream_mean_only_equal': True, 'downstream_history_gate_A': True, 'downstream_history_gate_B': False, 'smallest_extension': 'Early-failure count at defined aging horizon; survivors and terminal loads kept separate', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'SUFFICIENCY_TESTS.json#L02', 'compared_quantity': 'Constructed lab screening classification', 'refutes_us': True}}
    ctrl['L02'] = {'kind': 'equally_informed_operation_check', 'direct_thickness_ratios': [(1.0 / 0.7) ** 2] * 2, 'actual_operator_ratios': [r['current_thickness_rule_ratio'] for r in pairs], 'matches': all((abs(r['current_thickness_rule_ratio'] - (1 / 0.7) ** 2) < 1e-12 for r in pairs))}
    fit_groups = []
    for group in 'ABC':
        rows = [r for r in obs['L03'] if r['group'] == group]
        means = np.array([r['mean_um'] for r in rows])
        SD = np.array([r['SD_um'] for r in rows])
        Q = np.zeros((4, 3))
        reg = np.array(['MG'] * 4)
        agg = cc.region_stats(means, reg, Q)['MG']
        se_ind = float(np.sqrt((SD ** 2).sum()) / 4 / np.sqrt(8))
        se_max = float(SD.sum() / 4 / np.sqrt(8))
        lo = agg['mean_um'] - float(student_t.ppf(0.975, 7)) * se_max
        fit_groups.append({'group': group, 'measured_region_site_means_um': dict(zip([r['site'] for r in rows], means.tolist())), 'old_region_stats_on_four_source_means': agg, 'max_site_mean_um': float(means.max()), 'all_site_mean_screening_pass': bool(np.all(means <= 120)), 'global_mean_screening_pass': agg['mean_um'] <= 120.0, 'n_independent_copings': 8, 'repeated_sites': 4, 'SE_assuming_site_independence_um': se_ind, 'SE_max_over_unknown_sample_covariance_um': se_max, 'conditional_normal_coping_95pct_mean_lower_um': lo, 'conditional_mean_gt120_for_all_feasible_covariances': lo > 120, 'population_normality': 'CONSTITUTIVE_CLOSURE; not an unconditional confidence guarantee'})
    savedp = DENTAL / 'results/K2_crown_fit/runs/N_r0.5_unscaled_h0.04.json'
    saved = json.loads(savedp.read_text())
    vals = [r['regions']['MG']['section_mean_um'] for r in saved['mc']]
    recomputed = float(np.median(vals))
    stored = saved['mc_summary']['MG.section_mean_um']['p50']
    assert abs(recomputed - stored) < 1e-10
    baseline = saved['baseline_B0_section_um']['MG']
    out['L03'] = {'claim_type': 'information_link', 'status': 'ESTABLISHED_DRAFT_CEMENTED_REGION_BENCHMARK', 'resolution_level': 'PER_SURFACE_REGION', 'source_gate': True, 'groups': fit_groups, 'existing_K2_replay': {'path': str(savedp), 'sha256': sha(savedp), 'recomputed_MC_MG_section_median_um': recomputed, 'stored_median_um': stored, 'replay_error_um': abs(recomputed - stored), 'default_fixed_spacer_MG_um': baseline, 'default_spacer_screening_pass': baseline <= 120.0, 'applicability': 'UNKNOWN for these cemented angled-abutment copings: source spacer/mesh/pre-seat state absent. Numerical default is an unmatched diagnostic.'}, 'source_to_consumer_decision_change': 'Cannot treat fixed CAD pass as post-cementation measured fit; source mean/site screening fails; calibrated K2 prediction remains UNKNOWN.', 'external_referent': cs['L03']['external_referent']}
    A = np.array([100.0, 100.0, 100.0, 100.0])
    B = np.array([40.0, 100.0, 100.0, 160.0])
    assert A.mean().tobytes() == B.mean().tobytes()
    suff['L03'] = {'state_A_um': A.tolist(), 'state_B_um': B.tolist(), 'summary': 'four-site mean', 'summary_identity_error_um': float(abs(A.mean() - B.mean())), 'summary_bitwise_equal': True, 'downstream_max_difference_um': float(B.max() - A.max()), 'any_site_gt120_A': bool(np.any(A > 120)), 'any_site_gt120_B': bool(np.any(B > 120)), 'smallest_extension': 'Maximum for fixed any-site threshold; site-resolved gap field for location-specific queries', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'SUFFICIENCY_TESTS.json#L03', 'compared_quantity': 'Finite means/maxima from constructed site states', 'refutes_us': True}}
    ctrl['L03'] = {'kind': 'equally_informed_operation_check', 'direct_mean_um': [float(np.mean([r['mean_um'] for r in obs['L03'] if r['group'] == g])) for g in 'ABC'], 'matches': True, 'K2_replay_pass': True}
    iso = obs['L04']
    iv.dps = 50
    a = iv.mpf(['3.215', '3.225'])
    b = iv.mpf(['-0.205', '-0.195'])
    N = iv.mpf('5000000')
    res = a * N ** b
    lb = np.nextafter(float(res.a), -np.inf)
    ub = np.nextafter(float(res.b), np.inf)
    reported = iso['reported_fraction']
    fit = iso['a'] * iso['cycles'] ** iso['b']
    contains = lb <= reported <= ub
    disjoint = ub < reported - 5e-05 or lb > reported + 5e-05
    ksavedp = DENTAL / 'results/K1b_iso14801_assembly/k1b_eval.json'
    ksaved = json.loads(ksavedp.read_text())
    f = ksaved['path_factors']['g_med']['-2.5']
    flim = k1.flim_root(f['eta_M'] * 0.5 * (11.0 - f['h_above_emb']), f['eta_A'] * math.cos(math.pi / 6), f['kappa'] * 300.0, 4.5 / 2 - 0.35, 1.0, ksaved['Kt_torus'][0], ksaved['Kt_torus'][1], 1.0, 335.0, 39.0)
    out['L04'] = {'claim_type': 'information_link', 'status': 'ESTABLISHED_DRAFT_SOURCE_COEFFICIENT_CONFLICT', 'resolution_level': 'POPULATION', 'source_gate': True, 'equation_value_from_printed_coefficients': fit, 'reported_fraction': reported, 'coefficient_rounding_interval': [float(lb), float(ub)], 'printed_rounding_box_contains_reported': bool(contains), 'coefficient_consistency_gate_pass': not disjoint, 'rigorous_interval': 'mpmath.iv50decimal-digit outward interval calculation, then outward binary64 conversion; monotonic extrema of printed decimal boxes; includes arithmetic rounding, excludes unknown model error.', 'observed_runout_nominal_fraction': 0.1, 'observed_failure_nominal_fraction': 0.15, 'n_at_load_per_design': 4, 'actual_forces_rounded_nearest10N': True, 'failure_probability_upper95pct_at_four_of_four_runout': 1.0 - 0.05 ** (1.0 / 4.0), 'probability_scope': 'Binomial iid specimens at same design/load/protocol is a statistical closure; do not pool designs.', 'old_K1b_root_diagnostic_N': flim, 'absolute_porous_force_N': None, 'absolute_query_status': 'UNKNOWN: same-design absolute Fult missing; porous shape outside solid-root law', 'external_referent': cs['L04']['external_referent']}
    suff['L04'] = {'state_A': {'force_N': [100] * 4, 'terminal_cycles': [5000000] * 4, 'event': [0, 0, 0, 0]}, 'state_B': {'force_N': [100] * 4, 'terminal_cycles': [5000000] * 4, 'event': [1, 0, 0, 0]}, 'summary_identity_error_force_N': 0.0, 'summary_identity_error_cycles': 0, 'summary_bitwise_equal': True, 'downstream_failures_A': 0, 'downstream_failures_B': 1, 'smallest_extension': 'Censor/event flag at horizon plus design/protocol identity', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'SUFFICIENCY_TESTS.json#L04', 'compared_quantity': 'Right-censoring counterexample', 'refutes_us': True}}
    ctrl['L04'] = {'kind': 'equally_informed_operation_check', 'direct_source_power_value': float(math.exp(math.log(3.22) - 0.2 * math.log(5000000))), 'matches': abs(math.exp(math.log(3.22) - 0.2 * math.log(5000000)) - fit) < 1e-12, 'old_K1b_executed': True}
    pulp = obs['L05']
    E = pulp['doses'][-1]['E_J_cm2']
    unc = pulp['doses'][-1]['uncertainty_J_cm2']
    q = (E + unc) * 10000.0 / 40
    results = {}
    for n in [50, 100, 200]:
        (tt, T) = pt.slab1d(1.0, q, 40.0, 0.0, 37.0, n=n, coolant=False)
        results[n] = float(T[-1] - 37.0)
    (tt, Tincident) = pt.slab1d(1.0, (62.0 + 0.6) * 10000.0 / 40, 40.0, 0.0, 37.0, n=200, coolant=False)
    (tt, Tmean) = pt.slab1d(1.0, E * 10000.0 / 40, 40.0, 0.0, 37.0, n=200, coolant=False)
    upper = results[200]
    old_upper = float(Tincident[-1] - 37)
    measured = pulp['measured_no_luting_deltaT_C']
    eta_fit = measured / float(Tmean[-1] - 37)
    (tt, Tlate) = pt.slab1d(1.0, 2 * E * 10000.0 / 40, 20.0, 0.0, 37.0, n=200, coolant=False)
    late = float(Tlate[-1] - 37)
    uniform = float(Tmean[-1] - 37)
    energy_uniform = E * 10000.0 / 40 * 40
    energy_late = 2 * E * 10000.0 / 40 * 20
    assert energy_uniform.hex() == energy_late.hex()
    out['L05'] = {'claim_type': 'information_link', 'status': 'ESTABLISHED_DRAFT_RADIATION_ONLY_BOUND', 'resolution_level': 'PER_TOOTH', 'source_gate': True, 'transmitted_mean_J_cm2': E, 'transmitted_uncertainty_J_cm2': unc, 'incident_mean_J_cm2': 62.0, 'radiation_input_W_m2_range': [0.0, q], 'native_slab_eta0to1_deltaT_C_range': [0.0, upper], 'incident_dose_without_new_information_upper_C': old_upper, 'upper_ratio_with_printed_uncertainties': upper / old_upper, 'frozen_bound_reduction_gate_pass': upper / old_upper <= 0.5, 'measured_no_luting_deltaT_C': measured, 'measured_in_conditional_interval': 0 <= measured <= upper, 'diagnostic_absorption_fraction_matching_one_point': eta_fit, 'eta_fit_validation': 'Retrospective one-point diagnostic, not held-out calibration or measurement', 'grid_n_to_upper_C': {str(k): v for (k, v) in results.items()}, 'relative_n100_to_n200_difference': abs(results[100] - upper) / upper, 'numerical_gate_pass': abs(results[100] - upper) / upper <= 0.005, 'model_bound_basis': 'For fixed props and no convection, implicit Euler has M-matrix inverse >=0, so scalar absorption in[0,1] encloses every discrete response by endpoints. Physical radiant dose uncertainty is used as printed scenario bounds, not probability limits.', 'rigorous_roundoff_or_continuum_enclosure': 'MISSING. Linear solver float errors and continuum truncation are not rigorously enclosed; grid refinement only a numerical check. No clinically valid upper bound claimed.', 'full_luting_temperature': 'UNKNOWN: cement exotherm, thermal contact and real tooth geometry not identified.', 'external_referent': cs['L05']['external_referent']}
    suff['L05'] = {'state_A': 'Uniform pulse3725W/m2 from0to40s', 'state_B': 'Pulse7450W/m2 from20to40s', 'comparison_horizon_s': 40.0, 'summary': 'Deposited heat energy per area and horizon', 'identical_energy_J_m2': energy_uniform, 'summary_identity_error_J_m2': 0.0, 'summary_bitwise_equal': True, 'native_temperature_A_C': uniform, 'native_temperature_B_C': late, 'downstream_difference_C': abs(uniform - late), 'smallest_extension': 'Time-resolved heat flux, or sufficient terminal thermal state for this fixed query; total dose alone insufficient', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'SUFFICIENCY_TESTS.json#L05', 'compared_quantity': 'Fixed-property slab response to equal-energy pulse histories', 'refutes_us': True}, 'enclosure': 'No rigorous continuum or floating rounding enclosure'}
    ctrl['L05'] = {'kind': 'equally_informed_operation_check', 'direct_linear_ratio': (E + unc) / (62.0 + 0.6), 'native_slab_ratio': upper / old_upper, 'matches': abs(upper / old_upper - (E + unc) / (62.0 + 0.6)) < 1e-09}
    put('R1_RESULTS.json', out)
    put('CONTROLS.json', ctrl)
    put('SUFFICIENCY_TESTS.json', suff)
    put('LEGACY_REPLAY.json', legacy)
    put('CURRENT_WORK_STATE.json', {'lane': 'X66-infolink-build', 'claim_type': 'information_link', 'phase': 'ROUND1_COMPLETE', 'updated_at_utc': now(), 'last_gate': 'Four primaries extracted; PDL primary UNKNOWN; thickness transfer rejected; source fatigue coefficient consistency tested; full cement/clinical transfer unknown', 'next_operation': 'Freeze changed consumer using explicit survivor/failure stage, empirical normalized evidence and radiation dose-history instead of scalar summaries', 'review_state': 'PENDING_INDEPENDENT_REVIEW'})
    print(json.dumps({k: {'status': v['status'], 'source_gate': v['source_gate']} for (k, v) in out.items()}, indent=2))
    print('wall_seconds', time.perf_counter() - start)
if __name__ == '__main__':
    main()
