import json, time, resource
from decimal import Decimal, localcontext
from pathlib import Path
import numpy as np
from r1 import assert_freeze, dump

def decimal_pathway(s):
    with localcontext() as ctx:
        ctx.prec = 50
        n = s['N_completed']
        t = s['time_first_visit']
        p = (Decimal(n['P']) * Decimal(str(t['P']['mean'])) + Decimal(n['fallback_RCT']) * Decimal(str(t['fallback_RCT']['mean']))) / Decimal(s['N_assigned']['P_plan'])
        return float(p)

def main():
    assert_freeze('R5')
    start = time.perf_counter()
    s = json.loads(Path('raw/R5_TRIAL_INPUTS.json').read_text())
    n = s['N_completed']
    t = s['time_first_visit']
    (np_, nf) = (n['P'], n['fallback_RCT'])
    na = s['N_assigned']['P_plan']
    mu = (np_ * t['P']['mean'] + nf * t['fallback_RCT']['mean']) / na
    radius = (np_ * t['P']['mean_half_unit'] + nf * t['fallback_RCT']['mean_half_unit']) / na
    direct = t['direct_RCT']['mean']
    dr = t['direct_RCT']['mean_half_unit']
    saving = [direct - dr - (mu + radius), direct + dr - (mu - radius)]
    ss = sum(((n[key] - 1) * t[key]['sd'] ** 2 + n[key] * (t[key]['mean'] - mu) ** 2 for key in ['P', 'fallback_RCT']))
    combined_sd = float(np.sqrt(ss / (na - 1)))
    curve_checks = []
    grid = np.array(s['curve_grid_days'], float)
    weights = np.r_[np.diff(grid)[0] / 2, (np.diff(grid)[:-1] + np.diff(grid)[1:]) / 2, np.diff(grid)[-1] / 2]
    for (key, vals) in s['VAS_curves'].items():
        area = float(weights @ np.array(vals))
        rad = float(weights @ np.array(s['curve_half_units'][key]))
        reported = s['reported_AUC'][key]
        rr = s['AUC_half_units'][key]
        curve_checks.append({'group': key, 'trapezoid_from_reported_mean_curve': area, 'curve_rounding_interval': [area - rad, area + rad], 'reported_AUC_rounding_interval': [reported - rr, reported + rr], 'complete_panel_compatible': abs(area - reported) <= rad + rr})
    implied_complete_panel_AUC = (np_ * s['reported_AUC']['P'] + nf * s['reported_AUC']['fallback_RCT']) / na
    result = {'round': 'R5', 'claim_type': 'information_link', 'fallback_fraction': nf / na, 'assigned_path_first_visit_mean_conditional_min': mu, 'source_rounding_interval_min': [mu - radius, mu + radius], 'assigned_path_first_visit_sd_conditional_min': combined_sd, 'first_visit_saving_conditional_min': saving, 'completed_procedure_saving_min': direct - t['P']['mean'], 'omitted_fallback_time_penalty_min': mu - t['P']['mean'], 'Decimal_control_error': abs(mu - decimal_pathway(s)), 'source_curve_checks': curve_checks, 'complete_panel_assigned_P_AUC_hypothetical': implied_complete_panel_AUC, 'pain_effect_status': 'UNKNOWN: complete-panel reconstruction fails for source RCT AUC and per-time sample denominators absent. Hypothetical AUC is not an empirical arm effect.', 'causal_ITT_status': 'UNKNOWN: two excluded control-arm patients and outcome denominators unresolved. Assignment-path descriptive first-time result conditional on timing N.', 'full_treatment_time_status': 'UNKNOWN: second RCT appointment time absent.', 'gates': {'count_mass': 'PASS' if np_ + nf == na and n['direct_RCT'] + n['excluded_RCT'] == s['N_assigned']['RCT'] else 'FAIL', 'source_rounding_first_visit_saving_ge_10min': 'CONDITIONAL_PASS' if min(saving) >= 10 else 'FAIL', 'source_complete_panel_AUC': 'PASS' if all((r['complete_panel_compatible'] for r in curve_checks)) else 'FAIL', 'causal_full_course_choice': 'UNKNOWN'}, 'cost': {'seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'external_subject_measurements': 'Published data, no new measurements', 'prior_all_rounds_inherited': True}}
    dump('raw/R5_RESULTS.json', result)
    dump('CURRENT_WORK_STATE.json', {'tag': 'XBREAK-hunt-4', 'stage': 'R5_COMPLETE', 'latest_gate': result['gates'], 'next_operation': 'Freeze R6: target preoperative risk of intraoperative RCT fallback using published CBCT x final-state counts. Preserve missing CBCT and leave-one-out evaluation.'})
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
