from model import *

def run(rows):
    start = time.monotonic()
    elig = [r for r in rows if r['substrate'] == 'titanium' and r['height_mm'] is not None and r['independent_first_state']]
    m = fit(elig, 'reduced')
    loo = loso(elig, 'reduced')
    p = json.loads((ROOT / 'PREREG_RETENTION_R2.json').read_text())['gates']
    gate = 'UNKNOWN' if len({r['study'] for r in elig}) < p['min_studies'] or not loo['all_folds_fitted'] else 'PASS' if loo['log_rmse'] <= p['loso_log_rmse_max'] and loo['95pct_coverage'] >= p['loso_95pct_coverage_min'] else 'FAIL_TRANSFER'
    a = next((r for r in rows if r['study'] == 'PMC10154144' and r['group'] == 'H6-P0'))
    b = next((r for r in rows if r['study'] == 'PMC10154144' and r['group'] == 'H6-P4'))
    keys = ['height_mm', 'total_convergence_deg', 'cement', 'thermal_cycles', 'substrate', 'crown_material']
    sa = [a[k] for k in keys]
    sb = [b[k] for k in keys]
    identity = float(max((abs(float(a[k]) - float(b[k])) for k in ['height_mm', 'total_convergence_deg', 'thermal_cycles'])))
    lr = math.log(b['mean_force_N'] / a['mean_force_N'])
    se = math.sqrt((a['sd_force_N'] / a['mean_force_N']) ** 2 / a['n'] + (b['sd_force_N'] / b['mean_force_N']) ** 2 / b['n'])
    dof = (se ** 2) ** 2 / (((a['sd_force_N'] / a['mean_force_N']) ** 2 / a['n']) ** 2 / (a['n'] - 1) + ((b['sd_force_N'] / b['mean_force_N']) ** 2 / b['n']) ** 2 / (b['n'] - 1))
    from scipy.stats import t
    z = t.ppf(0.975, dof)
    pair = {'summary_keys': keys, 'state_A_summary': sa, 'state_B_summary': sb, 'summary_byte_identity': json.dumps(sa) == json.dumps(sb), 'numeric_identity_error': identity, 'state_A': a['row_id'], 'state_B': b['row_id'], 'changed_feature': [a['surface_feature'], b['surface_feature']], 'force_A_N': a['mean_force_N'], 'force_B_N': b['mean_force_N'], 'difference_N': b['mean_force_N'] - a['mean_force_N'], 'relative_difference': b['mean_force_N'] / a['mean_force_N'] - 1, 'ratio_95pct_CI_approx': [math.exp(lr - z * se), math.exp(lr + z * se)], 'sufficient': False, 'minimal_extension_in_this_pair': 'Dimple count (0 vs4), or the retained local interface-feature geometry; global force law still UNKNOWN.', 'external_referent': {'kind': 'independent_measurement', 'locator': 'https://doi.org/10.4047/jap.2023.15.2.63#T1', 'compared_quantity': 'H6-P0 vs H6-P4 pull-off force_N at identical recorded height/TOC/cement/thermal protocol', 'refutes_us': True}, 'resolution': 'POPULATION means of PER_TOOTH tests; feature lives at PER_SURFACE_REGION'}
    aging = []
    for a in rows:
        for b in rows:
            if a['study'] == b['study'] and a['thermal_cycles'] == 0 and (b['thermal_cycles'] is not None) and (b['thermal_cycles'] > 0) and all((a[k] == b[k] for k in ['height_mm', 'cement', 'crown_material', 'substrate', 'surface_feature', 'surface_treatment'])):
                aging.append([a['row_id'], b['row_id']])
    out = {'round': 'R2', 'claim_type': 'information_link', 'eligible_groups': len(elig), 'eligible_studies': len({r['study'] for r in elig}), 'excluded_groups': len(rows) - len(elig), 'model': m, 'loso': loo, 'gate': gate, 'sufficiency_test': pair, 'aging': {'matched_aged_unaged_pairs': aging, 'aged_unaged_ratio': 'UNKNOWN', 'reason': 'Image-only Table2 in PMC10315082 not locally available; all numeric primary force groups lack matched before/after protocol.'}, 'X1b_absolute_force': 'UNKNOWN_OUT_OF_SUBSTRATE_SUPPORT', 'wall_s': time.monotonic() - start, 'external_referent': pair['external_referent']}
    (ROOT / 'rounds/RESULTS_R2.json').write_text(json.dumps(out, indent=2) + '\n')
    return out
if __name__ == '__main__':
    o = run(json.loads((ROOT / 'raw/RETENTION_GROUPS.json').read_text()))
    print({k: o[k] for k in ['eligible_groups', 'eligible_studies', 'gate', 'model']})
    print('LOSO', o['loso']['log_rmse'], o['loso']['95pct_coverage'])
    print('pair', o['sufficiency_test']['relative_difference'])
