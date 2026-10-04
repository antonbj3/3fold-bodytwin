"""R2: identified mean contrasts, sample agreement bounds, refused target transport."""
from collections import Counter
from decimal import Decimal
from pathlib import Path
import csv, hashlib, json, time
import numpy as np
from extract_sources import ROOT, save, extract_cunali, x13_rows
D = Decimal

def contrast(r):
    (rm, cm, rs, cs) = [D(r[k]) for k in ['replica_mean_um', 'ct_mean_um', 'replica_sd_um', 'ct_sd_um']]
    bias = rm - cm
    sd_lo = max(D(0), abs(rs - cs) - D('.010'))
    sd_hi = rs + cs + D('.010')
    return dict(**r, bias_replica_minus_ct_um=str(bias), bias_rounding_interval_um=[str(bias - D('.010')), str(bias + D('.010'))], sample_difference_sd_identified_set_um=[str(abs(rs - cs)), str(rs + cs)], sample_difference_sd_rounding_enclosure_um=[str(sd_lo), str(sd_hi)], nominal_normal_loa_lower_endpoint_identified_set_um=[str(bias - D('1.96') * (rs + cs)), str(bias - D('1.96') * abs(rs - cs))], nominal_normal_loa_upper_endpoint_identified_set_um=[str(bias + D('1.96') * abs(rs - cs)), str(bias + D('1.96') * (rs + cs))], estimated_actual_loa=None, agreement_status='UNKNOWN_WITHIN_REGION_COVARIANCE_AND_DIFFERENCE_DISTRIBUTION', external_referent=dict(kind='independent_measurement', locator=r['source_locator'], compared_quantity='same-cohort regional method means / marginal SDs', refutes_us=True))

def ecological_check(rows):
    refs = {'Amann': dict(slope=0.4736, intercept=36.206, r2=0.97432, pearson=0.98708), 'Dentsply_Sirona': dict(slope=0.4695, intercept=44.553, r2=0.86219, pearson=0.92854)}
    out = []
    for (sys, ref) in refs.items():
        rr = [r for r in rows if r['system'] == sys]
        x = np.array([float(r['replica_mean_um']) for r in rr])
        y = np.array([float(r['ct_mean_um']) for r in rr])
        coeff = np.linalg.lstsq(np.column_stack([x, np.ones(4)]), y, rcond=None)[0]
        slope = float(np.sum((x - x.mean()) * (y - y.mean())) / np.sum((x - x.mean()) ** 2))
        intercept = float(y.mean() - slope * x.mean())
        rho = float(np.corrcoef(x, y)[0, 1])
        r2 = rho * rho
        errors = dict(slope=abs(slope - ref['slope']), intercept=abs(intercept - ref['intercept']), r2=abs(r2 - ref['r2']))
        gate = errors['slope'] <= 0.0001 and errors['intercept'] <= 0.005 and (errors['r2'] <= 5e-05)
        out.append(dict(system=sys, n_regional_means=4, resolution='PER_SURFACE_REGION', computed_slope=slope, computed_intercept_um=intercept, computed_r2=r2, computed_pearson=rho, published=ref, errors=errors, external_coefficient_gate='PASS' if gate else 'FAIL', independent_OLS_max_coefficient_error=float(np.max(np.abs(coeff - [slope, intercept]))), source_locator='https://doi.org/10.1590/0103-6440201601531 Figure4 p470', within_region_specimen_correlation=None, interpretation='Figure coefficients reproduced using four regional means, not specimen pairs. No affine prediction bound for individual crowns or target labs is available.'))
    return out

def witness():
    r = np.array([64.0, 64.0, 128.0, 128.0])
    c = np.array([72.0, 72.0, 136.0, 136.0])
    states = [c, c[::-1].copy()]
    summary = lambda cc: np.array([np.mean(r), np.std(r, ddof=1), np.mean(cc), np.std(cc, ddof=1)])
    s = [summary(cc) for cc in states]
    assert np.array_equal(s[0], s[1])
    vals = []
    for cc in states:
        d = r - cc
        sd = float(np.std(d, ddof=1))
        bias = float(np.mean(d))
        pred = r - bias
        vals.append(dict(replica_um=r.tolist(), ct_um=cc.tolist(), bias_um=bias, paired_sd_um=sd, nominal_normal_loa_um=[bias - 1.96 * sd, bias + 1.96 * sd], incorrect_120um_classifications=int(np.sum((pred <= 120) != (cc <= 120)))))
    return dict(kind='our_own_fixture', role='exact constructive identifiability proof, no empirical calibration', resolution='PER_POINT', summary=['replica mean', 'replica sample SD', 'CT mean', 'CT sample SD'], identical_marginal_value_multisets=True, summary_values_um=s[0].tolist(), identity_max_abs_error_um=float(np.max(np.abs(s[0] - s[1]))), states=vals, downstream_paired_sd_difference_um=vals[1]['paired_sd_um'] - vals[0]['paired_sd_um'], downstream_nominal_loa_width_difference_um=3.92 * (vals[1]['paired_sd_um'] - vals[0]['paired_sd_um']), downstream_classification_error_difference=vals[1]['incorrect_120um_classifications'] - vals[0]['incorrect_120um_classifications'], minimum_extension_for_nominal_normal_loa='Within-region sample cross-covariance, equivalently sample SD of the paired differences. Normality is still a closure.', minimum_extension_for_actual_specimen_threshold_error='Actual paired differences with specimen identity, region and state; tails/exceedance counts cannot follow from a covariance alone.')

def interval_class(lo, hi, lower, upper):
    if lo >= lower and hi <= upper:
        return 'WITHIN'
    if hi < lower or lo > upper:
        return 'OUTSIDE'
    return 'UNKNOWN_INTERVAL_CROSSES_BOUNDARY'

def query_x13(rows, donors, persist=True):
    out = []
    for r in rows:
        reasons = []
        if r['method'] != 'silicone_replica':
            reasons.append('different_measurement_operator')
        if r['material'] != 'zirconia':
            reasons.append('different_material')
        if r['restoration'] not in ['crown', 'coping']:
            reasons.append('different_restoration_geometry')
        if r['state'] != 'replica':
            reasons.append('different_observation_state')
        if r['region'] not in ['marginal', 'axial', 'occlusal', 'transition']:
            reasons.append('different_or_unresolved_quantity')
        q = dict(row_id=r['row_id'], study=r['study'], arm=r['arm'], region=r['region'], material=r['material'], method=r['method'], resolution='PER_SURFACE_REGION', reported_mean_um=r['measured_mean_um'], source_locator=('https://doi.org/' + r['doi'] if r.get('doi') else 'https://pmc.ncbi.nlm.nih.gov/articles/' + r['study'] + '/') + '#' + r['source_table'], source_cell=dict(row=r['source_row'], column=r['source_column'], text=r['source_cell']), authorized_correction=False, authorized_decision='UNKNOWN_NO_TARGET_PAIRED_VALIDATION', authorized_rejection_reason='Missing matching target lab/die/force/state registration and paired validation', conditional_scope_rejections=reasons, conditional_scenarios=[])
        if not reasons:
            ds = [d for d in donors if d['region'] == r['region']]
            for d in ds:
                x = D(str(r['measured_mean_um']))
                bias = D(d['bias_replica_minus_ct_um'])
                corrected = x - bias
                lo = corrected - D('.015')
                hi = corrected + D('.015')
                threshold_cases = []
                if r['region'] == 'marginal':
                    for threshold in [100, 120]:
                        original = interval_class(x - D('.005'), x + D('.005'), D(0), D(threshold))
                        decision = interval_class(lo, hi, D(0), D(threshold))
                        threshold_cases.append(dict(upper_um=threshold, lower_um=0, threshold_kind='HISTORICAL_RESEARCH_COMPARATOR', threshold_locator='DOI10.1590/0103-6440201601531 Introduction and Discussion p471 refs7/8', raw_decision=original, conditional_decision=decision, flip=original != decision and (not original.startswith('UNKNOWN')) and (not decision.startswith('UNKNOWN'))))
                elif r['region'] in ['axial', 'occlusal']:
                    original = interval_class(x - D('.005'), x + D('.005'), D(50), D(100))
                    decision = interval_class(lo, hi, D(50), D(100))
                    threshold_cases.append(dict(upper_um=100, lower_um=50, threshold_kind='X13_ENGINEERING_SCENARIO', threshold_resolution='PHENOMENOLOGICAL', threshold_locator=str(ROOT.parent / 'LANE_X13_CEMENT_GAP/PREREG_R1.json'), raw_decision=original, conditional_decision=decision, flip=original != decision and (not original.startswith('UNKNOWN')) and (not decision.startswith('UNKNOWN'))))
                q['conditional_scenarios'].append(dict(donor_system=d['system'], donor_source=d['source_locator'], assumed_bias_um=str(bias), corrected_mean_um=str(corrected), source_rounding_only_interval_um=[str(lo), str(hi)], target_prediction_interval_um=None, scenario_only='Assumes source region mean offset transports unchanged; this is not empirically established and there is no rigorous target error enclosure.', decisions=threshold_cases))
        out.append(q)
    if persist:
        save('X13_DECISIONS_R2.json', out)
        with (ROOT / 'X13_DECISIONS.csv').open('w') as f:
            fields = ['row_id', 'study', 'region', 'material', 'method', 'reported_mean_um', 'authorized_correction', 'conditional_scope_rejections']
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for q in out:
                writer.writerow({k: q[k] for k in fields})
    eligible = [q for q in out if not q['conditional_scope_rejections']]
    summaries = []
    for sys in ['Amann', 'Dentsply_Sirona']:
        for (reg, upper, lower) in [('marginal', 100, 0), ('marginal', 120, 0), ('axial', 100, 50), ('occlusal', 100, 50)]:
            cases = []
            for q in eligible:
                if q['region'] != reg:
                    continue
                sc = next((s for s in q['conditional_scenarios'] if s['donor_system'] == sys))
                d = next((d for d in sc['decisions'] if d['upper_um'] == upper and d['lower_um'] == lower))
                cases.append((q, sc, d))
            summaries.append(dict(donor_system=sys, region=reg, lower_um=lower, upper_um=upper, resolution='POPULATION', eligible_cells=len(cases), conditional_flips=sum((d['flip'] for (q, s, d) in cases)), confirmed_flips=0, flipped_cells=[dict(row_id=q['row_id'], raw_mean_um=q['reported_mean_um'], scenario_mean_um=s['corrected_mean_um'], raw=d['raw_decision'], scenario=d['conditional_decision'], locator=q['source_locator']) for (q, s, d) in cases if d['flip']]))
    return dict(total_cells=len(rows), resolution='POPULATION', authorized_cells=0, confirmed_flips=0, actual_physical_changed_decisions=None, actual_physical_changed_decisions_status='UNKNOWN_NO_TARGET_PAIRED_VALIDATION', authorized_correction_unknown_cells=len(rows), conditional_eligible_cells=len(eligible), conditional_rejected_cells=len(out) - len(eligible), conditional_rejection_fraction=(len(out) - len(eligible)) / len(out), rejection_reasons_multilabel=dict(Counter((z for q in out for z in q['conditional_scope_rejections']))), conditional_scenarios=summaries)

def run():
    started = time.perf_counter()
    rows = extract_cunali()
    donors = [contrast(r) for r in rows]
    checks = ecological_check(rows)
    suff = witness()
    targets = query_x13(x13_rows(), donors)
    save('REGIONAL_METHOD_CONTRASTS.json', donors)
    universal_bias_intersection = [max((float(d['bias_rounding_interval_um'][0]) for d in donors)), min((float(d['bias_rounding_interval_um'][1]) for d in donors))]
    result = dict(round='R2', claim_type='information_link', outcome='REGIONAL_CONTRAST_IDENTIFIED_SPECIMEN_AGREEMENT_AND_TARGET_TRANSFER_UNKNOWN', external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.1590/0103-6440201601531 Table1 and Figure4', compared_quantity='same-cohort regional method contrast; four-region OLS diagnostic', refutes_us=True), source_regions=len(rows), source_resolution='PER_SURFACE_REGION', source_n_per_system=10, shared_mean_bias_possible=universal_bias_intersection[0] <= universal_bias_intersection[1], source_mean_bias_intersection_um=universal_bias_intersection, universal_offset_claim='No single value reproduces all source region mean differences to reporting precision. This does not by itself reject a population common-effect statistical model.', ecological_correlation=checks, sufficiency_test=suff, x13=targets, equally_informed_control='Independent direct OLS arithmetic agrees with lstsq; exact covariance identities and permutation witness agree. No algorithm superiority claimed.', empirical_prediction_enclosure=None, seconds=time.perf_counter() - started)
    save('RESULTS_R2.json', result)
    print(json.dumps(dict(round='R2', region_pairs=len(rows), ecological_gate=[r['external_coefficient_gate'] for r in checks], conditional_cells=targets['conditional_eligible_cells'], confirmed_flips=0, scenarios=[(s['donor_system'], s['region'], s['lower_um'], s['upper_um'], s['conditional_flips']) for s in targets['conditional_scenarios']]), ensure_ascii=False))
    return result
if __name__ == '__main__':
    run()
