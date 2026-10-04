"""Measured competence handover; no invented law from temperature to healing weeks."""
import math
import re
import time
import xml.etree.ElementTree as ET
from common import ROOT, CORPUS, load, dump, sha, state, verify_freeze
GROUPS = ['chisel', '50rpm', '500rpm', '1000rpm']

def pairs(text):
    return [tuple(map(float, p)) for p in re.findall('(\\d+(?:\\.\\d+)?)\\s*%?\\s*±\\s*(\\d+(?:\\.\\d+)?)', text)]

def extract():
    path = CORPUS / 'PMC10631188.xml'
    root = ET.parse(path).getroot()
    table = []
    for tbid in ['Tab1', 'Tab2']:
        tb = root.find('.//table-wrap[@id="' + tbid + '"]')
        vals = []
        for tr in tb.findall('.//tbody/tr'):
            cells = [' '.join(e.itertext()) for e in tr]
            if len(cells) > 1 and pairs(cells[1]):
                vals.append(pairs(cells[1])[0])
        if len(vals) != 8:
            raise ValueError('Unexpected temperature table ' + tbid)
        table.append(vals)
    p21 = ' '.join(root.find('.//sec[@id="Sec21"]/p').itertext())
    p22 = ' '.join(root.find('.//sec[@id="Sec22"]/p').itertext())
    p23 = ' '.join(root.find('.//sec[@id="Sec23"]/p').itertext())
    a = pairs(p21)
    b = pairs(p22)
    c = pairs(p23)
    if len(a) != 8 or len(b) != 8 or len(c) != 4:
        raise ValueError('Source outcome paragraph structure changed ' + str([len(a), len(b), len(c)]))
    rows = []
    for (i, group) in enumerate(GROUPS):
        rows.append(dict(group=group, n_reported=4, source_sha256=sha(path), population='canine mandibular harvested bone chips', resolution='POPULATION', finest_shared_region='PER_SURFACE_REGION: donor site/process group; individual specimen pairing unavailable', process=dict(rpm=None if group == 'chisel' else int(group[:-3]), irrigation=False, drill_diameter_mm=1, depth_mm=4, reported_motor_torque_Ncm=35, force_N='NOT_REPORTED', feed_mm_s='NOT_REPORTED', temperature_window_s=10, full_temperature_history='NOT_AVAILABLE', actual_heat_duration='NOT_ESTABLISHED_BY_T10'), surface_T10_C={'mean': table[0][i][0], 'SD': table[0][i][1], 'locator': 'doi:10.1186/s12903-023-03611-w; Table1; cortical central surface'}, intracortical_T10_C={'mean': table[1][i][0], 'SD': table[1][i][1], 'locator': 'doi:10.1186/s12903-023-03611-w; Table2;1mm apical,2mm depth'}, viability_day3_pct={'mean': a[i][0], 'SD': a[i][1], 'locator': 'doi:10.1186/s12903-023-03611-w;Sec21/p1;Fig5e'}, cell_yield_day14_x1e4_per_g={'mean': a[i + 4][0], 'SD': a[i + 4][1], 'locator': 'doi:10.1186/s12903-023-03611-w;Sec21/p1;Fig5f'}, ALP_after_week1_king_unit_per_gprot={'mean': b[i + 4][0], 'SD': b[i + 4][1], 'locator': 'doi:10.1186/s12903-023-03611-w;Sec22/p1;Fig6g'}, induced_mineralization_day21_OD560={'mean': c[i][0], 'SD': c[i][1], 'locator': 'doi:10.1186/s12903-023-03611-w;Sec23/p1;FigS1e'}, tissue='harvested chips, not retained peri-implant wall', assay_time_basis=dict(viability='reported assay day3; precise offset from surgery not established', cell_yield='two weeks of bone-chip outgrowth culture', ALP='after one week of osteogenic induction of cultured cells; not postoperative day7', mineralization='after21days of osteogenic induction; not postoperative day21'), later_healing_delay_weeks='UNKNOWN', ISQ_over_weeks='UNKNOWN'))
    literal = dict(surface=[(32.26, 0.32), (32.83, 0.78), (36.64, 0.44), (49.62, 1.68)], intracortical=[(33.5, 0.26), (35.15, 0.21), (37.85, 0.37), (50.33, 0.15)], viable=[(97.65, 0.68), (95.64, 1.44), (94.13, 1.56), (56.54, 2.75)], yield14=[(213.0, 11.27), (215.8, 12.03), (207.6, 6.65), (59.6, 4.93)], ALP=[(53.0, 1.04), (55.33, 1.68), (77.35, 4.56), (15.61, 1.91)], CaOD=[(2.89, 0.05), (2.91, 0.08), (6.42, 0.04), (1.69, 0.02)])
    keys = ['surface_T10_C', 'intracortical_T10_C', 'viability_day3_pct', 'cell_yield_day14_x1e4_per_g', 'ALP_after_week1_king_unit_per_gprot', 'induced_mineralization_day21_OD560']

    def parity(test):
        for (short, key) in zip(literal, keys):
            for (i, r) in enumerate(test):
                if max((abs(r[key][part] - literal[short][i][j]) for (j, part) in enumerate(['mean', 'SD']))) > 1e-08:
                    return False
        return True
    assert parity(rows)
    dump('raw/canine_group_measurements.json', rows)
    dump('raw/canine_independent_transcription.json', literal)
    return (rows, parity)

def mean_ci(mean, sd, n=4):
    margin = 3.182446305284263 * sd / math.sqrt(n)
    return [mean - margin, mean + margin]

def guard(record, tissue, day):
    if tissue != 'harvested bone chips' or day not in [3, 7, 14, 21]:
        return 'UNKNOWN'
    return 'OBSERVED_GROUP_SUMMARY'

def main():
    start = time.perf_counter()
    verify_freeze('PREREG_R3.json')
    verify_freeze('FROZEN_PREDICTIONS_R3.json')
    (rows, parity) = extract()
    (a, b, c, d) = rows
    ports = []
    for r in rows:
        v = r['viability_day3_pct']
        y = r['cell_yield_day14_x1e4_per_g']
        alp = r['ALP_after_week1_king_unit_per_gprot']
        ports.append(dict(group=r['group'], timescale='HANDOVER', source_sha256=r['source_sha256'], resolution='POPULATION summary; physical support PER_SURFACE_REGION donor site', fast_process_seconds=10, fast_time_note='measured T10 window; full duration unreported', state=dict(immediate_viable_fraction='UNKNOWN', viable_fraction_day3=v['mean'] / 100, viable_fraction_day3_SD=v['SD'] / 100, viable_fraction_day3_mean95CI=[x / 100 for x in mean_ci(v['mean'], v['SD'])], cell_yield_day14_per_g=y['mean'] * 10000.0, cell_yield_day14_mean95CI_per_g=[x * 10000.0 for x in mean_ci(y['mean'], y['SD'])], competence_ALP_after_week1=alp['mean'], competence_units='king unit/gprot', mineralization_after21days=r['induced_mineralization_day21_OD560']['mean']), biological_closure_debt=['Immediate retained-wall viability', 'damage-zone removal/revascularization rate', 'bone apposition and contact stiffness vs time', 'ISQ geometry and interface transfer law'], assay_time_basis=r['assay_time_basis'], human_recommendation='UNKNOWN', implant_ISQ='UNKNOWN', healing_delay_weeks='UNKNOWN'))
    v500 = c['viability_day3_pct']
    v1000 = d['viability_day3_pct']
    difference = v500['mean'] - v1000['mean']
    low = mean_ci(v500['mean'], v500['SD'])[0] - mean_ci(v1000['mean'], v1000['SD'])[1]
    high = mean_ci(v500['mean'], v500['SD'])[1] - mean_ci(v1000['mean'], v1000['SD'])[0]
    reversals = []
    for i in range(4):
        for j in range(i + 1, 4):
            vi = rows[i]['viability_day3_pct']
            vj = rows[j]['viability_day3_pct']
            ai = rows[i]['ALP_after_week1_king_unit_per_gprot']
            aj = rows[j]['ALP_after_week1_king_unit_per_gprot']
            if (vi['mean'] - vj['mean']) * (ai['mean'] - aj['mean']) < 0:
                (ci, cj) = (mean_ci(vi['mean'], vi['SD']), mean_ci(vj['mean'], vj['SD']))
                reversals.append(dict(groups=[GROUPS[i], GROUPS[j]], viability_diff_pp=vi['mean'] - vj['mean'], ALP_diff=ai['mean'] - aj['mean'], viability_intervals_overlap=max(ci[0], cj[0]) <= min(ci[1], cj[1]), status='mean-rank challenge, not a causal effect proof'))
    faults = {}
    for (name, key, factor, offset) in [('temperature_plus10', 'surface_T10_C', 1, 10), ('viability_percent_as_fraction', 'viability_day3_pct', 0.01, 0), ('yield_units_times10000', 'cell_yield_day14_x1e4_per_g', 10000.0, 0)]:
        bad = [{k: dict(v) if isinstance(v, dict) else v for (k, v) in r.items()} for r in rows]
        bad[0][key]['mean'] = bad[0][key]['mean'] * factor + offset
        faults[name + '_rejected'] = not parity(bad)
    faults['retained_wall_transfer_rejected'] = guard(rows[2], 'retained periimplant wall', 14) == 'UNKNOWN'
    faults['day14_as_day0_rejected'] = guard(rows[2], 'harvested bone chips', 0) == 'UNKNOWN'
    faults['human_ISQ_rejected'] = all((p['implant_ISQ'] == 'UNKNOWN' for p in ports))
    assert all(faults.values()), faults
    dump('BIOLOGICAL_HANDOVER_PORTS.json', ports)
    out = dict(round='R3', claim_type='information_link', outcome='Measured short-process to biological competence link recovered in canine bone-chip data; wall healing delay and ISQ remain unresolved', external_referent=load('PREREG_R3.json')['external_referent'], measured_contrasts=dict(viability_500minus1000_percentage_points=difference, conservative_difference_of_marginal95CI_pp=[low, high], yield_day14_500_over1000=c['cell_yield_day14_x1e4_per_g']['mean'] / d['cell_yield_day14_x1e4_per_g']['mean'], ALP_week1_500_over50=c['ALP_after_week1_king_unit_per_gprot']['mean'] / b['ALP_after_week1_king_unit_per_gprot']['mean'], resolution='POPULATION', interpretation='reported group summaries, not individual paired predictions'), gates={'G1_source_parity': parity(rows), 'G2_viability_contrast_supported': difference > 0 and low > 0, 'G3_single_channel_sufficiency': 'FAIL for mean-rank sufficiency; viability CI overlap limits inference', 'G4_faults': all(faults.values()), 'G5_no_clinical_extrapolation': all((p['implant_ISQ'] == 'UNKNOWN' for p in ports))}, rank_reversals=reversals, fault_controls=faults, comparison_current_practice=dict(within_10s_window_flags=[False] * 4, complete_47C60s_flags=['UNKNOWN'] * 4, reason='T10 is not full thermal history; measured interval cannot establish60s exposure', data_link='Viability/ALP/yield differ across process groups even when the flag supplies no later outcome'), attrition={'process_outcome_group_links': {'candidates': 4, 'kept': 4, 'rejected': 0}, 'individual_voxel_or_specimen_links': {'candidates': 4, 'kept': 0, 'rejected': 4, 'fraction': 1.0, 'reason': 'only group means; individual pairing and spatial chip coordinates unavailable'}, 'periimplant_wall_or_ISQ_transfers': {'candidates': 4, 'kept': 0, 'rejected': 4, 'fraction': 1.0, 'reason': 'wrong tissue support and missing outcome'}}, timing_limit='Assay days are culture/induction times, not measured postoperative healing days; ALP oneweekinduction and mineralization21daysinduction.', uncertainty=dict(source='SD,n4 reported', derived='Student-t marginal mean95CI under normal independent replicate assumption', difference_bounds='difference of marginal intervals; NOT a paired95CI', covariance='UNKNOWN', causal_thermal_effect='UNKNOWN: rpm, force, heat, chip structure co-vary'), full_cost=dict(fit_s=0, compute_validation_s=time.perf_counter() - start, queries=4, preparation='two source tables and three outcome paragraphs', discovery='outcome-seen reconstruction disclosed', fallback='wall/ISQ/healing delay UNKNOWN'))
    frozen = load('FROZEN_PREDICTIONS_R3.json')['predictions']
    contrasts = out['measured_contrasts']
    mappings = [('500rpm_minus1000rpm_day3_viability_pp', 'viability_500minus1000_percentage_points'), ('500rpm_over1000rpm_day14_cell_yield', 'yield_day14_500_over1000'), ('500rpm_over50rpm_week1_ALP', 'ALP_week1_500_over50')]
    frozen_checks = [dict(prediction_name=a, expected=frozen[a], measured=contrasts[b], absolute_error=abs(frozen[a] - contrasts[b]), tolerance=1e-08, passed=abs(frozen[a] - contrasts[b]) <= 1e-08) for (a, b) in mappings]
    frozen_check = dict(round='R3', source_parity_tolerance_used=1e-08, checks=frozen_checks, all_passed=all((x['passed'] for x in frozen_checks)), diagnosis='Frozen39.10pp is50rpm minus1000rpm; label says500rpm. Actual37.59pp. Original preserved; failure stays FAIL.')
    dump('FROZEN_PREDICTION_CHECK_R3.json', frozen_check)
    out['frozen_prediction_check'] = frozen_check
    dump('results_R3.json', out)
    text = f"# R3 handoff\nCanine bone-chip study DOI10.1186/s12903-023-03611-w now connects process groups to day3 viable fractions,\nday14 cell yield and week1 osteogenic competence; assays after21days retained separately.4/4 group links kept.\n500rpm vs1000rpm has {difference:.2f}percentage-point viability contrast and {out['measured_contrasts']['yield_day14_500_over1000']:.6f}x cell yield.\nMarginal t-interval difference bounds are [{low:.6f},{high:.6f}]pp under n4 independent normal-replicate assumptions;\nthey are NOT a paired95%CI. No independent replication or lab measurement was performed.\nSingle viability channel has opposite mean ordering to ALP; all reversed viability pairs have overlapping intervals.\nFrozen39.10pp contrast FAILS by1.51pp due to wrong group label; measured500-minus1000 contrast is37.59pp.\n10s temperature window does not establish complete exposure. No CEM43 dose or healing-week prediction inferred from T10.\nAll4 individual point/specimen links and all4 peri-implant-wall/ISQ transfers rejected. Source numbers were seen before\npreregistration; this is disclosed source reconstruction. Current controls and injected errors pass.\nNext construction: same-specimen wall-temperature history, depth-resolved viability/microcracks and serial stability;\nuse the supplied minimum measurement contract and bench template. Need empirical healing closure; no invented rate.\n"
    (ROOT / 'HANDOFF_R3.md').write_text(text)
    (ROOT / 'HANDOFF.md').write_text(text)
    state('R3_DECIDED', out['gates'], 'Assemble one-command demo,figure,measurement contract and pending graph feedback')
    print(out['outcome'])
if __name__ == '__main__':
    main()
