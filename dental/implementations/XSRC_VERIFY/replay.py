"""Offline audit replay. No network and no writes to recipient lanes.
Manually reviewed source contracts remain pending independent scientific review.
"""
from dental_release.paths import expand as _release_expand
import copy, csv, datetime, hashlib, json, math, resource, time, re
import xml.etree.ElementTree as E
from collections import Counter
from decimal import Decimal as Dec
from fractions import Fraction as F
from pathlib import Path
from scipy.stats import beta, binomtest
from scipy.optimize import brentq
R = Path(__file__).resolve().parent

def read(name):
    return json.loads((R / name).read_text())

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write(name, data):
    p = R / name
    p.parent.mkdir(exist_ok=True, parents=True)
    p.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def leaves(x):
    if isinstance(x, dict):
        for v in x.values():
            yield from leaves(v)
    elif isinstance(x, list):
        for v in x:
            yield from leaves(v)
    else:
        yield str(x)

def source_text(p):
    if p.suffix == '.xml':
        tree = E.parse(p).getroot()
        return '\n'.join((' '.join(''.join(x.itertext()).split()) for tag in ['p', 'td', 'th', 'article-title'] for x in tree.findall('.//' + tag)))
    if p.suffix == '.json':
        return p.read_text()
    if p.suffix == '.pdf':
        return p.with_suffix('.txt').read_text()
    return p.read_text()

def valid_numeric(row):
    if sha(R / row['source_file']) != row['source_sha256']:
        return (False, 'HASH_MISMATCH')
    txt = source_text(R / row['source_file'])
    if row['source_file'].endswith('.pdf'):
        return (True, 'PRIMARY_DIAGRAM_MANUAL_READ_PLUS_HASH')
    printed_words = {'P04_Ramani': ('ninety-three', '93'), 'P05_Hoang': ('twenty-three', '23'), 'P06_Leong': ('six articles', '6'), 'P14_Raghav': ('twenty-one', '21'), 'P17_Li_review': ('seventeen', '17'), 'P19_Zanini': ('fifty-three', '53')}
    if row['id'] in printed_words:
        (phrase, count) = printed_words[row['id']]
        assert phrase in txt.lower(), (row['id'], phrase)
        txt += '\n' + count
    missing = []
    for s in leaves(row['values']):
        value = Dec(s)
        normalized = re.sub('[-−–]\\s+(?=\\d)', '-', txt.replace('−', '-'))
        tokens = re.findall('(?<![\\d.])(?:\\d+\\.\\d+|\\.\\d+|\\d+)(?![\\d.])', normalized)
        tokens += re.findall('(?<![\\d.])-(?:\\d+\\.\\d+|\\.\\d+|\\d+)(?![\\d.])', normalized)
        if not any((Dec(t) == value for t in tokens)):
            missing.append(s)
    return (not missing, 'DECIMALS_PRESENT' if not missing else 'MISSING:' + ','.join(missing))

def contract_valid(row, reference):
    fields = ['values', 'unit', 'population', 'protocol', 'quantity', 'table_locator', 'observation_resolution', 'statistic_resolution', 'time_scale', 'value_metadata', 'transfer_status', 'source_sha256']
    return all((row[k] == reference[k] for k in fields))

def exact_ci(k, n):
    return [0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1)), 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))]

def xml_rows(p, tid):
    root = E.parse(R / p).getroot()
    tab = next((t for t in root.findall('.//table-wrap') if t.get('id') == tid))
    return [[' '.join(''.join(c.itertext()).split()) for c in row] for row in tab.findall('.//tr')]
start = time.perf_counter()
lock = read('CURATED_SOURCE_LOCK.json')
for (name, digest) in lock['sha256'].items():
    assert sha(R / name) == digest, name
for item in read('INPUT_INVENTORY_R4.json'):
    if item.get('frozen_report_snapshot'):
        assert sha(R / item['frozen_report_snapshot']) == item['sha256']
rescore_lock = read('CROWN_RESCORE_INPUT_LOCK.json')
assert sha(R / rescore_lock['snapshot']) == rescore_lock['sha256']
assert sha(R / 'SOURCE_ADMISSION.json') == read('SOURCE_ADMISSION_LOCK.json')['sha256']
for row in read('RECIPIENT_MANIFEST.json'):
    assert sha(R / 'recipient_snapshots' / Path(row['snapshot']).name) == row['sha256'], row['snapshot']
observations = read('VERIFIED_OBSERVATIONS.json')
audit = read('SOURCE_AUDIT.json')
admission = read('SOURCE_ADMISSION.json')
checks = []
for row in observations:
    (ok, detail) = valid_numeric(row)
    checks.append({'id': row['id'], 'pass': ok, 'detail': detail})
write('raw/NUMERIC_PRESENCE_CHECKS.json', checks)
assert all((c['pass'] for c in checks)), [c for c in checks if not c['pass']]
from scoped_source_gate import check as scoped_check
c04 = next((x for x in observations if x['id'] == 'C04_regions'))
assert scoped_check(c04, R), 'C04 named primary table cells differ'
c04_bad = copy.deepcopy(c04)
(c04_bad['values']['PIC_occlusal'], c04_bad['values']['PIC_axial']) = (c04_bad['values']['PIC_axial'], c04_bad['values']['PIC_occlusal'])
assert not scoped_check(c04_bad, R), 'C04 swapped-region fault accepted'
write('raw/SCOPED_SOURCE_GATE.json', dict(valid_row_accepted=True, known_numbers_wrong_region_rejected=True, conclusion_changed=False, scope='C04 Table2 only'))
negative = []
for row in observations:
    assert row['values'].keys() == row['value_metadata'].keys()
    for meta in row['value_metadata'].values():
        assert set(meta) >= {'unit', 'observation_resolution', 'statistic_resolution', 'time_scale'}
        assert meta['time_scale'] in {'SIMULTANEOUS', 'HANDOVER'}
    for (field, bad) in [('unit', 'INJECTED_WRONG_UNIT'), ('population', 'INJECTED_WRONG_POPULATION'), ('protocol', 'INJECTED_WRONG_PROTOCOL'), ('quantity', 'INJECTED_WRONG_ENDPOINT'), ('values', {'injected': '999999999'}), ('transfer_status', 'MATCHED_CLINICAL_CALIBRATION'), ('statistic_resolution', 'INJECTED_WRONG_STATISTIC'), ('time_scale', 'INJECTED_WRONG_TIME_EDGE'), ('value_metadata', {'injected': {'unit': 'INJECTED_WRONG_UNIT'}}), ('table_locator', 'INJECTED_WRONG_TABLE')]:
        mutant = copy.deepcopy(row)
        mutant[field] = bad
        rejected = not contract_valid(mutant, row)
        assert rejected
        negative.append({'id': row['id'], 'field': field, 'injected_error_rejected': rejected})
write('raw/CONTRACT_INJECTIONS.json', negative)
assert next((x for x in observations if x['id'] == 'F05_maxshare'))['value_metadata']['controls_teeth_in_occlusion_range']['unit'] == 'teeth'
assert all((m['unit'] == ('degree' if 'angle' in k else 'mm') for x in observations if x['id'] in {'G01_Varga', 'G02_Wu'} for (k, m) in x['value_metadata'].items()))
assert all((x['time_scale'] == 'SIMULTANEOUS' for x in observations if x['id'] in {'P10_Donnermeyer', 'P11_Jang', 'P12_Allihaibi', 'P13_Piecha', 'P14_Raghav', 'P15_Hazard', 'P16_Ma', 'P20_Mayo'}))
force_tab = xml_rows('sources/PMC13168838.xml', 'joor70173-tbl-0003')
women = [r for r in force_tab if r and r[0] == 'Women']
assert women[0][1] == '179' and women[1][1] == '107', women
assert women[0][4].startswith('818.6') and women[1][4].startswith('443.7'), women
force_ratio = F('818.6') / F('443.7')
region_sum = sum(map(Dec, ['14.3', '50.1', '43.9']))
assert region_sum == Dec('108.3')
reader_tab = xml_rows('sources/PMC12872107.xml', 't2')
assert reader_tab[0][0] == 'Periapical Radiography' and 'CBCT' in reader_tab[0][1]
assert [int(re.match('\\d+', x)[0]) for x in reader_tab[2][1:]] == [199, 54, 253]
assert [int(re.match('\\d+', x)[0]) for x in reader_tab[3][1:]] == [204, 143, 347]
reader_delta = F(403, 600) - F(253, 600)
reader_mcnemar = float(binomtest(54, 258, 0.5).pvalue)
table_controls = [{'control': 'Women without dentures, Table 3', 'expected': women[0][4].split()[0], 'injected_value': '443.7', 'injection_rejected': not women[0][4].startswith('443.7')}, {'control': 'Paired reader table PR-correct row', 'expected': [199, 54, 253], 'injected_value': [199, 204, 403], 'injection_rejected': [199, 204, 403] != [int(re.match('\\d+', x)[0]) for x in reader_tab[2][1:]]}, {'control': 'Regional relative force conservation', 'expected_total_percent': '100', 'observed_total_percent': str(region_sum), 'injected_value': '100', 'injection_rejected': Dec('100') != region_sum}]
assert all((x['injection_rejected'] for x in table_controls))
icc_table = xml_rows('sources/PMC8700671.xml', 'diagnostics-11-02292-t002')
icc_raw = next((r for r in icc_table if any((c.startswith('B to A point') for c in r))))
icc_index = next((i for (i, c) in enumerate(icc_raw) if c.startswith('B to A point')))
icc_row = icc_raw[icc_index:]
assert list(map(Dec, icc_row[7:10])) == [Dec('.653'), Dec('.136'), Dec('.866')], icc_row
diameter_table = xml_rows('sources/PMC6190233.xml', 'micromachines-08-00221-t004')
diameter_row = next((r for r in diameter_table if r and r[0] == 'Tool diameter'))
assert diameter_row[1] == '802.2 μm', diameter_row
table_controls.extend([{'control': 'Cephalometry Table2 3D ICC estimate, not CI lower bound', 'expected': icc_row[7], 'injected_value': '.136', 'injection_rejected': Dec(icc_row[7]) != Dec('.136')}, {'control': 'Micro-milling Table4 tool diameter with printed unit', 'expected': diameter_row[1], 'injected_value': '.802 μm', 'injection_rejected': diameter_row[1] != '.802 μm'}])
assert all((x['injection_rejected'] for x in table_controls))
write('raw/TABLE_CONTROLS.json', table_controls)
pulp = read('recipient_snapshots/LANE_X88_PULPOTOMY_DECISION_results.json')
late = F(3, 47)
intra = F(25, 86)
late_result = {'later_RCT': {'k': 3, 'n': 47, 'fraction': float(late), 'ci95': exact_ci(3, 47), 'resolution': 'POPULATION', 'endpoint': 'Later RCT among trauma-treated teeth', 'table': 'Tables2/4, DOI10.1016/j.jds.2025.10.031'}, 'intraoperative_comparator': {'k': 25, 'n': 86, 'fraction': float(intra), 'ci95': exact_ci(25, 86), 'resolution': 'POPULATION', 'endpoint': 'Intraoperative conversion among planned caries pulpotomies', 'locator': 'DOI10.1111/iej.14144 Table4+randomizedarmflow'}, 'difference_pp': float((late - intra) * 100), 'pooling_allowed': False, 'reason': 'Different intervention time, indication and denominator; not a validation/refutation of the25/86rate.', 'table2_PDL_descriptive': {'widened': [2, 4], 'normal': [1, 43], 'risk_ratio': float(F(2, 4) / F(1, 43)), 'widened_ci95': exact_ci(2, 4), 'normal_ci95': exact_ci(1, 43), 'model_calibration': 'REJECTED_TABLE_CONFLICT', 'conflict': 'Table3 reverses group labels;Table4 operator/restoration disagreesTable2; conditional counts only.'}}
late_result['CI_assumption'] = 'Clopper-Pearson intervals under a binomial sampling model; unequal follow-up, treatment selection and table conflicts are not covered. No physical transport or floating enclosure is claimed.'
profiles = {'freehand': {'e': (1.84, 1.05), 'a': (2.3, 1.02), 't': (6.98, 3.78)}, 'pilot_guided': {'e': (1.69, 1.06), 'a': (1.95, 1.28), 't': (6.01, 3.44)}, 'fully_guided': {'e': (1.34, 0.58), 'a': (1.62, 0.61), 't': (3.05, 1.42)}, 'dynamic_navigation': {'e': (1.36, 0.65), 'a': (1.48, 0.65), 't': (3.71, 1.32)}}
guide_rows = []
guide_maxerr = 0.0
for (label, p) in profiles.items():
    for coverage in [0.9, 0.95]:
        a = (1 - coverage) / 6
        k = math.sqrt((1 - a) / a)
        e = p['e'][0] + p['e'][1] * k
        ap = p['a'][0] + p['a'][1] * k
        angle = p['t'][0] + p['t'][1] * k
        budget = max(e, ap) + 4 * math.sin(min(math.radians(angle), math.pi) / 2)
        roots = [mu + brentq(lambda z: s * s / (s * s + z * z) - a, 0, s * 1000) for (mu, s) in p.values()]
        control = max(roots[:2]) + 4 * math.sin(min(math.radians(roots[2]), math.pi) / 2)
        err = abs(budget - control)
        guide_maxerr = max(err, guide_maxerr)
        assert err < 1e-10
        guide_rows.append({'guide': label, 'target': coverage, 'guide_budget_mm': budget, 'independent_control_error_mm': err, 'resolution': 'PHENOMENOLOGICAL', 'assumption': 'samplemoments treated as exactpopulationmoments; physicaljointcoverage UNKNOWN', 'injected_plus1mm_rejected': abs(budget + 1 - control) > 1e-10})
guide_original = read('recipient_snapshots/LANE_X87_NERVE_MARGIN_TABLE_results.json')
guide_total_rows = []
for original in guide_original['main_table']:
    fresh = next((x for x in guide_rows if x['guide'] == original['guide'] and x['target'] == original['target']))
    error = abs(fresh['guide_budget_mm'] - original['guide_budget_mm'])
    assert error < 1e-10
    if original['revision_budget_mm'] == 'INF':
        total = 'INF'
        total_error = None
        assert original['margin_mm'] == 'INF'
    else:
        total = fresh['guide_budget_mm'] + original['revision_budget_mm'] + original['residual_clearance_mm']
        total_error = abs(total - original['margin_mm'])
        assert total_error < 1e-10
    guide_total_rows.append({'guide': original['guide'], 'segment': original['segment'], 'target': original['target'], 'guide_budget_mm': fresh['guide_budget_mm'], 'inherited_revision_budget_mm': original['revision_budget_mm'], 'inherited_residual_mm': original['residual_clearance_mm'], 'total_margin_mm': total, 'original_margin_error_mm': total_error, 'resolution': 'PHENOMENOLOGICAL', 'physical_margin_mm': 'UNKNOWN', 'revision_geometry_replayed': False})
wall_rows = []
if (R / 'recipient_snapshots/CROWN_RESCORE.json').exists():
    data = read('recipient_snapshots/CROWN_RESCORE.json')
    candidates = [x for x in data['rows'] if x.get('tag') != 'scalar']
    assert len(candidates) == 72
    assert sum((bool(x.get('shape_envelope_pass', False)) for x in candidates)) == 0
    for t in [0.4, 0.5, 0.8, 1.0]:
        rows = [x for x in data['rows'] if x.get('status') == 'RESCORED' and x.get('tag') != 'scalar']
        known = [x for x in rows if x.get('wall', {}).get('sampled_minimum_mm') is not None]
        violated = sum((x['wall']['sampled_minimum_mm'] < t for x in known))
        lowerok = sum((x['wall'].get('vertex_cover_lower_mm', -1) >= t for x in known))
        wall_rows.append({'scenario_wall_mm': t, 'candidate_count': len(candidates), 'scored': len(known), 'unscored': len(candidates) - len(known), 'sampled_violations': violated, 'exact_arithmetic_cover_above': lowerok, 'unresolved': len(known) - violated - lowerok, 'rigorous_float_certificate': 'MISSING', 'resolution': 'PER_TOOTH', 'shape_plus_function_pass': 0, 'material_assignment': 'UNASSIGNED_SCENARIO;not a new clinical material validation'})
d1 = [0.125] * 100
d2 = [0.0] * 90 + [1.25] * 10
mu1 = sum(d1) / 100
mu2 = sum(d2) / 100
assert mu1 == mu2 == 0.125
q1 = sorted(d1)[94]
q2 = sorted(d2)[94]
g1 = [1.0, -1.0]
g2 = [-1.0, 1.0]
rms1 = math.sqrt(sum((x * x for x in g1)) / 2)
rms2 = math.sqrt(sum((x * x for x in g2)) / 2)
assert rms1 == rms2 and sorted(g1) == sorted(g2)
suff = {'distance_mean': {'summary_mm': mu1, 'identity_error_mm': abs(mu1 - mu2), 'p95_A_mm': q1, 'p95_B_mm': q2, 'downstream_difference_mm': q2 - q1, 'minimal_extension_for_this_query': 'Observed95thquantile atsame surface/registration; mean alone insufficient.', 'resolution': 'PER_SURFACE_REGION'}, 'regional_signed_loss': {'summary_RMS_mm': rms1, 'same_signed_multiset': True, 'identity_error_mm': abs(rms1 - rms2), 'nerve_region_loss_A_mm': g1[0], 'nerve_region_loss_B_mm': g2[0], 'downstream_difference_mm': 2.0, 'minimal_extension_for_this_pair': 'Keep signed loss at the queried location; pooledRMS and even multiset insufficient.', 'resolution': 'PER_POINT'}, 'external_referent': {'kind': 'our_own_fixture', 'locator': 'replay.py summarywitnessblock', 'compared_quantity': 'Exact finite-vector summary identity and downstream discrepancy', 'refutes_us': True}, 'claim': 'Mathematical insufficiency example; not external dental evidence.'}
mill = next((x for x in observations if x['id'] == 'M01_RMS'))
mill_recalc = [{'region': k, 'RMS_um': float(F(v[0])), 'ratio_to_20um': float(F(v[0]) / 20), 'resolution': 'PER_SURFACE_REGION', 'source': 'DOI10.3390/healthcare9080983 Table1', 'not_tool_deflection': True} for (k, v) in mill['values'].items()]
recomputed = {'force_population_correction': {'wrong_denture_median_N': '443.7', 'correct_non_denture_median_N': '818.6', 'median_ratio': float(force_ratio), 'median_difference_N': float(F('818.6') - F('443.7')), 'scope': 'POPULATION unpairedoldergroups;no forcevalidation', 'conclusion': 'WithdrawBunnyS4non-dentureclaim andBunnyT1T2nooverlap inference.'}, 'force_region_mass': {'printed_total_percent': str(region_sum), 'excess_pp': str(region_sum - 100), 'accepted_partition': False, 'resolution': 'POPULATION'}, 'reader_accuracy': {'PRcorrect': [253, 600], 'PR_CBCTcorrect': [403, 600], 'absolute_change_pp': float(reader_delta * 100), 'discordant_gain_loss': [204, 54], 'McNemar_exact_p': reader_mcnemar, 'resolution': 'POPULATION ofreader-casepairs,not600independentpatients', 'pulp_conversion_transport': 'REJECTED'}, 'pulpotomy': late_result, 'guide_replay': guide_rows, 'guide_control_maximum_error_mm': guide_maxerr, 'crown_wall_scenarios': wall_rows, 'milling': mill_recalc, 'sufficiency': suff, 'linear_sensitivity': 'None reported. Cantelli is an exact inequality under exactmoment assumption; numerical floatingenclosure and physical samplinguncertainty are MISSING.'}
write('raw/RECOMPUTED.json', recomputed)
recomputed['guide_total_replay'] = guide_total_rows
reader = recomputed['reader_accuracy']
reader['McNemar_iid_nominal_p'] = reader.pop('McNemar_exact_p')
reader['uncertainty'] = 'The nominal binomial/McNemar calculation assumes independent discordant pairs. These repeated reader/case observations are clustered. A robust uncertainty interval is UNKNOWN without reader/case-level responses; the descriptive 25 percentage-point change needs no independence assumption.'
write('raw/RECOMPUTED.json', recomputed)
new_scenarios = []
for (sid, groups) in [('G07_2024clinical', ['free', 'static', 'dynamic']), ('G12_positional', ['guide', 'free'])]:
    source = next((x for x in observations if x['id'] == sid))
    for group in groups:
        for target in [0.9, 0.95]:
            aa = (1 - target) / 6
            kk = math.sqrt((1 - aa) / aa)
            vv = [tuple(map(float, source['values'][group + '_' + field])) for field in ['entry', 'apex', 'angle']]
            z = [mu + kk * sd for (mu, sd) in vv]
            budget = max(z[:2]) + 4 * math.sin(min(math.radians(z[2]), math.pi) / 2)
            independent = [mu + brentq(lambda q: s * s / (s * s + q * q) - aa, 0, s * 1000) for (mu, s) in vv]
            control = max(independent[:2]) + 4 * math.sin(min(math.radians(independent[2]), math.pi) / 2)
            assert abs(budget - control) < 1e-10
            new_scenarios.append({'source_id': sid, 'group': group, 'target': target, 'conditional_budget_mm': budget, 'control_error_mm': abs(budget - control), 'injected_plus1mm_rejected': abs(budget + 1 - control) > 1e-10, 'resolution': 'PHENOMENOLOGICAL', 'physical_coverage': 'UNKNOWN', 'assumptions': 'Sample mean/SD treated as exact population moments; unchanged 2mm radius. New population, not replacement for mandibular recipient.', 'float_enclosure': 'MISSING'})
mechanics = {'beam_equations': 'delta=FL^3/(3EI); I=pi*d^4/64; k=3EI/L^3', 'doubling_d_stiffness_factor_exact': str(F(2) ** 4), 'doubling_d_deflection_factor_exact': str(F(1, 2) ** 4), 'equivalent_0p8D_I_factor_exact': str(F(4, 5) ** 4), 'equivalent_0p8D_deflection_factor_exact': str(F(5, 4) ** 4), 'bunny_k_inverse_d_cubed_rejected': True, 'resolution': 'PHENOMENOLOGICAL', 'empirical_calibration': 'UNKNOWN; these are exact algebra under beam closure, not empirical tool properties', 'injected_doubling_stiffness_8_rejected': F(2) ** 4 != 8}
recomputed['new_guide_source_scenarios'] = new_scenarios
recomputed['milling_beam_consistency'] = mechanics
recomputed['milling_primary_deflection'] = {'source_id': 'M07_Mamedov_deflection', 'X_measured_tip_range_um': [-38.99, 39.34], 'Y_measured_tip_range_um': [-41, 41.81], 'observed_peaks_exceed_declared20um_scenario': True, 'resolution': 'PER_POINT', 'recipient20um_scenario_refuted': False, 'reason': 'Different1.5mm WC tool and Al7050 cut; no matched force/geometry or1N displacement. Laser1mm from tip then extrapolation.', 'matched_force_displacement_stiffness': 'UNKNOWN', 'locator': 'DOI10.1016/j.procir.2013.06.080 printedp150 Section3.2/Fig5'}
write('raw/RECOMPUTED.json', recomputed)
impacts = [{'consumer': 'PROOF_LANE_CROWN_FIX_PREP_R2', 'changed_conclusion': False, 'changed_interval': False, 'answer': 'ManufacturerML/HT andadhesiveLS2 criteria confirmed; no new seated margin/film measurement. .5mmKATANAshoulder is not specified in primary diagram. Recipientalreadylabelsit scenario. Originalpartialgeometry/wholecrownfailure remains.', 'source_ids': ['C04_regions', 'C07_Kale', 'C08_Dauti', 'C09_Sarhan', 'C12_KATANA', 'C13_Ivoclar'], 'next_operation': 'Matched manufactured/seated regionalfilm andobservedfinishline; no pooledmean substitution.'}, {'consumer': 'PROOF_LANE_FULL_CROWN_DIAG', 'changed_conclusion': False, 'changed_interval': False, 'answer': 'No sourcecalibrates wholeouter-surfacearea p95 clinicalfloor. Confirmstechnician6pointobservations andlabialmeans mismatch. Conditionalwallscenario rerunsavedrows, shapeplusfunction remains0; generic.5mm ismaterialscenario only.', 'source_ids': ['S01_Abduo', 'S02_Nagata', 'S04_Oxilia', 'S07_Gurian', 'C12_KATANA', 'C13_Ivoclar'], 'next_operation': 'Two acceptedCADdesigns ofsamepreparation, savedregion/signedgapvectors andsamep95pipeline.'}, {'consumer': _release_expand('X94'), 'changed_conclusion': False, 'changed_interval': False, 'answer': 'X94usesexactRubioTable1 already; no update to itsquantilefailures. BunnyS4population andT1T2refutationwithdrawn. S5regional108.3% rejected asjointpartition.', 'source_ids': ['F04_Prescale', 'F05_maxshare', 'F06_regions_conflict', 'F07_reference', 'F08_Hattori'], 'next_operation': 'Sameindividual/sameMIP registeredabsolute andrelative force; no olderpopulationmedian ascalibration.'}, {'consumer': _release_expand('X87'), 'changed_conclusion': False, 'changed_interval': False, 'answer': 'Bunnyguideharvestabsent. IndependentlycheckedrecipientVarga/Wumoments yield sameCantellibudgets, no clinicalmarginidentified. No canalinterannotator primaryharvest to add.', 'source_ids': ['G01_Varga', 'G02_Wu'], 'next_operation': 'Registered local signed gaploss andindependentwallcontour; sampleSDnot quantile.'}, {'consumer': _release_expand('X96'), 'changed_conclusion': False, 'changed_interval': False, 'answer': 'Unsignedimplantplacementmoments unchanged. Theydo not supplysigneddrill-deptherror orpairedinstrumentshape. Manufacturerlength/table results remainunchanged andphysicalrisk UNKNOWN.', 'source_ids': ['G01_Varga', 'G02_Wu'], 'next_operation': 'Sameboringtool/profile/datum andregisteredachievedpose; new signedgaplossmeasurement.'}, {'consumer': _release_expand('X88'), 'changed_conclusion': False, 'changed_interval': False, 'answer': 'No additionalmatchedintraoperativeconversionwithCBCTcategories. Later3/47traumaRCTseries isseparate; PDLHRcalibrationblockedbytableconflicts. Readerstudyaddsdiagnostic25ppcontrast, notconversionrisk. Existing25/86categoryintervals preserved.', 'source_ids': ['P01_Chew', 'P03_Gurler', 'P07_Sulaiman', 'P08_Li', 'P13_Piecha'], 'next_operation': 'Prospectiveplanned/treatmentreceived/latefailure denominators andpairedCBCT/bleedingprotocol, frozenpredictions.'}, {'consumer': _release_expand('X89'), 'changed_conclusion': False, 'changed_interval': False, 'answer': 'Bunnymillingharvestabsent. Recipient36.9-50.5umRMS confirmed;2.5/1.0/.6mmburPMMAprocess, notisolate mountedtooldeflection.1Nbeamlimit remainsPHENOMENOLOGICAL.', 'source_ids': ['M01_RMS'], 'next_operation': 'Measured mountedtoolforce-displacement envelope plus matchedCAMcut andregionalscan.'}]
write('CONSUMER_IMPACTS.json', impacts)
for item in impacts:
    if item['consumer'] == _release_expand('X87'):
        item.update(changed_information_coverage=True, answer='Actual guide harvest verified. Mandibular Varga/Wu replay remains unchanged. New SMCD measures centreline disagreement, not local wall loss; ICC .136 is CI lower bound, estimate .653. New primary source scenarios computed separately. No clinical joint clearance identified.', source_ids=['G01_Varga', 'G02_Wu', 'G07_2024clinical', 'G12_positional', 'G19_canal_centerline', 'G20_cadaver_index', 'G22_ceph_ICC', 'G24_Grunheid'], next_operation='Registered local signed wall-contour differences plus achieved tool pose; do not insert mean CI, patient-average SD or centreline SMCD as individual wall-error quantiles.')
    if item['consumer'] == _release_expand('X96'):
        item.update(changed_information_coverage=True, answer='Actual guide harvest verified; global unsigned placement/depth statistics and averaged patient poses do not supply signed drill overshoot. Native geometric/table conclusion and interval unchanged.', source_ids=['G01_Varga', 'G07_2024clinical', 'G12_positional', 'G15_Eldabe_patient', 'G19_canal_centerline'], next_operation='Pair exact tool cutting envelope/datum and measured achieved pose with local canal-wall contour; keep signed depth loss.')
    if item['consumer'] == _release_expand('X89'):
        item.update(changed_information_coverage=True, answer='Primary laser-deflection benchmark recovered: WC1.5mm tool/Al7050, X[-38.99,39.34],Y[-41,41.81]um. This changes benchmark availability, not the recipient force-displacement calibration. Runout/wear/PICN/RMS are distinct observables. Original150/256 conditional witnesses and UNKNOWN physical calibration unchanged.', source_ids=['M01_RMS', 'M02_equivalent_diameter', 'M03_runout_fit', 'M04_wear_force', 'M05_bur_wear', 'M06_PICN', 'M07_Mamedov_deflection', 'M08_crown_gaps'], next_operation='Use recovered primary protocol to design a matched mounted-force/displacement measurement; freeze geometry/force prediction before measurement. Do not invert separate force and displacement maxima into stiffness.')
write('CONSUMER_IMPACTS.json', impacts)
draft = ''.join((json.dumps(row, ensure_ascii=False, sort_keys=True) + '\n' for row in observations))
draftpath = R / 'raw/SOURCE_DRAFT.jsonl'
if draftpath.exists():
    assert draftpath.read_text() == draft, 'Existing draft differs; refusing overwrite.'
else:
    with draftpath.open('x') as f:
        f.write(draft)
write('raw/DRAFT_RECEIPT.json', {'path': str(draftpath), 'sha256': sha(draftpath), 'rows': len(observations), 'publication': 'DRAFT_ONLY_PENDING_REVIEW'})
job_counts = {}
for job in ['CROWN-MARGIN-WALL', 'TOOTH-FORCE-CLOSED-BITE', 'GUIDE-ACCURACY', 'PULPOTOMY-CONVERSION', 'MILLING-DEFLECTION', 'CONTRALATERAL-SYMMETRY']:
    aa = [x for x in audit if x['job'] == job]
    oo = [x for x in observations if x['job'] == job]
    ad = [x for x in admission if x['job'] == job]
    blocked = sum((not x['numeric_bundle_exported'] for x in ad))
    job_counts[job] = {'source_items_reviewed': len(aa), 'no_numerical_bundle_export_items': blocked, 'no_numerical_bundle_export_fraction': blocked / len(aa) if aa else None, 'dropout_reasons': [{'label': x['label'], 'reason': x['non_admission_reason']} for x in ad if not x['numeric_bundle_exported']], 'observation_bundles': len(oo), 'new_matched_physical_calibration_bundles': 0, 'specification_bundles': sum((x['transfer_status'] == 'SPECIFICATION_ONLY' for x in oo)), 'transfer_reasons': dict(Counter((x['transfer_status'] for x in oo))), 'bunny_report_delivered': next((x['exists'] for x in read('INPUT_INVENTORY_R4.json') if x['job'] == job))}
results = {'claim_type': 'information_link', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'main_answer': 'Scoped source observations now retain locator, unit, population and protocol for seven consumers. Harvest errors are corrected; no new matched physical calibration was identified.', 'source_jobs': job_counts, 'source_item_count': len(audit), 'observation_bundle_count': len(observations), 'numerical_field_contract_count': sum((len(x['value_metadata']) for x in observations)), 'no_numerical_bundle_export_items': sum((not x['numeric_bundle_exported'] for x in admission)), 'consumer_impacts': impacts, 'recomputed': recomputed, 'verification': {'literal_number_presence_pass': sum((x['pass'] for x in checks)), 'literal_number_presence_total': len(checks), 'contract_injections_rejected': len(negative), 'semantic_source_verification': 'Manual primary reading; replay cannot discover all scientific source errors. Injection checks test frozen contract integrity.', 'independent_guide_control_error_mm': guide_maxerr}, 'external_referent': {'kind': 'independent_measurement', 'locator': 'VERIFIED_OBSERVATIONS.json individual DOI/PMID/table locators; SOURCE_AUDIT.json rejected/corrected identities', 'compared_quantity': 'Printed numerical observations and their observation contracts', 'refutes_us': True}, 'missing_delivery': ['BT-DW48-HUNT-SRC-GUIDE-ACCURACY', 'BT-DW48-HUNT-SRC-MILLING-DEFLECTION'], 'not_claimed': ['physical calibration of recipient lanes', 'clinical recommendations', 'algorithm advantage', 'all harvests delivered', 'whole crown clinical p95 floor'], 'cost': {'replay_wall_seconds': time.perf_counter() - start, 'maxrss_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1, 'GPU': False, 'fit_seconds': 0, 'source_fetch_queries': len(read('FETCH_LOG.json')), 'source_query_errors': sum(('error' in x for x in read('FETCH_LOG.json'))), 'source_fulltext_errors': sum(('fulltext_error' in rec for x in read('FETCH_LOG.json') for rec in x.get('records', []))), 'discovery_and_manual_reading_seconds': 'UNKNOWN', 'physical_measurements': 0, 'fallback': 'Unavailable primary tables/quantities remain unadmitted; sources and conflicts are retained.'}, 'artifact_hashes': {'draft_jsonl': sha(draftpath), 'source_audit': sha(R / 'SOURCE_AUDIT.json'), 'observations': sha(R / 'VERIFIED_OBSERVATIONS.json'), 'prereg1': sha(R / 'PREREG_R1.json'), 'prereg2': sha(R / 'PREREG_R2.json'), 'prereg3': sha(R / 'PREREG_R3.json')}, 'large_arrays': []}
results['missing_delivery'] = []
results['not_claimed'].remove('all harvests delivered')
results['artifact_hashes']['prereg4'] = sha(R / 'PREREG_R4.json')
results['round'] = 'XSRC-verify-R4'
write('results.json', results)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
(fig, ax) = plt.subplots(1, 3, figsize=(13, 4))
ax[0].bar(['Denture women\nn = 107', 'Non-denture women\nn = 179'], [443.7, 818.6], color=['#a95050', '#35768b'])
ax[0].set_ylabel('Total arch median force (N)')
ax[0].set_title('Population label correction\nYoshikawa, Table 3')
vals = [float(intra), float(late)]
cis = [exact_ci(25, 86), exact_ci(3, 47)]
ax[1].bar(['Intraoperative\ncaries: 25/86', 'Later RCT\ntrauma: 3/47'], [v * 100 for v in vals], color=['#35768b', '#c69448'])
ax[1].errorbar([0, 1], [v * 100 for v in vals], yerr=[[100 * (v - c[0]) for (v, c) in zip(vals, cis)], [100 * (c[1] - v) for (v, c) in zip(vals, cis)]], fmt='none', color='black', capsize=4)
ax[1].set_ylabel('Observed proportion (%, exact 95% CI)')
ax[1].set_title('Different endpoints; no pooling\nPatel, Table 4 / Li, Tables 2, 4')
ax[2].bar([x['region'] for x in mill_recalc], [x['RMS_um'] for x in mill_recalc], color='#35768b', yerr=[float(mill['values'][x['region']][1]) for x in mill_recalc], capsize=4)
ax[2].axhline(20, color='#a95050', linestyle='--', label='Declared 20 um scenario')
ax[2].legend(fontsize=8)
ax[2].set_ylabel('Milled PMMA process RMS (um, +/- SD)')
ax[2].set_title('Measured manufacturing outcome\nHealthcare 2021, Table 1')
for a in ax:
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)
    a.tick_params(axis='x', labelsize=8)
fig.suptitle('Verified sources retain population, time and observation operator', fontsize=12)
fig.tight_layout()
fig.savefig(R / 'FIGURE.png', dpi=180)
fig.savefig(R / 'FIGURE.pdf')
plt.close(fig)
(fig, ax) = plt.subplots(1, 2, figsize=(10, 4))
for (i, (lo, hi)) in enumerate([(-38.99, 39.34), (-41, 41.81)]):
    ax[0].plot([lo, hi], [i, i], color='#35768b', linewidth=6, solid_capstyle='butt')
    ax[0].scatter([lo, hi], [i, i], color='#35768b')
    ax[0].text((lo + hi) / 2, i + 0.1, f'{lo:g} to {hi:g}', ha='center', fontsize=8)
ax[0].set_yticks([0, 1], ['X direction', 'Y direction'])
ax[0].set_ylim(-0.6, 1.6)
ax[0].set_xlim(-50, 67)
ax[0].axvline(0, color='gray', linewidth=0.7)
ax[0].set_xlabel('Observed tip-deflection range (um)')
ax[0].set_title('Recovered primary experiment\nMamedov 2013, printed p150 / Fig5')
ax[0].text(0.02, 0.04, 'WC tool 1.5 mm / Al7050; laser 1 mm from tip\n20 krpm; ap 300 um; feed 20 um/tooth', transform=ax[0].transAxes, fontsize=8)
ax[1].errorbar([0.653], [0], xerr=[[0.653 - 0.136], [0.866 - 0.653]], fmt='o', capsize=6, color='#35768b')
ax[1].scatter([0.136], [0.35], marker='x', s=70, color='#a95050')
ax[1].set_yticks([0, 0.35], ['ICC and 95% CI', 'Harvest misread as estimate'])
ax[1].set_ylim(-0.6, 0.8)
ax[1].set_xlim(0, 1)
ax[1].set_xlabel('ICC (dimensionless)')
ax[1].set_title('Estimate versus lower CI column\nDiagnostics 2021, Table2')
for a in ax:
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)
fig.tight_layout()
fig.savefig(R / 'FIGURE_R4.png', dpi=180)
fig.savefig(R / 'FIGURE_R4.pdf')
plt.close(fig)
results['cost']['replay_wall_seconds'] = time.perf_counter() - start
results['cost']['maxrss_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
write('results.json', results)
write('CURRENT_WORK_STATE.json', {'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'phase': 'ROUND_R3_REPLAY_COMPLETE', 'last_gate': f"51 observation contracts and {results['numerical_field_contract_count']} numerical field contracts checked; seven recipient decisions; {len(negative)} injected errors rejected; physical calibration remains UNKNOWN", 'next_operation': 'Independent review, two missing The swarm reports, matched physical measurements', 'review_state': 'PENDING_INDEPENDENT_REVIEW'})
write('CURRENT_WORK_STATE.json', {'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'phase': 'ROUND_R4_REPLAY_COMPLETE', 'last_gate': f"All six reports: {len(observations)} bundles, {results['numerical_field_contract_count']} fields, {len(audit)} named items; {len(negative)} integrity injections rejected; physical calibration remains UNKNOWN", 'next_operation': 'Independent source-contract review, then matched local wall/force-displacement measurements', 'review_state': 'PENDING_INDEPENDENT_REVIEW'})
actual_cost = results['cost']
write('raw/LATEST_REPLAY_COST.json', actual_cost)
cost_path = R / 'raw/RELEASE_COST.json'
if not cost_path.exists():
    write('raw/RELEASE_COST.json', actual_cost)
results['cost'] = read('raw/RELEASE_COST.json')
results['cost']['additional_crossref_identity_queries'] = 4 + len(read('raw/R4_CROSSREF_IDENTITY_CHECKS.json'))
results['cost']['web_tool_query_count'] = 'UNKNOWN'
results['artifact_hashes']['source_admission'] = sha(R / 'SOURCE_ADMISSION.json')
results['dropout'] = {'report_delivery': {'missing': 0, 'expected': 6, 'fraction': 0}, 'numeric_bundle_export': {'not_exported': sum((not x['numeric_bundle_exported'] for x in admission)), 'source_items': len(admission), 'fraction': sum((not x['numeric_bundle_exported'] for x in admission)) / len(admission), 'definition': 'A named audit item without an exported verified numeric bundle; specifications from two extra manufacturer documents are counted separately.'}, 'new_matched_physical_calibration': {'new_bundles': 0, 'screened_observation_bundles': len(observations) - 2, 'specifications_excluded': 2, 'definition': 'No new source measurement matched a missing physical consumer observable. Context observations and confirmed existing inputs remain useful.'}}
write('results.json', results)
print(json.dumps({'observations': len(observations), 'source_items': len(audit), 'injections_rejected': len(negative), 'guide_error_mm': guide_maxerr, 'wall_rows': wall_rows, 'draft_rows': len(observations)}, ensure_ascii=False))
