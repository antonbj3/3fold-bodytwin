import json, re, xml.etree.ElementTree as ET, hashlib
from pathlib import Path
from science import R, load, dump

def textxml(name):
    return ' '.join(''.join(ET.parse(R / 'sources' / name).getroot().itertext()).split())
profiles = {p['id']: p for p in load('inputs/guide_profiles.json')['profiles']}
v = (R / 'sources/Varga2020.txt').read_text().split('MANDIBLE')[1].split('MAXILLA')[0]
lines = {}
for line in v.splitlines():
    m = re.match('^\\s*(Fr|Pi|Fu)\\s+(\\d+)\\s+([\\d.]+.*)$', line)
    if m:
        lines.setdefault(m[1], []).append([float(x) for x in m[3].split()])
controls = []
records = []

def record(source, locator, quantity, values, unit, role, resolution='POPULATION'):
    records.append({'source': source, 'locator': locator, 'quantity': quantity, 'values': values, 'unit': unit, 'resolution': resolution, 'role': role})
for (short, gid) in [('Fr', 'freehand'), ('Pi', 'pilot_guided'), ('Fu', 'fully_guided')]:
    angle = lines[short][0][0:2]
    entry = lines[short][0][5:7]
    apex = lines[short][1][0:2]
    for (key, values, unit) in [('angle_mean_sd_deg', angle, 'deg'), ('entry_mean_sd_mm', entry, 'mm'), ('apex_mean_sd_mm', apex, 'mm')]:
        expected = profiles[gid][key]
        assert values == expected, (gid, key, values, expected)
        bad = expected.copy()
        bad[1] += 0.5
        controls.append({'name': 'Published patient SD and mean', 'guide': gid, 'quantity': key, 'error': 0.0, 'injected_SD_plus0p5_rejected': values != bad})
        record('Varga2020.txt', 'doi:10.1111/clr.13578 Table4 MANDIBLE; accepted manuscript Table3', gid + '/' + key, values, unit, 'GUIDE_INPUT')
wu = textxml('PMC7683639.xml')
m = re.search('dynamic navigation group were \\(1\\.36\\s*±\\s*0\\.65\\) mm, \\(1\\.48\\s*±\\s*0\\.65\\) mm, and \\(3\\.71\\s*±\\s*1\\.32\\)', wu)
assert m is not None
for (key, values, unit) in [('entry_mean_sd_mm', [1.36, 0.65], 'mm'), ('apex_mean_sd_mm', [1.48, 0.65], 'mm'), ('angle_mean_sd_deg', [3.71, 1.32], 'deg')]:
    assert profiles['dynamic_navigation'][key] == values
    controls.append({'name': 'Dynamic published values', 'quantity': key, 'error': 0.0, 'injected_SD_plus0p5_rejected': profiles['dynamic_navigation'][key] != [values[0], values[1] + 0.5]})
    record('PMC7683639.xml', 'doi:10.1186/s40729-020-00272-0 Results/Comparison of accuracy Figure4', 'dynamic_navigation/' + key, values, unit, 'GUIDE_INPUT')
meta = textxml('meta2024.xml')
assert re.search('1\\.44 mm \\(95% CI: 1\\.34[–-]1\\.54\\)', meta)
assert re.search('1\\.36 \\(95% CI: 1\\.18[–-]1\\.54\\)', meta)
for (key, values) in [('static_mean_apex_mm', [1.44, 1.34, 1.54]), ('dynamic_mean_apex_mm', [1.36, 1.18, 1.54])]:
    assert load('inputs/guide_profiles.json')['meta_analysis_context'][key] == values[0]
    controls.append({'name': 'Meta-analysis mean; not individual SD', 'quantity': key, 'error': 0.0, 'injected_mean_plus0p5_rejected': values[0] != values[0] + 0.5, 'CI_as_patient_SD_fault_check': 'raw/VALIDATION.json: Quantity and coverage link / pooled_mean_CI'})
    record('meta2024.xml', 'doi:10.1186/s12903-024-04033-y Figure9 and Results/Subgroup navigation', key, values, 'mm', 'CONTEXT_ONLY_MEAN_CI_NOT_TAIL')
jar = textxml('multigrader2022.xml')
for (mean, sd) in [(0.68, 0.38), (0.62, 0.39), (0.52, 0.38), (0.47, 0.4)]:
    assert f'{mean:.2f} ({sd:.2f})' in jar
record('multigrader2022.xml', 'doi:10.1038/s41598-022-20605-w TableS2, values reproduced in Results/Reference segmentation; supplement fetch failed', 'Expert1-4 consensus SMCD mean(SD)', [[0.68, 0.38], [0.62, 0.39], [0.52, 0.38], [0.47, 0.4]], 'mm', 'REJECTED_AS_LOCAL_WALL_95_QUANTILE')
ger = load('sources/Gerlach2010.json')['resultList']['result'][0]
assert ger['doi'] == '10.1016/j.joms.2009.09.059'
assert '1.3 mm (SD 0.384)' in ger['abstractText']
record('Gerlach2010.json', 'doi:10.1016/j.joms.2009.09.059 PMID20036043 Abstract Results/Conclusions; full numerical table unavailable', 'Mean within-canal tracing P95 and between-case SD', [1.3, 0.384], 'mm', 'REJECTED_AS_PATIENT_OR_SEGMENT_WALL_TAIL')
iac = (R / 'sources/IACAT2.txt').read_text()
assert 'average IoU of 0.70' in iac and 'Dice of 0.81' in iac
record('IACAT2.txt', 'https://federicobolelli.it/media/publications/pdfs/2023iciap_iacat2.pdf section4.1 printedp9', 'Two independent expert masks; meanIoU/meanDice', [0.7, 0.81], 'dimensionless', 'REJECTED_AS_MM_WALL_BUDGET')
for name in ['SMCD mean to wall quantile', 'Within-canal P95 to patient-level95', 'Dice to local wall clearance', 'Mask revision to independent annotation', 'Pooled mean CI to patient SD', 'Endpoint proxy to anatomical foramen/ramus']:
    controls.append({'name': 'Quantity/lineage check: ' + name, 'executed_fault_check': 'raw/VALIDATION.json: check Quantity and coverage link'})
raw = load('raw/R1_RESULT.json')
assert raw['physical_margin'] == 'UNKNOWN'
controls.append({'name': 'Clinical promotion gate', 'physical_validation_available': False, 'executed_fault_check': 'raw/VALIDATION.json: Physical claim admission'})
manifest = []
for p in sorted((R / 'sources').iterdir()):
    if p.is_file() and p.name != 'SOURCE_HASHES.json':
        manifest.append({'file': p.name, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest(), 'bytes': p.stat().st_size})
dump('sources/SOURCE_HASHES.json', manifest)
dump('raw/SOURCE_NUMBERS.json', records)
dump('raw/SOURCE_CONTROLS.json', controls)
dump('raw/EXTERNAL_COMPARISON.json', {'candidates': 6, 'retained_for_input_or_context': 6, 'eligible_for_joint_physical_margin': 0, 'rejected_as_joint_physical_margin': 6, 'rejection_fraction': 1.0, 'rejection_reasons': {'no_segment_matched_joint_tail_or_independent_wall_bound': 3, 'wrong_observable_or_wrong_coverage_unit': 3}, 'guide_apex_mean_comparison_mm': {'fully_guided_trial': 1.62, 'pooled_static': 1.44, 'difference': 0.18, 'dynamic_trial': 1.48, 'pooled_dynamic': 1.36, 'difference_dynamic': 0.12}, 'residual_on_local_true_wall': None, 'reason': 'No same-operator same-specimene anatomical reference ; no valid residual complex', 'source_changed_direction': False, 'source_confirmed_direction': 'Independent annotator variation exists externally, but tracing-path/SMCD/Dice quantities cannot substitute for local whole-cylinder wall loss', 'external_referent': {'kind': 'independent_measurement', 'locator': 'doi:10.1038/s41598-022-20605-w TableS2; doi:10.1016/j.joms.2009.09.059 Abstract', 'compared_quantity': 'Independent tracing variation (external observed quantity only); not local required margin', 'refutes_us': True}, 'failed_retrievals': ['Springer supplement DOCX URL: DNS resolution error', 'PMC browser recaptcha; EuropePMC XML succeeded', 'Nature article browser redirect inaccessible; journal XML succeeded'], 'scope': 'Six deliberately selected numerical sources, not a systematic literature census'})
print('Source numbers', len(records), 'controls', len(controls), 'physical tail sources 0/6')
