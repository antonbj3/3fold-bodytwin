"""Source-conditioned dental reference inputs. Local data only; no graph writes."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
from lxml import etree as ET
import xml.etree.ElementTree as SE
import json, csv, hashlib, re, copy, datetime, time, resource, sys, importlib.util, struct, statistics, io, contextlib
from decimal import Decimal
import numpy as np
H = Path(__file__).resolve().parent
R = H.parent.parent
P = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
D = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/New Volume/coldstore/render_match_benchmarks/ti64-lpbf-psp'))
OUT = H / 'raw'
OUT.mkdir(exist_ok=True)
T0 = time.perf_counter()
sources = {}
faults = []
facts = {}
links = []
suff = []
checks = []

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(n, x):
    (H / n).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def src(p):
    p = Path(p)
    sources[str(p)] = {'path': str(p), 'sha256': sha(p), 'bytes': p.stat().st_size}
    return sources[str(p)]['sha256']

def txt(e):
    return re.sub('\\s+', ' ', ' '.join(e.itertext())).strip()

def nums(s):
    return [float(x.replace('−', '-')) for x in re.findall('[−+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?', s)]

def xloc(table, r, c):
    return f'(//table-wrap[@id="{table}"]//tr)[{r}]/*[{c}]'

def state(phase, last, nextop):
    save('CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-infolink4', phase=phase, last_gate=last, next_operation=nextop, updated_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))

def frozen(n):
    assert sha(H / n) == (H / (n + '.sha256')).read_text().strip(), f'Frozen drift: {n}'
for n in ['PREREG_INFOLINK4_R1.json', 'PREREG_INFOLINK4_R2.json', 'FROZEN_PREDICTIONS.json', 'LINK_CONTRACTS.json', 'LINK_CONTRACTS_V2.json', 'SOURCE_CONFLICTS.json']:
    frozen(n)
contracts = json.loads((H / 'LINK_CONTRACTS_V2.json').read_text())

def manual_lookup(path, spec):
    rt = SE.parse(path).getroot()
    if 'xpath' in spec:
        e = rt.findall('./body/sec')[2].findall('./p')[20]
        return (nums(txt(e))[spec['number_index']], None)
    table = next((t for t in rt.iter('table-wrap') if t.get('id') == spec['table']))
    trs = list(table.iter('tr'))
    e = list(trs[spec['row'] - 1])[spec['column'] - 1]
    ns = nums(txt(e))
    if 'sd_row' in spec:
        sd = nums(txt(list(trs[spec['sd_row'] - 1])[spec['sd_column'] - 1]))[0]
    else:
        sd = ns[1]
    return (ns[0], sd)

def validate_record(r, facit, contract, sourcehash):
    assert abs(r['original_value'] - facit['value']) <= 1e-10, 'value mismatch'
    assert r['reported_sd'] == facit['sd'] or (r['reported_sd'] is not None and facit['sd'] is not None and (abs(r['reported_sd'] - facit['sd']) <= 1e-10)), 'SD mismatch'
    for (key, expected) in [('unit', facit['unit']), ('condition_key', facit['key']), ('consumer', contract['consumer']), ('time_scale', contract['time_scale']), ('resolution_level', contract['resolution_level']), ('source_sha256', sourcehash), ('statistic', facit['statistic']), ('n', facit['n'])]:
        assert r[key] == expected, key + ' mismatch'
    assert r['xpath'] == facit.get('xpath', xloc(facit.get('table', ''), facit.get('row', 0), facit.get('column', 0))), 'locator mismatch'
    return True

def fault_test(name, mutator, validator, obj):
    x = copy.deepcopy(obj)
    mutator(x)
    try:
        validator(x)
    except (AssertionError, ValueError, KeyError, IndexError) as e:
        faults.append({'control': name, 'rejected': True, 'reason': str(e)})
        return
    faults.append({'control': name, 'rejected': False})
    raise AssertionError('Mutation accepted: ' + name)
pairs = {'L01': (0, 4), 'L03': (0, 1), 'L04': (0, 3), 'L05': (0, 2), 'L06': (0, 4), 'L07': (1, 2), 'L09': (0, 2), 'L10': (0, 2), 'L11': (0, 3), 'L15': (0, 1)}

def summary_counterexample(c, observations):
    (i, j) = pairs[c['id']]
    A = [o['original_value'] for o in observations]
    B = A.copy()
    (B[i], B[j]) = (B[j], B[i])
    SA = sorted(A)
    SB = sorted(B)
    err = max((abs(a - b) for (a, b) in zip(SA, SB)))
    assert err == 0 and struct.pack(f'{len(SA)}d', *SA) == struct.pack(f'{len(SB)}d', *SB)
    delta = float(Decimal(str(A[i])) - Decimal(str(B[i])))
    assert delta != 0
    changed = [k for k in observations[i]['condition_key'] if observations[i]['condition_key'][k] != observations[j]['condition_key'][k]]
    z = {'id': c['id'], 'state_type': 'CONSTRUCTED_INFORMATION_STATES_NOT_TWO_MEASURED_PHYSICAL_SAMPLES', 'state_A_values': A, 'state_B_values': B, 'summary': 'entire sorted unlabelled value multiset (stronger than mean or variance)', 'summary_identity_error': err, 'bitwise_equal': True, 'summary_resolution': 'POPULATION', 'downstream_query_key': observations[i]['condition_key'], 'query_A': A[i], 'query_B': B[i], 'downstream_difference_A_minus_B': delta, 'unit': observations[i]['unit'], 'minimal_extension_for_this_pair': changed, 'full_contract_extension': 'joint condition key attached to each source observation', 'interpretation': 'Loss of condition identity, not a newly measured biological effect; direct keyed lookup preserves it.'}

    def validate(z):
        assert sorted(z['state_A_values']) == sorted(z['state_B_values'])
        assert z['query_A'] == z['state_A_values'][i] and z['query_B'] == z['state_B_values'][i]
        assert z['downstream_difference_A_minus_B'] == float(Decimal(str(z['query_A'])) - Decimal(str(z['query_B'])))
    fault_test(c['id'] + '/sufficiency_query_wrong', lambda q: q.__setitem__('query_B', q['query_B'] + 1), validate, z)
    suff.append(z)
state('EXTRACT_PRIMARY', 'Frozen retrospective contracts loaded', 'Extract 85 source-conditioned observations; check with second XML parser')
for c in contracts:
    p = P / (c['source'] + '.xml')
    ss = src(p)
    root = ET.parse(str(p))
    obs = []
    me = root.xpath(c['method_xpath'])[0]
    method = txt(me)
    for term in c['method_terms']:
        assert term in method, (c['id'], 'missing method term', term)
    unittext = txt(root.xpath(f'''//table-wrap[@id="{c['unit_table']}"]''')[0]) if c['unit_table'] else method + ' ' + txt(root.xpath('/article/body')[0])
    assert c['unit_term'] in unittext, (c['id'], 'missing unit witness')
    for (k, spec) in enumerate(c['records']):
        xp = spec.get('xpath', xloc(spec.get('table', ''), spec.get('row', 0), spec.get('column', 0)))
        e = root.xpath(xp)
        assert len(e) == 1, (c['id'], xp, len(e))
        text = txt(e[0])
        ns = nums(text)
        val = ns[spec['number_index']] if 'number_index' in spec else ns[0]
        sd = None if spec['sd'] is None else nums(txt(root.xpath(xloc(spec['table'], spec['sd_row'], spec['sd_column']))[0]))[0] if 'sd_row' in spec else ns[1]
        r = dict(link_id=c['id'], observation_index=k, source_path=str(p), xpath=xp, source_locator=str(p) + '#' + xp, source_sha256=ss, original_text=text, original_value=val, reported_sd=sd, unit=spec['unit'], condition_key=spec['key'], consumer=c['consumer'], time_scale=c['time_scale'], resolution_level=c['resolution_level'], statistic=spec['statistic'], n=spec['n'], uncertainty_kind='reported SD; not SEM or confidence interval' if sd is not None else 'SD UNKNOWN in text paragraph', leaf_status='EXTERNALLY_MEASURED')
        if 'sd_row' in spec:
            r['sd_xpath'] = xloc(spec['table'], spec['sd_row'], spec['sd_column'])
        validate_record(r, spec, c, ss)
        (vv, dd) = manual_lookup(p, spec)
        assert abs(vv - val) <= 1e-10 and dd == sd
        checks.append({'id': c['id'], 'observation_index': k, 'two_parser_agreement': True, 'manual_facit_agreement': True})
        obs.append(r)
    f = c['records'][0]
    r = obs[0]
    validate = lambda z: validate_record(z, f, c, ss)
    mutations = {'wrong_value': ('original_value', r['original_value'] + 1), 'wrong_SD': ('reported_sd', 999.0), 'wrong_unit': ('unit', 'WRONG_UNIT'), 'wrong_condition': ('condition_key', {'wrong': 'assignment'}), 'wrong_consumer': ('consumer', 'DENT-IF-IA-PRELOAD' if c['consumer'] != 'DENT-IF-IA-PRELOAD' else 'K34'), 'wrong_timescale': ('time_scale', 'SIMULTANEOUS' if c['time_scale'] == 'HANDOVER' else 'HANDOVER'), 'wrong_resolution': ('resolution_level', 'PER_POINT'), 'wrong_hash': ('source_sha256', '0' * 64), 'wrong_locator': ('xpath', '/article/front/article-title'), 'wrong_statistic': ('statistic', 'individual measurement'), 'wrong_n': ('n', 999)}
    for (name, (key, value)) in mutations.items():
        fault_test(c['id'] + '/' + name, lambda q, k=key, v=value: q.__setitem__(k, v), validate, r)
    facts[c['id']] = obs
    summary_counterexample(c, obs)
    title = txt(root.xpath('/article/front/article-meta/title-group/article-title')[0])
    doi = root.xpath('/article/front/article-meta/article-id[@pub-id-type="doi"]/text()')[0]
    lic = txt(root.xpath('/article/front/article-meta/permissions')[0])
    methods_locator = str(p) + '#' + c['method_xpath']
    first = obs[0]
    link = {'id': 'DENT-PROOF_LANE4-INFOLINK-' + c['id'], 'type': 'INFORMATION-LINK', 'status': 'OPEN', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'claim_type': 'information_link', 'claim': c['title'], 'source_id': c['source'], 'consumer': c['consumer'], 'original_value': first['original_value'], 'unit': first['unit'], 'statistic': first['statistic'], 'source_locator': first['source_locator'], 'resolution_level': c['resolution_level'], 'time_scale': c['time_scale'], 'depends_on': [], 'evidence': [str(p), 'results/PROOF_LANE_INFOLINK4/results.json'], 'relations': [{'type': 'input_for', 'source': c['source'], 'target': c['consumer'], 'target_kind': c['target_kind'], 'role': c['quantity'], 'claim_type': 'information_link', 'resolution_level': c['resolution_level'], 'time_scale': c['time_scale'], 'source_locator': first['source_locator'], 'semantics_status': 'PENDING_INDEPENDENT_REVIEW', 'note': c['limitations']}], 'source_contract': {'doi': doi, 'article_title': title, 'source_sha256': ss, 'license': lic, 'methods_locator': methods_locator, 'method_witness_terms': c['method_terms'], 'scope': c['scope'], 'observations': obs, 'native_resolution': 'published group statistics with joint condition keys; no individual data or covariance invented', 'limitations': c['limitations']}, 'uncertainty': 'Source SD retained per condition; sampling dependence, raw-data quality and instrument bias remain unresolved. No new confidence interval or equivalence claim.', 'current_practice_without_information': {'baseline': c['current_practice_without_information'], 'baseline_query_result': 'UNKNOWN_SOURCE_CONDITION', 'matched_source_mean_queries_resolved': len(obs), 'control_status': 'executed absent-input lookup; not a clinical-practice prevalence claim'}, 'sufficiency_test': 'results/PROOF_LANE_INFOLINK4/SUFFICIENCY_TESTS.json#' + c['id'], 'outcome': 'ESTABLISHED_DRAFT', 'physical_predictive_improvement': 'NOT_TESTED'}
    links.append(link)
    save('raw/' + c['source'] + '_observations.json', obs)
state('LEGACY_PORE_LINK', 'All publication values matched manual facit and independent parser', 'Read old cell and re-extract per-pore threshold observations, retaining refuted narrow hypothesis')
oldcell = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/projects/cad-to-simulation-G/scripts/physics_exp/g1335_tail_family_cert_real_ti64_lpbf_mixture_caveat_forced_on_real_keyhole_contaminated_specimens.py'))
oldresult = oldcell.parents[2] / 'evidence' / (oldcell.stem + '.json')
module = oldcell.parents[2] / 'src/cad_to_sim/tail_family_cert.py'
for p in [oldcell, oldresult, module]:
    src(p)
spec = importlib.util.spec_from_file_location('legacy_tail_family', module)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
legacy_out = io.StringIO()
with contextlib.redirect_stdout(legacy_out):
    assert m.selftest() == 0
(H / 'LEGACY_SELFTEST.stdout').write_text(legacy_out.getvalue())
old = json.loads(oldresult.read_text())
po = []
pmetrics = {}
pore_rows = []
for pid in [37, 38, 40, 36, 35, 11]:
    p = D / 'pores' / f'sample {pid}.csv'
    sh = src(p)
    with p.open(encoding='utf-8-sig') as f:
        lines = list(csv.reader(f))
    assert lines[1][2] == 'EqDiameter (mm)' and lines[1][4] == 'Sphericity'
    rows = [{'csv_line': i, 'pore_id': int(v[0]), 'diameter_mm': float(v[2]), 'sphericity': float(v[4])} for (i, v) in enumerate(lines[2:], 3)]
    assert all((x['diameter_mm'] > 0 for x in rows))
    ds = np.asarray([x['diameter_mm'] for x in rows])
    over = [x for x in rows if x['diameter_mm'] >= 0.5]
    fam = m.tail_family(ds)
    with p.open(encoding='utf-8-sig') as f:
        next(f)
        reader = csv.DictReader(f)
        independent = [x for x in reader if Decimal(x['EqDiameter (mm)']) >= Decimal('0.5')]
    assert len(independent) == len(over)
    expected = {37: (10649, 36), 38: (2475, 35)}
    if pid in expected:
        assert (len(rows), len(over)) == expected[pid]
    else:
        assert len(over) == 0
    if pid in [37, 38]:
        e = old['numbers']['contaminated_specimens'][str(pid)]
        assert fam['tail_family'] == e['kind'] and abs(fam['cv'] - e['cv']) <= 0.0005
        pmetrics[str(pid)] = {'n_pores': len(rows), 'count_d_ge_0_5mm': len(over), 'fraction': len(over) / len(rows), 'legacy_family': fam, 'cv_branch_hypothesis': 'REFUTED; actual CV exceeds .5, so this is not a hidden single_scale call', 'resolution_level': 'POPULATION', 'native_resolution': 'within one XCT pore catalog, not independent material samples'}
    else:
        pmetrics[str(pid)] = {'n_pores': len(rows), 'count_d_ge_0_5mm': len(over), 'resolution_level': 'POPULATION'}
    pore_rows.extend((dict(process_set=pid, source_path=str(p), **x) for x in rows))
    for r in over:
        po.append(dict(link_id='L12', source_path=str(p), source_locator=str(p) + f":{r['csv_line']}:EqDiameter (mm)", csv_line=r['csv_line'], pore_id=r['pore_id'], original_value=r['diameter_mm'], unit='mm', resolution_level='PER_POINT', native_resolution='one segmented pore, identified within XCT catalog', time_scale='SIMULTANEOUS', source_sha256=sh, process_set=pid, sphericity=r['sphericity']))
    z = {'count': len(over)}
    fault_test(f'L12/process_set_{pid}/wrong_count', lambda q: q.__setitem__('count', q['count'] + 1), lambda q: (_ for _ in ()).throw(AssertionError('count differs from Decimal CSV count')) if q['count'] != len(independent) else True, z)
    save(f'raw/pores_{pid}.json', rows)
A = np.array([1, 4, 4], dtype=float) / 16
B = np.array([2, 2, 5], dtype=float) / 16
SA = [len(A), float(A.mean()), float(np.mean(A * A))]
SB = [len(B), float(B.mean()), float(np.mean(B * B))]
assert SA == SB and struct.pack('3d', *SA) == struct.pack('3d', *SB)
cut = 9 / 32
qa = int(np.sum(A > cut))
qb = int(np.sum(B > cut))
suff.append({'id': 'L12', 'state_type': 'CONSTRUCTED_POSITIVE_PORE_DIAMETERS_NOT_NEW_MEASUREMENTS', 'A_mm': A.tolist(), 'B_mm': B.tolist(), 'summary_N_mean_second_moment': SA, 'summary_identity_error': 0.0, 'bitwise_equal': True, 'threshold_mm': cut, 'query': 'count diameter > threshold', 'query_A': qa, 'query_B': qb, 'downstream_difference_A_minus_B': qa - qb, 'unit': 'pore_count', 'minimal_extension_for_this_query': 'threshold exceedance count; maximum suffices only for the Boolean any-exceedance query', 'scope': 'No fatigue or defect mechanism inferred from this information counterexample.'})
assert qa != qb
facts['L12'] = po
save('raw/PORE_METRICS.json', pmetrics)
source_node = next((json.loads(s)['id'] for s in (R / 'notes/expansion/old_cells_20261002.jsonl').read_text().splitlines() if oldcell.name in s))
f = po[0]
link = {'id': 'DENT-PROOF_LANE4-INFOLINK-L12', 'type': 'INFORMATION-LINK', 'status': 'OPEN', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'claim_type': 'information_link', 'claim': 'Original pore catalogs retain upper-tail counts missing from a single family label; the old narrow CV-branch claim remains refuted.', 'source_id': source_node, 'source_node_id': source_node, 'consumer': 'DENT-MAT-TI64-LPBF', 'original_value': f['original_value'], 'unit': 'mm', 'source_locator': f['source_locator'], 'resolution_level': 'PER_POINT', 'time_scale': 'SIMULTANEOUS', 'depends_on': [], 'evidence': [str(oldcell), str(oldresult), str(D / 'pores/sample 37.csv'), str(D / 'pores/sample 38.csv'), 'results/PROOF_LANE_INFOLINK4/results.json'], 'relations': [{'type': 'input_for', 'source': source_node, 'target': 'DENT-MAT-TI64-LPBF', 'target_kind': 'draft_node', 'role': 'measured_pore_upper_tail_for_inspection', 'claim_type': 'information_link', 'resolution_level': 'PER_POINT', 'native_resolution': 'pore record + catalog identity; consumer aggregates', 'time_scale': 'SIMULTANEOUS', 'source_locator': f['source_locator'], 'semantics_status': 'PENDING_INDEPENDENT_REVIEW', 'note': 'Measured diameters and threshold counts only. No mechanism assignment, fatigue law, or same-coupon mechanical pairing.'}], 'source_contract': {'old_cell_locator': str(oldcell) + ':15-18;171-181', 'observations': po, 'checked_clean_process_sets': [40, 36, 35, 11], 'raw_registry': 'results/PROOF_LANE_INFOLINK4/raw/PORE_METRICS.json', 'local_primary_data': 'Luo et al. Ti64 PSP pore catalogs; EqDiameter (mm) header row2', 'publication_doi': '10.1016/j.addma.2022.102915', 'deposit': 'Zenodo 6587905 (inherited locator; local workbook metadata checked in rounds2/3)', 'license': 'UNKNOWN in available local original; originals not redistributed', 'limits': 'Large diameter alone does not identify keyhole versus lack-of-fusion. BIC mixture result not rerun, not evidence of physical mechanism. 0.5 mm is the inherited descriptive threshold, not a dental acceptance criterion. Six targeted catalogs, not all 42 process sets.'}, 'derived_result': pmetrics, 'current_practice_without_information': {'baseline': 'Legacy tail-family label without per-pore threshold query; process average/single family cannot return measured exceedance count', 'baseline_query_result': 'UNKNOWN_COUNT', 'added_information': 'catalog-keyed segmented pore rows', 'not_an_algorithm_competition': True}, 'sufficiency_test': 'results/PROOF_LANE_INFOLINK4/SUFFICIENCY_TESTS.json#L12', 'outcome': 'ESTABLISHED_DRAFT', 'uncertainty': 'Segmentation/detection uncertainty not measured. No failure probability assigned.', 'physical_predictive_improvement': 'NOT_TESTED'}
links.append(link)

def pore_check(z):
    assert z['original_value'] == f['original_value'] and z['unit'] == 'mm' and (z['process_set'] == f['process_set']) and (z['source_locator'] == f['source_locator'])
for (name, key, val) in [('wrong_value', 'original_value', f['original_value'] + 1), ('wrong_unit', 'unit', 'um'), ('wrong_catalog', 'process_set', 0), ('wrong_locator', 'source_locator', 'none')]:
    fault_test('L12/' + name, lambda q, k=key, v=val: q.__setitem__(k, v), pore_check, f)
baseline = {}
candidate = {}
for c in contracts:
    i = c['id']
    baseline[i] = ['UNKNOWN_SOURCE_CONDITION'] * len(facts[i])
    candidate[i] = [r['original_value'] for r in facts[i]]
assert all((all((v == 'UNKNOWN_SOURCE_CONDITION' for v in x)) for x in baseline.values()))
manifest = list(sources.values())
manifest.sort(key=lambda q: q['path'])
if (H / 'SOURCE_MANIFEST.json').exists():
    assert json.loads((H / 'SOURCE_MANIFEST.json').read_text()) == manifest, 'Source bytes changed'
else:
    save('SOURCE_MANIFEST.json', manifest)
rej = json.loads((H / 'REJECTIONS.json').read_text())
nall = len(links) + len(rej)
attr = {'assessed_candidates': nall, 'established': len(links), 'rejected': len(rej), 'rejected_fraction': len(rej) / nall, 'unique_candidates': 15, 'unique_rejected': 4, 'unique_rejected_fraction': 4 / 15, 'reason_counts': {k: sum((x['reason_code'] == k for x in rej)) for k in sorted(set((x['reason_code'] for x in rej)))}, 'denominator_note': 'Explicit candidate links only; 196 old index entries are not all adjudicated. 3 old-cell candidates deeply read; 1 retained,2 lack local primary facit.', 'record_attrition': {'curated_publication_observations': 85, 'accepted': 85, 'rejected': 0, 'scope': 'These selected cells only; not every numeric field or study in corpus.'}, 'pore_record_attrition': {'selected_catalogs': 6, 'read_rows': len(pore_rows), 'invalid_rows': 0, 'tail_records_retained': len(po), 'below_threshold_records': len(pore_rows) - len(po), 'below_threshold_is_not_rejection': True}}
assert 10 <= len(links) <= 15 and len(po) == 71
save('FACIT.json', facts)
save('OBSERVATION_REGISTRY.json', facts)
save('VALIDATION.json', {'all_passed': True, 'numeric_observations': len(checks), 'extraction_checks': checks, 'fault_injections': len(faults), 'all_faults_rejected': all((x['rejected'] for x in faults)), 'physical_model_validated': False})
save('FAULT_INJECTION.json', faults)
save('SUFFICIENCY_TESTS.json', suff)
save('ATTRITION.json', attr)
save('CONTROLS.json', {'claim_type': 'information_link', 'baseline_without_source_conditions': baseline, 'candidate_keyed_means': candidate, 'strong_equally_informed_control_role': 'independent source parser agrees; extraction validation only', 'queries_resolved': sum((len(x) for x in candidate.values())), 'algorithm_gain_claimed': False, 'prediction_accuracy_gain': 'NOT_MEASURED'})
(H / 'information_links.preview.jsonl').write_text(''.join((json.dumps(x, ensure_ascii=False, allow_nan=False) + '\n' for x in links)))
tq = time.perf_counter()
qdict = {json.dumps(x['condition_key'], sort_keys=True): x['original_value'] for x in facts['L01']}
construction = time.perf_counter() - tq
keys = list(qdict)
tq = time.perf_counter()
for _ in range(1000):
    for k in keys:
        _ = qdict[k]
qtime = time.perf_counter() - tq
summary = {'lane': 'PROOF_LANE-infolink4', 'claim_type': 'information_link', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'established_links': len(links), 'attrition': attr, 'numerical_observations': 85, 'per_pore_tail_records': len(po), 'external_lineages': len(contracts) + 1, 'fault_injections_rejected': len(faults), 'sufficiency_tests': len(suff), 'max_summary_identity_error': max((s['summary_identity_error'] for s in suff)), 'scientific_claim': 'Condition-keyed source inputs established; no new physical accuracy or clinical validity claim.', 'negative_results': ['Priority links already present.', 'L02 wrong endpoint; replaced by tensile L15 without changing thresholds.', 'L08 already used and independently checked.', 'Two periodontal/retention cell claims lack located local primary facit.', 'Legacy CV<.5 hypothesis remains refuted.', 'Source conflicts remain attached.'], 'physical_predictive_improvement': 'NOT_TESTED', 'affine_sensitivity': 'NONE_REPORTED; exact lookup/permutation/count operations only', 'cost': {'demo_wall_seconds_before_plot': time.perf_counter() - T0, 'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads_limit': 4, 'source_bytes_read_and_hashed': sum((s['bytes'] for s in manifest)), 'lookup_construction_seconds': construction, 'lookup_30000_questions_seconds': qtime, 'fit_seconds': 0, 'external_measurements': 0, 'user_questions': 0, 'network_requests': 0, 'preparation_discovery_and_reasoning': 'Not separately metered; start/final timestamps reported in workflow state. No invented model token estimate.', 'fallbacks': 'Explicit exclusions, no full dataset copies, no GPU, no archive extraction.', 'legacy_scope': 'unmodified tail_family and selftest rerun on six catalogs; no Gaussian-mixture or bootstrap rerun'}, 'prereg_sha256': sha(H / 'PREREG_INFOLINK4_R1.json'), 'prereg_revision_sha256': sha(H / 'PREREG_INFOLINK4_R2.json'), 'source_manifest_sha256': sha(H / 'SOURCE_MANIFEST.json'), 'large_arrays': [], 'generated_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
save('results.json', summary)
state('SOURCE_LINKS_VALIDATED', '11 links, source checks and exact insufficiency tests passed', 'Generate reviewable report/figure, install one new expansion file and run working build')
print(json.dumps({k: summary[k] for k in ['established_links', 'numerical_observations', 'per_pore_tail_records', 'fault_injections_rejected', 'sufficiency_tests']}, indent=2))
