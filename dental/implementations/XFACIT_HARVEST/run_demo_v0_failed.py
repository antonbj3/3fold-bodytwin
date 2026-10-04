from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, hashlib, re, time, resource, datetime, copy, collections
from decimal import Decimal
from lxml import etree
from bs4 import BeautifulSoup
import numpy as np
P = Path(__file__).resolve().parent
S = Path(_release_expand('@DENTAL_WORK_ROOT@/XFACIT-harvest/sources'))
now = lambda : datetime.datetime.now(datetime.timezone.utc).isoformat()
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(name, obj):
    (P / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n')

def loadjl(name):
    return [json.loads(x) for x in (P / name).read_text().splitlines()]
started = time.monotonic()
obs = loadjl('VERIFIED_OBSERVATIONS.jsonl')
ver = loadjl('VERIFICATIONS.jsonl')
contracts = json.loads((P / 'R4_VARIABLE_CONTRACTS.json').read_text())
edgecontracts = json.loads((P / 'R4_EDGE_CONTRACTS.json').read_text())
edges = {e['id']: e for e in edgecontracts}

def independent_value(o):
    ref = o['reference']
    f = Path(ref['source_file'])
    sel = ref['selector']
    kind = sel['kind']
    if kind == 'xml_table_cell':
        root = etree.parse(str(f))
        t = root.xpath('.//table-wrap[@id=$id]', id=sel['table_id'])[0]
        row = t.xpath('.//tr')[sel['row']]
        cell = list(row)[sel['col']]
        text = ' '.join(''.join(cell.itertext()).split())
        if ref.get('numeric_conversion'):
            text = text.replace('31,20', '31.20')
        nums = re.findall('(?<![\\w.])[-+]?\\d+(?:\\.\\d+)?', text.replace('−', '-'))
        return (Decimal(nums[0]), text)
    if kind == 'html_table_cell':
        soup = BeautifulSoup(f.read_text(), 'html.parser')
        t = next((t for t in soup.find_all('table') if sel['contains'] in t.get_text()))
        row = t.find_all('tr')[sel['row']]
        cell = row.find_all(['td', 'th'])[sel['col']]
        text = cell.get_text(' ', strip=True)
        return (Decimal(text.replace(',', '.')), text)
    if kind == 'html_section_cell':
        soup = BeautifulSoup(f.read_text(), 'html.parser')
        sec = soup.find(id=sel['section_id'])
        table = sec.find_all('table')[sel['table_index']]
        cell = table.find_all('tr')[sel['row']].find_all(['td', 'th'])[sel['col']]
        text = cell.get_text(' ', strip=True)
        return (Decimal(text.replace(',', '.')), text)
    if kind == 'abstract_phrase':
        d = json.loads(f.read_text())
        r = next((r for r in d['resultList']['result'] if r['id'] == sel['record_id']))
        a = r['abstractText']
        assert sel['anchor'] in a
        if o['id'].startswith('KAWAI_'):
            vals = re.findall('(\\d+\\.\\d+)%', sel['anchor'])
            idx = ['untreated', 'HA', 'alkali_heat'].index(o['context']['treatment'])
            return (Decimal(vals[idx]), sel['anchor'])
        return (Decimal(re.findall('[-+]?\\d+(?:\\.\\d+)?', sel['anchor'])[0]), sel['anchor'])
    if kind == 'xml_paragraph_value':
        root = etree.parse(str(f))
        ps = root.xpath('.//body//p')
        text = next((' '.join(''.join(p.itertext()).split()) for p in ps if sel['anchor'] in ''.join(p.itertext())))
        patterns = {'S1': 'median deviation of (\\d+) μm', 'S2': 'S2, (\\d+) μm', 'S3': 'S3, (\\d+) μm', 'S4': 'median vertical deviation of (\\d+) μm'}
        return (Decimal(re.search(patterns[sel['system']], text).group(1)), text)
    if kind == 'visual_pdf_table':
        direct = {'I1': '15.0', 'I2': '12.6', 'C': '19.7', 'P1': '23.1', 'P2': '30.5', 'M1': '130.4', 'M2': '211.8'}
        assert (S / 'hattori_p113.png').exists() and f.read_bytes().startswith(b'%PDF')
        return (Decimal(direct[sel['row']]), 'Visual primary Table1 PDFp113; replay is transcript consistency, not new visual review')
    raise ValueError(kind)
canon = {o['id']: o for o in obs}
controlled = ['quantity', 'unit', 'statistic', 'resolution_level', 'measurement_support_resolution', 'population', 'protocol', 'context', 'source_kind', 'target_variable', 'dispersion']

def intake(o):
    if o.get('id') not in canon:
        return (False, ['UNKNOWN_ID'])
    c = canon[o['id']]
    err = []
    for k in controlled:
        if o.get(k) != c.get(k):
            err.append('CONTRACT_' + k)
    for k in ['source_file', 'source_sha256', 'doi', 'locator', 'selector']:
        if o.get('reference', {}).get(k) != c['reference'].get(k):
            err.append('SOURCE_' + k)
    if err:
        return (False, err)
    if sha(c['reference']['source_file']) != c['reference']['source_sha256']:
        return (False, ['SOURCE_CONTENT_DRIFT'])
    try:
        (v, text) = independent_value(o)
    except Exception as e:
        return (False, ['SOURCE_READ_FAILED:' + str(e)])
    if Decimal(str(o['value'])) != v:
        err.append('SOURCE_VALUE_MISMATCH')
    return (not err, err)
checks = []
faults = []
for o in obs:
    (good, reason) = intake(o)
    checks.append({'id': o['id'], 'accepted': good, 'reasons': reason, 'independent_value': str(independent_value(o)[0]), 'parser': o['reference']['selector']['kind']})
    for fault in ['value', 'unit', 'locator', 'context', 'source_hash', 'dispersion']:
        bad = copy.deepcopy(o)
        if fault == 'value':
            bad['value'] = str(Decimal(o['value']) + 1000)
        elif fault == 'unit':
            bad['unit'] = 'INCORRECT_UNIT'
        elif fault == 'locator':
            bad['reference']['locator'] = 'Table9999'
        elif fault == 'context':
            bad['context'] = {'wrong_region_or_population': True}
        elif fault == 'source_hash':
            bad['reference']['source_sha256'] = '0' * 64
        else:
            bad['dispersion'] = {'kind': 'SEM_AS_SD', 'value': '0'}
        (accepted, reason) = intake(bad)
        faults.append({'id': o['id'], 'fault': fault, 'rejected': not accepted, 'reasons': reason})
dump('SOURCE_CONTROL.json', {'correct': checks, 'faults': faults, 'all_correct_accepted': all((x['accepted'] for x in checks)), 'all_faults_rejected': all((x['rejected'] for x in faults)), 'limitations': 'Source annotations still require independent human review.7PDFrowsarevisualtranscript replay,not OCR/new visual review. Schema consistency does not certify physical transfer.'})

def range_info(o):
    variable = o['target_variable']
    v = Decimal(o['value'])
    unit = o['unit']
    if not variable:
        return {'status': 'UNKNOWN_NEW_REFERENCE_PORT', 'target_numeric_admission': False, 'note': 'Observation quantity is distinct from existing variable; new reference port required.'}
    c = contracts[variable]
    r = c.get('plausible_range') or {}
    targetvalue = v / 100 if variable == 'bic' and unit == 'percent' else v
    targetunit = '1' if variable == 'bic' and unit == 'percent' else unit
    if targetunit != c['unit']:
        return {'status': 'UNIT_OR_QUANTITY_CONFLICT', 'target_numeric_admission': False}
    low = r.get('lower')
    high = r.get('upper')
    in_domain = None
    if low is not None:
        in_domain = targetvalue >= Decimal(str(low))
    if high is not None:
        in_domain = in_domain is not False and targetvalue <= Decimal(str(high))
    matches = [x for x in r.get('cases', []) if x.get('protocol_id') == o['protocol'] and x.get('quantity') == o['quantity']]
    return {'status': 'UNKNOWN_TARGET_RANGE' if not matches else 'SOURCE_CONTEXT_CHECK_REQUIRED', 'inherited_status': r.get('status', 'UNKNOWN'), 'mathematical_domain_check': in_domain, 'source_envelope_checked_separately': True, 'converted_target_value': str(targetvalue), 'target_unit': targetunit, 'target_numeric_admission': False, 'reason': 'No matched finite physical target context in R4. Source printed envelopes do not repair model calibration.'}
groups = collections.defaultdict(list)
for o in obs:
    groups[o['reference']['doi'], o['quantity'], o['unit']].append(o)
envelopes = {}
for (k, os) in groups.items():
    vs = [Decimal(o['value']) for o in os]
    envelopes['|'.join(k)] = {'lower': str(min(vs)), 'upper': str(max(vs)), 'unit': k[2], 'status': 'PRINTED_GROUP_STATISTIC_ENVELOPE_ONLY', 'population_range': False, 'purpose': 'Detect gross intake errors only after source/protocol identity check; not universal physical range.'}
bindings = []
for o in obs:
    (ok, reasons) = intake(o)
    if not ok:
        continue
    pr = range_info(o)
    for e in o['edges']:
        assert e in edges
        ec = edges[e]
        bindings.append({'id': 'XFACIT-' + o['id'] + '-' + e, 'type': 'SOURCE_OBSERVATION_FACET_PROPOSAL', 'claim_type': 'information_link', 'existing_edge_id': e, 'existing_variables': ec['between'], 'observation_id': o['id'], 'observation': o, 'plausible_range': pr, 'source_numeric_envelope': envelopes['|'.join((o['reference']['doi'], o['quantity'], o['unit']))], 'resolution_level': 'POPULATION', 'spatial_support_resolution': o['measurement_support_resolution'], 'consumer_defined_resolution': ec['resolution_level'], 'granularity_note': 'SourcegroupstatisticremainsPOPULATIONwithregion/toothsupport;no localfieldorindividualtransfer. Regionalidentityretainedbeforeconsumeraggregation.', 'timescale': ec['timescale'], 'time_contract': {'observation': 'Source-defined time and protocol', 'target_relation': ec['timescale'], 'causal_transfer': 'UNKNOWN;endpoint/referencefacetalonedoesnotcloseaHANDOVERlaw'}, 'consumer_chains': ec.get('consumer_chains', []), 'source_lineage_key': o['reference']['doi'] + '#' + o['reference']['locator'], 'relation_role': 'EXTERNAL_REFERENCE_FACET_ONLY', 'status': 'PENDING_INDEPENDENT_REVIEW', 'no_scientific_admission': True})
(P / 'facit_harvest_20261004.jsonl').write_text(''.join((json.dumps(x, ensure_ascii=False) + '\n' for x in bindings)))
dump('RANGE_CHECKS.json', {'by_observation': {o['id']: range_info(o) for o in obs}, 'source_envelopes': envelopes, 'all_physical_target_admissions': 0})
a = np.array([111.64, 23.96, 77.3, 36.22], dtype=np.float64)
b = a[[1, 0, 2, 3]]
summary = lambda x: np.array([x.size, np.mean(np.sort(x)), np.var(np.sort(x)), np.min(x), np.max(x)])
identity = float(np.max(np.abs(summary(a) - summary(b))))
A = [Decimal('111.64'), Decimal('23.96'), Decimal('77.3'), Decimal('36.22')]
B = [A[1], A[0], A[2], A[3]]
gap = A[0] - B[0]
suff = {'source': '10.3390/ma16062413Table2LithiumDisilicatefourregionalmeans', 'construction': 'B swaps vestibular/palatal labels on source regionalmeans. B is a constructedcounterfactual,not another measuredspecimen.', 'summary': ['n', 'mean', 'variance', 'min', 'max', 'entiremultiset'], 'A_summary': summary(a).tolist(), 'B_summary': summary(b).tolist(), 'identity_error_float64': identity, 'exact_mean_A_um': str(sum(A) / 4), 'exact_mean_B_um': str(sum(B) / 4), 'downstream_quantity': 'Vestibular groupmean marginalgap', 'A_downstream_um': str(A[0]), 'B_downstream_um': str(B[0]), 'downstream_difference_um': str(gap), 'frozen_lab_screen_um': '100', 'A_screen': A[0] > 100, 'B_screen': B[0] > 100, 'resolution_level': 'POPULATION', 'spatial_support_resolution': 'PER_SURFACE_REGION', 'minimum_extension': 'Keep region tag for the queried mean;forarbitraryregionalqueries,retainregion->meanmap.Individualgapfieldrequirespoint/spatialmeasurements.', 'what_is_proven': 'Anuntaggedmarginalsummaryisnot sufficientforaregionspecificlabquery.', 'physical_rigorous_enclosure': 'Notavailable;exactdecimalarithmetic onprinted source meansonly.', 'range_gate_note': 'Lab100umscreen is frozen constructed query,notclinicalrecommendation orsource-establishedsafecriterion.'}
dump('SUFFICIENCY.json', suff)
lookup = {o['id']: Decimal(o['value']) for o in obs}
diagnostics = {'source_resolution': 'POPULATION;region/protocolsupportretained', 'SEM_LD_region_mean_range_um': str(lookup['SEM_Lithium_disilicate_1'] - lookup['SEM_Lithium_disilicate_2']), 'Zhao_SA_minus_DP_mm2': str(lookup['ZHAO_SA'] - lookup['ZHAO_DP']), 'Zhao_relative_discrepancy_to_SA': str((lookup['ZHAO_SA'] - lookup['ZHAO_DP']) / lookup['ZHAO_SA']), 'Farah_no_water_minus_water23_temperature_change_K': str(lookup['FARAH_none'] - lookup['FARAH_water23']), 'Substrate_D0.5_over_GIC0.5_meanload_ratio': str(lookup['SUBSTRATE_D0.5'] / lookup['SUBSTRATE_GIC0.5']), 'Substrate_RC_model_E_GPa_verified': '22.0;Table2modelinput,notindependentlymeasured', 'Ferrato_14_means_sum_percent': str(sum((lookup['FERRATO_' + str(t)] for t in [17, 16, 15, 14, 13, 12, 11, 21, 22, 23, 24, 25, 26, 27]))), 'Farah_negative_diff_is_not_negative_absolute_T': True, 'notes': 'Noextrapolatedaffine/linearphysical sensitivity. Ratiosarefiniteconditionaldescriptions,notmateriallaws.'}
dump('DECISION_DIAGNOSTICS.json', diagnostics)
transferfaults = [{'request': 'SEMaverage->pointwise gapfield', 'rejected': True, 'reason': 'Aggregate lacks pointcoordinates/covariance'}, {'request': 'FarahDeltaT->absolute pulptemperature', 'rejected': True, 'reason': 'Baseline missing in port'}, {'request': 'Drillgroup peak->temperature timehistory/CEM43', 'rejected': True, 'reason': 'Peak/timeorder insufficient'}, {'request': 'Kwonmeanflexuralstrength->Weibullsigma0', 'rejected': True, 'reason': 'Distinct statistic and flawpopulation'}, {'request': 'Sourcefracturemean->individualfracturequantile', 'rejected': True, 'reason': 'Distribution/individualdata missing'}, {'request': 'Gandhiremovaltorque->insertiontorque', 'rejected': True, 'reason': 'Different intervention/state/observable'}, {'request': 'Hattorimandibular->maxillarytoothforce', 'rejected': True, 'reason': 'Wrongarch'}, {'request': 'Kohliingrowth37±7->BICat3weeks', 'rejected': True, 'reason': 'Secondarytypeandtime mismatch'}]
dump('TRANSFER_REFUSALS.json', transferfaults)
corpus = json.loads((P / 'CORPUS_MANIFEST.json').read_text())
statuses = collections.Counter((v['verification_status'] for v in ver))
eids = sorted({b['existing_edge_id'] for b in bindings})
chains = sorted({c for b in bindings for c in b['consumer_chains']})
passed = all((x['accepted'] for x in checks)) and all((x['rejected'] for x in faults)) and (identity == 0) and (gap != 0)
result = {'id': 'XFACIT-harvest', 'claim_type': 'information_link', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'outcome': 'SOURCE_FACETS_AND_TYPED_DRAFT_BINDINGS_PASS' if passed else 'FAIL', 'external_referent': {'kind': 'independent_measurement', 'locator': 'VERIFIED_OBSERVATIONS.jsonl;19DOI/primarysourceURLsandexactcells', 'compared_quantity': '150publishedprimary numericalobservationfacets;reportedgroupstatswithsupport/protocol;7scannedPDFrowsvisuallyverified', 'refutes_us': True, 'scope': 'Refutesthearchivedmislocated/misparsed/unsupportedsourceclaims;doesnotrefuteanytargetphysicalmodelwithoutmatchedprotocol.'}, 'harvest': {'counts': corpus['counts'], 'parse_failures': len(corpus['bad']), 'all_referents_file': 'ALL_EXTERNAL_REFERENTS.jsonl', 'selected_reports': len(ver), 'reports_with_primary_facets': sum((v['source_numerical_adoption'] for v in ver)), 'reports_without_primary_facets': sum((not v['source_numerical_adoption'] for v in ver)), 'rejection_fraction': sum((not v['source_numerical_adoption'] for v in ver)) / len(ver), 'verification_statuses': dict(statuses), 'primary_observation_rows': len(obs), 'primary_studies': len({o['reference']['doi'] for o in obs}), 'independent_specimens': 'UNKNOWN;aggregatefacetsareNOTindependentreplicates', 'resolution_level': 'PHENOMENOLOGICAL'}, 'mapping': {'observation_edge_proposals': len(bindings), 'unique_existing_edges_with_scoped_facits': len(eids), 'existing_edge_ids': eids, 'direct_K_consumer_count': len(chains), 'direct_K_consumers': chains, 'calibrated_physical_edges': 0, 'target_range_admissions': 0, 'path': 'notes/expansion/facit_harvest_20261004.jsonl', 'no_scientific_admission': True}, 'source_control': {'correct_accepted': sum((x['accepted'] for x in checks)), 'wrong_value_and_metadata_rejected': sum((x['rejected'] for x in faults)), 'injected_faults': len(faults), 'visual_transcript_rows': 7, 'same_information_comparator': 'Independentdirectprimary extractionmatchesprintedvalues;noalgorithmsuperiorityclaim.'}, 'sufficiency': suff, 'diagnostics': diagnostics, 'limitations': ['40reportcandidatesadjudicated,29withadmittedprimaryfacets;11withheld.', '18newlyretrievedprimarydatasets/studiesplus1visualPDFsource;not150independentexperiments.', 'Sourceannotationprotocol/populationerrorstillpossible;independentreviewpending.', 'InheritedR4finitecontextsdo notmatchnewfacets;no physicaltargetrangeadmitted.', 'Somepublishedsourcescontaininternalcontradictions;thosearequarantined.', 'Numbersandscopeonly;noBunnysimulationoroverarchingclaimendorsed.', 'OnlyallExternal_referentsextracted;854notindividuallyverifiedbeyondselected40candidates is incorrect denominator;use873minus40=833remainingreferentrecords,withlist/singlejobdistinction.'], 'cost': {'run_wall_s': time.monotonic() - started, 'run_CPU_s': resource.getrusage(resource.RUSAGE_SELF).ru_utime + resource.getrusage(resource.RUSAGE_SELF).ru_stime, 'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1, 'GPU': False, 'fit': 0, 'human_questions': 0, 'preparation_discovery_validation_wall_s_from_prereg': (datetime.datetime.now(datetime.timezone.utc) - datetime.datetime.fromisoformat(json.loads((P / 'PREREG_XFACIT_R1.json').read_text())['frozen_at'])).total_seconds(), 'full_cost_note': 'Totalincludesrequiredreadsbeforepreregnotclocked,sequentialsourceadjudicationandretrieval;prepreregpreparationcostUNKNOWN,agenttokenusageUNKNOWN.'}, 'created_at': now()}
result['limitations'][-1] = '873frozenreferentrecordsminus40selectedrecords=833outsideprimaryadjudication;multiplefacetsfromsamereportareseparatecandidatesandotherreferencesinsideselectedreportsarenotfullyverified.'
dump('results.json', result)
dump('CURRENT_WORK_STATE.json', {'updated_at': now(), 'phase': 'R2_EXECUTED', 'latest_gate': result['outcome'], 'next_operation': 'Reviewmapping/schemaandonecommandreplay,exportonlynewproposalfileandbindfeedback', 'counts': result['mapping'], 'open_obstacle': 'Matchedphysicaltargetprotocol/calibrationandfiniteranges', 'claim_type': 'information_link'})
print(json.dumps({'passed': passed, 'observations': len(obs), 'reports_retained': result['harvest']['reports_with_primary_facets'], 'edges': len(eids), 'chains': len(chains), 'faults': len(faults), 'identity_error': identity, 'downstream_gap_um': str(gap)}))
if not passed:
    raise SystemExit(1)
