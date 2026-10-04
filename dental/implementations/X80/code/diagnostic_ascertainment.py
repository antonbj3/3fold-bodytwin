from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, hashlib, re, datetime
from fractions import Fraction
P = Path(__file__).resolve().parents[1]
contract = json.loads((P / 'PREREG_R6_DIAGNOSTIC_ASCERTAINMENT.json').read_text())
assert contract['claim_type'] == 'information_link'
meta = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/meta'))
ids = {'18705814', '22925444', '17261997'}
records = {}
scanned = 0
manifest = []
for fn in ['dental_implant.jsonl', 'dental_crown.jsonl']:
    path = meta / fn
    b = path.read_bytes()
    manifest.append(dict(path=str(path), sha256=hashlib.sha256(b).hexdigest(), bytes=len(b)))
    for (line_number, line) in enumerate(b.splitlines(), 1):
        scanned += 1
        r = json.loads(line)
        if str(r.get('id')) not in ids:
            continue
        id = str(r['id'])
        item = {'pmid': id, 'doi': r.get('doi'), 'title': r['title'], 'locator': 'line' + str(line_number), 'source_path': str(path), 'line_sha256': hashlib.sha256(line).hexdigest(), 'abstractText': r.get('abstractText', '')}
        if id in records:
            records[id]['duplicate_origins'].append({k: item[k] for k in ['source_path', 'locator', 'line_sha256']})
        else:
            records[id] = dict(item, duplicate_origins=[])
s = records['18705814']['abstractText']
assert 'Thirty-five subjects out of 1500' in s and 'Sixteen' in s and ('19' in s) and ('Thirty-five controls' in s) and ('Nine out of the 1500' in s)
known = Fraction(9, 1500)
selected = Fraction(9, 35)
missing = Fraction(1500 - 35, 1500)
upper = Fraction(9 + 1500 - 35, 1500)
sourcefacts = dict(cohort=1500, tested_selected_in_cohort=35, selected_ACRG=16, selected_PFG=19, additional_controls=35, positive_selected=9, positive_ACRG=8, positive_PFG=1, positive_controls=0, cohort_identified_positive_fraction=float(known), selected_tested_positive_fraction=float(selected), untested_fraction=float(missing), same_test_positive_rate_partial_identification=[float(known), float(upper)], disease_prevalence='UNKNOWN', test_sensitivity='UNKNOWN', test_specificity='UNKNOWN', patient_PPV='UNKNOWN', resolution='POPULATION', source_level='PRIMARY_ABSTRACT; full primary table not present locally', denominator_caveat='Random controls are sampled at allergy centre and not assumed to verify the remaining1465 cohort subjects. Bounds are same-test positivity, not true disease prevalence.')
a = {'TP': 4, 'FP': 16, 'FN': 1, 'TN': 79}
b = {'TP': 18, 'FP': 2, 'FN': 2, 'TN': 78}
summary_a = a['TP'] + a['FP']
summary_b = b['TP'] + b['FP']
ppv_a = Fraction(a['TP'], summary_a)
ppv_b = Fraction(b['TP'], summary_b)
assert sum(a.values()) == sum(b.values()) == 100 and summary_a == summary_b

def valid_observation_claim(claim):
    return claim.get('selected_denominator') == 35 and claim.get('cohort_denominator') == 1500 and (claim.get('disease_prevalence') == 'UNKNOWN') and (claim.get('reference_endpoint') == 'UNKNOWN')
good_claim = {'selected_denominator': 35, 'cohort_denominator': 1500, 'disease_prevalence': 'UNKNOWN', 'reference_endpoint': 'UNKNOWN'}
assert valid_observation_claim(good_claim)
challenges = []
for (name, field, bad_value) in [('selected denominator', 'selected_denominator', 1500), ('disease prevalence from identified positives', 'disease_prevalence', float(known)), ('implant failure as diagnostic gold standard', 'reference_endpoint', 'implant_loss')]:
    bad_claim = dict(good_claim)
    bad_claim[field] = bad_value
    rejected = not valid_observation_claim(bad_claim)
    assert rejected
    challenges.append({'name': name, 'valid_claim_accepted': True, 'injected_claim': bad_claim, 'rejected': rejected})
public = [dict(pmid=r['pmid'], doi=r['doi'], title=r['title'], source_path=r['source_path'], locator=r['locator'], line_sha256=r['line_sha256'], duplicate_origins=r['duplicate_origins']) for r in records.values()]
(P / 'raw/diagnostic_source_manifest.json').write_text(json.dumps({'files': manifest, 'records': public}, indent=2) + '\n')
res = {'round': 'R6', 'claim_type': 'information_link', 'outcome': 'TEST_ASCERTAINMENT_RESOLVED_DISEASE_PREVALENCE_UNKNOWN', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'external_referent': contract['external_referent'], 'source_facts': sourcefacts, 'sufficiency': {'summary': '20 test positives among100; no reference labels', 'identity_error': 0, 'integer_identical': True, 'states': [a, b], 'downstream_PPV': [float(ppv_a), float(ppv_b)], 'downstream_difference': float(ppv_b - ppv_a), 'minimum_extension': 'Independent disease-reference confusion matrix, with ascertainment/selection for target population', 'resolution': 'PHENOMENOLOGICAL', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'code/diagnostic_ascertainment.py', 'compared_quantity': 'PPV under two hypothetical reference-label assignments', 'refutes_us': True}}, 'control_challenges': challenges, 'other_primary_record': {'pmid': '22925444', 'doi': records['22925444']['doi'], 'quantity': 'TNFα and IL1β after titanium stimulation are associated with implant loss in retrospective109subjects', 'inference': 'Different endpoint from hypersensitivity; no Ti-test diagnostic confusion matrix in local abstract. Not a peri-implantitis causal experiment.'}, 'missing_primary_records': ['17261997'], 'screening': {'target_PMIDs': 3, 'unique_retained': len(records), 'missing': 3 - len(records), 'missing_fraction': (3 - len(records)) / 3, 'metadata_lines_scanned': scanned, 'duplicate_records_collapsed': sum((len(r['duplicate_origins']) for r in records.values()))}, 'clinical_recommendation': None, 'control': 'Exact rational fraction calculation and primary-abstract integer counts agree; every injected denominator/operator mismatch rejected.'}
(P / 'rounds/R6/results.json').write_text(json.dumps(res, indent=2) + '\n')
(P / 'rounds/R6/HANDOFF.md').write_text('The local primary abstract shows selective testing:35/1500 subjects plus35controls. The identified-positive share0.6% and selected positive share25.7% answer different questions. Unverified same-test positivity has broad count bounds; true disease prevalence and PPV remain UNKNOWN. MELISA cohort PMID17261997 is not in the targeted local metadata files and was not numerically incorporated.\n\nNext construction: obtain reference-diagnosis labels and the sampling strata before a patient-conditioned material comparison.\n')
r = json.loads((P / 'results.json').read_text())
r['rounds']['R6'] = res
r['preregistrations']['PREREG_R6_DIAGNOSTIC_ASCERTAINMENT.json'] = hashlib.sha256((P / 'PREREG_R6_DIAGNOSTIC_ASCERTAINMENT.json').read_bytes()).hexdigest()
r['diagnostic_sources_sha256'] = hashlib.sha256((P / 'raw/diagnostic_source_manifest.json').read_bytes()).hexdigest()
(P / 'results.json').write_text(json.dumps(r, indent=2, ensure_ascii=False) + '\n')
print(json.dumps({'R6': res['outcome'], 'identified_positive_pct': 100 * float(known), 'selected_positive_pct': 100 * float(selected), 'unverified_fraction': float(missing), 'same_test_rate_bounds': sourcefacts['same_test_positive_rate_partial_identification'], 'disease_prevalence': 'UNKNOWN'}, indent=2))
