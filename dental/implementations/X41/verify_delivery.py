"""Fail closed on missing required outputs, wrong claims, or untested controls."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import datetime
import hashlib
import json
import re
ROOT = Path(__file__).resolve().parent

def read(p):
    return json.loads((ROOT / p).read_text())

def main():
    required = ['MANUSCRIPT.md', 'README_DEMO.md', 'CLAIM_CHECKLIST.md', 'results.json', 'RESULTS_R2.json', 'SOURCE_MANIFEST.json', 'SOURCE_MANIFEST_REFERENCES.json', 'SOURCE_MANIFEST_GUIDE_AND_DATASET.json', 'FROZEN_PREDICTIONS.json', 'NUMERIC_CLAIMS.json', 'MANUSCRIPT_NUMBER_LEDGER.json', 'EXTERNAL_FACIT_CHECK.json', 'ENVIRONMENT.json', 'FIGURE_ENVIRONMENT.json', 'raw/CONTROLS.json', 'raw/GUIDE_PARAMETER_REPLAY.json', 'EXAMPLE_ASSESSMENT_WITH_BOUNDS.json', 'DATASET_CHECK.json', _release_expand('MEASUREMENT_PLAN')]
    required += [f'figures/FIGURE_{i}_{name}.{ext}' for (i, name) in [(1, 'DECISION_ERRORS'), (2, 'RESOLUTION_FLOOR')] for ext in ['png', 'pdf', 'svg']]
    required += [f'tables/TABLE_{i}_{name}.{ext}' for (i, name) in [(1, 'CONNECTED_ERRORS'), (2, 'DECIDABILITY'), (3, 'GUIDE_SCENARIOS'), (4, 'IDENTIFIABILITY')] for ext in ['csv', 'md']]
    for p in required:
        if not (ROOT / p).is_file() or (ROOT / p).stat().st_size == 0:
            raise ValueError('Required artifact missing: ' + p)
    for p in ROOT.glob('PREREG*.json'):
        if hashlib.sha256(p.read_bytes()).hexdigest() != p.with_suffix('.sha256').read_text().strip():
            raise ValueError('Prereg drift: ' + p.name)
    for p in ROOT.glob('FROZEN_PREDICTIONS*.json'):
        if hashlib.sha256(p.read_bytes()).hexdigest() != p.with_suffix('.sha256').read_text().strip():
            raise ValueError('Prediction drift: ' + p.name)

    def valid_control_list(rows):
        return len(rows) > 0 and all((isinstance(r, dict) and r.get('pass') is True and (r.get('injected_rejected') is True) for r in rows))
    controls = read('raw/CONTROLS.json')
    guides = read('raw/GUIDE_PARAMETER_REPLAY.json')
    if not valid_control_list(controls) or not valid_control_list(guides['controls']):
        raise ValueError('Missing/failed numerical controls')
    if valid_control_list([]) or valid_control_list([{'pass': True}]):
        raise ValueError('Vacuous control acceptance')
    rows = read('RESULTS_R2.json')['checks']
    if not rows or not all((r.get('valid_pass') is True and r.get('fault_rejected') is True for r in rows)):
        raise ValueError('R2 controls missing/fail')
    r = read('results.json')
    c = r['conditional_canal']
    if len(r['counterexamples']) != 6 or c['class1'] != 3648 or c['class2'] != 484 or (c['physical_certificates'] != 0):
        raise ValueError('Summary changed')
    if read('EXAMPLE_ASSESSMENT_WITH_BOUNDS.json')['bounded_reference_decision']['decision'] != 'UNKNOWN_MISSING_BOUND_INFORMATION':
        raise ValueError('Missing reference bound silently certified')
    m = (ROOT / 'MANUSCRIPT.md').read_text()
    if '{{' in m or not all((x in m for x in ['88.29', '87.17', '89.39', '99.9828', '99.9965', 'PHENOMENOLOGICAL', '1 mm residual'])):
        raise ValueError('Manuscript numeric/scope render incomplete')
    cited = set(re.findall('https://doi\\.org/([^\\s)]+)', m))
    verified = read('REFERENCES_FINAL.json')['references']
    if cited != {x['doi'] for x in verified} or not all((x['identity_verified'] for x in verified)):
        raise ValueError('DOI citation coverage mismatch')
    claim_items = re.findall('^\\| (\\d+) \\|', (ROOT / 'CLAIM_CHECKLIST.md').read_text(), flags=re.M)
    if claim_items != [str(i) for i in range(1, 45)]:
        raise ValueError('CLAIM numbering incomplete')
    ds = read('references/TF2_dataset.json')
    if ds['license'] != 'CC-BY-SA 4.0' or read('DATASET_CHECK.json')['label_members_in_archive'] != 480:
        raise ValueError('Dataset terms/count mismatch')
    dscheck = read('DATASET_CHECK.json')
    if hashlib.sha256((ROOT / dscheck['member_locator']).read_bytes()).hexdigest() != dscheck['member_index_sha256'] or len(read(dscheck['member_locator'])['labels']) != 480:
        raise ValueError('Dataset member index drift')
    for name in ['SOURCE_MANIFEST.json', 'SOURCE_MANIFEST_REFERENCES.json', 'SOURCE_MANIFEST_GUIDE_AND_DATASET.json']:
        d = read(name)
        rows = d['files'] if isinstance(d, dict) else d
        for e in rows:
            if hashlib.sha256((ROOT / e['local']).read_bytes()).hexdigest() != e['sha256']:
                raise ValueError('Source drift: ' + e['local'])
    validation = {'status': 'PASS_COMPUTATIONAL_DELIVERY', 'checked_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'required_artifacts': len(required), 'numeric_controls': len(controls), 'R2_controls': len(read('RESULTS_R2.json')['checks']), 'guide_controls': len(guides['controls']), 'vacuous_controls_rejected': True, 'verified_cited_DOIs': len(cited), 'CLAIM_items_mapped': 44, 'physical_facit': 'UNKNOWN', 'scientific_admission': False, 'clinical_recommendation': False, 'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    (ROOT / 'VALIDATION.json').write_text(json.dumps(validation, indent=2) + '\n')
    r['reference_interval_construction'] = {k: v for (k, v) in read('RESULTS_R2.json').items() if k != 'checks'}
    r['guide_parameter_replay'] = {'controls': len(guides['controls']), 'residual_clearance_mm': guides['residual_clearance_mm'], 'nominal_comparator_mm': guides['nominal_comparator_mm'], 'physical_certificates': 0, 'locator': 'raw/GUIDE_PARAMETER_REPLAY.json'}
    r['manuscript_status'] = 'RESEARCH_DRAFT_COMPUTATION_REPRODUCIBLE_PHYSICAL_AND_REPORTING_GAPS_RETAINED'
    r['delivery_validation'] = {k: v for (k, v) in validation.items() if k != 'checked_utc'}
    (ROOT / 'results.json').write_text(json.dumps(r, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    (ROOT / 'CURRENT_WORK_STATE.json').write_text(json.dumps({'lane': 'X41-manuscript-decision', 'status': 'MANUSCRIPT_AND_TOOL_DELIVERED_PENDING_REVIEW', 'updated_utc': validation['checked_utc'], 'latest_gate': validation, 'next_operation': 'Independent manuscript review; freeze an actual segmentation model before matched anatomical acquisition', 'physical_status': 'UNKNOWN', 'human_input_required': False, 'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW'}, indent=2) + '\n')
    print(json.dumps(validation))
if __name__ == '__main__':
    main()
