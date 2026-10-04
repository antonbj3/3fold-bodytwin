"""Check manuscript outputs against evidence and exercise real wrong-value failures."""
import copy
import datetime
import hashlib
import importlib.metadata
import json
import re
import sys
from pathlib import Path
P = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    r = json.loads((P / 'results.json').read_text())
    v = json.loads((P / 'VALIDATION.json').read_text())
    manuscript = (P / 'MANUSCRIPT.md').read_text()
    controls = []

    def exercise(name, predicate, good, bad):
        passes = bool(predicate(good))
        rejects = not bool(predicate(bad))
        controls.append(dict(name=name, passes=passes, injected_wrong_value_rejected=rejects))
        if not passes or not rejects:
            raise AssertionError(name)
    core = '160/192 (0.833; 95% bootstrap interval 0.781–0.885)'
    exercise('Abstract primary estimate', lambda s: core in s, manuscript, manuscript.replace(core, core.replace('0.833', '0.933')))
    core = 'The paired difference was 0.078 (0.016–0.141).'
    exercise('Abstract paired contrast', lambda s: core in s, manuscript, manuscript.replace(core, core.replace('0.078', '0.178')))
    core = '7/72 (9.7%)'
    exercise('Abstract corrected report denominator', lambda s: core in s, manuscript, manuscript.replace(core, '7/989 (9.7%)'))
    for name in ['table1_flow', 'table2_arch_profile', 'table3a_opposition', 'table3b_corrected_crossbite', 'table4_report_conflicts']:
        text = (P / 'tables' / (name + '.md')).read_text().strip()
        exercise('Exact table insertion ' + name, lambda s: text in s, manuscript, manuscript.replace(text, text.replace('|', '!', 1)))
    checklist = json.loads((P / 'CHECKLIST.json').read_text())
    exercise('Checklist completeness of mapping', lambda rows: len(rows) == 78 and len({(q['scheme'], q['item']) for q in rows}) == 78, checklist['items'], checklist['items'][:-1])
    exercise('No compliance promotion', lambda q: q['formal_compliance_claim'] is False, checklist, dict(checklist, formal_compliance_claim=True))
    exercise('No clinical/publication promotion', lambda q: q['publication_ready'] is False, r, dict(r, publication_ready=True))
    exercise('No unverified license promotion', lambda s: 'unconfirmed designation' in s and 'exact Bite2Text grant was not verified' in s, manuscript, manuscript.replace('unconfirmed designation', 'verified dataset license'))
    for figure in ['figure1_case_flow', 'figure2_agreement_and_conflict']:
        for ext in ['png', 'pdf', 'svg']:
            f = P / 'figures' / (figure + '.' + ext)
            assert f.exists() and f.stat().st_size > 1000
    refs = json.loads((P / 'sources/DOI_VERIFICATION.json').read_text())['records']
    for ref in refs:
        doi = ref.get('doi', ref.get('supplied_doi'))
        if doi:
            exercise('Reference DOI ' + ref['name'], lambda s: doi in s, manuscript, manuscript.replace(doi, '10.0000/wrong'))
    for key in ['X7_R4', 'X21_R2', 'X21_R5']:
        exercise('Retain failed gate ' + key, lambda q: q[key] == 'FAIL', r['original_gates'], dict(r['original_gates'], **{key: 'PASS'}))
    versions = {'python': sys.version.split()[0]}
    for pkg in ['numpy', 'scipy', 'scikit-learn', 'trimesh']:
        try:
            versions[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            versions[pkg] = 'NOT_INSTALLED'
    versions['matplotlib_provider'] = '../LANE_X7_BITE2TEXT/code/vendor_plotting'
    versions['threads_max'] = 4
    (P / 'raw/SOFTWARE_VERSIONS.json').write_text(json.dumps(versions, indent=2) + '\n')
    v['package_checks'] = controls
    v['package_pass'] = all((q['passes'] and q['injected_wrong_value_rejected'] for q in controls))
    (P / 'VALIDATION.json').write_text(json.dumps(v, indent=2) + '\n')
    (P / 'GRAPH_FEEDBACK.json').write_text(json.dumps(dict(schema='X39-manuscript-feedback-v1', lane='X39-manuscript-bite', review_state='PENDING_INDEPENDENT_REVIEW', target_id=None, proposed_target_id='DENT-VAL-BITE2TEXT-REPORT-AGREEMENT', coverage_status='MISSING_MATCHING_CANONICAL_TARGET; exact packet request failed', result_path=str(P / 'results.json'), result_sha256=sha(P / 'results.json'), manuscript_path=str(P / 'MANUSCRIPT.md'), manuscript_sha256=sha(P / 'MANUSCRIPT.md'), claim_type=['information_link', 'capability'], measured_quantity='Existing report category agreement and within-case report conflict', units='dimensionless proportions', resolution_level='PER_ARCH -> POPULATION', time_scale='SIMULTANEOUS', population_regime='989 local report-eligible cases; 198 reused internal-test cases; field-specific subsets; no external clinical testing', uncertainty='Existing fixed-model patient-case bootstrap and scoped Wilson intervals; no anatomy/contact uncertainty', frozen_gate=dict(reporting_contract='PREREG_REPORTING.json', sha256=sha(P / 'PREREG_REPORTING.json'), numeric_tolerance=1e-12), baseline='Scan-free training-majority proxy; paired existing X7 system for descriptive comparison; matched HGB/conformal equality retained without superiority claim', outcome='MANUSCRIPT_DRAFT_COMPLETE; category/source limits retained; publication-ready FALSE', negative_result=True, negative_scope='X21 whole-method and X7 joint-informativeness failures preserved; clinical diagnosis/contact/landmark validity not established', next_operation='Independent clinician adjudication and new cohort / target point and loaded-contact measurements from already frozen parent protocols', proposed_edges=[dict(from_artifact='../LANE_X7_BITE2TEXT/results.json', to_proposal='DENT-VAL-BITE2TEXT-REPORT-AGREEMENT', quantity='Arch category / readable report version', resolution_level='PER_ARCH', time_scale='SIMULTANEOUS'), dict(from_artifact='../LANE_X21_CONTACT_MAP/raw/COMPARISON_ROWS_R2.json', to_proposal='DENT-VAL-BITE2TEXT-REPORT-AGREEMENT', quantity='First-readable report category linked to geometric assessment', resolution_level='PER_ARCH', time_scale='SIMULTANEOUS')], dispatch_status='NOT_DISPATCHED_NO_SUITABLE_TARGET', scientific_admission=False), indent=2) + '\n')
    artifacts = []
    for path in sorted(P.rglob('*')):
        if path.is_file() and path.name not in ['ARTIFACT_MANIFEST.json', 'CURRENT_WORK_STATE.json'] and ('__pycache__' not in path.parts) and ('logs' not in path.parts):
            artifacts.append(dict(path=str(path.relative_to(P)), bytes=path.stat().st_size, sha256=sha(path)))
    (P / 'ARTIFACT_MANIFEST.json').write_text(json.dumps(dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), files=artifacts), indent=2) + '\n')
    (P / 'CURRENT_WORK_STATE.json').write_text(json.dumps(dict(lane='X39-manuscript-bite', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='MANUSCRIPT_COMPLETE_PENDING_INDEPENDENT_REVIEW', latest_gate=f"{len(v['checks'])} evidence checks + {len(controls)} package checks passed; all wrong values rejected", next_operation='Independent clinician/reference review and matched new-cohort/landmark acquisition; submission declarations still missing', primary_artifact='MANUSCRIPT.md', publication_ready=False), indent=2) + '\n')
    print(f"Package PASS: {len(v['checks'])} evidence and {len(controls)} output checks; all injected wrong values rejected.")
if __name__ == '__main__':
    main()
