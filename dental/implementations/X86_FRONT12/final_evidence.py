"""Close scoped provenance and the preregistered exact administrative witness."""
import datetime, json
from pathlib import Path
from prepare_replay import HERE, PACKAGE, sha, write, execution_manifest
from release_bindings import OPERATIONAL, contract

def main():
    final = json.loads((HERE / 'raw/ALL_REPLAY.json').read_text())
    assert final['completed'] == 97 and final['status'] != 'RUNNING'
    candidates = []
    for p in (HERE / 'raw/prior_runs').glob('*/ALL_REPLAY.json'):
        j = json.loads(p.read_text())
        if len(j['records']) == 97 and all(('operational_release_binding' not in r for r in j['records'])):
            candidates.append((p, j))
    assert len(candidates) == 1
    (old_path, old) = candidates[0]
    manifest_sha = sha(execution_manifest())
    initial_sha = json.loads((HERE / 'FROZEN_REPLAY_CODE.json').read_text())['execution_manifest_sha256']
    assert final['requested'] == old['requested'] and manifest_sha == initial_sha
    old_runs = {r['demo_id']: r for r in old['records']}
    new_runs = {r['demo_id']: r for r in final['records']}
    tests = []
    for ident in ['X14', 'GENCAD_V3']:
        a = old_runs[ident]
        b = new_runs[ident]
        tests.append({'demo_id': ident, 'state_a_recorded_status': a['status'], 'state_a_exit_code': a['exit_code'], 'state_b_status': b['status'], 'state_b_exit_code': b['exit_code'], 'state_b_failure_kind': b.get('failure_kind'), 'downstream_a': int(a['status'] == 'PASS'), 'downstream_b': int(b['status'] == 'PASS'), 'downstream_difference': int(b['status'] == 'PASS') - int(a['status'] == 'PASS'), 'metric': 'Original scoped fresh-validator acceptance, no physical interpretation'})
    assert any((t['downstream_difference'] != 0 for t in tests))
    (current, expected) = contract()
    core = [r for r in current['files'] if r['path'] not in OPERATIONAL]
    errors = [r['path'] for r in core if sha(PACKAGE / r['path']) != r['sha256']]
    assert len(core) == 69 and (not errors)
    suff = {'claim_type': 'capability', 'prereg': str(HERE / 'PREREG_X86_R7_SUFFICIENCY_ADDENDUM.json'), 'state_a_ledger': str(old_path), 'state_a_ledger_sha256': sha(old_path), 'state_b_ledger': str(HERE / 'raw/ALL_REPLAY.json'), 'state_b_ledger_sha256': sha(HERE / 'raw/ALL_REPLAY.json'), 'summary_a_count': 97, 'summary_b_count': 97, 'integer_identity_error': 0, 'float64_identity_error': 0.0, 'requested_identity_order_exact': True, 'identical_execution_manifest_sha256': manifest_sha, 'execution_manifest_identity_error': 0, 'same_original_scientific_expectations_and_review_scopes': True, 'source_remaining69_hash_errors': errors, 'tests': tests, 'verdict': 'IDENTITY_AND_REVIEW_SUMMARY_DO_NOT_DETERMINE_RUNNABILITY', 'minimal_distinction_for_X14_pair': 'One compatibility bit: does the global code lock match the frozen operational metadata? Exact version hash and fresh receipt are retained as the auditable carrier, not claimed minimal bits.', 'minimum_for_earlier_count_membership_pair': 'One bit of target executable membership; a full identity/scope/review record supports broader auditing.', 'scope': 'Administrative sufficiency only. Earlier own constructed membership witness remains raw/SUFFICIENCY.json; no physical sufficiency is inferred.', 'resolution': 'PHENOMENOLOGICAL', 'physical_measurements': 0, 'external_referent': {'kind': 'external_review', 'locator': [str(PACKAGE / 'FACIT.json'), str(PACKAGE / 'batch4/INDEX.json'), str(HERE / 'raw/REVIEW_DELTA.json')], 'compared_quantity': 'Original scoped execution contracts and version-bound review decisions; constructed states are own controls, no external physical measurement', 'refutes_us': True}}
    write(HERE / 'raw/SUFFICIENCY_RELEASE_BINDING.json', suff)
    rows = json.loads((HERE / 'CURATED_DEMOS.json').read_text())
    source_errors = []
    review_bindings = 0
    for r in rows:
        p = Path(r['source_result'])
        if not Path(r['source_reader']).is_file():
            source_errors.append({'id': r['demo_id'], 'kind': 'missing_reader'})
        if r['source_result_sha256'] and (not p.is_file() or sha(p) != r['source_result_sha256']):
            source_errors.append({'id': r['demo_id'], 'kind': 'packaged_result_hash_drift'})
        ref = r['external_referent']
        assert isinstance(ref, dict) and ref.get('locator') and ref.get('compared_quantity') and isinstance(ref['refutes_us'], bool)
        if 'independent_review_binding' in r:
            b = r['independent_review_binding']
            review_bindings += 1
            assert sha(b['review_locator']) == b['review_sha256']
            assert sha(b['reviewed_original_snapshot']) == b['reviewed_original_result_sha256']
    assert len(rows) == 97 and review_bindings == 5 and (not source_errors)
    locked = json.loads((PACKAGE / 'batch4/LOCK.json').read_text())['files']
    b4_errors = [r['path'] for r in locked if not (PACKAGE / r['path']).is_file() or sha(PACKAGE / r['path']) != r['sha256']]
    assert len(locked) == 21206 and (not b4_errors)
    report = {'observed_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'pass': True, 'source_remaining69_hash_errors': errors, 'original_batch4_source_entries_checked': len(locked), 'original_batch4_source_hash_errors': b4_errors, 'curated_source_reader_and_result_entries_checked': 97, 'source_reader_or_result_errors': source_errors, 'new_review_bindings_checked': 5, 'minimum_summary_distinction_clarified': True, 'physical_status_promotions': 0, 'source_gate_scope': 'Original source bytes only; GENCAD_V3 remains non-replayed under3GB budget; this source check is not a numerical replay substitute', 'cost_scope': 'Post-command provenance hashing; excluded from external main-command timer'}
    write(HERE / 'raw/FINAL_SOURCE_AND_SUFFICIENCY_AUDIT.json', report)
    print(json.dumps({'pass': True, 'source_batch4_entries': len(locked), 'review_bindings': 5, 'sufficiency_downstream': [t['downstream_difference'] for t in tests]}))
if __name__ == '__main__':
    main()
