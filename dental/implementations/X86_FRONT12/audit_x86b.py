"""Verify delivered reader, actual receipts, review bindings and retained files."""
import hashlib
import json
import os
import re
import time
from pathlib import Path
from prepare_replay import HERE, sha, write
from receipt_contract import validator_pass
from recover_receipt_bytes import recover

def read(p):
    return json.loads(Path(p).read_text())

def main():
    started = time.monotonic()
    frozen = read(HERE / 'FROZEN_FRONT_X86B.json')
    ledger = read(HERE / 'raw/x86b/FRONT_REPLAY.json')
    result = read(HERE / 'results.json')
    front = (HERE / 'FRONT.md').read_text()
    all_demos = (HERE / 'ALL_DEMOS.md').read_text()
    ids = re.findall('^\\| `\\./run_demos.sh ([^`]+)`', front, re.M)
    assert ids == frozen['front_ids'] and ledger['requested'] == ids and (ledger['completed'] == 12)
    assert ledger['pass_count'] == 12 and ledger['failure_count'] == 0
    assert result['scientific_status_promotions'] == result['physical_measurements'] == 0
    assert result['status_counts'] == {'PRESENT': 7, 'PARTIAL': 39, 'REFUTED': 2, 'RESEARCH': 3, 'MISSING': 1}
    assert len(re.findall('^\\| `', all_demos, re.M)) == 101
    assert result['command_wall_seconds'] == ledger['wall_seconds']
    assert f"{ledger['wall_seconds']:.3f} s" in front and 'IN_PROGRESS' not in front
    assert 'R7 exact summary FAIL' in front and 'target_id=null' in front
    assert 'No paired improvement against R1–R3' in front and 'REJECT / NOT_READY' in front
    assert 'no new rerun' in front and result['eligible_hash_only_count'] == 88
    code = read(HERE / 'FROZEN_CODE_X86B.json')
    assert ledger['runtime_code_sha256'] == sha(HERE / 'FROZEN_CODE_X86B.json')
    for (name, want) in code['files'].items():
        assert sha(HERE / name) == want, 'Executed code drift ' + name
    snapshots = []
    for r in ledger['records']:
        assert r['exit_code'] == 0 and validator_pass(r['validator_receipt'], r['demo_id'])
        assert r['mismatch_count'] == 0 and r['mismatches'] == []
        source = Path(r['validator_path'])
        live_matches = sha(source) == r['validator_sha256']
        raw = source.read_bytes() if live_matches else recover(r['validator_receipt'], r['validator_sha256'])
        assert raw is not None, 'Missing exact original receipt bytes ' + r['demo_id']
        target = HERE / 'raw/x86b/original_receipts' / f"{r['demo_id']}.json"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
        assert sha(target) == r['validator_sha256']
        snapshots.append({'demo_id': r['demo_id'], 'path': str(target), 'sha256': r['validator_sha256'], 'original_latest_locator': str(source), 'latest_still_matches': live_matches, 'mode': 'copy_exact_original_bytes' if live_matches else 'recover_only_if_independently_captured_original_SHA256_matches'})
        assert Path(r['log']).is_file()
        assert f"{r['wall_seconds']:.3f} s" in front
    for hashes in [ledger['hashes_before'], ledger['hashes_after']]:
        assert hashes['checked'] == hashes['pass_count'] == 101 and hashes['mismatch_count'] == 0
        for r in hashes['records']:
            assert sha(r['result_path']) == r['expected_sha256'] == r['observed_sha256']
            assert all((c['pass'] and sha(c['path']) == c['expected'] for c in r['binding_checks']))
    reviews = []
    for row in frozen['rows'][-3:] + frozen['rejections']:
        b = row['independent_review_binding']
        review = read(b['review_locator'])
        assert sha(b['review_locator']) == b['review_sha256']
        assert sha(b['reviewed_original_snapshot']) == b['reviewed_original_result_sha256']
        assert review['result_sha256'] == b['reviewed_original_result_sha256']
        assert review['demo_decision'] == row['delivery_review_status']
        reviews.append({'demo_id': row['demo_id'], 'review': b['review_locator'], 'sha256': b['review_sha256'], 'decision': review['decision'], 'target_id': review.get('target_id'), 'demo_decision': review['demo_decision']})
    hunt4 = next((r for r in reviews if r['demo_id'] == 'LANE_XBREAK_HUNT_4'))
    assert hunt4['target_id'] is None
    suff = ledger['faults']['sufficiency']
    assert suff['summary_a_sha256'] == suff['summary_b_sha256']
    assert suff['summary_byte_identity_error'] == suff['float64_identity_error'] == 0
    assert suff['downstream_difference'] == 1
    assert len(ledger['faults']['receipt_faults']) == 12 and len(ledger['faults']['hash_faults']) == 101
    assert all((all((r[k] for k in ['bad_exit_rejected', 'bad_receipt_rejected', 'wrong_identity_rejected'])) for r in ledger['faults']['receipt_faults']))
    arrays = []
    run = Path(ledger['run_directory'])
    for (base, dirs, files) in os.walk(run, followlinks=False):
        for name in files:
            p = Path(base) / name
            if p.is_symlink():
                continue
            n = p.stat().st_size
            if n > 50000000:
                arrays.append({'path': str(p), 'bytes': n, 'sha256': sha(p)})
    write(HERE / 'raw/x86b/LARGE_ARTIFACTS.json', arrays)
    write(HERE / 'raw/x86b/ORIGINAL_RECEIPT_BYTES.json', snapshots)
    broken = []
    for filename in ['FRONT.md', 'ALL_DEMOS.md', 'README_DEMO.md']:
        text = (HERE / filename).read_text()
        for target in re.findall('\\]\\(([^)]+)\\)', text):
            if target.startswith(('https://', 'http://', 'app://')):
                continue
            target = target.split(':')[0]
            if not (HERE / target).exists():
                broken.append({'file': filename, 'target': target})
    assert not broken, broken
    audit = {'pass': True, 'front_ids_exact': ids, 'fresh_passes': 12, 'scientific_mismatches': 0, 'result_hashes_verified': 101, 'new_review_bindings': reviews, 'unchanged_scientific_status_counts': True, 'no_physical_status_promotions': True, 'executed_code_hashes_verified': len(code['files']), 'fresh_receipt_hashes_verified': 12, 'receipt_mutants_rejected': 36, 'result_hash_mutants_rejected': 101, 'summary_identity_error': 0, 'downstream_difference': 1, 'broken_reader_links': broken, 'large_artifacts': arrays, 'original_receipt_snapshots': snapshots, 'overwritten_latest_aliases': sum((not r['latest_still_matches'] for r in snapshots)), 'audit_wall_seconds': time.monotonic() - started, 'cost_scope': 'Post-replay delivered-reader/byte audit; excluded from bounded12 command timer', 'independent_scientific_review': False}
    write(HERE / 'FINAL_AUDIT_X86B.json', audit)
    print(json.dumps({k: audit[k] for k in ['pass', 'fresh_passes', 'result_hashes_verified', 'receipt_mutants_rejected', 'result_hash_mutants_rejected', 'broken_reader_links', 'audit_wall_seconds']}))
if __name__ == '__main__':
    main()
