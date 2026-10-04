"""Audit the finished frozen receipts, original source bytes and archived locators."""
import json, os, datetime
from pathlib import Path
from prepare_replay import HERE, DATA, PACKAGE, sha, write, execution_manifest
from validate_delivery import complete_coverage, passed_run

def main():
    ledger = json.loads((HERE / 'raw/ALL_REPLAY.json').read_text())
    expected = [r['id'] for r in json.loads(execution_manifest().read_text())['executions']]
    assert ledger['status'] != 'RUNNING' and complete_coverage(ledger['records'], expected)
    frozen = json.loads((HERE / 'FROZEN_REPLAY_CODE_FINAL_V2.json').read_text())
    code_drift = [name for (name, h) in frozen['code_sha256'].items() if sha(HERE / name) != h]
    assert not code_drift
    source_metadata = []
    for r in json.loads((HERE / 'CLEAN_COPY_MANIFEST.json').read_text())['copied_files']:
        p = PACKAGE / r['relative_path']
        h = sha(p)
        if h != r['sha256']:
            source_metadata.append({'path': str(p), 'expected': r['sha256'], 'actual': h})
    seen = {}
    checked = 0
    bytes_checked = 0
    errors = []
    symlinks = 0
    manifests = list((HERE / 'raw/prior_runs').glob('*/output_manifests/*.json')) + list((HERE / 'raw/output_manifests').glob('*.json'))
    for m in manifests:
        for r in json.loads(m.read_text()):
            if not r.get('locator'):
                symlinks += 1
                continue
            p = Path(r['locator'])
            key = str(p)
            if key not in seen:
                if not p.is_file():
                    seen[key] = None
                else:
                    seen[key] = sha(p)
                    checked += 1
                    bytes_checked += p.stat().st_size
            if seen[key] != r['sha256']:
                errors.append({'manifest': str(m), 'locator': key, 'expected': r['sha256'], 'actual': seen[key]})
    invalid_positive = [r['demo_id'] for r in ledger['records'] if r['status'] == 'PASS' and (not passed_run(r))]
    hashes = []
    for p in HERE.glob('*.sha256'):
        content = p.read_text().split()
        name = content[1]
        candidate = HERE / name
        if p.name == 'PREREG_X86_REPLAY_R2_v1_WRONG_COUNT.sha256':
            candidate = HERE / 'PREREG_X86_REPLAY_R2_v1_WRONG_COUNT.json'
        if not candidate.exists():
            errors.append({'hash_record': str(p), 'missing': name})
        elif sha(candidate) != content[0]:
            errors.append({'hash_record': str(p), 'drift': name})
        hashes.append(name)
    report = {'completed_executables': len(expected), 'identity_coverage': True, 'valid_original_positive_receipts': not invalid_positive, 'invalid_positives': invalid_positive, 'frozen_runtime_code_drift': code_drift, 'source_metadata_drift_against_initial_snapshot': source_metadata, 'archive_locators_checked': checked, 'archive_bytes_hashed': bytes_checked, 'shared_symlink_records': symlinks, 'archive_or_freeze_hash_errors': errors, 'frozen_hash_records_checked': hashes, 'observed_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'source_live_metadata_unchanged': not source_metadata, 'source_change_scope': 'Live batch26 dispatcher/catalog arrived during the frozen batch<=25 test; recorded separately in LIVE_SHARED_SOURCE_RELEASE_DELTA.json. Our execution source mounts remain read-only.', 'pass': not errors and (not code_drift) and (not invalid_positive), 'scope': 'Administrative artifact and original-validator integrity; does not validate physical predictions or historic review semantics'}
    write(HERE / 'FINAL_ARTIFACT_AUDIT.json', report)
    print(json.dumps({k: report[k] for k in ['pass', 'completed_executables', 'archive_locators_checked', 'archive_bytes_hashed', 'archive_or_freeze_hash_errors', 'source_metadata_drift_against_initial_snapshot']}))
    if not report['pass']:
        raise SystemExit(1)
if __name__ == '__main__':
    main()
