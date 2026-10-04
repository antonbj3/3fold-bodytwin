"""Verify files actually consumed at their current package paths; all large artifacts hash-bound."""
import argparse, json
from pathlib import Path
from identity import ROOT, sha, write

def verify(package):
    manifest = json.loads((package / 'DELIVERY_MANIFEST.json').read_text())
    bad = []
    for row in manifest['files']:
        p = Path(row['path'])
        p = p if p.is_absolute() else package / p
        if not p.is_file() or sha(p) != row['sha256']:
            bad.append(str(p))
    return bad

def seal():
    files = []
    selected = []
    for folder in ['code', 'sources', 'figures', 'raw', 'tables']:
        selected.extend((p for p in (ROOT / folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and (p.suffix not in ['.log'])))
    selected.extend((p for p in ROOT.iterdir() if p.is_file() and p.name not in ['DELIVERY_MANIFEST.json', 'CURRENT_WORK_STATE.json', 'LAST_RUN.json', 'GRAPH_FEEDBACK.json', 'GRAPH_FEEDBACK_RECEIPT.json', 'results.sha256'] and (p.suffix in ['.json', '.sha256', '.sh', '.md'])))
    for p in sorted(set(selected)):
        if p.name in ['LOCAL_GUARD_CONTROLS.json', 'REPLAY_COMPARISON.json', 'FULL_PIPELINE_COST.json']:
            continue
        files.append({'path': str(p.relative_to(ROOT)), 'sha256': sha(p), 'bytes': p.stat().st_size})
    for name in ['R1_MANIFEST.json', 'R2_MANIFEST.json', 'R3_MANIFEST.json', 'EXPORT_RECEIPTS.json']:
        for item in json.loads((ROOT / 'raw' / name).read_text()):
            p = Path(item['path'])
            files.append({'path': str(p), 'sha256': sha(p), 'bytes': p.stat().st_size})
            side = p.with_suffix('.frame.json')
            if side.exists():
                files.append({'path': str(side), 'sha256': sha(side), 'bytes': side.stat().st_size})
    for item in json.loads((ROOT / 'results.json').read_text()).get('own_data_artifact_manifest', []):
        p = Path(item['path'])
        files.append({'path': str(p), 'sha256': sha(p), 'bytes': p.stat().st_size})
    for item in json.loads((ROOT / 'FROZEN_COHORT.json').read_text())['cohort']:
        from identity import source_path
        p = source_path(item['source_path'])
        files.append({'path': str(p), 'sha256': sha(p), 'bytes': p.stat().st_size})
    write(ROOT / 'DELIVERY_MANIFEST.json', {'schema': 'x95-current-package-local-and-external-artifact-manifest-v1', 'files': files})
    (ROOT / 'results.sha256').write_text(sha(ROOT / 'results.json') + '\n')

def self_test():
    assert not verify(ROOT)
    path = ROOT / 'raw/R1_ROWS.json'
    original = path.read_bytes()
    rows = json.loads(original)
    i = next((i for (i, r) in enumerate(rows) if r.get('status') == 'MEASURED_DIGITAL_LABEL_GEOMETRY'))
    rows[i]['cut_center_min_mm'] = 99.0
    try:
        write(path, rows)
        bad = verify(ROOT)
        rejected = str(path) in bad
    finally:
        path.write_bytes(original)
    results = {'valid_current_package_pass': not verify(ROOT), '99mm_corrupted_local_result_rejected': rejected, 'restored_current_package_pass': not verify(ROOT), 'guard_reads_current_work_copy': True, 'X12_source_status': 'NOT_READY_UNCHANGED'}
    write(ROOT / 'raw/LOCAL_GUARD_CONTROLS.json', results)
    assert all((results[k] for k in ['valid_current_package_pass', '99mm_corrupted_local_result_rejected', 'restored_current_package_pass']))
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--seal', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    a = parser.parse_args()
    if a.seal:
        seal()
    if a.self_test:
        self_test()
    bad = verify(ROOT)
    print(json.dumps({'pass': not bad, 'failed': bad}, indent=2))
    raise SystemExit(bool(bad))
