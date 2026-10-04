"""Bind every retained attempt to its canonical historical paths without rewriting it."""
import json
from pathlib import Path
from prepare_replay import HERE, sha, write
rows = []
paths = sorted((HERE / 'raw/prior_runs').glob('*/ALL_REPLAY.json'))
if (HERE / 'raw/ALL_REPLAY.json').exists():
    paths.append(HERE / 'raw/ALL_REPLAY.json')
for p in paths:
    j = json.loads(p.read_text())
    base = p.parent
    records = []
    for r in j['records']:
        ident = r['demo_id']
        m = base / 'output_manifests' / f'{ident}.json'
        receipt = base / 'receipts' / f'{ident}.json'
        log = base / 'replay_logs' / f'{ident}.log' if base != HERE / 'raw' else HERE / 'replay_logs' / f'{ident}.log'
        records.append({'id': ident, 'recorded_status': r['status'], 'exit_code': r['exit_code'], 'wall_seconds': r['wall_seconds'], 'canonical_receipt': str(receipt), 'receipt_sha256': sha(receipt) if receipt.exists() else None, 'canonical_output_manifest': str(m), 'output_manifest_sha256': sha(m) if m.exists() else None, 'canonical_stdout_stderr': str(log), 'log_sha256': sha(log) if log.exists() else None})
    rows.append({'ledger': str(p), 'ledger_sha256': sha(p), 'requested': j['requested'], 'status': j['status'], 'wall_seconds': j.get('command_total_wall_seconds', j.get('wall_seconds', j.get('wall_seconds_so_far'))), 'recorded_pass_count': j.get('pass_count'), 'records': records, 'path_rule': 'Old embedded raw/log absolute paths refer to their own attempt. Canonical paths here retain the byte-exact original records; physical data-artifact locators remain unchanged.'})
write(HERE / 'raw/ATTEMPT_CANONICAL_INDEX.json', {'attempts': rows, 'original_ledgers_rewritten': False})
print(len(rows), 'attempts indexed')
