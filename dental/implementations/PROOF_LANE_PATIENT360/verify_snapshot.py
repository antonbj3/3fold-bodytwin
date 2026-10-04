"""Check the fixed local dependency/code snapshot before and after replay."""
from pathlib import Path
import hashlib, json, sys
R = Path(__file__).resolve().parent
manifest = json.loads((R / 'REPLAY_LOCK.json').read_text())
for row in manifest['files']:
    p = R / row['path']
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest() != row['sha256']:
        raise SystemExit('Snapshot mismatch: ' + row['path'])
print('Fixed snapshot verified:', len(manifest['files']), 'files')
