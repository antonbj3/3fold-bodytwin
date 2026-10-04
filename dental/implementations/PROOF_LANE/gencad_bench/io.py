import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'review_data'
DENTAL = next((p for p in ROOT.parents if p.name == 'dental' and (p / 'results').is_dir()), None)
if DENTAL is None:
    from dental_release.paths import MissingInput
    raise MissingInput('Shared dental inputs require the external source workspace; see docs/DATA.md')

def canonical(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()

def digest(x):
    return hashlib.sha256(canonical(x)).hexdigest()

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda : f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def dump(path, x):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    tmp.replace(path)

def now():
    return datetime.now(timezone.utc).isoformat()

def state(phase, gate, next_operation):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane=ROOT.name, phase=phase, latest_gate=gate, next_operation=next_operation, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW'))

def freeze(path, payload):
    path = Path(path)
    if path.exists():
        old = json.loads(path.read_text())
        if old['payload_sha256'] != digest(old['payload']):
            raise ValueError('Stored frozen payload is corrupt: ' + str(path))
        if old['payload_sha256'] != digest(payload):
            raise ValueError('Frozen payload changed: ' + str(path))
        return old
    out = dict(frozen_utc=now(), payload_sha256=digest(payload), payload=payload)
    dump(path, out)
    return out

def result(status, scope, **kw):
    if status not in {'PASS', 'FAIL', 'UNKNOWN', 'INVALID'}:
        raise ValueError(status)
    return dict(status=status, scope=scope, **kw)
