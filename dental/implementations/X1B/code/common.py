from pathlib import Path
import json, hashlib, datetime, os
R = Path(__file__).resolve().parents[1]

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def read(p):
    return json.loads((R / p).read_text())

def dump(p, d):
    p = R / p
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(d, indent=2, allow_nan=False, default=lambda x: x.item() if hasattr(x, 'item') else x.tolist()) + '\n')

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def frozen(p, d):
    p = R / p
    b = (json.dumps(d, indent=2, allow_nan=False) + '\n').encode()
    if p.exists():
        if json.loads(p.read_bytes()) != d:
            raise ValueError('Frozen artifact already exists and differs: ' + p.name)
    else:
        p.write_bytes(b)
    p.with_suffix('.sha256').write_text(sha(p) + '\n')

def state(milestone, gate, next_operation):
    dump('CURRENT_WORK_STATE.json', dict(round_tag='X1b-crown-loop', status='ACTIVE', updated_utc=now(), milestone=milestone, latest_gate=gate, next_operation=next_operation, physical_measurements_performed=False, review_state='PENDING_INDEPENDENT_REVIEW'))
