import os, sys, json, hashlib, datetime, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
PAYLOAD = ROOT / 'payload'
FAMILIES = ['molar_crown', 'premolar_crown', 'anterior_crown', 'bridge3', 'implant_crown', 'inlay', 'onlay', 'veneer', 'lattice_onlay']
LEVELS = ['easy', 'normal', 'hard', 'boundary']

def clean(x):
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, np.generic):
        return clean(x.item())
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (tuple, list)):
        return [clean(v) for v in x]
    if isinstance(x, float) and (not np.isfinite(x)):
        return None
    return x

def canonical(x):
    return json.dumps(clean(x), sort_keys=True, separators=(',', ':'), allow_nan=False).encode()

def digest(x):
    return hashlib.sha256(canonical(x)).hexdigest()

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def read(p):
    return json.loads(Path(p).read_text())

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')
    tmp.replace(p)

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def state(phase, gate, next_operation, **kw):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-gencad-v3', phase=phase, latest_gate=gate, next_operation=next_operation, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW', **kw))

def freeze(p, data):
    if Path(p).exists():
        obj = read(p)
        if obj['payload_sha256'] != digest(data):
            raise ValueError('immutable freeze differs: ' + str(p))
        return obj
    obj = dict(frozen_utc=now(), payload_sha256=digest(data), payload=data)
    dump(p, obj)
    return obj

def budget():
    base = PAYLOAD.resolve().parent
    size = sum((p.stat().st_size for p in base.rglob('*') if p.is_file()))
    if size > 3000000000:
        raise RuntimeError('lane exceeds 3 GB: ' + str(size))
    return size
