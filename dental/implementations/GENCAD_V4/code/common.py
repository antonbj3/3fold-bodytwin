import os, sys, json, hashlib, datetime, time
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'payload'

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def clean(x):
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, np.generic):
        return clean(x.item())
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, float) and (not np.isfinite(x)):
        raise ValueError('Nonfinite report value')
    return x

def read(p):
    p = Path(p)
    if (ROOT / 'RELEASE_LOCK.json').exists():
        from integrity import checked_json_bytes
        try:
            p.absolute().relative_to(ROOT)
        except ValueError:
            return json.loads(p.read_bytes())
        return json.loads(checked_json_bytes(p))
    return json.loads(p.read_bytes())

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    q = p.with_suffix(p.suffix + '.tmp')
    q.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')
    q.replace(p)

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def freeze(p, x):
    p = Path(p)
    if p.exists():
        raise RuntimeError('Refusing to replace frozen file ' + str(p))
    dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')

def state(phase, gate, next_operation, **extra):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-gencad-v4', phase=phase, latest_gate=gate, next_operation=next_operation, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW', **extra))

def inventory(base):
    return {str(p.relative_to(base)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(Path(base).rglob('*')) if p.is_file() and '__pycache__' not in str(p)}
