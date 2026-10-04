from dental_release.paths import expand as _release_expand
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[key] = '4'
import json, hashlib, datetime, sys
sys.dont_write_bytecode = True
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
DENTAL = ROOT.parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X53-milling-tools'))

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
        return None
    return x

def dump(path, x):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')
    tmp.replace(path)

def read(path):
    return json.loads(Path(path).read_text())

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def digest(x):
    return hashlib.sha256(json.dumps(clean(x), sort_keys=True, separators=(',', ':')).encode()).hexdigest()

def freeze(path, payload):
    path = Path(path)
    obj = dict(timestamp_utc=now(), payload=payload, payload_sha256=digest(payload))
    if path.exists():
        old = read(path)
        if old['payload_sha256'] != digest(old['payload']):
            raise ValueError('freeze corrupt: ' + str(path))
        return old
    dump(path, obj)
    path.with_suffix('.sha256').write_text(sha(path) + '\n')
    return obj

def state(phase, gate, next_operation):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X53-milling-tools', claim_type='capability', phase=phase, latest_gate=gate, next_operation=next_operation, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW'))

def budget():
    n = sum((p.stat().st_size for root in [ROOT, DATA] for p in root.rglob('*') if p.is_file() and (not p.is_symlink())))
    if n > 3000000000:
        raise RuntimeError('3 GB lane limit')
    return n
