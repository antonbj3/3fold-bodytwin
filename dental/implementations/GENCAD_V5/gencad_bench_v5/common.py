from dental_release.paths import expand as _release_expand
import os, json, hashlib, datetime, time, sys
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
V4 = ROOT.parent / _release_expand('GENCAD_V4')
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V5'))
PYTHON = _release_expand('@DENTAL_PYTHON@')

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
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
        return None
    return x

def read(p):
    return json.loads(Path(p).read_text())

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def freeze(p, x):
    p = Path(p)
    if p.exists():
        raise ValueError('Frozen file exists: ' + str(p))
    dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')

def state(phase, gate, next_operation):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-gencad-v5', claim_type='capability', phase=phase, latest_gate=gate, next_operation=next_operation, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW'))

def npz(p):
    with np.load(p, allow_pickle=False) as a:
        return dict(a)
