from dental_release.paths import expand as _release_expand
import os, sys, json, hashlib, datetime, resource
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', MPLBACKEND='Agg')
sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_AS, (3500 * 1024 ** 2, 3500 * 1024 ** 2))
import numpy as np
SOURCE_ROOT = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('PROOF_LANE_INVERSE_RUN_ROOT', str(SOURCE_ROOT)))
BASE = SOURCE_ROOT.parent
DATA = Path(os.environ.get('PROOF_LANE_INVERSE_RUN_DATA', _release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_INVERSE_CROWN')))
DATA.mkdir(exist_ok=True)
X60 = BASE / _release_expand('X60')
X1 = BASE / _release_expand('X1B')

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def read(p):
    return json.loads(Path(p).read_text())

def clean(x):
    if isinstance(x, (np.ndarray, np.generic)):
        return clean(x.tolist())
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    return x

def save(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def freeze(p, x):
    p = Path(p)
    if p.exists():
        raise ValueError('Refuse frozen overwrite ' + str(p))
    save(p, dict(frozen_utc=now(), **x))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')

def state(phase, gate, nxt):
    save(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-inverse-crown', phase=phase, latest_gate=gate, next_operation=nxt, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW'))
