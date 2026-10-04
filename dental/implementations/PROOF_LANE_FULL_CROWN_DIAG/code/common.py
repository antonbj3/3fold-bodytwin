from dental_release.paths import expand as _release_expand
from pathlib import Path
import os, sys, json, hashlib, datetime, time
os.environ.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMEXPR_NUM_THREADS='4', MPLBACKEND='Agg')
sys.dont_write_bytecode = True
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
BASE = Path(os.environ.get('DENTAL_REVIEW_SOURCE_RESULTS', _release_expand('@DENTAL_INPUT_ROOT@/artifacts')))
DATA = ROOT / 'review_data'
DATA.mkdir(parents=True, exist_ok=True)
PUBLIC = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_FULL_CROWN/public'))
R3 = BASE / 'PROOF_LANE_FULL_CROWN_R3'
V4 = BASE / _release_expand('GENCAD_V4')
PYTHON = _release_expand('@DENTAL_PYTHON@')

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def read(p):
    return json.loads(Path(p).read_text())

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
        return None
    return x

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(exist_ok=True, parents=True)
    t = p.with_suffix(p.suffix + '.tmp')
    t.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')
    t.replace(p)

def freeze(p, x):
    p = Path(p)
    if p.exists():
        raise ValueError('Already frozen: ' + str(p))
    dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')

def npz(p):
    with np.load(p, allow_pickle=False) as z:
        return dict(z)

def state(phase, gate, nextop):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-full-crown-diag', claim_type=['capability', 'information_link'], phase=phase, latest_gate=gate, next_operation=nextop, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW'))
