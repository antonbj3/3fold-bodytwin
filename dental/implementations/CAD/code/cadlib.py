from dental_release.paths import expand as _release_expand
import os, sys, json, hashlib, datetime, importlib.util, time
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMEXPR_NUM_THREADS='4', MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1')
sys.dont_write_bytecode = True
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_CAD_10X'))
DATA.mkdir(parents=True, exist_ok=True)

def read(p):
    return json.loads(Path(p).read_text())

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def clean(x):
    if isinstance(x, (np.ndarray, np.generic)):
        return clean(x.tolist())
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (tuple, list)):
        return [clean(v) for v in x]
    if isinstance(x, float) and (not np.isfinite(x)):
        return None
    return x

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
        raise ValueError('Already frozen: ' + str(p))
    dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix('.sha256').write_text(sha(p) + '\n')

def module(name, path):
    s = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m
    s.loader.exec_module(m)
    return m
    def_unused = None

def state(phase, gate, nextop):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-cad-10x', phase=phase, latest_gate=gate, next_operation=nextop, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW', data_root=DATA, resource_contract={'threads_max': 4, 'intermediates_limit_bytes': 3000000000}))

def budget():
    n = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    if n > 3000000000:
        raise RuntimeError('3GB cap')
    return n
