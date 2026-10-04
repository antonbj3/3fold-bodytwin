from dental_release.paths import expand as _release_expand
import os, sys, json, hashlib, datetime, importlib.util, resource
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMEXPR_NUM_THREADS='4', IGL_NUM_THREADS='1', MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1')
sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_AS, (3500 * 1024 ** 2, 3500 * 1024 ** 2))
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
PKG = BASE / 'DEMO48_PACKAGE'
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_COMBINE_CROWN_CHAIN'))

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
    import numpy as np
    if isinstance(x, (np.ndarray, np.generic)):
        return clean(x.tolist())
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, float) and (not __import__('math').isfinite(x)):
        return None
    return x

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    t = p.with_suffix(p.suffix + '.tmp')
    t.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')
    t.replace(p)

def freeze(name, x):
    p = ROOT / name
    if p.exists():
        raise ValueError('Already frozen ' + name)
    dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')

def verify_freeze(p):
    p = Path(p)
    if sha(p) != p.with_suffix('.sha256').read_text().split()[0]:
        raise ValueError('Freeze drift ' + str(p))

def module(name, p):
    s = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m
    s.loader.exec_module(m)
    return m

def scoped_modules(path, names):
    saved = {k: sys.modules.get(k) for k in ['common', 'geometry', 'score', 'collision', 'run_r1']}
    old = list(sys.path)
    try:
        sys.path.insert(0, str(path))
        for k in saved:
            sys.modules.pop(k, None)
        return [module(name, Path(path) / (name + '.py')) for name in names]
    finally:
        sys.path[:] = old
        for (k, v) in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v

def state(phase, gate, nextop):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane=_release_expand('CROWN_CHAIN'), claim_type='capability', phase=phase, latest_gate=gate, next_operation=nextop, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW'))
