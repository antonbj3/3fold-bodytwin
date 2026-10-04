from dental_release.paths import expand as _release_expand
from pathlib import Path
import os, sys, json, hashlib, datetime, time, importlib.util, resource
os.environ.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMEXPR_NUM_THREADS='4', MPLBACKEND='Agg')
sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_AS, (3500 * 1024 ** 2, 3500 * 1024 ** 2))
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_FULL_CROWN_R4'))

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
    if isinstance(x, (np.ndarray, np.generic)):
        return clean(x.tolist())
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
    p.parent.mkdir(parents=True, exist_ok=True)
    q = p.with_suffix(p.suffix + '.tmp')
    q.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')
    q.replace(p)

def freeze(p, x):
    p = Path(p)
    if p.exists():
        raise ValueError('Already frozen ' + str(p))
    dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')

def npz(p):
    with np.load(p, allow_pickle=False) as z:
        return dict(z)

def state(phase, gate, nextop):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-full-crown-r4', claim_type='capability', phase=phase, latest_gate=gate, next_operation=nextop, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW'))

def module(name, p):
    s = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m
    s.loader.exec_module(m)
    return m
REFERENT = dict(kind='published_dataset', locator=['https://osf.io/xctdy/', str(BASE / 'PROOF_LANE_FULL_CROWN_DIAG/PREREG_ANNOTATED.json')], compared_quantity='Original-annotated natural crown triangles; withheld exterior above virtual cervical observation. Same-subject lower arch is measured geometry but bite registration and loaded contact are UNKNOWN.', refutes_us=True)
