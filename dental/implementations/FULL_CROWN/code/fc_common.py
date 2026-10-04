from dental_release.paths import expand as _release_expand
import os, sys, json, hashlib, datetime, time
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMEXPR_NUM_THREADS='4', PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
sys.dont_write_bytecode = True
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_FULL_CROWN'))
V4 = ROOT.parent / _release_expand('GENCAD_V4')
V5 = ROOT.parent / _release_expand('GENCAD_V5')
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
    if isinstance(x, (tuple, list)):
        return [clean(v) for v in x]
    if isinstance(x, float) and (not np.isfinite(x)):
        return None
    return x

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(exist_ok=True, parents=True)
    s = p.with_suffix(p.suffix + '.tmp')
    s.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')
    s.replace(p)

def freeze(p, x):
    p = Path(p)
    if p.exists():
        raise ValueError('already frozen ' + str(p))
    dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')

def npz(p):
    with np.load(p, allow_pickle=False) as a:
        return dict(a)

def state(phase, gate, nextop):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-full-crown', claim_type='capability', phase=phase, latest_gate=gate, next_operation=nextop, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW'))

def functional():
    sys.path.insert(0, str(V5 / 'gencad_bench_v5'))
    import functional as f
    return f
REFERENT = dict(kind='published_dataset', locator=['https://ditto.ing.unimore.it/bite2text/', 'https://ditto.ing.unimore.it/bits2bites/', str(V4 / 'payload/whole_private')], compared_quantity='Native measured tooth surface and supplied registered neighbour/antagonist geometry; inferred FDI and virtual preparation are not measurements', refutes_us=True)
