from dental_release.paths import expand as _release_expand
from pathlib import Path
import os, sys, json, hashlib, datetime, importlib.util
os.environ.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMEXPR_NUM_THREADS='4', MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1')
sys.dont_write_bytecode = True
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
OLD = BASE / 'PROOF_LANE_FULL_CROWN_R2'
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_FULL_CROWN_R3'))
OLD_DATA = DATA.parent / 'PROOF_LANE_FULL_CROWN_R2'
PUBLIC = DATA.parent / 'PROOF_LANE_FULL_CROWN/public'
V4 = BASE / _release_expand('GENCAD_V4')
V6 = BASE / _release_expand('GENCAD_V6')
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
    q = p.with_suffix(p.suffix + '.tmp')
    q.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')
    q.replace(p)

def freeze(p, x):
    p = Path(p)
    if p.exists():
        raise ValueError('frozen file exists ' + str(p))
    dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')

def npz(p):
    with np.load(p, allow_pickle=False) as z:
        return dict(z)

def module(name, p):
    s = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m
    s.loader.exec_module(m)
    return m

def state(phase, gate, nextop):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-full-crown-r3', claim_type='capability', phase=phase, latest_gate=gate, next_operation=nextop, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW'))
REFERENT = dict(kind='published_dataset', locator=['https://ditto.ing.unimore.it/bits2bites/', 'https://ditto.ing.unimore.it/bite2text/', str(V4 / 'payload/whole_private')], compared_quantity='Whole observed tooth surface and geometric opposing gap in supplied pose; X11 inferred labels. New contact prescription comes from the same scan: contact comparison is reconstruction of prescribed data, NOT independent contact prediction. Axial source surfaces not supplied to contact generator.', refutes_us=True)
