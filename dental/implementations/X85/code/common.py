from dental_release.paths import expand as _release_expand
import os
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1')
import sys, json, hashlib, datetime, importlib.util
from pathlib import Path
sys.dont_write_bytecode = True
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
DENT = ROOT.parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X85'))
X82 = ROOT.parent / _release_expand('X82')
V6 = ROOT.parent / _release_expand('GENCAD_V6')
V4 = ROOT.parent / _release_expand('GENCAD_V4')

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

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
    return x

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(clean(x), indent=2, allow_nan=False, ensure_ascii=False) + '\n')
    tmp.replace(p)

def read(p):
    return json.loads(Path(p).read_text())

def freeze(p, x):
    p = Path(p)
    if p.exists():
        raise ValueError('Frozen artifact already exists: ' + str(p))
    dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix(p.suffix + '.sha256').write_text(sha(p) + '  ' + p.name + '\n')

def state(phase, gate, next_operation):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X85-occlusal-adjust', claim_type='capability', updated_utc=now(), phase=phase, latest_gate=gate, next_operation=next_operation, review_state='PENDING_INDEPENDENT_REVIEW', physical_measurement='NOT_RUN', data_root=DATA))

def module(name, p):
    s = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m
    s.loader.exec_module(m)
    return m

def v6_contact():
    saved = sys.modules.get('common')
    paths = list(sys.path)
    try:
        vc = module('x85_v6_common_readonly', V6 / 'code/common.py')
        sys.modules['common'] = vc
        sys.path.insert(0, str(V6 / 'code'))
        return module('x85_v6_contact_readonly', V6 / 'code/contact.py')
    finally:
        if saved is None:
            sys.modules.pop('common', None)
        else:
            sys.modules['common'] = saved
        sys.path[:] = paths

def x82():
    paths = list(sys.path)
    try:
        sys.path.insert(0, str(X82 / 'code'))
        return module('x85_x82_height_readonly', X82 / 'code/height_force.py')
    finally:
        sys.path[:] = paths
