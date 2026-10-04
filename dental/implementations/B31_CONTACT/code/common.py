from dental_release.paths import expand as _release_expand
import os
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1')
import sys, json, hashlib, datetime, importlib.util
from pathlib import Path
sys.dont_write_bytecode = True
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
DENT = ROOT.parents[1]
UP = Path(os.environ.get('X95_UPSTREAM_ROOT', _release_expand('@DENTAL_INPUT_ROOT@/artifacts')))
DATA = Path(os.environ.get('X95_DATA_ROOT', _release_expand('@DENTAL_WORK_ROOT@/X95-spatial-contact-repair')))
PYTHON = _release_expand('@DENTAL_PYTHON@')

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
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, float) and (not np.isfinite(x)):
        return None
    return x

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    t = p.with_suffix(p.suffix + '.tmp')
    t.write_text(json.dumps(clean(x), indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    t.replace(p)

def read(p):
    return json.loads(Path(p).read_text())

def freeze(p, x):
    p = Path(p)
    if p.exists():
        raise ValueError('Frozen file already exists: ' + str(p))
    dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix(p.suffix + '.sha256').write_text(sha(p) + '  ' + p.name + '\n')

def state(phase, latest, next_operation):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X95-spatial-contact-repair', claim_type='capability', phase=phase, updated_utc=now(), latest_gate=latest, next_operation=next_operation, review_state='PENDING_INDEPENDENT_REVIEW', physical_measurement='NOT_RUN', data_root=DATA, intermediate_bytes_limit=3000000000, threads=1))

def module(name, p):
    s = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m
    s.loader.exec_module(m)
    return m

def parents():
    saved = sys.modules.get('common')
    paths = list(sys.path)
    try:
        vc = module('x95_v6_common', UP / 'PROOF_LANE_GENCAD_V6/code/common.py')
        sys.modules['common'] = vc
        sys.path.insert(0, str(UP / 'PROOF_LANE_GENCAD_V6/code'))
        contact = module('x95_contact', UP / 'PROOF_LANE_GENCAD_V6/code/contact.py')
        reuse = module('x95_reuse', UP / 'PROOF_LANE_GENCAD_V6/code/reuse.py')
    finally:
        if saved is None:
            sys.modules.pop('common', None)
        else:
            sys.modules['common'] = saved
        sys.path[:] = paths
    return (contact, reuse)

def check_lock():
    for (p, d) in read(ROOT / 'INPUT_LOCK_R1.json')['files'].items():
        if sha(p) != d['sha256']:
            raise ValueError('Input drift: ' + p)
    paths = [ROOT / 'INPUT_LOCK_R1.json']
    for pattern in ['PREREG_*.json', 'DECOMPOSITION_*.json', 'FROZEN*.json']:
        paths.extend(ROOT.glob(pattern))
    for p in paths:
        if p.exists() and sha(p) != p.with_suffix(p.suffix + '.sha256').read_text().split()[0]:
            raise ValueError('Contract drift: ' + str(p))
