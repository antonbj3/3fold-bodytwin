from dental_release.paths import expand as _release_expand
from pathlib import Path
import datetime, hashlib, importlib.util, json, os, sys, time
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT.parent
PARENT = RESULTS / 'PROOF_LANE_PATIENT360'
PARENT_RUN = PARENT / 'runs/20261003T074621625651Z'
PARENT_DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_patient360')) / PARENT_RUN.name
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_patient360_r2'))
SOURCE_FILES = {}

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda : f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (tuple, list)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, np.generic):
        return clean(x.item())
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, float) and (not np.isfinite(x)):
        return None
    return x

def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(value), indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def load(path):
    return json.loads(Path(path).read_text())

def digest(value):
    return hashlib.sha256(json.dumps(clean(value), sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()

def source(path):
    path = Path(path).resolve()
    row = dict(path=str(path), sha256=sha(path), bytes=path.stat().st_size)
    SOURCE_FILES[str(path)] = row
    return row

def imported(name, path):
    source(path)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

def freeze(path, value):
    path = Path(path)
    if path.exists():
        raise ValueError('IMMUTABLE_FREEZE_EXISTS:' + str(path))
    dump(path, dict(frozen_utc=now(), **value))
    path.with_suffix('.sha256').write_text(sha(path) + '\n')

def state(status, gate, next_operation, **kw):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-patient360-r2', claim_type='capability', updated_utc=now(), status=status, latest_gate=gate, next_operation=next_operation, review_state='PENDING_INDEPENDENT_REVIEW', **kw))

def budget():
    n = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file())) if DATA.exists() else 0
    if n > 3000000000:
        raise RuntimeError('LANE_3GB_LIMIT')
    return n

def check(name, normal_pass, fault_rejected, evidence=None):
    return dict(name=name, normal_pass=bool(normal_pass), injected_fault_rejected=bool(fault_rejected), pass_=bool(normal_pass and fault_rejected), evidence=evidence)

def rejects(fn, *args, **kw):
    try:
        fn(*args, **kw)
    except (ValueError, AssertionError):
        return True
    return False
