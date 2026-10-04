from pathlib import Path
import hashlib, json, datetime, sys, importlib.util
import os
R = Path(__file__).resolve().parents[1]
W = R / 'dependencies'
RESULTS = W / 'results'
DATA = R / 'data'

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def clean(x):
    import numpy as np
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, np.generic):
        return clean(x.item())
    if isinstance(x, float):
        return x if __import__('math').isfinite(x) else None
    if isinstance(x, Path):
        return str(x)
    return x

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(clean(x), ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def load(p):
    return json.loads(Path(p).read_text())

def module(name, path):
    path = Path(path)
    sys.path.insert(0, str(path.parent))
    s = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m
    s.loader.exec_module(m)
    return m

def artifact(p):
    return dict(path=str(p), sha256=sha(p), bytes=Path(p).stat().st_size)

def freeze(p, x):
    p = Path(p)
    if p.exists():
        raise RuntimeError('Immutable freeze already exists ' + str(p))
    dump(p, {'frozen_utc': now(), 'prereg_sha256': sha(R / 'PREREG_R1.json'), **x})
    p.with_suffix('.sha256').write_text(sha(p) + '\n')

def state(stage, gate, nextop):
    dump(R / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-patient360', claim_type='capability', updated_utc=now(), status=stage, latest_gate=gate, next_operation=nextop, review_state='PENDING_INDEPENDENT_REVIEW'))

def budget():
    n = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file())) if DATA.exists() else 0
    if n > 3000000000:
        raise RuntimeError('3 GB total lane budget exceeded')
    return n
