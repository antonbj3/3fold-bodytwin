from dental_release.paths import expand as _release_expand
import datetime, hashlib, json, os, sys
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
ROOT = Path(os.environ.get('X60_RUN_ROOT', str(Path(__file__).resolve().parents[1])))
DATA = Path(os.environ.get('X60_RUN_DATA', _release_expand('@DENTAL_WORK_ROOT@/X60')))
V4 = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/GENCAD_V4'))
X1 = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/X1B'))
X23 = Path(_release_expand('@DENTAL_INPUT_ROOT@/artifacts/LANE_X23_CROWN_WEIBULL_FLOOR'))

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
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
        raise ValueError('nonfinite report')
    return x

def read(p):
    return json.loads(Path(p).read_text())

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def freeze(p, x):
    p = Path(p)
    if p.exists():
        old = read(p)
        if {k: v for (k, v) in old.items() if k != 'frozen_utc'} != clean(x):
            raise ValueError('frozen drift ' + str(p))
    else:
        dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')

def state(phase, gate, next_operation):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X60-crown-optimizer', updated_utc=now(), phase=phase, latest_gate=gate, next_operation=next_operation, review_state='PENDING_INDEPENDENT_REVIEW', physical_measurement=False))

def inventory(base):
    return {str(p.relative_to(base)): dict(sha256=sha(p), bytes=p.stat().st_size) for p in sorted(Path(base).rglob('*')) if p.is_file() and '__pycache__' not in str(p)}

def budget():
    n = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    if n > 3000000000:
        raise RuntimeError('3GB budget exceeded')
    return n
