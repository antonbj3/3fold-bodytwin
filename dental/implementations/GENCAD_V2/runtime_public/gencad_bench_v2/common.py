from dental_release.paths import expand as _release_expand
import hashlib, json, os, time, datetime
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
DENTAL = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace'))
DATA = Path(os.environ.get('GENCAD_DATA', _release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V2')))
FAMILIES = ['molar_crown', 'premolar_crown', 'anterior_crown', 'bridge3', 'implant_crown', 'inlay', 'onlay', 'veneer', 'lattice_onlay']
LEVELS = ['easy', 'normal', 'hard', 'boundary']

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

def canonical(x):
    return json.dumps(clean(x), sort_keys=True, separators=(',', ':'), allow_nan=False).encode()

def digest(x):
    return hashlib.sha256(canonical(x)).hexdigest()

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')
    tmp.replace(p)

def read(p):
    return json.loads(Path(p).read_text())

def freeze(p, payload):
    p = Path(p)
    if p.exists():
        old = read(p)
        if old['payload_sha256'] != digest(old['payload']) or old['payload_sha256'] != digest(payload):
            raise ValueError('Frozen artifact mismatch: ' + str(p))
        return old
    out = dict(frozen_utc=now(), payload_sha256=digest(payload), payload=payload)
    dump(p, out)
    return out

def state(phase, gate, next_op, **kw):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-gencad-v2', phase=phase, latest_gate=gate, next_operation=next_op, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW', **kw))

def budget():
    size = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    if size > 3000000000:
        raise RuntimeError('Lane disk cap exceeded')
    return size
