from dental_release.paths import expand as _release_expand
from pathlib import Path
import os, json, hashlib, datetime, time, resource
import numpy as np
R = Path(os.environ.get('X95_OUTPUT_ROOT', Path(__file__).resolve().parents[1]))
SRC = Path(__file__).resolve().parents[1]
D = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace'))
X1 = D / 'results/LANE_X1B_CROWN_LOOP'
X60 = D / 'results/LANE_X60_CROWN_OPTIMIZER'
DATA = Path(os.environ.get('X95_DATA_ROOT', _release_expand('@DENTAL_WORK_ROOT@/X95-native-support-shape')))

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def read(p):
    return json.loads(Path(p).read_text())

def clean(v):
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, dict):
        return {k: clean(x) for (k, x) in v.items()}
    if isinstance(v, (list, tuple)):
        return [clean(x) for x in v]
    return v

def dump(p, obj):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(clean(obj), indent=2, allow_nan=False, ensure_ascii=False) + '\n')

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def freeze(p, obj):
    p = Path(p)
    if p.exists():
        raise ValueError('Cannot replace frozen artifact ' + str(p))
    dump(p, obj)
    Path(str(p) + '.sha256').write_text(sha(p) + '\n')

def state(m, gate, n):
    dump(R / 'CURRENT_WORK_STATE.json', dict(status='ACTIVE', updated_utc=now(), milestone=m, latest_gate=gate, next_operation=n, review_state='PENDING_INDEPENDENT_REVIEW', physical_measurements_performed=False))

def verify_inputs():
    for row in read(SRC / 'INPUT_LOCK.json')['files']:
        if sha(row['path']) != row['sha256']:
            raise ValueError('Input hash drift: ' + row['path'])

def cost(t0):
    return dict(wall_seconds=time.perf_counter() - t0, process_CPU_seconds=time.process_time(), peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
