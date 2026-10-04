"""Read-only reuse of measured source import and continuous triangle extrema."""
from dental_release.paths import expand as _release_expand
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import sys, json, hashlib, datetime, time
from pathlib import Path
import numpy as np
H = Path(__file__).resolve().parents[1]
D = Path(_release_expand('@DENTAL_WORK_ROOT@/X18b_prep_surface'))
X18 = H.parent / _release_expand('X18')
P = H.parent / 'LANE_NEXT_P_OCCLUSION_VALIDATION'
sys.path.insert(0, str(X18 / 'code'))
from geometry import pair, query_height, write_stl
from experiment import contact, score
sys.path.insert(0, str(P / 'code'))
import continuous_gap as cg

def sha(p):
    q = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            q.update(b)
    return q.hexdigest()

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, allow_nan=False, default=lambda a: a.tolist() if isinstance(a, np.ndarray) else a.item() if isinstance(a, np.generic) else str(a)) + '\n')

def state(status, gate, next_operation, **kw):
    dump(H / 'CURRENT_WORK_STATE.json', dict(lane='X18b-prep-surface', status=status, latest_gate=gate, next_operation=next_operation, updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), **kw))

def artifact(p):
    return dict(path=str(p), sha256=sha(p), bytes=Path(p).stat().st_size)

def data_budget():
    n = sum((p.stat().st_size for p in D.rglob('*') if p.is_file())) if D.exists() else 0
    if n > 3000000000:
        raise RuntimeError('Lane intermediate size limit exceeded')
    return n

def verify_freeze(tag):
    fr = json.loads((H / f'FROZEN_PREDICTIONS_{tag}.json').read_text())
    assert sha(H / f'PREREG_{tag}.json') == fr['prereg_sha256']
    for f in fr['files']:
        assert sha(f['path']) == f['sha256'], f['path']
    for (p, s) in fr['code_sha256'].items():
        assert sha(p) == s, p
    return fr
