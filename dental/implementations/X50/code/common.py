"""Local research only. No network, no source writes, at most two BLAS threads."""
from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '2'
import datetime, hashlib, json
from pathlib import Path
P = Path(__file__).resolve().parents[1]
X7 = P.parent / _release_expand('X7')
X21 = P.parent / _release_expand('X21')
D = Path(_release_expand('@DENTAL_WORK_ROOT@/X50-bite-screening'))
FINDINGS = ['crossbite', 'open_bite', 'deep_bite', 'angle_II', 'angle_III', 'scissor_bite']
NAMES = ['Korsbett', 'Anteriort open bite', 'Deep bite', 'Angle II (molar)', 'Angle III (molar)', 'Saxbett']

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def read(p):
    return json.loads(Path(p).read_text())

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def write(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def freeze(p, x):
    p = Path(p)
    if p.exists():
        raise RuntimeError('Refusing to replace frozen artifact ' + str(p))
    write(p, x)
    p.with_name(p.name + '.sha256').write_text(sha(p) + '\n')

def state(status, gate, next_operation):
    write(P / 'CURRENT_WORK_STATE.json', dict(lane='X50-bite-screening', updated_utc=now(), status=status, latest_gate=gate, next_operation=next_operation, review_state='PENDING_INDEPENDENT_REVIEW'))
