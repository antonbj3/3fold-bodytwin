from dental_release.paths import expand as _release_expand
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '2'
from pathlib import Path
import json, hashlib, datetime, time, resource
ROOT = Path(__file__).resolve().parents[1]
DENT = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace'))
DATA = ROOT / 'models'
ZIP = Path(_release_expand('@DENTAL_DATA_ROOT@/geometry/Bite2Text/Bite2Text.zip'))
X11 = Path(_release_expand('@DENTAL_WORK_ROOT@/X11/targets/Bite2Text'))
X7 = DENT / 'results/LANE_X7_BITE2TEXT'
X21 = DENT / 'results/LANE_X21_CONTACT_MAP'
X76 = DENT / 'results/LANE_X76_DATASET_COMPLETION'

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def digest(b):
    return hashlib.sha256(b).hexdigest()

def relocated(value):
    if isinstance(value, dict):
        return {k: relocated(v) for (k, v) in value.items()}
    if isinstance(value, list):
        return [relocated(v) for v in value]
    if isinstance(value, str):
        for (old, new) in [(_release_expand('@DENTAL_IMPLEMENTATIONS@/X79'), str(ROOT)), (_release_expand('@DENTAL_WORK_ROOT@/X79'), str(DATA))]:
            if value == old or value.startswith(old + '/'):
                return new + value[len(old):]
    return value

def read(p):
    p = Path(relocated(str(p)))
    return relocated(json.loads((p if p.is_absolute() else ROOT / p).read_text()))

def write(p, x):
    p = Path(p)
    p = p if p.is_absolute() else ROOT / p
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def freeze(p, x):
    p = ROOT / p
    if p.exists():
        assert sha(p) == p.with_suffix('.sha256').read_text().strip(), str(p) + ' drift'
        return
    write(p, x)
    p.with_suffix('.sha256').write_text(sha(p) + '\n')

def state(stage, gate, next_operation, **extra):
    write('CURRENT_WORK_STATE.json', dict(lane='X79-report-findings', claim_type=['information_link', 'capability'], updated_utc=now(), stage=stage, latest_gate=gate, next_operation=next_operation, **extra))

def cost(start, cpu):
    return dict(wall_s=time.perf_counter() - start, cpu_s=time.process_time() - cpu, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, threads=2, GPU=False)
