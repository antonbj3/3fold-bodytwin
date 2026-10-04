from dental_release.paths import expand as _release_expand
import os
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', MPLBACKEND='Agg')
from pathlib import Path
import datetime, hashlib, json, time, resource
resource.setrlimit(resource.RLIMIT_AS, (3500 * 1024 ** 2, 3500 * 1024 ** 2))
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get('X95_DATA_DIR', _release_expand('@DENTAL_WORK_ROOT@/X95_crown_prep_envelope')))
DATA.mkdir(parents=True, exist_ok=True)

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def read(p):
    return json.loads(Path(p).read_text())

def clean(x):
    if isinstance(x, (np.ndarray, np.generic)):
        return clean(x.tolist())
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, float) and (not np.isfinite(x)):
        return None
    return x

def dump(p, d):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(clean(d), indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    tmp.replace(p)

def state(phase, gate, nextop):
    dump(ROOT / 'CURRENT_WORK_STATE.json', {'updated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'phase': phase, 'latest_gate': gate, 'next_operation': nextop, 'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW'})

def load_np(p):
    with np.load(p, allow_pickle=False) as z:
        return dict(z)

def usage():
    return {'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'threads_configured': 1, 'data_bytes': sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))}
