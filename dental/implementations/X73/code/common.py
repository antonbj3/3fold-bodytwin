from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import datetime, hashlib, io, json, subprocess, sys, time, zipfile
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
DENT = ROOT.parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X73-canal-wall-full'))
X58 = DENT / 'results/LANE_X58_CANAL_WALL_SPREAD'
X8 = DENT / 'results/LANE_X8_GUIDE_NERVE_RISK'
TF1 = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/8838D60F38D5FBDE/ToothFairy_Dataset.zip'))
MAX = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/8838D60F38D5FBDE/IAN_Maxillo_dataset.zip'))
TF2 = Path(_release_expand('@DENTAL_CORPUS_ROOT@/geometry/ToothFairy2/ToothFairy2_Dataset.zip'))
SP = 0.3

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(4 << 20), b''):
            h.update(b)
    return h.hexdigest()

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    q = p.with_suffix(p.suffix + '.tmp')
    q.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    q.replace(p)

def state(stage, gate, next_operation, **kw):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X73-canal-wall-full', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), stage=stage, latest_gate=gate, next_operation=next_operation, **kw))

def check_frozen():
    for p in ROOT.glob('PREREG_*.json'):
        assert sha(p) == p.with_suffix('.sha256').read_text().split()[0], str(p)

def diskcheck():
    DATA.mkdir(parents=True, exist_ok=True)
    n = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    assert n < 3000000000, '3 GB lane budget exceeded'
    return n

def member(z, n):
    if z.getinfo(n).compress_type == 9:
        return subprocess.check_output(['7z', 'e', '-so', '-mmt=1', str(z.filename), n], stderr=subprocess.DEVNULL)
    return z.read(n)

def npy(z, n):
    b = member(z, n)
    return (np.load(io.BytesIO(b), allow_pickle=False), hashlib.sha256(b).hexdigest())

def mha(z, n):
    b = member(z, n)
    end = b.index(b'\n', b.index(b'ElementDataFile')) + 1
    h = dict((s.split(' = ', 1) for s in b[:end].decode().strip().split('\n') if ' = ' in s))
    assert h.get('CompressedData', 'False') == 'False' and h['ElementDataFile'].strip() == 'LOCAL'
    dt = {'MET_UCHAR': 'u1', 'MET_CHAR': 'i1', 'MET_DOUBLE': 'f8', 'MET_FLOAT': 'f4', 'MET_SHORT': 'i2', 'MET_INT': 'i4', 'MET_UINT': 'u4', 'MET_USHORT': 'u2'}[h['ElementType']]
    assert h.get('BinaryDataByteOrderMSB', 'False') == 'False'
    shape = tuple(map(int, h['DimSize'].split()[::-1]))
    a = np.frombuffer(b, dtype='<' + dt, count=int(np.prod(shape)), offset=end).reshape(shape)
    return (a, h, hashlib.sha256(b).hexdigest())

def arrsha(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def classify(lo, hi, t=2.0):
    assert np.isfinite(lo) and np.isfinite(hi) and (0 <= lo <= hi)
    return 'BELOW' if hi < t else 'ABOVE' if lo >= t else 'UNKNOWN'

def summary(a):
    a = np.asarray(a)
    a = a[np.isfinite(a)]
    return dict(n=int(a.size), median_mm=float(np.median(a)) if a.size else None, p95_mm=float(np.quantile(a, 0.95)) if a.size else None, max_mm=float(a.max()) if a.size else None)

def writecsv(p, rows):
    import csv
    with open(p, 'w') as f:
        if not rows:
            return
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
