from dental_release.paths import expand as _release_expand
import datetime, hashlib, json, os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DENT = ROOT.parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X78-same-patient-tissue'))
X77 = DENT / 'results/LANE_X77_SAME_PATIENT_FLOW'
X12 = DENT / 'results/LANE_X12_PULPY3D'
X73 = DENT / 'results/LANE_X73_CANAL_WALL_FULL'
X69 = DENT / 'results/LANE_X69_IOS_CBCT_PAIR'

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(4 << 20), b''):
            h.update(b)
    return h.hexdigest()

def read(p):
    return json.loads(Path(p).read_text())

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def freeze(name, x):
    p = ROOT / name
    if p.exists():
        assert sha(p) == p.with_suffix(p.suffix + '.sha256').read_text().split()[0], name
        return
    dump(p, x)
    p.with_suffix(p.suffix + '.sha256').write_text(sha(p) + '\n')

def state(stage, gate, next_operation):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X78-same-patient-tissue', updated_utc=now(), stage=stage, latest_gate=gate, next_operation=next_operation, threads_max=4, gpu=False, intermediate_limit_bytes=3000000000))

def diskcheck():
    DATA.mkdir(parents=True, exist_ok=True)
    n = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    assert n < 3000000000
    return n
