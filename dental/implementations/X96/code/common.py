from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import json, hashlib, datetime, sys, time, resource
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
DENT = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace'))
DATA = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/storage/tmp_dental_audit/XREV-batch32/runtime/X96'))
X8 = DENT / 'results/LANE_X8_GUIDE_NERVE_RISK'
X73 = DENT / 'results/LANE_X73_CANAL_WALL_FULL'
CERT = Path(_release_expand('@DENTAL_WORK_ROOT@/X8-guide-nerve-risk/R5_FULL_GEOMETRY.jsonl'))
sys.path.insert(0, str(X8))
from full_geometry import voxel_cylinder_bracket
sys.path.insert(0, str(DENT / 'cells/geometry'))
import tf2_io

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(4 << 20), b''):
            h.update(b)
    return h.hexdigest()

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    t = p.with_suffix(p.suffix + '.tmp')
    t.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    t.replace(p)

def load(p):
    return json.loads(Path(p).read_text())

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def state(stage, gate, next_operation, **kw):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X96-drill-overshoot', updated_utc=now(), stage=stage, latest_gate=gate, next_operation=next_operation, **kw))

def freeze(name, x):
    p = ROOT / name
    if p.exists():
        assert load(p) == x, 'Frozen prereg changed'
    else:
        dump(p, x)
    p.with_suffix(p.suffix + '.sha256').write_text(sha(p) + '\n')

def check():
    for p in ROOT.glob('PREREG_*.json'):
        assert sha(p) == p.with_suffix(p.suffix + '.sha256').read_text().strip()
    pr = ROOT / 'PREREG_R1.json'
    if pr.exists():
        for (p, h) in load(pr)['inputs'].items():
            assert sha(p) == h, 'Frozen source changed: ' + p
    if (ROOT / 'PREREG_R2.json').exists():
        assert load(ROOT / 'PROTOCOLS.json')['protocols'] == load(ROOT / 'PREREG_R2.json')['selection']['protocols'], 'Protocol/datum binding changed'

def diskcheck():
    DATA.mkdir(parents=True, exist_ok=True)
    n = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    assert n < 3000000000
    return n

def certs():
    return [json.loads(l) for l in CERT.read_text().splitlines()]

def sources():
    return {str(p): sha(p) for p in [CERT, X73 / 'raw/LINEAGE_SITES.json', X73 / 'raw/IMAGE_GROUPS.json', X8 / 'full_geometry.py', X8 / 'GUIDE_PROFILES.json']}
