from dental_release.paths import expand as _release_expand
import os, sys
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '4'
sys.dont_write_bytecode = True
import json, hashlib, datetime, importlib.util
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X72'))
sys.path.insert(0, str(ROOT / 'vendor'))

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def enc(x):
    if isinstance(x, (Q, Decimal)):
        return str(x)
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {str(k): enc(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [enc(v) for v in x]
    if hasattr(x, 'tolist'):
        return enc(x.tolist())
    return x

def write(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(enc(x), indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    tmp.replace(p)

def freeze(p, x):
    p = Path(p)
    if not p.exists():
        write(p, x)
        p.with_suffix(p.suffix + '.sha256').write_text(sha(p) + '\n')
    assert sha(p) == p.with_suffix(p.suffix + '.sha256').read_text().strip()

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def state(stage, gate, next_op):
    write(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X72-falt-fixes', status=stage, updated_utc=now(), latest_gate=gate, next_operation=next_op, data_directory=DATA, intermediate_limit_bytes=3000000000))

def record(p):
    p = Path(p).resolve()
    if p.is_relative_to(ROOT) or p.is_relative_to(DATA):
        return p
    manifest = ROOT / 'SOURCE_MANIFEST.json'
    d = json.loads(manifest.read_text()) if manifest.exists() else {}
    s = sha(p)
    if str(p) in d:
        assert d[str(p)]['sha256'] == s, 'SOURCE_DRIFT ' + str(p)
    else:
        d[str(p)] = dict(sha256=s, bytes=p.stat().st_size)
        write(manifest, d)
    return p

def read(p, exact=False):
    p = record(p)
    return json.loads(p.read_text(), parse_float=Decimal) if exact else json.loads(p.read_text())

def module(name, p):
    p = record(p)
    spec = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m

def verify_sources():
    d = json.loads((ROOT / 'SOURCE_MANIFEST.json').read_text())
    return dict(files=len(d), changed=[p for (p, r) in d.items() if sha(p) != r['sha256']])
