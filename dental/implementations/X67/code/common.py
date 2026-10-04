from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, hashlib, datetime
ROOT = Path(__file__).resolve().parents[1]
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
DENTAL = ROOT.parents[1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def dump(n, x):
    p = ROOT / n
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def load(n):
    return json.loads((ROOT / n).read_text())

def freeze(n, x):
    if (ROOT / n).exists():
        raise RuntimeError('Refuse to replace frozen ' + n)
    x = dict(x, created_at_utc=now())
    dump(n, x)
    (ROOT / Path(n).with_suffix('.sha256')).write_text(sha(ROOT / n) + '  ' + n + '\n')

def verify(n):
    assert sha(ROOT / n) == (ROOT / Path(n).with_suffix('.sha256')).read_text().split()[0], n

def state(stage, gate, next_operation):
    old = load('CURRENT_WORK_STATE.json')
    old.update(status=stage, latest_gate=gate, next_operation=next_operation, timestamp_utc=now())
    dump('CURRENT_WORK_STATE.json', old)
