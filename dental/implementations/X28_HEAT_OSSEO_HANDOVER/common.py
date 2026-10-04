"""Shared provenance and fixed thermal arithmetic; no fitted tissue-death law."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
ROOT = Path(__file__).resolve().parent
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
DENTAL = ROOT.parents[1]

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def now():
    return datetime.now(timezone.utc).isoformat()

def dump(name, value):
    path = ROOT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def load(name):
    return json.loads((ROOT / name).read_text())

def verify_freeze(name):
    expected = (ROOT / name.replace('.json', '.sha256')).read_text().split()[0]
    if sha(ROOT / name) != expected:
        raise ValueError('Frozen file changed: ' + name)

def freeze(name, obj):
    if (ROOT / name).exists():
        raise ValueError('Refuse to overwrite frozen artifact: ' + name)
    dump(name, obj)
    (ROOT / name.replace('.json', '.sha256')).write_text(sha(ROOT / name) + '  ' + name + '\n')

def cem_rate(temp_c):
    r = 0.5 if temp_c >= 43 else 0.25
    return r ** (43 - temp_c)

def state(stage, gate, next_operation):
    dump('CURRENT_WORK_STATE.json', dict(lane='X28-heat-osseo', updated_at_utc=now(), stage=stage, latest_gate=gate, next_operation=next_operation, claim_type='information_link', review_state='PENDING_INDEPENDENT_REVIEW'))
