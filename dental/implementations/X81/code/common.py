import hashlib, json, datetime
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def read(key):
    rec = json.loads((ROOT / 'INPUT_LOCK.json').read_text())[key]
    p = ROOT / rec['snapshot']
    if hashlib.sha256(p.read_bytes()).hexdigest() != rec['sha256']:
        raise ValueError('input hash mismatch: ' + key)
    obj = json.loads(p.read_text())
    mapping = json.loads((ROOT / 'raw/RELOCATED_ARRAYS.json').read_text()) if (ROOT / 'raw/RELOCATED_ARRAYS.json').exists() else {}

    def relocate(value):
        if isinstance(value, dict):
            return {k: relocate(v) for (k, v) in value.items()}
        if isinstance(value, list):
            return [relocate(v) for v in value]
        if isinstance(value, str) and value in mapping:
            return str(ROOT / mapping[value])
        return value
    return relocate(obj)

def save(name, obj):
    p = ROOT / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def verify_freeze(name):
    p = ROOT / name
    assert hashlib.sha256(p.read_bytes()).hexdigest() == (ROOT / (name + '.sha256')).read_text().strip(), name
    obj = json.loads(p.read_text())
    mapping = json.loads((ROOT / 'raw/RELOCATED_ARRAYS.json').read_text()) if (ROOT / 'raw/RELOCATED_ARRAYS.json').exists() else {}

    def relocate(value):
        if isinstance(value, dict):
            return {k: relocate(v) for (k, v) in value.items()}
        if isinstance(value, list):
            return [relocate(v) for v in value]
        if isinstance(value, str) and value in mapping:
            return str(ROOT / mapping[value])
        return value
    return relocate(obj)

def state(stage, gate, next_operation):
    save('CURRENT_WORK_STATE.json', {'lane': 'X81-decide-within-noise', 'claim_type': 'capability', 'status': stage, 'updated_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'latest_gate': gate, 'next_operation': next_operation, 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'resources': {'threads_max': 4, 'intermediates_limit_bytes': 3000000000, 'gpu': False}})

def scalar_decision(lo, hi, tau):
    if lo is None or hi is None:
        return 'UNKNOWN'
    if lo > hi:
        raise ValueError('unordered interval')
    return 'PASS' if hi <= tau else 'FAIL' if lo > tau else 'UNRESOLVED'
