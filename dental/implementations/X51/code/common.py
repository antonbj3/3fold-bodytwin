import datetime, hashlib, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for (k, v) in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if hasattr(o, 'tolist'):
        return clean(o.tolist())
    return o

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(name):
    return json.loads((ROOT / name).read_text())

def save(name, obj):
    p = ROOT / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(clean(obj), indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def check_lock(name):
    p = ROOT / name
    if sha(p) != p.with_suffix('.sha256').read_text().split()[0]:
        raise ValueError('Frozen contract changed: ' + name)

def state(stage, gate, next_operation):
    save('CURRENT_WORK_STATE.json', dict(lane='X51-lab-acceptance', claim_type='capability', stage=stage, latest_gate=gate, next_operation=next_operation, updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat()))
