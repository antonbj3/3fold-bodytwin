from pathlib import Path
import datetime, hashlib, json, os
ROOT = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def dump(p, value):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def check_frozen(name):
    p = ROOT / name
    assert sha(p) == (ROOT / (name + '.sha256')).read_text().strip(), name
    return json.loads(p.read_text())

def inputs():
    manifest = check_frozen('SOURCE_MANIFEST.json')
    for (key, row) in manifest.items():
        assert sha(ROOT / row['path']) == row['sha256'], key
    return manifest

def state(**kwargs):
    p = ROOT / 'CURRENT_WORK_STATE.json'
    s = json.loads(p.read_text())
    s.update(updated_utc=stamp(), **kwargs)
    dump(p, s)
