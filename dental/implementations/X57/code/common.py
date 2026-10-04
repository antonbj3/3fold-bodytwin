from dental_release.paths import expand as _release_expand
import datetime
import hashlib
import json
import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = Path(os.environ.get('X57_OUTPUT', str(ROOT))).resolve()
ARCHIVE = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/datasets/dental_3fold_extra/MMDental/MMDental.zip'))
MEMBER = 'MMDental/medical_records.csv'

def sha(data):
    return hashlib.sha256(data).hexdigest()

def read(name):
    p = OUTPUT / name
    return json.loads((p if p.exists() else ROOT / name).read_text())

def write(name, obj):
    path = OUTPUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def freeze(name):
    if OUTPUT != ROOT:
        raise ValueError('Cannot freeze prereg during replay')
    p = ROOT / name
    obj = read(name)
    obj['frozen_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    write(name, obj)
    (ROOT / name.replace('.json', '.sha256')).write_text(sha(p.read_bytes()) + '  ' + name + '\n')

def check_frozen(name):
    expected = (ROOT / name.replace('.json', '.sha256')).read_text().split()[0]
    if sha((ROOT / name).read_bytes()) != expected:
        raise ValueError('Frozen file changed: ' + name)

def checkpoint(phase, gate, next_operation):
    state = read('CURRENT_WORK_STATE.json')
    state.update(updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), phase=phase, latest_gate=gate, next_operation=next_operation)
    write('CURRENT_WORK_STATE.json', state)
