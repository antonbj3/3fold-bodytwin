import hashlib, json, os, time
from pathlib import Path
from datetime import datetime, timezone
P = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False, ensure_ascii=False) + '\n')
    tmp.replace(path)

def utc():
    return datetime.now(timezone.utc).isoformat()

def check(round):
    f = P / ('PREREG_' + round + '.json')
    assert sha(f) == f.with_suffix('.sha256').read_text().strip(), 'Frozen prereg hash drift'
    return sha(f)

def state(phase, gate, next_operation, directory=None):
    dump((Path(directory) if directory is not None else P) / 'CURRENT_WORK_STATE.json', dict(lane='X95-seating-identifiable-assay', phase=phase, latest_gate=gate, next_operation=next_operation, updated_utc=utc(), physical_status='UNKNOWN_MISSING_SAME_SPECIMEN_OBSERVATIONS', review_state='PENDING_INDEPENDENT_REVIEW'))
