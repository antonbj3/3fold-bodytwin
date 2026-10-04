from dental_release.paths import expand as _release_expand
from pathlib import Path
import datetime, hashlib, json
R = Path(__file__).resolve().parents[1]
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
X1B = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/X1B'))

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def dump(p, x):
    p = R / p
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, allow_nan=False, default=lambda v: v.item() if hasattr(v, 'item') else v.tolist()) + '\n')

def read(p):
    return json.loads((R / p).read_text())

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def state(phase, gate, next_operation):
    dump('CURRENT_WORK_STATE.json', dict(lane='X29-bridge-connector', updated_utc=now(), phase=phase, latest_gate=gate, next_operation=next_operation, claim_type='information_link', review_state='PENDING_INDEPENDENT_REVIEW', intermediate_budget_bytes=3000000000, physical_measurements_performed=False))

def freeze(p, x):
    p = R / p
    b = (json.dumps(x, indent=2, allow_nan=False) + '\n').encode()
    if p.exists():
        if p.read_bytes() != b:
            raise ValueError('Frozen file differs: ' + str(p))
    else:
        p.write_bytes(b)
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')
