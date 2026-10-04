import json, hashlib, datetime
from pathlib import Path
H = Path(__file__).resolve().parent

def read(n):
    return json.loads((H / n).read_text())

def write(n, x):
    (H / n).write_text(json.dumps(x, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def exact_equal(a, b):
    return json.dumps(a, sort_keys=True, ensure_ascii=False, allow_nan=False) == json.dumps(b, sort_keys=True, ensure_ascii=False, allow_nan=False)

def entries():
    out = read('INPUT_LOCK.json')['files']
    if (H / 'AUX_INPUT_LOCK.json').exists():
        out.update(read('AUX_INPUT_LOCK.json')['files'])
    return out

def entry(k):
    return entries()[k]

def source_path(k):
    r = entry(k)
    p = H / r['snapshot']
    if sha(p) != r['sha256']:
        raise ValueError('SOURCE_HASH:' + k)
    return p

def source(k):
    return json.loads(source_path(k).read_text())

def verify():
    for k in entries():
        source_path(k)
    for n in ['INPUT_LOCK.json', 'AUX_INPUT_LOCK.json', 'PREREG_R4A.json', 'DECOMPOSITION.json', 'FROZEN_PREDICTIONS.json', 'PREREG_R4B.json', 'RANGE_CONTRACTS.json']:
        if (H / (n + '.sha256')).exists() and sha(H / n) != (H / (n + '.sha256')).read_text().strip():
            raise ValueError('FROZEN_HASH:' + n)

def state(status, gate, next_op):
    write('CURRENT_WORK_STATE.json', {'lane': 'PROOF_LANE-constraint-net-r4', 'updated_utc': utc(), 'status': status, 'latest_gate': gate, 'next_operation': next_op, 'intermediate_limit_bytes': 3000000000})
