"""Exact finite acquisition semantics; no empirical probabilities or laboratory data."""
from pathlib import Path
import hashlib, json, datetime
HERE = Path(__file__).resolve().parent

def read(name):
    return json.loads((HERE / name).read_text())

def write(name, o):
    p = HERE / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(o, ensure_ascii=False, indent=2) + '\n')

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def checkpoint(status, gate, next_operation):
    state = dict(lane='X71-measurement-priority', status=status, updated_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), latest_gate=gate, next_operation=next_operation, physical_measurements_performed=0, threads=1, intermediate_limit_bytes=3000000000)
    write('CURRENT_WORK_STATE.json', state)
    with (HERE / 'MILESTONES.jsonl').open('a') as f:
        f.write(json.dumps(state, ensure_ascii=False) + '\n')

def verify_lock():
    d = read('INPUT_LOCK.json')
    for x in d['inputs']:
        if sha(HERE / x['path']) != x['sha256']:
            raise ValueError('INPUT_HASH_DRIFT:' + x['path'])
    if sha(HERE / 'SOURCE_MANIFEST.json') != d['manifest_sha256']:
        raise ValueError('SOURCE_MANIFEST_DRIFT')

def evaluate(ids, decisions, measures):
    chosen = [m for m in measures if m['id'] in ids]
    supplied = set((s for m in chosen for s in m['supplies']))
    ready = [d['id'] for d in decisions if d.get('count_in_value', True) and d['requires'] and (set(d['requires']) <= supplied) and (not d['other_blockers'])]
    touched = [d['id'] for d in decisions if d.get('count_in_value', True) and set(d['requires']) & supplied]
    eligible = [d for d in decisions if d.get('count_in_value', True)]
    chains = {d['chain'] for d in eligible}
    complete = [c for c in chains if all((d['id'] in ready for d in eligible if d['chain'] == c))]
    return dict(ids=sorted(ids), supplied=sorted(supplied), evaluation_ready=ready, touched=touched, complete_chains=sorted(complete), operator_hours=[sum((m['operator_hours_interval'][i] for m in chosen)) for i in [0, 1]], instrument_hours=[sum((m['instrument_hours_interval'][i] for m in chosen)) for i in [0, 1]], new_specimens=sum((m['new_specimens'] for m in chosen)), guaranteed_decidable_without_outcomes=0, already_physically_decided=0)

def admissible(ids, measures):
    if 'M04' in ids and 'M05' in ids:
        return False
    if 'M05' in ids and 'M06' in ids:
        return False
    for m in measures:
        if m['id'] in ids and (not set(m['prerequisites']) <= set(ids)):
            return False
    return True

def portfolios(decisions, measures):
    for mask in range(1 << len(measures)):
        ids = [m['id'] for (i, m) in enumerate(measures) if mask & 1 << i]
        if admissible(ids, measures):
            yield evaluate(ids, decisions, measures)
