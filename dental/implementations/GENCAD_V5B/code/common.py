from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, hashlib, datetime, sys, time, os
import numpy as np
P = Path(__file__).resolve().parents[1]
D = P.parent
V5 = D / _release_expand('GENCAD_V5')
V4 = D / _release_expand('GENCAD_V4')
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V5B'))
OLD_DATA = DATA.parent / _release_expand('GENCAD_V5')
sys.dont_write_bytecode = True

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def clean(x):
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, np.generic):
        return clean(x.item())
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (tuple, list)):
        return [clean(v) for v in x]
    if isinstance(x, float) and (not np.isfinite(x)):
        return None
    return x

def read(p):
    return json.loads(Path(p).read_text())

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def freeze(p, x):
    p = Path(p)
    x = clean(x)
    if p.exists():
        old = read(p)
        if {k: v for (k, v) in old.items() if k != 'frozen_utc'} != x:
            raise ValueError('Frozen artifact drift: ' + str(p))
    else:
        dump(p, dict(frozen_utc=now(), **x))
    p.with_suffix('.sha256').write_text(sha(p) + '\n')

def load_npz(p):
    with np.load(p, allow_pickle=False) as z:
        return dict(z)

def state(phase, gate, next_op):
    dump(P / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-gencad-v5b', status='ACTIVE', phase=phase, latest_gate=gate, next_operation=next_op, updated_utc=now(), claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW'))

def load_module(name, p):
    import importlib.util
    s = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m

def cohort():
    rows = read(V5 / 'raw/FUNCTIONAL_ROWS.json') + read(V5 / 'raw/R6_EXTERNAL_SUPPORT.json')['rows']
    out = []
    for r in rows:
        q = dict(r)
        track = q['track']
        name = q['participant']
        key = q['key']
        if track in ['R3', 'R5']:
            path = V4 / 'payload' / ('whole_predictions' if track == 'R3' else 'prep_predictions') / name / key / 'mesh.npz'
        elif track == 'ToothCraft':
            path = OLD_DATA / 'external_meshes' / key / 'mesh.npz'
        elif track == 'R6':
            path = OLD_DATA / 'external_support' / name.replace('ToothCraft_', '') / key / 'mesh.npz'
        else:
            raise ValueError(track)
        q.update(mesh_path=str(path), uid=track + '__' + name + '__' + key, prep_path=str(V4 / 'payload/whole_inputs' / key / 'preparation.npz'))
        out.append(q)
    return out
