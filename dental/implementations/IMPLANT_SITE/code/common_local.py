from dental_release.paths import expand as _release_expand
import os
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[k] = '1'
import json, pathlib, hashlib, datetime, csv, importlib.util, sys
sys.dont_write_bytecode = True
import numpy as np
P = pathlib.Path(__file__).resolve().parents[1]
D = P.parent
DATA = pathlib.Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE-combine-implant-site'))

def load(p):
    return json.loads(pathlib.Path(p).read_text())

def lines(p):
    return [json.loads(l) for l in pathlib.Path(p).read_text().splitlines() if l.strip()]

def dump(p, x):
    p = pathlib.Path(p)
    p.parent.mkdir(exist_ok=True, parents=True)
    p.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def module(name, path):
    path = pathlib.Path(path)
    sys.path.insert(0, str(path.parent))
    s = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m
X8 = D / _release_expand('X8')
X73 = D / _release_expand('X73')
X75 = D / _release_expand('X75')
X87 = D / _release_expand('X87')

def csvout(p, rows):
    with pathlib.Path(p).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def state(stage, gate, next_operation, **kw):
    dump(P / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-combine-implant-site', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), stage=stage, latest_gate=gate, next_operation=next_operation, **kw))

def check_frozen():
    for p in P.glob('PREREG_*.json'):
        assert sha(p) == p.with_suffix(p.suffix + '.sha256').read_text().strip(), str(p)

def register_inputs(paths):
    dest = P / 'INPUT_MANIFEST.json'
    existing = load(dest) if dest.exists() else []
    by = {r['path']: r for r in existing}
    for p in paths:
        p = pathlib.Path(p).resolve()
        v = dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size)
        if str(p) in by:
            assert by[str(p)] == v, 'Source changed: ' + str(p)
        by[str(p)] = v
    dump(dest, list(by.values()))

def verify_inputs():
    for m in load(P / 'INPUT_MANIFEST.json'):
        assert sha(m['path']) == m['sha256'], 'Source drift: ' + m['path']

def boot(rows, field, seed=20261004, B=2000):
    grouped = {}
    for r in rows:
        grouped.setdefault(r['image_group'], []).append(r[field])
    a = np.array([(sum(v), len(v)) for v in grouped.values()], float)
    if not len(a):
        return dict(n=0, events=0, fraction=None, ci95=[None, None], image_groups=0)
    rng = np.random.default_rng(seed)
    ix = rng.integers(0, len(a), (B, len(a)))
    s = a[ix].sum(1)
    dist = s[:, 0] / s[:, 1]
    return dict(n=int(a[:, 1].sum()), events=int(a[:, 0].sum()), fraction=float(a[:, 0].sum() / a[:, 1].sum()), ci95=np.quantile(dist, [0.025, 0.975]).tolist(), image_groups=len(a), CI_scope='Conditional percentile image-cluster bootstrap; convenience cohort, guide moments treated as fixed; physical error not included')
