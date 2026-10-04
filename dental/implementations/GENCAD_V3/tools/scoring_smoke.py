"""Three development cases only; exercise complete scorer before release sealing."""
import sys, shutil, time, subprocess, types
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from util import *
from sandbox import public_code, command
from integrity import safe_bytes, h, load_npz, strict_json
base = PAYLOAD / 'scoring_smoke'
inp = base / 'public'
out = base / 'output'
(inp / 'tasks').mkdir(parents=True, exist_ok=True)
(inp / 'scenes').mkdir(exist_ok=True)
out.mkdir(exist_ok=True)
paths = [p for p in sorted((PAYLOAD / 'public/tasks').glob('*.json')) if read(p)[0]['split'] == 'dev'][:3]
names = set(['SCORER_PARAMETERS.json', 'PREREG_R1.json', 'raw/COHORT_SUMMARY.json', 'payload/private/COHORT.json'])
for path in paths:
    names.add(str(path.relative_to(ROOT)))
    names.add('payload/private/references/' + path.stem + '.npz')
    shutil.copy2(path, inp / 'tasks' / path.name)
    for task in read(path):
        if task['status'] == 'READY':
            names.add('payload/public/' + task['geometry_file'])
            shutil.copy2(PAYLOAD / 'public' / task['geometry_file'], inp / task['geometry_file'])
files = {n: dict(sha256=sha(ROOT / n), bytes=(ROOT / n).stat().st_size) for n in sorted(names)}
freeze(ROOT / 'rounds/SCORING_SMOKE_INPUTS.json', files)

class Frozen:
    root = ROOT
    payload = PAYLOAD
    lock = {'closed_directories': ['payload/public', 'payload/private']}

    def __init__(self):
        self.files = files

    def location(self, n):
        return (PAYLOAD, n[8:]) if n.startswith('payload/') else (ROOT, n)

    def bytes(self, n):
        b = safe_bytes(*self.location(n))
        if h(b) != files[n]['sha256']:
            raise ValueError('smoke frozen drift')
        return b

    def json(self, n):
        return strict_json(self.bytes(n))

    def npz(self, n):
        return load_npz(self.bytes(n))

    def verify(self):
        for n in files:
            self.bytes(n)

    def tasks(self):
        for p in paths:
            n = str(p.relative_to(ROOT))
            yield (n, self.json(n))
public_code()
p = subprocess.run(command(out, inputs=inp), capture_output=True, text=True, timeout=300)
if p.returncode:
    raise RuntimeError(p.stderr)
anchor = digest(files)
sys.modules['release_anchor'] = types.SimpleNamespace(BENCHMARK_SHA256=anchor)
manifest = dict(benchmark_sha256=anchor, participants=['population', 'parametric', 'constraint_optimizer'], files={str(p.relative_to(out)): dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.rglob('*.npz')})
dump(out / 'FROZEN_PREDICTIONS.json', manifest)
freeze(ROOT / 'rounds/SCORING_SMOKE_PREDICTIONS.json', manifest)
from scorer import evaluate
result = evaluate(Frozen(), out, sha(out / 'FROZEN_PREDICTIONS.json'), base / 'evaluation')
dump(ROOT / 'raw/SCORER_DEVELOPMENT_SMOKE.json', result)
print(result['totals'])
shutil.rmtree(base)
