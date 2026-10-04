"""Development-only runtime probe; not a heldout score or prediction freeze."""
import sys, shutil, time, subprocess
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from util import *
from sandbox import public_code, command
base = PAYLOAD / 'engineering_smoke'
inp = base / 'public'
out = base / 'output'
(inp / 'tasks').mkdir(parents=True, exist_ok=True)
(inp / 'scenes').mkdir(exist_ok=True)
out.mkdir(exist_ok=True)
for path in [p for p in sorted((PAYLOAD / 'public/tasks').glob('*.json')) if read(p)[0]['split'] == 'dev'][:3]:
    shutil.copy2(path, inp / 'tasks' / path.name)
    for task in read(path):
        if task['status'] == 'READY':
            shutil.copy2(PAYLOAD / 'public' / task['geometry_file'], inp / task['geometry_file'])
public_code()
start = time.perf_counter()
p = subprocess.run(command(out, inputs=inp), capture_output=True, text=True, timeout=300)
dump(ROOT / 'raw/GENERATOR_SMOKE.json', dict(exit_code=p.returncode, seconds=time.perf_counter() - start, stdout=p.stdout, stderr=p.stderr, scope='three development cases, generation only; no target read or score'))
if p.returncode:
    raise RuntimeError(p.stderr)
shutil.rmtree(base)
