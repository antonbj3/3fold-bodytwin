from dental_release.paths import expand as _release_expand
import json, time, os, subprocess
from pathlib import Path
root = Path(__file__).resolve().parents[1]
pid = 751294
while True:
    try:
        rows = json.loads((root / 'raw/C_GENERATION.json').read_text())
    except (OSError, json.JSONDecodeError):
        rows = []
    if len(rows) >= 4:
        break
    if not Path('/proc/' + str(pid)).exists():
        break
    time.sleep(0.5)
try:
    os.kill(pid, 15)
except ProcessLookupError:
    pass
time.sleep(0.5)
env = os.environ.copy()
env.update(IGL_NUM_THREADS='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
with (root / 'raw/construct_c_bounded.log').open('w') as f:
    p = subprocess.run(['taskset', '-c', '0-3', _release_expand('@DENTAL_PYTHON@'), 'code/construct_c.py'], cwd=root, env=env, stdout=f, stderr=subprocess.STDOUT)
print('bounded C finished', p.returncode)
