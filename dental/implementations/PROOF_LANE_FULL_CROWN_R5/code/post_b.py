from dental_release.paths import expand as _release_expand
import time, subprocess, os
from pathlib import Path
r = Path(__file__).resolve().parents[1]
while not (r / 'FROZEN_PREDICTIONS_B.json').exists():
    time.sleep(0.5)
env = os.environ.copy()
env.update(IGL_NUM_THREADS='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
for (script, log) in [('score_fast.py', 'score_b.log'), ('x49_checks.py', 'x49_b.log')]:
    with (r / 'raw' / log).open('w') as f:
        p = subprocess.run(['taskset', '-c', '0-3', _release_expand('@DENTAL_PYTHON@'), 'code/' + script, 'B'], cwd=r, env=env, stdout=f, stderr=subprocess.STDOUT)
    if p.returncode:
        raise RuntimeError(script + ' failed ' + str(p.returncode))
print('B scoring and X49 done')
