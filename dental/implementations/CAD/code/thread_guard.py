import os
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'IGL_NUM_THREADS']:
    os.environ[k] = '1'
if hasattr(os, 'sched_getaffinity'):
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:4])
from scipy.optimize import linprog
from scipy.optimize._highspy._core import _Highs
_Highs.resetGlobalScheduler(True)
r = linprog([1.0], bounds=[(1.0, None)], options={'threads': 1})
if not r.success:
    raise RuntimeError('HiGHS scheduler initialization failed')
