"""Run the frozen probe code with an explicit one-thread HiGHS limit."""
import hashlib
import json
import platform
import runpy
import warnings
from pathlib import Path
import numpy
import scipy
import scipy.optimize
root = Path(__file__).resolve().parent
original = scipy.optimize.linprog
calls = []

def limited_linprog(*args, **kwargs):
    options = dict(kwargs.pop('options', {}) or {})
    options['threads'] = 1
    kwargs['options'] = options
    calls.append({'method': kwargs.get('method'), 'threads': 1})
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Unrecognized options detected.*', category=scipy.optimize.OptimizeWarning)
        return original(*args, **kwargs)
scipy.optimize.linprog = limited_linprog
runpy.run_path(str(root / 'review.py'), run_name='__main__')
(root / 'EXECUTION_ENVIRONMENT.json').write_text(json.dumps({'python': platform.python_version(), 'numpy': numpy.__version__, 'scipy': scipy.__version__, 'HiGHS_linprog_calls': calls, 'BLAS_OpenMP_limit': 1, 'review_sha256': hashlib.sha256((root / 'review.py').read_bytes()).hexdigest(), 'execution_wrapper_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'scope': 'Execution resource limit only; source functions, summaries, thresholds and metrics unchanged'}, indent=2) + '\n')
