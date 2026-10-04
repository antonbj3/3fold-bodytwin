"""Keep inherited LPs on the same one-thread HiGHS scheduler.

SciPy forwards the threads option to HiGHS. Mixing its default scheduler with
an explicit one-thread LP returned UNKNOWN_SOLVER in composition tests.
Only execution policy changes here; objectives and constraints are unchanged.
"""
import warnings
from scipy.optimize import linprog as scipy_linprog, OptimizeWarning

def linprog(*args, **kwargs):
    options = dict(kwargs.get('options') or {})
    options['threads'] = 1
    kwargs['options'] = options
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Unrecognized options detected: .*threads.*', category=OptimizeWarning)
        return scipy_linprog(*args, **kwargs)
