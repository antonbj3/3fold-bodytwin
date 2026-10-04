import math
from .io import result

def level1(checks, required):
    if not required:
        return dict(status='INVALID', score=0, eligible=False, reason='Empty requirement set')
    rows = [checks.get(k, dict(status='UNKNOWN', reason='Missing check')) for k in required]
    if any((x['status'] == 'INVALID' for x in rows)):
        return dict(status='INVALID', score=0, eligible=False)
    if any((x['status'] == 'FAIL' for x in rows)):
        return dict(status='FAIL', score=0, eligible=False)
    if any((x['status'] != 'PASS' for x in rows)):
        return dict(status='UNKNOWN', score=None, eligible=False)
    return dict(status='PASS', score=1, eligible=True)

def interval_score(lo, hi, y, alpha=0.05):
    if not all((math.isfinite(v) for v in [lo, hi, y, alpha])) or hi < lo or (not 0 < alpha < 1):
        raise ValueError('Invalid predictive interval')
    return hi - lo + 2 / alpha * max(lo - y, 0) + 2 / alpha * max(y - hi, 0)

def setup_match(task, row):
    required = ['material_product', 'specimen', 'thickness_mm', 'support', 'cement', 'load_angle_deg', 'ageing']
    return all((task.get(k) is not None and row.get(k) is not None and (task[k] == row[k]) for k in required))
