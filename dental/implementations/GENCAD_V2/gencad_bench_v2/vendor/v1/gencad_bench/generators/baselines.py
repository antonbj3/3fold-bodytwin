"""Equal information: each function receives the identical public task object."""
from copy import deepcopy
import numpy as np
from ..geometry import dome_profile, planes_from_cap
from ..checks.exact import q, norm2, vec, sqrt_bounds

def _design(task, outer, name):
    return dict(schema_version='0.1', kind='radial_crown_v1', task_id=task['task_id'], units='mm', frame=task['frame'], outer=outer, intaglio=deepcopy(task['public']['intaglio_template']), generator=name, predictions=None)

def mirror(task):
    return _design(task, deepcopy(task['public']['contralateral_profile']), 'contralateral_mirror')

def population(task):
    return _design(task, deepcopy(task['public']['population_profile']), 'population_mean')

def parametric(task):
    inner = deepcopy(task['public']['intaglio_template'])
    try:
        pl = planes_from_cap(inner)
        support = min((float(q(p[3]) / sqrt_bounds(norm2(vec(p[:3])))[1]) for p in pl))
    except ValueError:
        support = min(np.min(inner['radii_mm']), inner['apex_mm'])
    factor = 1 + (task['requirements']['wall_min_mm'] + 0.02) / support
    outer = dict(z_mm=(np.asarray(inner['z_mm']) * factor).round(6).tolist(), radii_mm=(np.asarray(inner['radii_mm']) * factor).round(6).tolist(), apex_mm=round(inner['apex_mm'] * factor, 6))
    return _design(task, outer, 'parametric_IFU')
