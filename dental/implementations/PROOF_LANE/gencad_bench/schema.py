import math
from .geometry import cap

def validate_task(t):
    required = ['schema_version', 'task_id', 'geometry_id', 'subject_group', 'split', 'units', 'frame', 'geometry_reference', 'preparation', 'material', 'requirements', 'truth_tiers', 'public']
    if any((k not in t for k in required)):
        raise ValueError('Missing task field')
    if t['schema_version'] != '0.1' or t['units'] != 'mm':
        raise ValueError('Unsupported version/units')
    if t['split'] not in ['training', 'development', 'test']:
        raise ValueError('Invalid split')
    for k in ['wall_min_mm', 'cement_min_mm', 'mill_radius_mm', 'mill_allowance_mm']:
        v = t['requirements'][k]
        if isinstance(v, bool) or not isinstance(v, (int, float)) or (not math.isfinite(v)) or (v < 0):
            raise ValueError('Invalid requirement ' + k)
    if t['requirements']['mill_radius_mm'] == 0:
        raise ValueError('Zero cutter radius')
    if 'cement_max_mm' in t['requirements']:
        v = t['requirements']['cement_max_mm']
        if isinstance(v, bool) or not isinstance(v, (int, float)) or (not math.isfinite(v)) or (v < 0):
            raise ValueError('Invalid maximum-film requirement')
    if not t['requirements']['required_checks']:
        raise ValueError('No required checks')
    return t

def validate_design(d, t):
    if d.get('task_id') != t['task_id'] or d.get('units') != t['units'] or d.get('frame') != t['frame']:
        raise ValueError('Design/task identity, unit or frame mismatch')
    if d.get('kind') == 'radial_crown_v1':
        cap(d['outer'])
        cap(d['intaglio'])
    elif d.get('kind') == 'triangle_mesh_v1':
        import numpy as np
        for key in ['outer_mesh', 'intaglio_mesh']:
            v = np.asarray(d[key]['vertices_mm'], float)
            f = np.asarray(d[key]['faces'])
            if v.ndim != 2 or v.shape[1] != 3 or (not len(v)) or (not np.isfinite(v).all()):
                raise ValueError('Invalid mesh vertices')
            if f.ndim != 2 or f.shape[1] != 3 or (not len(f)) or (not np.issubdtype(f.dtype, np.integer)) or (f.min() < 0) or (f.max() >= len(v)):
                raise ValueError('Invalid mesh faces')
    else:
        raise ValueError('Unsupported geometry representation')
    return d

def validate_split(tasks):
    seen = {}
    training = set()
    for t in tasks:
        s = t['subject_group']
        if s in seen and seen[s] != t['split']:
            raise ValueError('Subject split leakage')
        seen[s] = t['split']
        training.update(t['public']['population_training_groups'])
    if training & set(seen):
        raise ValueError('Template/test leakage')
    return True
