import time
import numpy as np
import trimesh
from .geometry import cap, shell, planes_from_cap, surface_samples
from .checks import cement, milling, wall, inherited, cement_max
from .io import result, digest
from .schema import validate_design
from .scoring import level1

def shape_diagnostics(design, reference):
    if design['kind'] == 'triangle_mesh_v1':
        ov = np.array(design['outer_mesh']['vertices_mm'])
        of = np.array(design['outer_mesh']['faces'])
    else:
        (ov, of) = cap(design['outer'])
    (rv, rf) = (reference['vertices'], reference['faces'])
    P = surface_samples(ov, of)
    Q = surface_samples(rv, rf)
    a = trimesh.Trimesh(ov, of, process=False)
    b = trimesh.Trimesh(rv, rf, process=False)
    da = trimesh.proximity.closest_point(b, P)[1]
    db = trimesh.proximity.closest_point(a, Q)[1]
    dist = np.r_[da, db]
    return dict(symmetric_sampled_RMS_mm=float(np.sqrt(np.mean(dist ** 2))), symmetric_sampled_p95_mm=float(np.quantile(dist, 0.95)), samples_per_surface=2048, reference='withheld native IOS tooth surface; heuristic cervical crop', clinical_quality='UNKNOWN', surface_sampling_error='UNKNOWN; not a Hausdorff certificate')

def evaluate(task, design, reference, cache):
    times = {}
    checks = {}
    tic = time.perf_counter()
    try:
        validate_design(design, task)
        if design['kind'] == 'triangle_mesh_v1':
            checks = {k: result('UNKNOWN', 'arbitrary mesh import', reason='No complete proof adapter for this representation') for k in task['requirements']['required_checks']}
            return dict(checks=checks, level1=level1(checks, task['requirements']['required_checks']), shape=shape_diagnostics(design, reference), timings_s=dict(mesh_import=time.perf_counter() - tic))
        mesh = shell(design['outer'], design['intaglio'])
        checks['geometry'] = result('PASS' if mesh.is_watertight and mesh.is_winding_consistent else 'INVALID', 'evaluator-built ring shell topology; geometric nesting tested separately by wall', watertight=bool(mesh.is_watertight), winding_consistent=bool(mesh.is_winding_consistent), n_faces=len(mesh.faces))
    except (KeyError, ValueError, TypeError) as e:
        checks['geometry'] = result('INVALID', 'submitted CAD representation', reason=str(e))
        return dict(checks=checks, level1=level1(checks, task['requirements']['required_checks']), shape=None, timings_s={})
    times['geometry'] = time.perf_counter() - tic
    (prepv, _) = cap(task['preparation']['profile'])
    (iv, _) = cap(design['intaglio'])
    try:
        planes = planes_from_cap(design['intaglio'])
    except ValueError as e:
        planes = None
        for name in ['cement', 'milling', 'insertion']:
            checks[name] = result('UNKNOWN', 'convex intaglio proof domain', reason=str(e))
    if planes is not None:
        operations = {'cement': lambda : cement.check(prepv.tolist(), planes, task['requirements']['cement_min_mm']), 'milling': lambda : milling.check(planes, iv.tolist(), task['requirements']['mill_radius_mm'], task['requirements']['mill_allowance_mm'], task['requirements']['allowed_exit_directions'][0]), 'insertion': lambda : inherited.insertion(planes)}
        for (name, op) in operations.items():
            relevant = {'cement': ['cement_min_mm'], 'milling': ['mill_radius_mm', 'mill_allowance_mm', 'allowed_exit_directions'], 'insertion': []}[name]
            key = digest([name, task['preparation'], design['intaglio'], {k: task['requirements'][k] for k in relevant}])
            tic = time.perf_counter()
            hit = key in cache
            if not hit:
                cache[key] = op()
            checks[name] = cache[key]
            times[name] = time.perf_counter() - tic
            times[name + '_cache_hit'] = hit
    tic = time.perf_counter()
    checks['wall'] = wall.check(design['outer'], design['intaglio'], task['requirements']['wall_min_mm'])
    times['wall'] = time.perf_counter() - tic
    checks['occlusion'] = result('UNKNOWN', 'full antagonist and neighbour nonpenetration', reason='No patient-linked registered antagonist and no trusted complete neighbour envelope in this task')
    if 'cement_max_mm' in task['requirements']:
        checks['cement_max'] = cement_max.check(task['preparation']['profile'], design['intaglio'], task['requirements']['cement_max_mm'])
    tic = time.perf_counter()
    shape = shape_diagnostics(design, reference)
    times['shape'] = time.perf_counter() - tic
    return dict(checks=checks, level1=level1(checks, task['requirements']['required_checks']), shape=shape, timings_s=times, shape_can_override_geometric_failure=False, fracture_prediction=result('UNKNOWN', 'published-measurement compatibility', reason='No matched physical setup for synthetic preparation'))
