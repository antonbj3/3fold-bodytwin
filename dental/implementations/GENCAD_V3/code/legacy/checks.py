import numpy as np
from .geometry import shell, area_weights, section_weights

def validate(task, design):
    if not isinstance(design, dict):
        raise ValueError('object required')
    if design.get('task_id') != task['task_id']:
        raise ValueError('wrong task')
    if design.get('status') == 'ABSTAIN':
        if not isinstance(design.get('reason'), str):
            raise ValueError('reason required')
        return None
    if design.get('units') != 'mm' or design.get('frame') != task['frame']:
        raise ValueError('units/frame')
    if set(design) - {'task_id', 'status', 'units', 'frame', 'outer_vertices', 'inner_vertices', 'faces', 'diagnostics'}:
        raise ValueError('unknown output field')
    xy = np.asarray(task['xy'], float)
    f = np.asarray(task['faces'], int)
    out = np.asarray(design.get('outer_vertices'), float)
    inn = np.asarray(design.get('inner_vertices'), float)
    fs = np.asarray(design.get('faces'))
    if out.shape != (len(xy), 3) or inn.shape != out.shape:
        raise ValueError('coordinate dimensions')
    if fs.shape != f.shape or fs.dtype.kind not in 'iu' or (not np.array_equal(fs, f)):
        raise ValueError('topology/indices')
    if not np.isfinite(out).all() or not np.isfinite(inn).all():
        raise ValueError('nonfinite')
    if not np.allclose(out[:, :2], xy, rtol=0, atol=1e-09) or not np.allclose(inn[:, :2], xy, rtol=0, atol=1e-09):
        raise ValueError('footprint changed')
    if np.max(np.abs(out[:, 2])) > 100 or np.max(np.abs(inn[:, 2])) > 100:
        raise ValueError('extent')
    return (out[:, 2], inn[:, 2])

def all_checks(task, design):
    try:
        parsed = validate(task, design)
    except (TypeError, ValueError, KeyError, OverflowError) as e:
        return dict(validity='INVALID', digital_score=0, reason=str(e), checks={})
    if parsed is None:
        return dict(validity='ABSTAIN', digital_score=None, checks={})
    (z, inner) = parsed
    xy = np.asarray(task['xy'])
    faces = np.asarray(task['faces'])
    prep = np.asarray(task['preparation_z'])
    tol = 1e-06
    req = task['requirements']
    checks = {}

    def put(k, value, passed, scope):
        checks[k] = dict(value=float(value), status='PASS' if passed else 'FAIL', scope=scope)
    put('nesting', np.min(z - inner), np.min(z - inner) > tol, 'outer graph strictly above intaglio at every affine node')
    gap = inner - prep
    put('film_min', gap.min(), gap.min() >= req['film_min_mm'] - tol, 'axial nominal film, not measured normal cement thickness')
    put('film_max', gap.max(), gap.max() <= req['film_max_mm'] + tol, 'axial nominal film')
    t = np.c_[xy, inner][faces]
    a = t[:, 1] - t[:, 0]
    b = t[:, 2] - t[:, 0]
    det = np.cross(a[:, :2], b[:, :2])
    gx = (a[:, 2] * b[:, 1] - b[:, 2] * a[:, 1]) / det
    gy = (a[:, 0] * b[:, 2] - b[:, 0] * a[:, 2]) / det
    L = float(np.max(np.hypot(gx, gy)))
    upper = float(np.min(z - inner))
    lower = upper / np.sqrt(1 + L * L)
    wall_status = 'PASS' if lower >= req['wall_mm'] - tol else 'FAIL' if upper < req['wall_mm'] - tol else 'UNKNOWN'
    checks['wall'] = dict(status=wall_status, lower_mm=lower, upper_mm=upper, scope='entire Lipschitz graph Euclidean separation bound; vertical rim excluded')
    A = task['A']
    b = np.asarray(task['obstacle_b'], float)
    if A.shape[0]:
        slack = b - A @ z
        put('antagonist', slack.min(), slack.min() >= -tol, 'all projected triangle intersection vertices, 1e-6 mm arithmetic tolerance; native scanner error unknown')
    else:
        checks['antagonist'] = dict(status='UNKNOWN', scope='No overlapping measured antagonist surface along task axis')
    checks['insertion'] = dict(status='PASS', scope='straight removal from horizontal preparation in fixed footprint surface model; axial retention and remote obstacles UNKNOWN')
    if task['family'] == 'bridge3':
        vals = [section_weights(xy, faces, x) @ (z - inner) for x in task['connector_x']]
        put('connector_area', min(vals), min(vals) >= req['connector_mm2'] - tol, 'actual shell cross-section, not reported area')
    if task['family'] == 'implant_crown':
        axis = np.asarray(task['implant_axis'], float)
        angle = float(np.degrees(np.arccos(np.clip(axis[2] / np.linalg.norm(axis), -1, 1))))
        put('channel_angle', angle, angle <= req['channel_max_deg'] + tol, 'geometric open channel vs specified implant-axis scenario; connection absent')
    if task['family'] == 'lattice_onlay':
        put('strut_width', task['minimum_strut_mm'], task['minimum_strut_mm'] >= req['strut_mm'] - tol, 'locked planar strut neck; fatigue and biocompatibility UNKNOWN')
    checks['milling'] = dict(status='PASS' if L < 1e-10 else 'UNKNOWN', scope='horizontal intaglio patch, unbounded approach halfspace, ideal ball only; boundary walls and shaft not certified')
    statuses = [v['status'] for (k, v) in checks.items() if k not in (['antagonist'] if task['family'] == 'veneer' else [])]
    verdict = 'FAIL' if 'FAIL' in statuses else 'UNKNOWN' if 'UNKNOWN' in statuses else 'PASS'
    return dict(validity=verdict, digital_score=0 if verdict == 'FAIL' else 1 if verdict == 'PASS' else None, complete_restoration_eligibility='UNKNOWN_MISSING_AXIAL_RETENTION_EMERGENCE_AND_PHYSICAL_VALIDATION', checks=checks)

def feasibility(task):
    """Independent necessary obstruction; LP additionally tests the locked affine model."""
    from scipy.optimize import linprog
    from scipy.sparse import csr_matrix, vstack
    xy = np.asarray(task['xy'])
    f = np.asarray(task['faces'])
    A = task['A']
    b = np.asarray(task['obstacle_b'])
    r = task['requirements']
    prep = np.asarray(task['preparation_z'])
    lower = prep + r['film_min_mm'] + r['wall_mm']
    constraints = []
    rhs = []
    if A.shape[0]:
        constraints.append(csr_matrix(A))
        rhs.extend(b)
    if task['family'] == 'bridge3':
        for x in task['connector_x']:
            w = section_weights(xy, f, x)
            constraints.append(csr_matrix(-w[None, :]))
            rhs.append(-r['connector_mm2'] - w @ (prep + r['film_min_mm']))
    if task['family'] == 'implant_crown':
        a = np.asarray(task['implant_axis'])
        ang = np.degrees(np.arccos(np.clip(a[2] / np.linalg.norm(a), -1, 1)))
        if ang > r['channel_max_deg'] + 1e-06:
            return dict(status='INFEASIBLE', witness='locked straight channel cannot meet implant access angle', margin_deg=float(r['channel_max_deg'] - ang))
    if task['family'] == 'lattice_onlay' and task['minimum_strut_mm'] < r['strut_mm'] - 1e-06:
        return dict(status='INFEASIBLE', witness='locked strut neck too narrow')
    if A.shape[0]:
        margins = b - A @ lower
        if margins.min() < -1e-06:
            return dict(status='INFEASIBLE', witness='nonnegative barycentric wall/film lower bound exceeds native antagonist ceiling', margin_mm=float(margins.min()), row=int(np.argmin(margins)))
    sol = linprog(np.zeros(len(lower)), A_ub=vstack(constraints) if constraints else None, b_ub=np.asarray(rhs) if rhs else None, bounds=[(v, None) for v in lower], method='highs', options={'threads': 1})
    if sol.success:
        return dict(status='FEASIBLE', witness='LP construction in restricted model', z=sol.x.tolist())
    return dict(status='INFEASIBLE' if sol.status == 2 else 'UNKNOWN', witness=sol.message)
