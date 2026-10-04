"""Five actual executable participants. Adapters name the exact inherited operation."""
import sys, time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, eye, hstack, vstack
from scipy.optimize import linprog
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator, RegularGridInterpolator
sys.path.insert(0, str(Path(__file__).resolve().parent / 'vendor'))
from obstacle import hessian, qp
import crown_fit_geometry as X1
from .geometry import section_weights

def submit(t, z, inner=None, diagnostics=None):
    if inner is None:
        inner = np.asarray(t['preparation_z']) + 0.08
    return dict(task_id=t['task_id'], status='DESIGN', units='mm', frame=t['frame'], outer_vertices=np.c_[t['xy'], z].tolist(), inner_vertices=np.c_[t['xy'], inner].tolist(), faces=np.asarray(t['faces']).tolist(), diagnostics=diagnostics or {})

def population(t):
    return submit(t, np.asarray(t['prior']))

def parametric(t):
    return submit(t, np.asarray(t['preparation_z']) + 0.08 + t['requirements']['wall_mm'] + 0.05)
_x1_cache = {}

def crown_loop(t):
    """Execute X1B's unchanged closing/redistance on a prior roof SDF; no FE is implied."""
    key = (t['geometry_file'], t['requirements']['tool_radius_mm'])
    if key in _x1_cache:
        return submit(t, _x1_cache[key], diagnostics=dict(source_operation='X1B crown_fit_geometry.closing', cache_hit=True))
    xy = np.asarray(t['xy'])
    z = np.asarray(t['prior'])
    rad = t['requirements']['tool_radius_mm']
    h = 0.5
    g = X1.Grid([xy[:, 0].min() - 2 * rad, xy[:, 1].min() - 2 * rad, z.min() - 2.5], [xy[:, 0].max() + 2 * rad, xy[:, 1].max() + 2 * rad, z.max() + 2.5], h)
    (x, y, zz) = g.mesh()
    sample = np.c_[x[:, :, 0].ravel(), y[:, :, 0].ravel()]
    f = LinearNDInterpolator(xy, z)
    roof = f(sample)
    bad = ~np.isfinite(roof)
    if bad.any():
        roof[bad] = NearestNDInterpolator(xy, z)(sample[bad])
    phi = (zz - roof.reshape(x.shape[:2])[:, :, None]).astype(np.float32)
    phi = X1.redistance_lowmem(phi, g, levels=(0.0, rad), workers=1)
    result = X1.closing(phi, g, rad, workers=1)
    axes = g.axes()
    heights = np.full(result.shape[:2], np.nan)
    for i in range(result.shape[0]):
        for j in range(result.shape[1]):
            hits = np.flatnonzero((result[i, j, :-1] <= 0) & (result[i, j, 1:] > 0))
            if len(hits):
                k = hits[-1]
                (a, b) = result[i, j, k:k + 2]
                heights[i, j] = axes[2][k] + h * float(-a / (b - a))
    if not np.isfinite(heights).all():
        raise ValueError('X1B closing lost roof crossing')
    out = RegularGridInterpolator(axes[:2], heights)(xy)
    _x1_cache[key] = out
    return submit(t, out, diagnostics=dict(source_operation='X1B crown_fit_geometry.closing', grid_h_mm=h, grid_shape=g.shape, cache_hit=False, physical_FE_run=False))

def antagonist(t):
    A = t['A']
    b = np.asarray(t['obstacle_b'])
    prior = np.asarray(t['prior'])
    if not A.shape[0]:
        return submit(t, prior, diagnostics=dict(source_operation='X18 obstacle.qp', no_obstacle=True))
    H = hessian(np.asarray(t['faces']), len(prior), 0.15)
    (z, info) = qp(A, b, H, prior + 0.2, maxiter=350, gtol=1e-08)
    violation = max(0, float(np.max(A @ z - b)))
    z -= violation + 1e-08
    return submit(t, z, diagnostics=dict(source_operation='X18 obstacle.qp plus global continuous-constraint repair', qp=info, repair_mm=violation))

def optimize(t):
    """Constrained L1 morphology objective; all active geometric constraints enter the solve."""
    xy = np.asarray(t['xy'])
    prior = np.asarray(t['prior'])
    n = len(prior)
    r = t['requirements']
    inner = np.asarray(t['preparation_z']) + r['film_min_mm']
    lower = inner + r['wall_mm']
    if t['family'] == 'implant_crown':
        a = np.asarray(t['implant_axis'])
        angle = np.degrees(np.arccos(np.clip(a[2] / np.linalg.norm(a), -1, 1)))
        if angle > r['channel_max_deg'] + 1e-06:
            return dict(task_id=t['task_id'], status='ABSTAIN', reason='fixed through-channel exceeds access-angle constraint')
    if t['family'] == 'lattice_onlay' and t['minimum_strut_mm'] < r['strut_mm'] - 1e-06:
        return dict(task_id=t['task_id'], status='ABSTAIN', reason='fixed lattice neck below requested strut width')
    I = eye(n, format='csr')
    A = t['A']
    rows = [hstack([A, csr_matrix((A.shape[0], n))]), hstack([I, -I]), hstack([-I, -I])]
    rhs = [np.asarray(t['obstacle_b']), prior, -prior]
    if t['family'] == 'bridge3':
        for x in t['connector_x']:
            w = section_weights(xy, np.asarray(t['faces']), x)
            rows.append(csr_matrix(np.r_[-w, np.zeros(n)][None, :]))
            rhs.append(np.array([-r['connector_mm2'] - w @ inner]))
    sol = linprog(np.r_[np.zeros(n), np.asarray(t['weights'])], A_ub=vstack(rows), b_ub=np.concatenate(rhs), bounds=[(v, None) for v in lower] + [(0, None)] * n, method='highs', options={'threads': 1})
    if not sol.success:
        return dict(task_id=t['task_id'], status='ABSTAIN', reason='LP ' + str(sol.status) + ': ' + sol.message)
    return submit(t, sol.x[:n], inner, dict(source_operation='new sparse constrained surface LP', objective=float(sol.fun), iterations=int(sol.nit)))
PARTICIPANTS = {'population': population, 'parametric': parametric, 'X1B_morphology': crown_loop, 'X18_antagonist': antagonist, 'constraint_optimizer': optimize}
