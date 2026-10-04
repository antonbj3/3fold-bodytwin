"""GenCAD dictionary -> dictionary plugin. Public geometry only."""
import os, time
import numpy as np
from scipy.optimize import minimize, LinearConstraint
from scipy.sparse import csr_matrix
from legacy.generators import submit
from legacy.geometry import section_weights
from fast_geometry import optimize as morphology
from mechanics import roof_model

def generate(t):
    start = time.perf_counter()
    budget = float(os.environ.get('X60_VOLUME_BUDGET', '1.0'))
    base = morphology(t)
    if base['status'] != 'DESIGN':
        return base
    z = np.asarray(base['outer_vertices'])[:, 2]
    inner = np.asarray(base['inner_vertices'])[:, 2]
    (model, M) = roof_model(t, z, inner)
    np_ = M.shape[1]
    zero = np.zeros(np_)
    initial = model.evaluate(zero)
    w = np.asarray(t['weights'])
    r = t['requirements']
    rows = [M]
    lo = [inner + r['wall_mm'] - z]
    hi = [np.full(len(z), np.inf)]
    if t['A'].shape[0]:
        rows.append(np.asarray(t['A'] @ M))
        lo.append(np.full(t['A'].shape[0], -np.inf))
        hi.append(np.asarray(t['obstacle_b']) - t['A'] @ z)
    rows.append(M)
    lo.append(np.full(len(z), -0.25))
    hi.append(np.full(len(z), 0.25))
    rows.append((w @ M)[None, :])
    lo.append(np.array([-np.inf]))
    hi.append(np.array([(budget - 1) * initial['volume_mm3']]))
    if t['family'] == 'bridge3':
        for section in t['connector_x']:
            ww = section_weights(t['xy'], t['faces'], section)
            rows.append((ww @ M)[None, :])
            lo.append(np.array([r['connector_mm2'] - ww @ (z - inner)]))
            hi.append(np.array([np.inf]))
    L = np.vstack(rows)
    lower = np.concatenate(lo)
    upper = np.concatenate(hi)
    usehi = np.isfinite(upper)
    uselo = np.isfinite(lower)
    G = np.vstack([L[usehi], -L[uselo]])
    b = np.r_[upper[usehi], -lower[uselo]]
    norm = np.maximum(np.linalg.norm(G, axis=1), 1e-12)
    G = G / norm[:, None]
    b = b / norm
    lc = LinearConstraint(G, -np.inf, b)
    history = []

    def fun(x):
        q = model.evaluate(x)
        return (q['J'], q['gradient'])

    def cb(x):
        q = model.evaluate(x)
        history.append(dict(J=q['J'], volume_mm3=q['volume_mm3']))
    sol = minimize(fun, zero, jac=True, method='SLSQP', constraints=[lc], bounds=[(-0.25, 0.25)] * np_, callback=cb, options=dict(maxiter=35, ftol=1e-08))
    final = model.evaluate(sol.x)
    residual = float(np.max(G @ sol.x - b))
    improvement = float(np.exp(initial['J'] - final['J']))
    return submit(t, z + M @ sol.x, inner, dict(operation='linear tetra elasticity + exact local shape adjoint + constrained SLSQP', optimizer_success=bool(sol.success), optimizer_message=str(sol.message), iterations=int(sol.nit), FEM_queries=model.calls, volume_budget=budget, base_volume_mm3=initial['volume_mm3'], volume_mm3=final['volume_mm3'], model_force_quantile_ratio=improvement, initial_log_hazard_over_m=initial['J'], final_log_hazard_over_m=final['J'], equilibrium_relative_residual=final['equilibrium_relative_residual'], constraint_scaled_violation=residual, coefficients_mm=sol.x.tolist(), seconds=time.perf_counter() - start, calibrated_force05_N=None, absolute_force_status='UNKNOWN_UNMATCHED_BATCH_SUPPORT_CONTACT_GEOMETRY', rigorous_design_box_enclosure='MISSING', local_history=history))
