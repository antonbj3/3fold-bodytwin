"""Regional actuation calibration and interval inverse. Physical validity is separate."""
from common import np, sha, clean
from scipy.optimize import linprog

def fit_response(probes):
    if probes.get('schema') != 'regional-actuation-probes-v1':
        raise ValueError('Unknown probe schema')
    if probes.get('units') != {'height': 'mm', 'force': 'N'}:
        raise ValueError('Units must be explicit')
    if probes.get('observation_operator') != 'calibrated_regional_axial_reactions':
        return dict(status='UNKNOWN', reason='Missing identified regional reaction observation operator')
    if probes.get('linear_support_closure') is not True or probes.get('matched_load_rate_history') is not True:
        return dict(status='UNKNOWN', reason='Missing linear closure or matched wrench/rate/history')
    if probes.get('height_error_mm') != 0.0:
        return dict(status='UNKNOWN', reason='This version requires exact simulated probe heights; real height uncertainty needs denominator enclosure')
    y0 = np.asarray(probes['baseline_force_N'], float)
    yp = np.asarray(probes['plus_force_N'], float)
    ym = np.asarray(probes['minus_force_N'], float)
    s = float(probes['step_mm'])
    eps = float(probes['raw_force_channel_bound_N'])
    if not np.isfinite(s) or s <= 0 or (not np.isfinite(eps)) or (eps < 0):
        raise ValueError('Invalid uncertainty or step')
    if len(y0.shape) != 1 or yp.shape != ym.shape or yp.shape[1] != len(y0):
        raise ValueError('Probe dimensions')
    if not np.isfinite(np.r_[y0, yp.ravel(), ym.ravel()]).all():
        raise ValueError('Nonfinite force')
    J = ((yp - ym) / (2 * s)).T
    eta = float(np.sqrt(len(y0)) * eps)
    beta = np.full(J.shape[1], eta / s)
    return dict(status='CONDITIONAL_CALIBRATION', f0=y0.tolist(), J_N_per_mm=J.tolist(), eta_N=eta, beta_N_per_mm=beta.tolist(), geometry_sha256=probes['geometry_sha256'], basis_sha256=probes['basis_sha256'], case=probes['case'], frame=probes['frame'], kind=probes['kind'], step_mm=s, nonlinear_remainder='ZERO_ONLY_UNDER_DECLARED_LINEAR_SUPPORT', resolution='PER_SURFACE_REGION', physical_validity='UNKNOWN' if probes['kind'] == 'SIMULATION' else 'CONDITIONAL_ON_EXTERNAL_CLOSURE_AND_OBSERVATION')

def predict(model, a):
    if model.get('status') != 'CONDITIONAL_CALIBRATION':
        return dict(status='UNKNOWN', reason=model.get('reason', 'No calibration'))
    a = np.asarray(a, float)
    J = np.asarray(model['J_N_per_mm'])
    f0 = np.asarray(model['f0'])
    if a.shape != (J.shape[1],) or not np.isfinite(a).all() or np.max(abs(a)) > 0.05 + 1e-12:
        raise ValueError('Actuation outside box')
    R = float(model['eta_N'] + np.asarray(model['beta_N_per_mm']) @ abs(a))
    f = f0 + J @ a
    intervals = np.c_[f - R, f + R]
    corridor = bool(np.min(f0 - model['eta_N']) > 0 and intervals[:, 0].min() > 0)
    return dict(status='CONDITIONAL_FORCE_INTERVAL' if corridor else 'UNKNOWN_CONTACT_CHANGE', force_N=f.tolist(), force_interval_N=intervals.tolist(), joint_l2_error_radius_N=R, positive_contact_path=corridor, enclosure='Analytic deterministic real-arithmetic bound under affine full contact; floating enclosure MISSING', resolution='PER_SURFACE_REGION', physical_validity=model['physical_validity'])

def inverse(model, B, weights, wall_allowance, band=(7.0, 9.0), additional_geometry=True):
    if model.get('status') != 'CONDITIONAL_CALIBRATION':
        return dict(status='UNKNOWN', reason=model.get('reason', 'No calibration'))
    J = np.asarray(model['J_N_per_mm'])
    f0 = np.asarray(model['f0'])
    eta = model['eta_N']
    beta = np.asarray(model['beta_N_per_mm'])
    m = J.shape[1]
    G = []
    rhs = []
    for i in range(2):
        G.extend([-J[i] + beta, J[i] + beta])
        rhs.extend([band[1] - eta - f0[i], f0[i] - eta - band[0]])
    G.extend(J + beta[None, :])
    rhs.extend(f0 - eta - 1e-06)
    if additional_geometry:
        G.extend(B)
        rhs.extend(wall_allowance)
    sol = linprog(np.asarray(weights) @ B, A_ub=np.asarray(G), b_ub=np.asarray(rhs), bounds=[(0, 0.05)] * m, method='highs', options={'threads': 1})
    if not sol.success:
        return dict(status='INFEASIBLE_IN_DECLARED_MODEL' if sol.status == 2 else 'UNKNOWN_SOLVER', message=sol.message, solver_status=sol.status)
    a = -sol.x
    ans = predict(model, a)
    good = ans['positive_contact_path'] and all((lo >= band[0] - 1e-07 and hi <= band[1] + 1e-07 for (lo, hi) in ans['force_interval_N'][:2]))
    return dict(status='CONDITIONAL_TARGET_DESIGN' if good else 'UNKNOWN', coefficients_mm=a.tolist(), at_selected_height=ans, target_region_band_N=list(band), objective_mean_removal_mm=float(sol.fun / np.sum(weights)), iterations=int(sol.nit))
