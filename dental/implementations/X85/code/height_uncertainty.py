from common import *
from regional import fit_response
import copy

def fit_bounded_heights(probes, final_coefficient_errors_mm=None):
    s = float(probes['step_mm'])
    d = float(probes['height_error_mm'])
    if not np.isfinite(d) or d < 0 or (not np.isfinite(s)) or (d >= s):
        raise ValueError('Need0<=height_error<step')
    p = copy.deepcopy(probes)
    p['height_error_mm'] = 0.0
    out = fit_response(p)
    if out['status'] != 'CONDITIONAL_CALIBRATION':
        return out
    J = np.array(out['J_N_per_mm'])
    norm = np.linalg.norm(J, axis=0)
    eta = out['eta_N']
    upper = (norm + eta / s) / (1 - d / s)
    beta = (eta + d * norm) / (s - d)
    u = np.zeros(J.shape[1]) if final_coefficient_errors_mm is None else np.asarray(final_coefficient_errors_mm, float)
    if u.shape != (J.shape[1],) or not np.isfinite(u).all() or np.min(u) < 0:
        raise ValueError('Invalid final coefficient errors')
    out.update(beta_N_per_mm=beta.tolist(), true_column_norm_upper_N_per_mm=upper.tolist(), probe_height_error_mm=d, final_coefficient_error_mm=u.tolist(), raw_baseline_error_N=eta, eta_N=float(eta + upper @ u), enclosure='Exact denominator norm bound plus finite actuation coefficients; linear full-contact closure required; roundoff MISSING')
    return out
