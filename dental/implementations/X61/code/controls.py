"""Independent quantitative controls, measured errors and deliberate corruptions."""
import copy, json, math
import numpy as np
from scipy.integrate import quad
from scipy.stats import norm, weibull_min
from scipy.special import logsumexp
from numpy.polynomial.hermite import hermgauss
from lab_alarm import *
from synthetic import series

def run():
    p = verified(R / 'inputs/DEMO_PROFILE.json')
    f = json.loads((R / 'inputs/X1B_FROZEN_PREDICTIONS.json').read_text())
    groups = model_groups(f)
    rows = series(6100, p, f)
    errs = []
    for m in (3.0, 6.0, 13.0):
        for x in (200.0, 1800.0, 3000.0, 6000.0):
            cdf = quad(lambda e: weibull_min.cdf(x * math.exp(-e), m, scale=3000) * norm.pdf(e, scale=0.01), -0.12, 0.12, epsabs=2e-13)[0]
            dens = quad(lambda e: weibull_min.pdf(x * math.exp(-e), m, scale=3000) * math.exp(-e) * norm.pdf(e, scale=0.01), -0.12, 0.12, epsabs=2e-13)[0]
            actual = math.exp(float(log_force_likelihood(x, np.array([m]), np.array([3000.0]), 0.01)[0])) / x
            errs.append(dict(shape=m, force_N=x, CDF_abs_error=abs(force_cdf(x, m, 3000, 0.01) - cdf), density_abs_error=abs(actual - dens), locator='SciPy weibull_min and scipy.integrate.quad, Gaussian log-force convolution'))
    post = JointPosterior(p, groups)
    for r in rows[:18]:
        post.update_force(r)
    z = post.initial.copy()
    (xx, ww) = hermgauss(32)
    ww /= math.sqrt(math.pi)
    for row in rows[:18]:
        key = (row['design_id'], int(row['load_angle_deg']))
        g = groups[key]
        fac = g['base_factor'] * np.where(post.form == 1, g['Q'] if g['interaction'] else 1.0, 1.0)
        fv = float(row['fracture_force_N'])
        true = fv * np.exp(-math.sqrt(2) * 0.01 * xx)
        dens = weibull_min.pdf(true[None, :], post.m[:, None], scale=(post.lam * fac)[:, None]) * true[None, :]
        z += np.log(np.maximum(dens @ ww, 1e-300))
    w = np.exp(z - logsumexp(z))
    err = float(np.max(abs(w - post.weights())))
    p2 = copy.deepcopy(p)
    p2['priors']['m_grid'][2] = 97
    p2['priors']['lambda_ratio_grid'][2] = 129
    a = JointPosterior(p, groups)
    b = JointPosterior(p2, groups)
    for r in rows:
        a.update_force(r)
        b.update_force(r)
    (qa, qb) = (a.report(), b.report())
    refine = {k: dict(coarse_mean=qa[k]['mean'], refined_mean=qb[k]['mean'], relative_difference=abs(qa[k]['mean'] / qb[k]['mean'] - 1)) for k in ('m', 'sigma0_force_equivalent')}
    out = dict(claim_type='capability', external_referent=dict(kind='published_code', locator='https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.weibull_min.html', compared_quantity='Weibull CDF/density and direct full density posterior', refutes_us=True), CDF_density_probes=errs, max_CDF_abs_error=max((e['CDF_abs_error'] for e in errs)), max_density_abs_error=max((e['density_abs_error'] for e in errs)), posterior_max_weight_abs_error=err, grid_refinement=refine, gates=dict(CDF_density=all((e['CDF_abs_error'] <= 1e-10 and e['density_abs_error'] <= 1e-10 for e in errs)), direct_full_posterior=err <= 1e-10), rigorous_numerical_enclosure='MISSING: finite probes measure quadrature/grid error, not a global bound. Alpha theorem assumes exact predictive CDF/categorical probabilities; no rigorous empirical or numerical alpha transfer claimed.')
    write(R / 'raw/CONTROL_PARITY.json', out)
    return out
if __name__ == '__main__':
    print(json.dumps(run()))
