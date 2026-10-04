"""Fixed-plan inference. Exact and approximate outputs use different contracts."""
import thread_budget
import argparse, math, json
import numpy as np
from scipy.stats import beta, binom, chi2, gamma, weibull_min
from scipy.special import gammaln, logsumexp
from scipy.optimize import brentq, minimize_scalar

def upper(n, k, alpha=0.05):
    if not isinstance(n, int) or not isinstance(k, int) or n < 1 or (k < 0) or (k > n) or (not 0 < alpha < 1):
        raise ValueError('Require integer0<=k<=n, n>0 and0<alpha<1')
    return 1.0 if k == n else float(beta.isf(alpha, k + 1, n - k))

def n_zero(p_bad=0.1, alpha=0.05, groups=1):
    if not 0 < p_bad < 1 or not 0 < alpha < 1 or (not isinstance(groups, int)) or (groups < 1):
        raise ValueError('Invalid risk, alpha or number of groups')
    n = math.ceil(math.log(alpha / groups) / math.log1p(-p_bad))
    while upper(n, 0, alpha / groups) > p_bad:
        n += 1
    return {'n': n, 'failures_allowed': 0, 'upper': upper(n, 0, alpha / groups), 'groups': groups, 'scope': 'conditional minimum if all n survive, not expected test count or80% power', 'resolution': 'POPULATION within frozen configuration/regime'}

def plan(p_bad=0.1, p_good=0.01, power=0.8, alpha=0.05, max_n=10000):
    if not 0 < p_good < p_bad < 1 or not 0 < power < 1 or (not 0 < alpha < 1):
        raise ValueError('Require0<p_good<p_bad<1 and valid alpha/power')
    for n in range(1, max_n + 1):
        c = int(binom.ppf(alpha, n, p_bad))
        if binom.cdf(c, n, p_bad) > alpha:
            c -= 1
        if c >= 0 and binom.cdf(c, n, p_good) >= power:
            return {'n': n, 'failures_allowed': c, 'false_acceptance_at_p_bad': float(binom.cdf(c, n, p_bad)), 'power_at_p_good': float(binom.cdf(c, n, p_good)), 'upper_at_acceptance_boundary': upper(n, c, alpha), 'p_bad': p_bad, 'p_good': p_good, 'alpha': alpha, 'power_required': power, 'fixed_plan': True, 'resolution': 'POPULATION within frozen specimen/process regime'}
    raise ValueError('No feasible plan below max_n')

def endpoint_rows(rows, threshold):
    """Failures <= threshold known bad; survival past threshold known good; retain unknowns."""
    bad = good = unknown = 0
    for r in rows:
        x = float(r['value'])
        event = r['failure']
        if not np.isfinite(x) or x <= 0 or event not in (0, 1):
            raise ValueError('Invalid strength/time/event')
        if event and x <= threshold:
            bad += 1
        elif x >= threshold:
            good += 1
        else:
            unknown += 1
    n = bad + good + unknown
    return {'n': n, 'bad': bad, 'good': good, 'unknown': unknown, 'sample_risk_interval': [bad / n, (bad + unknown) / n], 'upper95_conservative_all_unknown_bad': upper(n, bad + unknown), 'resolution': 'POPULATION', 'censor_independence_required': False}

def known_shape_lower(values, m, alpha=0.05, censoring='COMPLETE', r=None, shape_independently_known=False):
    if not shape_independently_known:
        raise ValueError('Estimated shape cannot enter an exact known-shape pivot')
    if not np.isfinite(m) or m <= 0 or (not 0 < alpha < 1):
        raise ValueError('Invalid shape/alpha')
    v = np.asarray(values, float)
    if v.ndim != 1 or not len(v) or np.any(~np.isfinite(v)) or np.any(v <= 0):
        raise ValueError('Positive finite sample required')
    if censoring == 'COMPLETE':
        nr = len(v)
        logT = logsumexp(m * np.log(v))
    elif censoring == 'TYPE_II':
        if r is None or not isinstance(r, int) or r < 1 or (r > len(v)):
            raise ValueError('Fixed TypeII r required')
        if np.any(np.diff(v) < 0) or np.any(v[r:] != v[r - 1]):
            raise ValueError('TypeII records: sorted r failures and every survivor stopped at rth failure')
        nr = r
        logT = logsumexp(m * np.log(v))
    else:
        raise ValueError('Exact pivot available only for COMPLETE or fixed TYPE_II censoring')
    eta = float(np.exp((math.log(2) + logT - math.log(chi2.isf(alpha, 2 * nr))) / m))
    return {'eta_lower': eta, 'mean_lower': eta * math.gamma(1 + 1 / m), 'q05_lower': eta * (-math.log(0.95)) ** (1 / m), 'independently_known_m': m, 'failure_count_pivot': nr, 'method': 'exact conditional on known m, iid and declared censoring scheme'}

def known_shape_plan(m, mean_ratio=1.2, power=0.8, alpha=0.05):
    if m <= 0 or mean_ratio <= 1 or (not 0 < power < 1) or (not 0 < alpha < 1):
        raise ValueError('Invalid planning assumptions')
    for n in range(1, 10001):
        p = float(gamma.sf(chi2.isf(alpha, 2 * n) / (2 * mean_ratio ** m), a=n))
        if p >= power:
            return {'m': m, 'true_mean_over_threshold': mean_ratio, 'n_statistical': n, 'n_with_ISO_bar_minimum': max(n, 10), 'power_statistical': p, 'power_at_ISO_bar_n': float(gamma.sf(chi2.isf(alpha, 2 * max(n, 10)) / (2 * mean_ratio ** m), a=max(n, 10))), 'closure': 'known shape from independent applicable evidence; COMPLETE iid bars', 'resolution': 'POPULATION'}
    return {'status': 'UNKNOWN_TOO_LARGE'}

def nll(logm, logeta, x, d):
    m = np.exp(logm)
    z = m * (np.log(x) - logeta)
    if np.max(z) > 700:
        return float('inf')
    return float(np.sum(np.exp(z)) - np.sum(d * (logm + (m - 1) * np.log(x) - m * logeta)))

def weibull_fit(x, event=None):
    x = np.asarray(x, float)
    d = np.ones(len(x)) if event is None else np.asarray(event, int)
    if len(x) == 0 or np.any(~np.isfinite(x)) or np.any(x <= 0) or (len(d) != len(x)) or np.any((d != 0) & (d != 1)):
        raise ValueError('Invalid Weibull data')
    r = int(d.sum())
    if r == 0:
        return {'status': 'UNKNOWN_NO_EVENTS', 'inference': 'profile MLE not identified'}

    def profile(a):
        m = np.exp(a)
        eta = (logsumexp(m * np.log(x)) - math.log(r)) / m
        return nll(a, eta, x, d)
    opt = minimize_scalar(profile, bounds=(math.log(0.15), math.log(100)), method='bounded', options={'xatol': 1e-12})
    m = float(np.exp(opt.x))
    logeta = float((logsumexp(m * np.log(x)) - math.log(r)) / m)
    if not opt.success or abs(opt.x - math.log(0.15)) < 1e-05 or abs(opt.x - math.log(100)) < 1e-05:
        return {'status': 'UNKNOWN_BOUNDARY_FIT'}
    return {'status': 'FITTED', 'm': m, 'eta': float(np.exp(logeta)), 'nll': float(opt.fun), 'n': len(x), 'failures': r, 'mean': float(np.exp(logeta) * math.gamma(1 + 1 / m)), 'q05': float(np.exp(logeta) * (-math.log(0.95)) ** (1 / m)), 'inference': 'model fit; no exact finite-sample confidence for estimated shape'}

def profile_lower(x, event=None, quantity='mean', alpha=0.05):
    x = np.asarray(x, float)
    d = np.ones(len(x)) if event is None else np.asarray(event, int)
    f = weibull_fit(x, d)
    if f['status'] != 'FITTED':
        return f
    c = lambda m: math.gamma(1 + 1 / m) if quantity == 'mean' else (-math.log(0.95)) ** (1 / m)

    def fixed(q):

        def obj(a):
            m = np.exp(a)
            return nll(a, math.log(q) - math.log(c(m)), x, d)
        opt = minimize_scalar(obj, bounds=(math.log(0.15), math.log(100)), method='bounded', options={'xatol': 1e-10})
        return opt.fun
    crit = float(chi2.isf(2 * alpha, 1))
    qhat = f[quantity]
    root = brentq(lambda lq: 2 * (fixed(math.exp(lq)) - f['nll']) - crit, math.log(qhat) - 12, math.log(qhat), xtol=1e-11)
    return {'estimate': qhat, 'lower': math.exp(root), 'm': f['m'], 'eta': f['eta'], 'quantity': quantity, 'confidence_nominal': 1 - alpha, 'method': 'asymptotic one-sided profile likelihood, chi2_1(1-2alpha)', 'finite_sample_guarantee': 'UNKNOWN', 'eligible_for_exact_lot_release': False, 'resolution': 'POPULATION'}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('plan')
    s.add_argument('--p-bad', type=float, default=0.1)
    s.add_argument('--p-good', type=float, default=0.01)
    s.add_argument('--power', type=float, default=0.8)
    s.add_argument('--alpha', type=float, default=0.05)
    s = sub.add_parser('upper')
    s.add_argument('--n', type=int, required=True)
    s.add_argument('--failures', type=int, required=True)
    s.add_argument('--alpha', type=float, default=0.05)
    s = sub.add_parser('weibull')
    s.add_argument('--input', required=True)
    s.add_argument('--quantity', choices=['mean', 'q05'], default='mean')
    a = p.parse_args()
    if a.cmd == 'plan':
        o = plan(a.p_bad, a.p_good, a.power, a.alpha)
    elif a.cmd == 'upper':
        o = {'upper': upper(a.n, a.failures, a.alpha), 'resolution': 'POPULATION', 'fixed_n_required': True}
    else:
        rows = json.load(open(a.input))
        o = profile_lower([r['value'] for r in rows], [r['failure'] for r in rows], a.quantity)
    print(json.dumps(o, indent=2, allow_nan=False))
if __name__ == '__main__':
    main()
