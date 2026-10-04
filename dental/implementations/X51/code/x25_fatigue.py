"""Small, explicit right-censored AFT implementation. No fitting package defaults."""
import numpy as np
from scipy.optimize import minimize, brentq
from scipy.special import log_ndtr, ndtr, ndtri
NSTAR = 5000000.0
FREF = 130.0

def arrays(rows, groups, pooled=False):
    g = np.zeros(len(rows), dtype=int) if pooled else np.array([groups.index(r['configuration']) for r in rows])
    x = np.log(np.array([r['max_load_N'] for r in rows]) / FREF)
    y = np.log(np.array([r['cycles'] for r in rows]) / NSTAR)
    d = np.array([r['failure'] for r in rows], dtype=float)
    return (g, x, y, d)

def objective(theta, arr, ng, family):
    (g, x, y, d) = arr
    (b, s) = np.exp(theta[-2:])
    mu = theta[g] - b * x
    z = (y - mu) / s
    if family == 'lognormal':
        logphi = -0.5 * z * z - 0.5 * np.log(2 * np.pi)
        logsf = log_ndtr(-z)
        h = np.exp(logphi - logsf)
        loss = d * (np.log(s) - logphi) - (1 - d) * logsf
        v = d * z + (1 - d) * h
        gm = -v / s
        gs = d * (1 - z * z) - (1 - d) * h * z
    else:
        ez = np.exp(np.clip(z, -700, 600))
        loss = d * (np.log(s) - z) + ez
        gm = (d - ez) / s
        gs = d + (d - ez) * z
    grad = np.r_[np.bincount(g, weights=gm, minlength=ng), np.sum(gm * (-b * x)), np.sum(gs)]
    return (float(np.sum(loss)), grad)

def fit(rows, groups, family='lognormal', pooled=False, initial=None, multistart=True):
    ng = 1 if pooled else len(groups)
    arr = arrays(rows, groups, pooled)
    (g, x, y, d) = arr
    X = np.c_[np.eye(ng)[g], -x]
    q = np.linalg.lstsq(X, y, rcond=None)[0]
    b0 = max(float(q[-1]), 0.5)
    a0 = [float(np.mean(y[g == k] + b0 * x[g == k])) for k in range(ng)]
    s0 = max(float(np.std(y - X @ q)), 0.3)
    base = np.r_[a0, np.log(b0), np.log(s0)]
    starts = [base] if initial is None else [np.array(initial)]
    if multistart:
        starts.extend([np.r_[a0, np.log(b0 * 0.5), np.log(s0 * 0.7)], np.r_[a0, np.log(b0 * 2), np.log(s0 * 1.5)]])
    bounds = [(-40.0, 40.0)] * ng + [(-5.0, 7.0), (-6.0, 5.0)]
    runs = []
    for start in starts:
        r = minimize(objective, start, args=(arr, ng, family), jac=True, method='L-BFGS-B', bounds=bounds, options={'maxiter': 1800, 'ftol': 1e-13, 'gtol': 1e-07, 'maxls': 60})
        (loss, grad) = objective(r.x, arr, ng, family)
        hit = any((min(t - lo, hi - t) < 0.0001 for (t, (lo, hi)) in zip(r.x, bounds)))
        runs.append({'theta': r.x.tolist(), 'nll_without_logtime': loss, 'gradient_max': float(np.max(np.abs(grad))), 'boundary_hit': bool(hit), 'success': bool(r.success), 'message': str(r.message)})
    best = min(runs, key=lambda r: r['nll_without_logtime'])
    spread = max((r['nll_without_logtime'] for r in runs)) - best['nll_without_logtime']
    valid = np.isfinite(best['nll_without_logtime']) and best['gradient_max'] < 0.001 and (not best['boundary_hit']) and (spread < 0.0001)
    return {'family': family, 'pooled': pooled, 'groups': groups, 'theta': best['theta'], 'nll': best['nll_without_logtime'] + float(np.sum(d * (y + np.log(NSTAR)))), 'gradient_max': best['gradient_max'], 'boundary_hit': best['boundary_hit'], 'start_nll_spread': float(spread), 'valid': bool(valid), 'starts': runs}

def naive_fit(rows, groups):
    """The explicit practice ablation: runout timestamp encoded as an event."""
    (g, x, y, d) = arrays(rows, groups)
    X = np.c_[np.eye(len(groups))[g], -x]
    q = np.linalg.lstsq(X, y, rcond=None)[0]
    s = np.sqrt(np.mean((y - X @ q) ** 2))
    return {'family': 'lognormal', 'pooled': False, 'groups': groups, 'theta': np.r_[q[:-1], np.log(max(q[-1], 1e-08)), np.log(max(s, 1e-08))].tolist(), 'raw_slope': float(q[-1]), 'valid': bool(q[-1] > 0 and s > 0), 'information_loss': 'right censoring timestamp as exact failure'}

def location(model, force, group):
    th = np.array(model['theta'])
    (b, s) = np.exp(th[-2:])
    k = 0 if model['pooled'] else model['groups'].index(group)
    return (th[k] - b * np.log(np.asarray(force) / FREF), s)

def risk(model, force, group, cycles=NSTAR):
    (mu, s) = location(model, force, group)
    z = (np.log(cycles / NSTAR) - mu) / s
    if model['family'] == 'lognormal':
        return ndtr(z)
    return -np.expm1(-np.exp(np.clip(z, -700, 600)))

def force_at_risk(model, p, group, cycles=NSTAR):
    th = np.array(model['theta'])
    (b, s) = np.exp(th[-2:])
    k = 0 if model['pooled'] else model['groups'].index(group)
    z = ndtri(p) if model['family'] == 'lognormal' else np.log(-np.log1p(-p))
    return float(FREF * np.exp((th[k] + s * z - np.log(cycles / NSTAR)) / b))

def observation_logscore(model, row):
    (mu, s) = location(model, row['max_load_N'], row['configuration'])
    z = (np.log(row['cycles'] / NSTAR) - mu) / s
    if model['family'] == 'lognormal':
        return float(-np.log(row['cycles']) - np.log(s) - 0.5 * np.log(2 * np.pi) - 0.5 * z * z if row['failure'] else log_ndtr(-z))
    ez = np.exp(np.clip(z, -700, 600))
    return float(-np.log(row['cycles']) - np.log(s) + z - ez if row['failure'] else -ez)

def endpoint_known(row, cycles=NSTAR):
    return bool(row['cycles'] >= cycles or row['failure'])

def endpoint(row, cycles=NSTAR):
    if not endpoint_known(row, cycles):
        return None
    return int(row['failure'] == 1 and row['cycles'] <= cycles)
