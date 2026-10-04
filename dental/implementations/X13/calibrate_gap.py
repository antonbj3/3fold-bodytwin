"""Regional study-intercept model and independent conventional GLS control."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '4'
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import solve_triangular
FEATURES = ['intercept', '(spacer-60)/40', 'zirconia', 'resin', 'PEEK_PEKK', 'printed', 'partial', 'cemented', 'CAD_exocad', 'CAD_3shape']
REGIONS = ('marginal', 'axial', 'occlusal')

def feature(d):
    cad = d['cad'].lower()
    return [1, (d['internal_spacer_um'] - 60) / 40, int(d['material'] == 'zirconia'), int(d['material'] == 'resin'), int(d['material'] == 'PEEK_PEKK'), int(d['process'] == 'printed'), int(d['restoration'] in ['endocrown', 'onlay']), int(d['state'] == 'cemented'), int('exocad' in cad), int('3shape' in cad)]

def arrays(rows):
    X = np.array([feature(d) for d in rows], float)
    y = np.array([d['measured_mean_um'] for d in rows])
    groups = [np.array([i for (i, d) in enumerate(rows) if d['study'] == s]) for s in sorted(set((d['study'] for d in rows)))]
    return (X, y, groups)

def fit(rows):
    (X, y, groups) = arrays(rows)
    P = np.diag([0] + [1 / 80 ** 2] * 9)

    def solve(logvar, return_details=False):
        (sigma, tau) = np.exp(logvar)
        s2 = sigma * sigma
        t2 = tau * tau
        Vi = np.zeros((len(y), len(y)))
        logdet = 0.0
        for ids in groups:
            n = len(ids)
            B = np.eye(n) / s2 - np.ones((n, n)) * t2 / (s2 * (s2 + n * t2))
            Vi[np.ix_(ids, ids)] = B
            logdet += n * np.log(s2) + np.log1p(n * t2 / s2)
        A = X.T @ Vi @ X + P
        beta = np.linalg.solve(A, X.T @ Vi @ y)
        res = y - X @ beta
        objective = logdet + res @ Vi @ res + beta @ P @ beta
        if return_details:
            return (beta, np.linalg.inv(A), sigma, tau, float(objective))
        return float(objective)
    opts = [minimize(solve, np.log(start), method='L-BFGS-B', bounds=[(np.log(2), np.log(1000))] * 2, options={'ftol': 1e-11, 'maxiter': 200}) for start in [(20, 20), (80, 80), (150, 50)]]
    best = min(opts, key=lambda x: x.fun)
    (beta, cov, sigma, tau, obj) = solve(best.x, True)
    return dict(beta=beta.tolist(), beta_cov=cov.tolist(), sigma_um=sigma, tau_um=tau, objective=obj, optimization_success=bool(best.success), starts=[dict(success=bool(o.success), objective=float(o.fun), message=str(o.message)) for o in opts], features=FEATURES, studies=sorted(set((d['study'] for d in rows))))

def predict(model, rows):
    X = np.array([feature(d) for d in rows], float)
    b = np.array(model['beta'])
    cov = np.array(model['beta_cov'])
    mu = X @ b
    sd = np.sqrt(model['sigma_um'] ** 2 + model['tau_um'] ** 2 + np.einsum('ij,jk,ik->i', X, cov, X))
    return (mu, sd)

def gls_control(rows, model, test):
    (X, y, groups) = arrays(rows)
    V = np.eye(len(y)) * model['sigma_um'] ** 2
    for ids in groups:
        V[np.ix_(ids, ids)] += model['tau_um'] ** 2
    L = np.linalg.cholesky(V)
    WX = solve_triangular(L, X, lower=True)
    Wy = solve_triangular(L, y, lower=True)
    prior = np.diag([0] + [1 / 80] * 9)
    beta = np.linalg.lstsq(np.vstack([WX, prior]), np.r_[Wy, np.zeros(10)], rcond=None)[0]
    return np.array([feature(d) for d in test]) @ beta

def simple_controls(rows, test):
    (X, y, groups) = arrays(rows)
    weights = np.array([1 / len(next((g for g in groups if i in g))) for i in range(len(y))])
    s = np.array([d['internal_spacer_um'] for d in rows])
    t = np.array([d['internal_spacer_um'] for d in test])
    mean = np.average(y, weights=weights)
    offset = np.average(y - s, weights=weights)
    A = np.stack([np.ones(len(y)), (s - 60) / 40], axis=1)
    b = np.linalg.lstsq(A * np.sqrt(weights[:, None]), y * np.sqrt(weights), rcond=None)[0]
    return {'mean': np.full(len(test), mean), 'offset': t + offset, 'affine': np.stack([np.ones(len(t)), (t - 60) / 40], axis=1) @ b}
