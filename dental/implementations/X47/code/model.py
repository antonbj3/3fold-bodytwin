"""Small marginal REML model with explicit rank/transport refusal.
No fitted linear sensitivity is claimed as a rigorous physical enclosure.
"""
import json, math, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.linalg import cho_factor, cho_solve, block_diag
ROOT = Path(__file__).resolve().parents[1]

def design(rows, kind):
    if kind == 'full':
        return np.array([[1, math.log(r['height_mm'] / 4), r['total_convergence_deg'] / 10, r['permanent_resin_cement'], math.log1p(r['thermal_cycles'] / 1000), r['substrate'] == 'titanium'] for r in rows], float)
    if kind == 'reduced':
        return np.array([[1, math.log(r['height_mm'] / 4), r['permanent_resin_cement']] for r in rows], float)
    raise ValueError(kind)

def covariance(rows, tau2):
    v = np.array([(r['sd_force_N'] / r['mean_force_N']) ** 2 / r['n'] for r in rows])
    studies = np.array([r['study'] for r in rows])
    return np.diag(v) + tau2 * (studies[:, None] == studies[None, :])

def solve_gls(X, y, V):
    f = cho_factor(V, lower=True)
    W = cho_solve(f, np.eye(len(y)))
    normal = X.T @ W @ X
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise ValueError('RANK_DEFICIENT')
    cov = np.linalg.inv(normal)
    beta = cov @ (X.T @ W @ y)
    return (beta, cov, f)

def fit(rows, kind):
    X = design(rows, kind)
    y = np.log([r['mean_force_N'] for r in rows])
    p = X.shape[1]
    if np.linalg.matrix_rank(X) < p:
        return {'status': 'UNKNOWN_RANK_DEFICIENT', 'rank': int(np.linalg.matrix_rank(X)), 'columns': p}

    def objective(tau2):
        V = covariance(rows, tau2)
        (b, c, f) = solve_gls(X, y, V)
        r = y - X @ b
        return float(2 * np.log(np.diag(f[0])).sum() - np.linalg.slogdet(c)[1] + r @ cho_solve(f, r) + (len(y) - p) * math.log(2 * math.pi))
    opt = minimize_scalar(objective, bounds=(0, 100), method='bounded', options={'xatol': 1e-10})
    candidates = [(0.0, objective(0)), (float(opt.x), float(opt.fun)), (100.0, objective(100))]
    tau2 = min(candidates, key=lambda t: t[1])[0]
    (b, c, f) = solve_gls(X, y, covariance(rows, tau2))
    V = covariance(rows, tau2)
    VX = np.linalg.solve(V, X)
    Vy = np.linalg.solve(V, y)
    bc = np.linalg.solve(X.T @ VX, X.T @ Vy)
    return {'status': 'FITTED', 'beta': b.tolist(), 'beta_cov': c.tolist(), 'tau2': tau2, 'optimizer_success': bool(opt.success), 'tau2_search_boundary': tau2 >= 99.99, 'gls_control_max_abs': float(np.max(np.abs(b - bc))), 'residual_log_rms': float(np.sqrt(np.mean((y - X @ b) ** 2)))}

def loso(rows, kind):
    out = []
    for study in sorted({r['study'] for r in rows}):
        tr = [r for r in rows if r['study'] != study]
        te = [r for r in rows if r['study'] == study]
        m = fit(tr, kind)
        entry = {'held_out_study': study, 'training_studies': len({r['study'] for r in tr}), 'model': m, 'predictions': []}
        if m['status'] == 'FITTED':
            for (r, x) in zip(te, design(te, kind)):
                mu = float(x @ np.array(m['beta']))
                var = float(x @ np.array(m['beta_cov']) @ x) + m['tau2'] + (r['sd_force_N'] / r['mean_force_N']) ** 2 / r['n']
                h = 1.96 * math.sqrt(var)
                entry['predictions'].append({'row_id': r['row_id'], 'observed_N': r['mean_force_N'], 'predicted_group_mean_N': math.exp(mu), '95pct_lower_N': math.exp(mu - h), '95pct_upper_N': math.exp(mu + h), 'log_error': math.log(r['mean_force_N']) - mu, 'covered': mu - h <= math.log(r['mean_force_N']) <= mu + h, 'interval_scope': 'Approximate new-study group-mean prediction interval; NOT an individual-crown interval.'})
        out.append(entry)
    predictions = [p for f in out for p in f['predictions']]
    rmse = math.sqrt(np.mean([r['log_error'] ** 2 for r in predictions])) if predictions else None
    coverage = float(np.mean([r['covered'] for r in predictions])) if predictions else None
    return {'folds': out, 'evaluated_groups': len(predictions), 'all_folds_fitted': all((x['model']['status'] == 'FITTED' for x in out)), 'log_rmse': rmse, '95pct_coverage': coverage}

def run_full(rows):
    start = time.monotonic()
    elig = [r for r in rows if r['independent_first_state'] and all((r[k] is not None for k in ['height_mm', 'total_convergence_deg', 'thermal_cycles']))]
    model = fit(elig, 'full')
    loo = loso(elig, 'full')
    gate = 'UNKNOWN_NOT_IDENTIFIABLE' if len({r['study'] for r in elig}) < 4 or model['status'] != 'FITTED' or (not loo['all_folds_fitted']) else 'PASS' if loo['log_rmse'] <= math.log(1.5) and loo['95pct_coverage'] >= 0.8 else 'FAIL_TRANSFER'
    out = {'round': 'R1', 'claim_type': 'information_link', 'eligible_groups': len(elig), 'total_extracted_groups': len(rows), 'eligible_studies': len({r['study'] for r in elig}), 'excluded_groups': len(rows) - len(elig), 'model': model, 'loso': loo, 'gate': gate, 'wall_s': time.monotonic() - start, 'external_referent': {'kind': 'independent_measurement', 'locator': 'raw/RETENTION_GROUPS.json: DOI + table_id + XML row/cell per group', 'compared_quantity': 'whole-crown pull-off mean force_N', 'refutes_us': True}}
    (ROOT / 'rounds/RESULTS_R1.json').write_text(json.dumps(out, indent=2) + '\n')
    return out
if __name__ == '__main__':
    rows = json.loads((ROOT / 'raw/RETENTION_GROUPS.json').read_text())
    out = run_full(rows)
    print({k: out[k] for k in ['eligible_groups', 'eligible_studies', 'model', 'gate']})
