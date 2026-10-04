"""Recompute GLS and heterogeneity from current observations, not saved PASS flags."""
import json, math
from pathlib import Path
import numpy as np
from statsmodels.regression.linear_model import GLS
P = Path(__file__).resolve().parent

def read(name):
    return json.loads((P / name).read_text())

def main():
    data = read('NORMALIZED_DATA.json')
    index = {r['row_id']: r for rows in data.values() for r in rows}
    checks = []
    for path in sorted(P.glob('MODEL_*.json')):
        m = json.loads(path.read_text())
        if 'beta' not in m:
            continue
        rows = [index[k] for k in m['row_ids']]
        terms = [c['term'] for c in m['coefficients']]
        X = np.array([[1.0] + [float(r[n.split('=')[0]] == n.split('=')[1]) if '=' in n else r[n] for n in terms[1:]] for r in rows])
        y = np.array([r['outcome'] for r in rows])
        variance = np.array([r['variance_no_n'] if m['no_n_division'] else r['variance'] for r in rows])
        S = np.diag(variance)
        for (i, a) in enumerate(rows):
            for (j, b) in enumerate(rows[:i]):
                if a['dataset'] == 'cement' and a['study'] == b['study'] and (a['arm'] == b['arm']):
                    S[i, j] = S[j, i] = m['rho_sampling'] * math.sqrt(variance[i] * variance[j])
        clusters = np.array([r['cluster'] for r in rows])
        B = (clusters[:, None] == clusters[None, :]).astype(float)
        V = S + m['tau2_study'] * B + m['omega2_cell'] * np.eye(len(rows))
        independent = GLS(y, X, sigma=V).fit()
        beta = np.array(m['beta'])
        beta_error = float(np.max(abs(independent.params - beta) / (1 + abs(independent.params))))
        cov_error = float(np.max(abs(independent.normalized_cov_params - np.array(m['coefficient_covariance']))) / (1 + np.max(abs(independent.normalized_cov_params))))
        coef_error = float(np.max(abs(beta - np.array([c['estimate'] for c in m['coefficients']]))))
        W = np.linalg.inv(S)
        Q = W - W @ X @ np.linalg.solve(X.T @ W @ X, X.T @ W)
        typical = (len(rows) - X.shape[1]) / np.trace(Q)
        total = m['tau2_study'] + m['omega2_cell']
        i2 = float(100 * total / (total + typical))
        wrong = beta.copy()
        wrong[0] += 1000
        fault_error = float(np.max(abs(independent.params - wrong) / (1 + abs(independent.params))))
        good = all((math.isfinite(x) for x in [beta_error, cov_error, coef_error, i2])) and max(beta_error, cov_error, coef_error) < 1e-05 and (abs(i2 - m['I2_total_percent']) < 1e-06)
        checks.append(dict(model=path.name, pass_gate=good, beta_error=beta_error, covariance_error=cov_error, coefficient_error=coef_error, recomputed_I2_percent=i2, injected_intercept_1000_rejected=fault_error > 1e-05))
    out = dict(status='PASS' if len(checks) == 6 and all((c['pass_gate'] and c['injected_intercept_1000_rejected'] for c in checks)) else 'FAIL', scope='Current-table GLS and generalized I2; variance-component optimization and clinical validity are separate.', checks=checks)
    (P / 'REVIEW_CONTROL_RESULT.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))
    if out['status'] != 'PASS':
        raise SystemExit(1)
if __name__ == '__main__':
    main()
