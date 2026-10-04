"""External published-code GLS check on actual retained observations and fitted V."""
from pathlib import Path
import json, numpy as np, statsmodels
from statsmodels.regression.linear_model import GLS
P = Path(__file__).resolve().parent
data = json.loads((P / 'NORMALIZED_DATA.json').read_text())
index = {r['row_id']: r for v in data.values() for r in v}
checks = []
for path in sorted(P.glob('MODEL_*.json')):
    m = json.loads(path.read_text())
    if 'beta' not in m:
        continue
    rows = [index[r] for r in m['row_ids']]
    names = [r['term'] for r in m['coefficients']]
    X = np.array([[1.0] + [r[n] if '=' not in n else float(r[n.split('=')[0]] == n.split('=')[1]) for n in names[1:]] for r in rows])
    y = np.array([r['outcome'] for r in rows])
    v = np.array([r['variance'] for r in rows])
    cl = np.array([r['cluster'] for r in rows])
    B = (cl[:, None] == cl[None, :]).astype(float)
    V = np.diag(v) + m['tau2_study'] * B + m['omega2_cell'] * np.eye(len(rows))
    result = GLS(y, X, sigma=V).fit()
    C = result.normalized_cov_params
    beta = np.array(m['beta'])
    errors = np.max(abs(result.params - beta) / (1 + abs(beta)))
    cov_error = np.max(abs(C - np.array(m['coefficient_covariance']))) / (1 + np.max(abs(C)))
    wrong_beta = beta.copy()
    wrong_beta[0] += 1
    rejected = np.max(abs(result.params - wrong_beta) / (1 + abs(beta))) > 1e-05
    checks.append(dict(model=path.name, relative_beta_error=float(errors), relative_known_covariance_error=float(cov_error), pass_gate=bool(max(errors, cov_error) < 1e-05), wrong_intercept_rejected=bool(rejected)))
out = dict(external_referent=dict(kind='published_code', locator='https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.GLS.html', compared_quantity='GLS beta and known-variance coefficient covariance at fitted REML V; does not independently validate variance-component optimization', refutes_us=False), statsmodels_version=statsmodels.__version__, checks=checks, all_pass=all((r['pass_gate'] and r['wrong_intercept_rejected'] for r in checks)))
(P / 'STATSMODELS_CONTROL.json').write_text(json.dumps(out, indent=2) + '\n')
print('External statsmodels GLS', len(checks), 'models', out['all_pass'])
assert out['all_pass']
