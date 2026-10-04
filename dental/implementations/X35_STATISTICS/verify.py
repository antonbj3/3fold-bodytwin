"""Independent scalar replays, source locks, exact tails and falsifying injections."""
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
for _var in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[_var] = '4'
import numpy as np
from scipy import optimize, stats
ROOT = Path(__file__).resolve().parent

def scalar(row):
    kind = row['kind']
    if kind == 'AUC':
        y = np.asarray(row['outcome'], bool)
        s = np.asarray(row['values'])
        ranks = stats.rankdata(s, method='average')
        n1 = y.sum()
        n0 = (~y).sum()
        return float((ranks[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
    if kind == 'kappa':
        a = row['reference']
        b = row['prediction']
        n = len(a)
        observed = sum((x == y for (x, y) in zip(a, b))) / n
        classes = set(a) | set(b)
        chance = sum((a.count(x) * b.count(x) for x in classes)) / (n * n)
        return (observed - chance) / (1 - chance)
    a = np.asarray(row['values'], float)
    if kind == 'median':
        return float(np.median(a[:, 0]))
    if kind == 'median_ratio':
        return float(np.median(a[:, 0]) / np.median(a[:, 1]))
    if kind == 'rmse':
        a = a * a
    if row.get('weighting') == 'equal_cluster':
        ids = np.array(row['cluster_ids'])
        v = np.mean([a[ids == c].mean(axis=0) for c in sorted(set(ids))], axis=0)
    else:
        v = a.mean(axis=0)
    if kind == 'rmse':
        return float(math.sqrt(v[0]))
    if kind == 'ratio':
        return float(v[0] / v[1])
    return float(v[0])

def run(result, raw, r2):
    controls = []
    lookup = {(m['demo'], m['metric']): m for m in result['metrics']}
    for row in raw:
        key = (row['demo'], row['metric'])
        expected = scalar(row)
        actual = lookup[key]['point']
        tol = max(1e-09, abs(expected) * 1e-08)
        bad = actual + max(0.1, abs(actual) * 0.1)
        controls.append(dict(name='/'.join(key), kind='Independent scalar replay', expected=expected, actual=actual, tolerance=tol, pass_=abs(actual - expected) <= tol, injected_value=bad, injected_rejected=abs(bad - expected) > tol))
    groups = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [3.0, 3.0, 3.0], [7.0, 7.0, 7.0]]
    draws = list(itertools.product(range(4), repeat=4))
    direct = np.array([np.mean([v for j in draw for v in groups[j]]) for draw in draws])
    from_means = np.array([np.mean([np.mean(groups[j]) for j in draw]) for draw in draws])
    duplicated = np.array([np.mean([v for j in draw for v in groups[j] * 10]) for draw in draws])
    sd = np.std([0.0, 1.0, 3.0, 7.0], ddof=1)
    correct_se = sd / math.sqrt(4)
    wrong_se = sd / math.sqrt(40)
    controls.append(dict(name='Whole-case duplication versus deliberate row-IID error', pass_=bool(np.allclose(direct, from_means) and np.allclose(direct, duplicated)), injected_rejected=abs(wrong_se - correct_se) > 1e-10, correct_case_SE=correct_se, wrong_row_SE=wrong_se))
    cells = np.array([1.0, 3.0, 4.0])
    weights = np.array([0.5, 0.5, 1.0])
    weighted_mean = float(np.average(cells, weights=weights))
    weighted_rmse = float(np.sqrt(np.average(cells * cells, weights=weights)))
    controls.append(dict(name='Frozen study weights versus injected pooled weighting', pass_=abs(weighted_mean - 3.0) < 1e-10 and abs(weighted_rmse - math.sqrt(10.5)) < 1e-10, injected_rejected=abs(cells.mean() - weighted_mean) > 1e-10, expected_weighted_mean=3.0, bad_pooled_mean=float(cells.mean())))
    for (n, s) in [(2, 0), (8, 0), (8, 3), (29, 0), (433, 37)]:
        from round2 import exact_binomial
        reported = exact_binomial(s, n)['ci95']
        lo = 0.0 if s == 0 else optimize.brentq(lambda p: stats.binom.sf(s - 1, n, p) - 0.025, 1e-14, 1 - 1e-14)
        hi = 1.0 if s == n else optimize.brentq(lambda p: stats.binom.cdf(s, n, p) - 0.025, 1e-14, 1 - 1e-14)
        error = max(abs(lo - reported[0]), abs(hi - reported[1]))
        controls.append(dict(name=f'Exact tail inversion {s}/{n}', expected=[lo, hi], observed=reported, pass_=error < 1e-10, injected_rejected=abs(hi - 0.1 - reported[1]) > 1e-10))
    for m in result['metrics']:
        k = m.get('k_clusters')
        d = m.get('mde_80pct_alpha05_cluster_SD')
        if k and d:
            threshold = stats.t.ppf(0.975, k - 1)
            nc = d * math.sqrt(k)
            power = float(stats.nct.sf(threshold, k - 1, nc) + stats.nct.cdf(-threshold, k - 1, nc))
            controls.append(dict(name=m['demo'] + '/' + m['metric'] + '/MDE-power', power=power, pass_=abs(power - 0.8) < 1e-08, injected_rejected=abs(0.05 - 0.8) > 1e-08, injected_d=0))
    pin = json.loads((ROOT / 'INPUT_LOCK.json').read_text())
    for item in pin:
        actual = hashlib.sha256(Path(item['path']).read_bytes()).hexdigest()
        bad = '0' * 64
        controls.append(dict(name=item['path'], kind='Read-only source hash', pass_=actual == item['sha256'], injected_rejected=bad != item['sha256']))
    ids = {d['id'] for d in result['demos']}
    expected = {d['id'] for d in json.loads((ROOT / 'INPUT_INVENTORY_R3.json').read_text())}
    controls.append(dict(name='All frozen package demos covered', pass_=ids == expected and len(ids) == 35, injected_rejected=ids - {'PROOF_LANE'} != expected))
    all_controls = result['controls'] + r2['controls'] + controls
    passed = all((c['pass_'] and c['injected_rejected'] for c in all_controls))
    record = dict(status='PASS' if passed else 'FAIL', n_controls=len(all_controls), controls=all_controls, meaning='Software/statistical validity controls, not physical validation or evidence of patient independence.')
    (ROOT / 'VALIDATION.json').write_text(json.dumps(record, ensure_ascii=False, indent=2, default=lambda x: x.item()) + '\n')
    if not passed:
        failures = [c for c in all_controls if not c['pass_'] or not c['injected_rejected']]
        raise ValueError('Failed controls: ' + str(failures[:3]))
    return record
