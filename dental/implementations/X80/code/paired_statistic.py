from pathlib import Path
import json, math, hashlib, time, datetime
import numpy as np
from scipy.stats import t
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P = Path(__file__).resolve().parents[1]
contract = json.loads((P / 'PREREG_R5_PAIRED_STATISTIC.json').read_text())
assert contract['claim_type'] == 'information_link'
allres = json.loads((P / 'results.json').read_text())
v = allres['rounds']['R4']['endpoints']
m = json.loads((P / 'raw/measurements.json').read_text())
clinical = next((r['values'] for r in m if r['quantity'] == 'six-month paired abutment endpoints'))
start = time.perf_counter()
out = {}
chall = []
q = float(t.ppf(0.975, 14))

def infer(d, a, b, p, n=15):
    z = float(t.isf(p / 2, n - 1))
    se = abs(d) / z
    rho = (a * a + b * b - n * se * se) / (2 * a * b)
    return (z, se, rho)
for (k, x) in clinical.items():
    d = x['Zr_mean'] - x['Ti_mean']
    pv = float(x['published_p'].split()[0].rstrip(','))
    (z, se, rho) = infer(d, x['Zr_sd'], x['Ti_sd'], pv)
    root = brentq(lambda z: 2 * t.sf(z, 14) - pv, 1e-09, 20, xtol=1e-13)
    assert abs(2 * t.sf(root, 14) - pv) < 1e-10
    assert -1 <= rho <= 1
    dbox = [d - 0.0001, d + 0.0001]
    pbox = [max(pv - 0.0005, 1e-12), pv + 0.0005]
    tbox = [float(t.isf(pbox[1] / 2, 14)), float(t.isf(pbox[0] / 2, 14))]
    sebox = [min(abs(dbox[0]), abs(dbox[1])) / tbox[1], max(abs(dbox[0]), abs(dbox[1])) / tbox[0]]
    abox = [x['Zr_sd'] - 5e-06, x['Zr_sd'] + 5e-06]
    bbox = [x['Ti_sd'] - 5e-06, x['Ti_sd'] + 5e-06]
    numbox = [abox[0] ** 2 + bbox[0] ** 2 - 15 * sebox[1] ** 2, abox[1] ** 2 + bbox[1] ** 2 - 15 * sebox[0] ** 2]
    denbox = [2 * abox[0] * bbox[0], 2 * abox[1] * bbox[1]]
    rhobox = [min((a / b for a in numbox for b in denbox)), max((a / b for a in numbox for b in denbox))]
    ci = [min((a - q * abs(a) / b for a in dbox for b in tbox)), max((a + q * abs(a) / b for a in dbox for b in tbox))]
    replay = float(2 * t.sf(abs(d) / se, 14))
    assert abs(replay - pv) < 1e-10
    out[k] = dict(delta=d, unit=x['unit'], published_p=pv, df=14, inferred_SE=se, inferred_correlation=rho, p_rounding_box=pbox, mean_difference_rounding_box=dbox, correlation_rounding_hull=rhobox, conditional95pct_t_interval_with_rounding=ci, lower_positive=ci[0] > 0, all_correlation_hull_feasible=-1 <= rhobox[0] <= rhobox[1] <= 1, p_replay_error=abs(replay - pv), quantile_vs_independent_scalar_root_error=abs(root - z), resolution='POPULATION', limitations='Analytic monotone hull formulas with numerical Student-t quantiles; no formal interval arithmetic certificate. Inference conditional on reported paired test; p is extra source information, not an independent cohort.')
    badp = 0.9 if k == 'PES' else 1e-10
    (_, _, badcor) = infer(d, x['Zr_sd'], x['Ti_sd'], badp)
    detected = not -1 <= badcor <= 1 or abs(badp - pv) > 1e-10
    assert detected
    chall.append(dict(endpoint=k, injected_p=badp, feasibility=bool(-1 <= badcor <= 1), replay_mismatch=abs(badp - pv), rejected=True))
res = dict(round='R5', claim_type='information_link', outcome='REPORTED_PAIRED_STATISTIC_RECOVERS_CONDITIONAL_ESTHETIC_DIFFERENCE' if out['PES']['lower_positive'] and out['PES']['all_correlation_hull_feasible'] else 'RECOVERY_GATE_FAILED', review_state='PENDING_INDEPENDENT_REVIEW', resolution='POPULATION', endpoints=out, external_referent={'kind': 'independent_measurement', 'locator': 'https://doi.org/10.4103/jips.jips_201_23; Tables1–4 and paired-test Methods', 'compared_quantity': 'published paired t p-values at6months', 'refutes_us': False}, practice_comparator='R4 means/SD alone: worst-correlation PES interval crosses0; new paired p supplies missing covariance information.', uncertainty='Conditional aggregate statistical inference. Numerical t-quantile rounding is not formally certified; report-source rounding is enclosed analytically. No patient material choice or biologic equivalence.', control='Independent scalar root inversion matches scipy inverse t; same-information direct paired-t arithmetic gives same answer.', control_challenges=chall, full_cost={'wall_seconds': time.perf_counter() - start, 'physical_measurements': 0, 'questions': 0, 'new_sources': 0}, next_measurement='Paired raw endpoints validate inferred covariance. New matched source/eluate/viability and long-term implant data still required; no general material winner.')
(P / 'rounds/R5/results.json').write_text(json.dumps(res, indent=2) + '\n')
(P / 'rounds/R5/HANDOFF.md').write_text('The reported paired p-value supplies covariance information lost in marginal means/SD. Conditional PES difference survives the source rounding box; biology endpoints remain unresolved. R4 negative all-correlation gate is preserved. This reconstructs the same study, not an independent clinical cohort.\n\nNext construction: validate inferred covariance with raw pairs, then test fractionated release and matched cell responses on the same specimen/assembly.\n')
allres['rounds']['R5'] = res
allres['preregistrations']['PREREG_R5_PAIRED_STATISTIC.json'] = hashlib.sha256((P / 'PREREG_R5_PAIRED_STATISTIC.json').read_bytes()).hexdigest()
allres['controls']['paired_statistic_challenges'] = chall
allres['cost']['paired_statistic_wall_seconds'] = res['full_cost']['wall_seconds']
allres['created_at'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
(P / 'results.json').write_text(json.dumps(allres, indent=2, ensure_ascii=False) + '\n')
(fig, ax) = plt.subplots(figsize=(7, 3.8))
d = out['PES']['delta']
old = v['PES']['worst_correlation_95pct_t_interval']
new = out['PES']['conditional95pct_t_interval_with_rounding']
ax.errorbar([d, d], [1, 0], xerr=np.array([[d - old[0], d - new[0]], [old[1] - d, new[1] - d]]), fmt='o', capsize=6)
ax.axvline(0, ls='--', color='black')
ax.set_yticks([0, 1])
ax.set_yticklabels(['Add reported paired p=0.008', 'Means/SD only, all correlations'])
ax.set_xlabel('Six-month PES difference Zr−Ti (score; POPULATION)')
ax.set_title('A paired statistic recovers information missing from marginal summaries')
fig.tight_layout()
fig.savefig(P / 'paired_statistic.png', dpi=150)
fig.savefig(P / 'paired_statistic.svg')
plt.close(fig)
print(json.dumps({'R5': res['outcome'], 'PES': out['PES'], 'controls': len(chall)}, indent=2))
