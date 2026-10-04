"""Frozen R2: case-events and paired case contrasts, not task-rate substitutes."""
import itertools
import json
import math
import os
from pathlib import Path
for _var in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[_var] = '4'
import numpy as np
from scipy import stats
ROOT = Path(__file__).resolve().parent

def exact_binomial(successes, n):
    if not 0 <= successes <= n or n < 1:
        raise ValueError('Invalid case-event counts')
    s = successes
    return dict(successes=s, n=n, point=s / n, ci95=[0.0 if s == 0 else float(stats.beta.ppf(0.025, s, n - s + 1)), 1.0 if s == n else float(stats.beta.ppf(0.975, s + 1, n - s))], one_sided95_lower=0.0 if s == 0 else float(stats.beta.ppf(0.05, s, n - s + 1)), one_sided95_upper=1.0 if s == n else float(stats.beta.ppf(0.95, s + 1, n - s)), assumption='IID sampled case events; population sampling and patient identity remain unverified')

def case_events(ids, values):
    groups = {}
    for (id, value) in zip(ids, values):
        groups.setdefault(id, []).append(float(value[0]))
    events = {id: all((v == 1 for v in rows)) for (id, rows) in groups.items()}
    return (events, exact_binomial(sum(events.values()), len(events)))

def contrast(ids, values, margin=None):
    groups = {}
    for (id, value) in zip(ids, values):
        groups.setdefault(id, []).append(float(value[0]))
    means = np.array([np.mean(groups[id]) for id in sorted(groups)])
    k = len(means)
    point = float(means.mean())
    se = float(means.std(ddof=1) / math.sqrt(k))
    result = dict(case_mean_differences=means.tolist(), k_cases=k, point_mm=point, case_t95=[point - stats.t.ppf(0.975, k - 1) * se, point + stats.t.ppf(0.975, k - 1) * se], case_t90=[point - stats.t.ppf(0.95, k - 1) * se, point + stats.t.ppf(0.95, k - 1) * se], assumptions='t/TOST: independent Gaussian case-mean differences. Sign flip: independent symmetric errors about null, not treatment randomization.', resolution_level='PER_ARCH -> POPULATION', time_scale='SIMULTANEOUS')
    flips = np.array(list(itertools.product([-1, 1], repeat=k)))
    distribution = (flips * means).mean(axis=1)
    result['exact_sign_flip_two_sided_p'] = float(np.mean(np.abs(distribution) >= abs(point) - 1e-12))
    result['sign_patterns'] = len(flips)
    result['exact_sign_flip_distribution_mm'] = distribution.tolist()
    if margin is not None:
        p_lower = float(stats.t.sf((point + margin) / se, k - 1)) if se else float(point <= -margin)
        p_upper = float(stats.t.cdf((point - margin) / se, k - 1)) if se else float(point >= margin)
        result.update(equivalence_margin_mm=margin, tost_p=max(p_lower, p_upper), equivalence_pass=max(p_lower, p_upper) < 0.05)
    return result

def run(raw):
    events = []
    contrasts = []
    controls = []
    for row in raw:
        if row['demo'] in {'PROOF_LANE', 'GENCAD_V2', 'GENCAD_V3'} and 'digital PASS' in row['metric']:
            (ev, interval) = case_events(row['cluster_ids'], row['values'])
            events.append(dict(demo=row['demo'], metric=row['metric'], analysis_prereg='PREREG_R3.json' if row['demo'] == 'GENCAD_V3' else 'PREREG_R2.json', event='All scheduled designs for this case return full digital PASS', case_events=ev, interval=interval, task_rate=float(np.mean(row['values'])), warning='Different estimand from task success; UNKNOWN/abstain are not PASS. No clinical viability claim.'))
            (_, duplicate) = case_events(row['cluster_ids'] * 2, row['values'] * 2)
            wrong = exact_binomial(1 if interval['successes'] == 0 else 0, interval['n'])
            controls.append(dict(name=row['demo'] + '/' + row['metric'] + '/case-event duplication', pass_=duplicate == interval, injected_rejected=wrong != interval, bad_success_count=wrong['successes']))
        if row['demo'] == 'X3' and 'paired p95 difference' in row['metric']:
            margin = 0.05 if 'collar_mesh' in row['metric'] else None
            c = contrast(row['cluster_ids'], row['values'], margin)
            contrasts.append(dict(demo='X3', metric=row['metric'], **c))
            dup = contrast(row['cluster_ids'] * 2, row['values'] * 2, margin)
            invariant = np.allclose(c['case_mean_differences'], dup['case_mean_differences'], atol=1e-12)
            if margin:
                bad = contrast(row['cluster_ids'], [[v[0] + 0.5] for v in row['values']], margin)
                controls.append(dict(name='X3 equivalence injection', pass_=invariant, injected_rejected=not bad['equivalence_pass'], injected_shift_mm=0.5, bad_tost_p=bad['tost_p']))
            else:
                bad = contrast(row['cluster_ids'], [[0.0] for v in row['values']])
                controls.append(dict(name='X3 sign-flip exact grid and shift', pass_=invariant and c['sign_patterns'] == 16, injected_rejected=bad['exact_sign_flip_two_sided_p'] == 1, injected_value_mm=0))
    for n in [2, 3, 8, 29]:
        upper = exact_binomial(0, n)['one_sided95_upper']
        expected = 1 - 0.05 ** (1 / n)
        controls.append(dict(name=f'Exact zero-event n={n}', expected=expected, observed=upper, pass_=abs(upper - expected) < 1e-10, injected_rejected=abs(exact_binomial(1, n)['one_sided95_upper'] - expected) > 1e-10))
    result = dict(preregs=['PREREG_R2.json', 'PREREG_R3.json'], claim_type='capability', case_events=events, paired_case_contrasts=contrasts, controls=controls, external_measurement='No new measurement. Original held-out scans remain external to X35; patient sampling and physical lab transport UNKNOWN.')
    original = {**result, 'prereg': 'PREREG_R2.json', 'case_events': [e for e in events if e['demo'] != 'GENCAD_V3'], 'controls': [c for c in controls if not c['name'].startswith('GENCAD_V3/')]}
    extension = dict(prereg='PREREG_R3.json', claim_type='capability', case_events=[e for e in events if e['demo'] == 'GENCAD_V3'], controls=[c for c in controls if c['name'].startswith('GENCAD_V3/')], scope='Versioned package expansion; same frozen statistical operations')
    (ROOT / 'RESULTS_R2.json').write_text(json.dumps(original, ensure_ascii=False, indent=2) + '\n')
    (ROOT / 'RESULTS_R3.json').write_text(json.dumps(extension, ensure_ascii=False, indent=2) + '\n')
    return result
