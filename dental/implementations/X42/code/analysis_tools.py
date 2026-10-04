import numpy as np, collections

def paired(rows, a, b, seed=6142, reps=2000):
    ix = {(r['case_key'], r['family'], r['level'], r['method']): r for r in rows}
    families = sorted({r['family'] for r in rows})
    out = []
    for fam in [*families, 'all']:
        grouped = collections.defaultdict(lambda : collections.defaultdict(list))
        counts = collections.Counter()
        for (key, ra) in ix.items():
            (ck, f, l, m) = key
            if m != a or (fam != 'all' and f != fam):
                continue
            rb = ix.get((ck, f, l, b))
            if rb is None:
                continue
            grouped[ck]['pass'].append(float(ra['verdict'] == 'PASS') - float(rb['verdict'] == 'PASS'))
            counts['paired_all_rows'] += 1
            if ra['verdict'] == 'PASS' and rb['verdict'] == 'PASS':
                for k in ['rmse_mm', 'contact_error_mm2']:
                    if ra.get(k) is not None and rb.get(k) is not None:
                        grouped[ck][k].append(ra[k] - rb[k])
                        counts['paired_' + k] += 1
        stat = {}
        for k in ['pass', 'rmse_mm', 'contact_error_mm2']:
            v = np.array([np.mean(d[k]) for d in grouped.values() if d[k]])
            if not len(v):
                stat[k] = None
                continue
            rng = np.random.default_rng(seed)
            s = np.array([v[rng.integers(len(v), size=len(v))].mean() for _ in range(reps)])
            stat[k] = {'mean_difference': float(v.mean()), 'case_bootstrap_95': np.quantile(s, [0.025, 0.975]).tolist(), 'case_clusters': len(v), 'resolution': 'POPULATION', 'direction': 'positive improves' if k == 'pass' else 'negative improves'}
        out.append({'family': fam, 'candidate': a, 'control': b, **dict(counts), 'paired': stat})
    return out

def decide(contrast):
    o = contrast['paired']
    p = o['pass']
    h = o['rmse_mm']
    c = o['contact_error_mm2']
    passgain = bool(p and p['mean_difference'] > 0)
    pareto = bool(h and c and (h['case_bootstrap_95'][1] <= 0) and (c['case_bootstrap_95'][1] <= 0) and (h['case_bootstrap_95'][1] < 0 or c['case_bootstrap_95'][1] < 0))
    return {'strict_pass_gain': passgain, 'height_contact_pareto': pareto, 'combined_gate': 'PASS' if passgain and pareto else 'FAIL'}
