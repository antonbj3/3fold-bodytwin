from common import *
import csv, gzip, collections, itertools, math
from scipy.stats import norm
Z_POWER = float((norm.ppf(0.975) + norm.ppf(0.8)) ** 2)

def load_rows():
    rows = list(csv.DictReader(gzip.open(V4 / 'raw/QUALITY_ROWS.csv.gz', 'rt')))
    for r in rows:
        for k in ['anatomy_rmse_mm', 'contact_symdiff_mm2', 'negative_gap_area_mm2', 'contact_support_mm2', 'predicted_contact_mm2']:
            r[k] = float(r[k]) if r.get(k) else None
    return rows

def entropy(p):
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p) if 0 < p < 1 else 0.0

def n_required(v):
    v = np.asarray(v, float)
    if len(v) < 3:
        return (None, 'INSUFFICIENT_CLUSTERS')
    mu = float(v.mean())
    var = float(v.var(ddof=1))
    if abs(mu) < 1e-14:
        return (None, 'ZERO_OBSERVED_EFFECT')
    if var < 1e-20:
        return (None, 'ZERO_EMPIRICAL_VARIANCE_NOT_A_POPULATION_CERTIFICATE')
    return (max(3, math.ceil(Z_POWER * var / mu ** 2)), 'APPROXIMATE')

def paired_plan(v):
    v = np.asarray(v, float)
    (n, st) = n_required(v)
    rng = np.random.default_rng(6105)
    boot = v[rng.integers(0, len(v), size=(1000, len(v)))] if len(v) else np.empty((0, 0))
    required = [n_required(x)[0] for x in boot]
    finite = [x for x in required if x is not None]
    mu = float(v.mean()) if len(v) else None
    sd = float(v.std(ddof=1)) if len(v) > 1 else None
    z = abs(mu) * np.sqrt(len(v)) / sd if sd and sd > 0 else None
    return dict(clusters=len(v), mean_difference=mu, sd_difference=sd, mean_CI95=np.quantile(boot.mean(1), [0.025, 0.975]) if len(v) else None, clusters_for_80pct=n, status=st, planning_power_at_current_n=float(norm.cdf(z - norm.ppf(0.975)) + norm.cdf(-z - norm.ppf(0.975))) if z is not None else None, bootstrap_finite_n_range95=np.quantile(finite, [0.025, 0.975]) if finite else None, bootstrap_unbounded_fraction=1 - len(finite) / 1000, scope='Retrospective normal-approximation paired cluster planning; not rigorous enclosure or held-out power validation')

def run():
    start = time.perf_counter()
    rows = load_rows()
    test = [r for r in rows if r['split'] == 'test']
    names = sorted({r['participant'] for r in test})
    by = collections.defaultdict(dict)
    for r in test:
        by[r['task_id']][r['participant']] = r
    keys = sorted(by)
    mat = np.array([[by[k][n]['L1'] == 'PASS' for n in names] for k in keys], float)
    items = []
    for (i, k) in enumerate(keys):
        p = mat[i].mean()
        total = np.delete(mat, i, axis=0).mean(0)
        corr = float(np.corrcoef(mat[i], total)[0, 1]) if np.std(mat[i]) and np.std(total) else None
        r = by[k][names[0]]
        items.append(dict(task_id=k, patient_group=r['patient_group'], family=r['family'], level=r['level'], difficulty_failure_fraction=1 - p, pass_probability=p, entropy_bits=entropy(p), discriminates=bool(0 < p < 1), corrected_item_total_correlation=corr, resolution='PER_TOOTH', status_counts=dict(collections.Counter((x['L1'] for x in by[k].values())))))
    families = []
    for fam in sorted({r['family'] for r in test}):
        ii = [i for (i, k) in enumerate(keys) if by[k][names[0]]['family'] == fam]
        mm = mat[ii]
        rates = mm.mean(0)
        info = [items[i] for i in ii]
        families.append(dict(family=fam, requested_items=len(ii), clusters=len({by[keys[i]][names[0]]['patient_group'] for i in ii}), participant_pass_rates=dict(zip(names, rates)), participant_range=float(np.ptp(rates)), discriminating_items=sum((x['discriminates'] for x in info)), constant_items=sum((not x['discriminates'] for x in info)), mean_entropy_bits=float(np.mean([x['entropy_bits'] for x in info])), designation='COVERAGE_SENTINEL_NO_PASS_DISCRIMINATION' if not any((x['discriminates'] for x in info)) else 'DISCRIMINATING_IN_THIS_PANEL', resolution='POPULATION'))
    pairs = []
    for fam in [r['family'] for r in families]:
        tasks = [v for v in by.values() if v[names[0]]['family'] == fam]
        for (a, b) in itertools.combinations(names, 2):
            for metric in ['pass', 'negative_gap_area_mm2', 'contact_symdiff_mm2']:
                cl = collections.defaultdict(list)
                excluded = 0
                for t in tasks:
                    (ra, rb) = (t[a], t[b])
                    if metric == 'pass':
                        d = float(ra['L1'] == 'PASS') - float(rb['L1'] == 'PASS')
                    elif ra['L1'] == rb['L1'] == 'PASS' and ra.get(metric) is not None and (rb.get(metric) is not None):
                        d = ra[metric] - rb[metric]
                    else:
                        excluded += 1
                        continue
                    cl[ra['patient_group']].append(d)
                vals = [np.mean(v) for v in cl.values()]
                pairs.append(dict(family=fam, a=a, b=b, metric=metric, excluded_task_pairs=excluded, requested_task_pairs=len(tasks), resolution='POPULATION', **paired_plan(vals)))
    out = dict(claim_type='capability', families=families, item_count=len(items), participants=names, patient_case_clusters=len({r['patient_group'] for r in test}), rows=len(test), pairs=pairs, dropout=dict(all_rows=len(rows), primary_rows=len(test), nonprimary_rows=len(rows) - len(test), primary_unknown_site_rows=sum((r['L1'] == 'UNKNOWN_SITE' for r in test)), policy='Non-PASS retained as zero pass; continuous metrics pairwise both-PASS; no post-hoc family removal'), seconds=time.perf_counter() - start)
    dump(ROOT / 'raw/ITEMS.json', items)
    dump(ROOT / 'raw/R1.json', out)
    (ROOT / 'history/HANDOFF_R1.md').write_text('R1 completed: descriptive item information and paired patient-case planning saved in raw/R1.json. Fixed thresholds and all nine families retained. Next construction changes information: external pretrained SDF generator. No general field ranking or IRT latent ability established.\n')
    state('R1_DECIDED', 'Item discrimination and cluster sample-size estimates computed', 'Run external pretrained model and score spatial function')
    return out
if __name__ == '__main__':
    run()
