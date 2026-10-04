import json, csv, time, math, resource
from fractions import Fraction as Q
import numpy as np
from geometry import R, write, sha, FDI14
KNOTS = np.array([0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1.0])
REF = {'R': np.array([25.2, 33.6, 38.3, 45.6, 49.6, 55.5, 60.2, 65.8, 77.0]) / 100, 'posterior': np.array([39.8, 54.5, 62.4, 75.8, 82.4, 87.3, 90.4, 94.0, 99.6]) / 100}

def w1_bounds(x, knots=KNOTS, vals=None, rank=1 / 178, rounding=0.0005):
    x = np.sort(x)
    n = len(x)
    cuts = np.unique(np.clip(np.r_[np.arange(n + 1) / n, knots - rank, knots + rank, 0, 1], 0, 1))
    lo = hi = 0.0
    for (a, b) in zip(cuts[:-1], cuts[1:]):
        u = (a + b) / 2
        v = x[min(int(u * n), n - 1)]
        left = np.where(knots <= u - rank)[0]
        right = np.where(knots >= u + rank)[0]
        ql = max(0, vals[left[-1]] - rounding) if len(left) else 0.0
        qr = min(1, vals[right[0]] + rounding) if len(right) else 1.0
        lo += (b - a) * max(ql - v, v - qr, 0)
        hi += (b - a) * max(abs(v - ql), abs(v - qr))
    return [lo, hi]

def quant(x):
    return np.quantile(x, KNOTS[1:-1], method='linear')

def metrics(x, tag):
    q = quant(x)
    target = REF[tag][1:-1]
    default = 0.5 if tag == 'R' else 8 / 14
    e = np.abs(q - target)
    base = float(np.max(abs(default - target)))
    return dict(quantile_levels=KNOTS[1:-1].tolist(), predicted_quantiles_share=q.tolist(), external_quantiles_share=target.tolist(), quantile_error_pp=(e * 100).tolist(), max_error_pp=float(e.max() * 100), default_max_error_pp=base * 100, improvement_pp=(base - e.max()) * 100, W1_partial_identification_bounds_share=w1_bounds(x, vals=REF[tag]), sample_quantile_gate='PASS' if e.max() <= 0.1 else 'FAIL', gain_gate='PASS' if base - e.max() >= 0.05 else 'FAIL', comparison='descriptive unmatched cohorts; not patient force validation')

def main():
    start = time.perf_counter()
    assert sha(R / 'PREREG_R2.json') == (R / 'PREREG_R2.json.sha256').read_text().split()[0]
    rows = list(csv.DictReader(open(R / 'raw/GEOMETRY_PROFILES.csv')))
    cases = {}
    for r in rows:
        if r['area_usable'] != 'True':
            continue
        cases.setdefault(r['case'], {}).setdefault(r['jaw'], {})[int(r['fdi'])] = float(r['area_share'])
    partition = json.load(open(R / 'raw/PATIENT_PARTITION.json'))['partition']
    groups = {**partition, 'all': sorted({r['case'] for r in rows})}
    sets = [set(v) for (k, v) in partition.items()]
    assert all((not a & b for (i, a) in enumerate(sets) for b in sets[i + 1:]))
    out = {}
    rng = np.random.default_rng(942)
    for (group, ids) in groups.items():
        used = sorted(set(ids) & set(cases))
        data = {}
        for j in ['upper', 'lower']:
            obs = {tag: np.array([sum((v for (t, v) in cases[c][j].items() if (t // 10 in (1, 4) if tag == 'R' else t % 10 >= 4))) for c in used]) for tag in REF}
            data[j] = {tag: metrics(x, tag) for (tag, x) in obs.items()}
            for (tag, x) in obs.items():
                bootstrap = np.array([quant(x[rng.integers(len(x), size=len(x))]) for _ in range(2000)])
                data[j][tag]['bootstrap95_quantiles_share'] = np.quantile(bootstrap, [0.025, 0.975], axis=0).T.tolist()
        out[group] = dict(n_partition=len(ids), n_usable=len(used), excluded_zero_area=len(ids) - len(used), groups_disjoint=group != 'all', metrics=data)
    primary = out['test']['metrics']['upper']
    gate = 'PASS' if all((v['sample_quantile_gate'] == 'PASS' and v['gain_gate'] == 'PASS' for v in primary.values())) else 'FAIL'
    frozen = json.load(open(R / 'FROZEN_PREDICTIONS.json'))
    old = json.load(open(R.parent / 'LANE_X2_OCCLUSION_B2B/raw/EXTERNAL_REFERENTS.json'))['studies']
    (ferr, hatt) = old[:2]
    diagnostics = []
    f = {t: v / 100 for (t, v) in zip(ferr['fdi'], ferr['mean_percent'])}
    fm = np.array([f[t] for t in FDI14['upper']])
    pm = np.array(frozen['profiles']['upper']['area_mean'])
    diagnostics.append(dict(source='Ferrato2017_Table3', raw_mass=float(fm.sum()), mass_deficit=float(1 - fm.sum()), raw_half_L1=float(np.abs(pm - fm).sum() / 2), categorical_TV='NOT_DEFINED: published14 means sum93%, not a probability distribution', normalization_gate='FAIL', conditional_shape_sensitivity_TV=float(np.abs(pm - fm / fm.sum()).sum() / 2), conditional_shape_status='Sensitivity only: normalized sum of means is not mean of conditional subject shares'))
    lower = np.array(frozen['profiles']['lower']['area_mean'])
    pooled = (lower[:7] + lower[7:]) / 2
    ref = np.array(hatt['single_side_mean_percent']) / 100
    diagnostics.append(dict(source='Hattori1996_Figure_open_bars', pooled_side_mass=float(ref.sum()), categorical_shape_TV=float(np.abs(2 * pooled - 2 * ref).sum() / 2), uniform7_shape_TV=float(np.abs(np.ones(7) / 7 - 2 * ref).sum() / 2), digitization='±1pp per single-side tooth, no raw individuals', symmetry_assumption='Pooled side comparison only, no predicted fixed left/right population target', physical_population_validation='UNKNOWN'))
    dkw = dict(external_n=178, external_epsilon=math.sqrt(math.log(2 / 0.025) / (2 * 178)), local_n=out['test']['n_usable'], local_epsilon=math.sqrt(math.log(2 / 0.025) / (2 * out['test']['n_usable'])), meaning='CDF confidence under iid cases and matched observable, not a remedy for cohort/protocol mismatch')
    mutation = {tag: metrics(np.zeros(out['test']['n_usable']), tag)['sample_quantile_gate'] == 'FAIL' for tag in REF}
    result = dict(primary_axis='upper', primary_group='test', sample_shape_gate=gate, physical_validation='UNKNOWN: M3 handling/selected Bite2Text malocclusion/maximum-force protocol and observation transfer unmatched', groups=out, source_quantile_table=dict(url='https://onlinelibrary.wiley.com/doi/full/10.1111/jopr.13838', doi='10.1111/jopr.13838', locator='Table1', n=178, knots=KNOTS.tolist(), values_share={k: v.tolist() for (k, v) in REF.items()}, rounding_share=0.0005, rank_slack=1 / 178, protocol='100um T-Scan III/Novus sensor; participant average of MIP frame over2–4 recordings; max-force single/multibite; customized maxillary chart; at least24 natural teeth; M3 contribution not specified'), categorical_diagnostics=diagnostics, DKW=dkw, injected_profiles_rejected=mutation, wall_s=time.perf_counter() - start, max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    write('raw/EXTERNAL_VALIDATION.json', result)
    write('CURRENT_WORK_STATE.json', dict(lane='X94-force-validation', status='DISTRIBUTION_TEST_COMPLETE', latest_gate=gate, next_operation='Changed construction: contact-graph transport of measured force constraints without support recovery', frozen_predictions_sha256=sha(R / 'FROZEN_PREDICTIONS.json')))
    print(json.dumps(dict(gate=gate, n=out['test']['n_usable'], primary=primary, categorical_diagnostics=diagnostics)))
if __name__ == '__main__':
    main()
