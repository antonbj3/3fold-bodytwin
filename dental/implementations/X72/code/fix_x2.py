from dental_release.paths import expand as _release_expand
from common import *
import numpy as np
import time

def classify(r):

    def cl(x):
        return 'III' if 'III' in x else 'II' if 'II' in x else 'I'
    a = cl(r['annotation']['Left Class'])
    b = cl(r['annotation']['Right Class'])
    return a if a == b else 'Mixed'

def family_gate(n, evaluated):
    assert n == evaluated == 18, 'WRONG_SIMULTANEOUS_FAMILY'
    return 0.05 / n

def main():
    t = time.perf_counter()
    alpha = family_gate(18, 18)
    rounds = []
    try:
        family_gate(12, 18)
        inj = False
    except AssertionError:
        inj = True
    for (rname, directory, outdir) in [('R1', 'cases', 'round1'), ('R2', 'cases_r2', 'round2'), ('R4', 'cases_3d', 'round4')]:
        rows = [read(p) for p in sorted((BASE / 'LANE_X2_OCCLUSION_B2B/raw' / directory).glob('*.json'))]
        assert len(rows) == 200
        old = read(BASE / _release_expand('X2') / outdir / 'results.json')
        old_angles = [c for c in old['contrasts'] if 'contrast' in c]
        assert len(old_angles) == 18
        corrected = []
        fronts = []
        maxerr = 0
        for (j, band) in enumerate([0.05, 0.1, 0.2]):
            eligible = [x for x in rows if x['metrics'][j]['mechanics'] is not None]
            rng = np.random.default_rng(6102)

            def bootstrap(x, y):
                if min(len(x), len(y)) < 10:
                    return None
                x = np.asarray(x)
                y = np.asarray(y)
                return x[rng.integers(0, len(x), (2000, len(x)))].mean(1) - y[rng.integers(0, len(y), (2000, len(y)))].mean(1)

            def anterior(w):
                return float((w[[0, 1, 4, 5]].sum() + w[[8, 9, 12, 13]].sum()) / 2)
            vals = {k: [] for k in ['Open Bite', 'Normal', 'Deep Bite', 'Inverted Bite']}
            for x in eligible:
                vals[x['annotation']['Anterior Bite']].append(anterior(np.array(x['metrics'][j]['area_shares_pp'])))
            boot = bootstrap(vals['Open Bite'], vals['Normal'])
            fronts.append(dict(band_mm=band, unchanged_interval=np.quantile(boot, [0.025, 0.975]).tolist() if boot is not None else None))
            for method in ['area', 'mechanics']:
                v = {k: [] for k in ['I', 'II', 'III']}
                for x in eligible:
                    c = classify(x)
                    if c in v:
                        m = x['metrics'][j]
                        v[c].append(anterior(np.array(m['area_shares_pp'] if method == 'area' else m['mechanics']['shares_pp'])))
                for (c1, c2) in [('II', 'I'), ('III', 'I'), ('III', 'II')]:
                    boot = bootstrap(v[c1], v[c2])
                    prev = next((c for c in old_angles if c['band_mm'] == band and c['method'] == method and (c['contrast'] == c1 + '-' + c2)))
                    if boot is None:
                        corrected.append(dict(band_mm=band, method=method, contrast=c1 + '-' + c2, status='UNKNOWN_SMALL_GROUP', n=[len(v[c1]), len(v[c2])], gate_before=prev['gate'], gate_after=False))
                        continue
                    before = np.quantile(boot, [0.05 / 24, 1 - 0.05 / 24])
                    after = np.quantile(boot, [alpha / 2, 1 - alpha / 2])
                    err = float(np.max(np.abs(before - np.array(prev['interval']))))
                    maxerr = max(maxerr, err)
                    mean = float(np.mean(v[c1]) - np.mean(v[c2]))
                    gate = abs(mean) >= 5 and (after[0] > 0 or after[1] < 0)
                    corrected.append(dict(band_mm=band, method=method, contrast=c1 + '-' + c2, mean_pp=mean, n=[len(v[c1]), len(v[c2])], before_interval_pp=before, after_interval_pp=after, old_interval_replay_error_pp=err, gate_before=prev['gate'], gate_after=bool(gate), significant_before=bool(before[0] > 0 or before[1] < 0), significant_after=bool(after[0] > 0 or after[1] < 0)))
        assert maxerr == 0, 'BOOTSTRAP_REPLAY_DRIFT'
        systematic = {}
        for method in ['area', 'mechanics']:
            for c in ['II-I', 'III-I', 'III-II']:
                group = [x for x in corrected if x['method'] == method and x['contrast'] == c]
                systematic[method + ':' + c] = dict(before=all((x['gate_before'] for x in group)), after=all((x['gate_after'] for x in group)))
        rounds.append(dict(round=rname, angles=corrected, fronts=fronts, systematic=systematic, old_replay_error_pp=maxerr, changed_gates=sum((x['gate_before'] != x['gate_after'] for x in corrected)), lost_significance=sum((x.get('significant_before', False) and (not x.get('significant_after', False)) for x in corrected)), dropout=dict(total=200, eligible_by_band=[sum((x['metrics'][j]['mechanics'] is not None for x in rows)) for j in range(3)], mixed_Angle=sum((classify(x) == 'Mixed' for x in rows)))))
    out = dict(claim_type='capability', resolution='POPULATION', rounds=rounds, alpha_per_contrast=alpha, family_size=18, nominal_family_coverage=0.95, coverage_status='BONFERRONI_BOOTSTRAP_NOMINAL; not rigorous finite-sample coverage;2000replicates imply unstable extreme quantiles', injections=dict(wrong12_family_rejected=inj), changed_conclusion=any((any((v['before'] != v['after'] for v in r['systematic'].values())) for r in rounds)), external_referent=dict(kind='published_dataset', locator='https://ditto.ing.unimore.it/bits2bites/; source raw case files in SOURCE_MANIFEST.json', compared_quantity='Angle-annotated anterior contact share; 18 pre-existing contrasts per round', refutes_us=True), cost_wall_s=time.perf_counter() - t)
    write(ROOT / 'raw/X2.json', out)
    state('X2_DONE', 'Exact original bootstrap replay; multiplicity18 evaluated', 'Prove corner licenses and replay contact port')
    print('X2', [(r['round'], r['changed_gates'], r['lost_significance']) for r in rounds])
if __name__ == '__main__':
    main()
