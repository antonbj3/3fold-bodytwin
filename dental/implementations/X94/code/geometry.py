from dental_release.paths import expand as _release_expand
import json, hashlib, time, resource, sys, csv
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
R = Path(__file__).resolve().parents[1]
X21 = R.parent / _release_expand('X21')
PARENT = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/FALT_TANDLAST'))
FDI = {j: [10 * q + t for q in ((1, 2) if j == 'upper' else (4, 3)) for t in range(1, 9)] for j in ['upper', 'lower']}
FDI14 = {j: [t for t in v if t % 10 != 8] for (j, v) in FDI.items()}

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(p, d):
    p = R / p
    p.parent.mkdir(exist_ok=True, parents=True)
    p.write_text(json.dumps(d, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + '\n')

def load():
    manifest = json.load(open(X21 / 'raw/INPUT_MANIFEST.json'))
    cases = []
    sources = {str(X21 / 'raw/INPUT_MANIFEST.json'): sha(X21 / 'raw/INPUT_MANIFEST.json')}
    for c in manifest['cases']:
        p = X21 / 'raw/cases' / f'{c}.json'
        s = json.load(open(p))
        sources[str(p)] = sha(p)
        p = X21 / 'raw/refined_certificates' / f'{c}.json'
        if not p.exists():
            p = X21 / 'raw/certificates' / f'{c}.json'
        b = json.load(open(p))
        sources[str(p)] = sha(p)
        e = [(int(r['upper_fdi']), int(r['lower_fdi'])) for r in b['pairs']]
        assert len(e) == len(set(e)) and e
        area = {(int(r['upper_fdi']), int(r['lower_fdi'])): Q(r['absolute_near_area_mm2']) for r in s['pairs']}
        assert set(e) == set(area)
        cases.append(dict(case=c, edges=e, gap_rows=b['pairs'], area=area, stats=s['stats'], cost=s['cost']))
    return (manifest, cases, sources)

def hull(e, teeth):
    I = [u in teeth or l in teeth for (u, l) in e]
    return [int(all(I)), int(any(I))]

def summary_labels():
    return ['I', 'C', 'P', 'M', 'anterior', 'posterior', 'R', 'L']

def teeth_for(jaw, tag, third=True):
    pool = FDI[jaw] if third else FDI14[jaw]
    return {t for t in pool if (t % 10 <= 2 if tag == 'I' else t % 10 == 3 if tag == 'C' else 4 <= t % 10 <= 5 if tag == 'P' else t % 10 >= 6 if tag == 'M' else t % 10 <= 3 if tag == 'anterior' else t % 10 >= 4 if tag == 'posterior' else t // 10 in (1, 4) if tag == 'R' else t // 10 in (2, 3))}

def aggregation(e, w, j):
    total = sum(w, Q(0))
    out = {t: Q(0) for t in FDI14[j]}
    if total == 0:
        return None
    ix = 0 if j == 'upper' else 1
    for (edge, v) in zip(e, w):
        out[edge[ix]] += v / total
    assert sum(out.values()) == 1
    return [float(out[t]) for t in FDI14[j]]

def main():
    start = time.perf_counter()
    assert sha(R / 'PREREG_R1.json') == (R / 'PREREG_R1.json.sha256').read_text().split()[0]
    (manifest, cases, sources) = load()
    rng = np.random.default_rng(94)
    proxy = {j: [] for j in FDI}
    degree = {j: [] for j in FDI}
    keep = []
    rows = []
    hulls = {}
    for third in [True, False]:
        for j in FDI:
            dims = [str(t) for t in (FDI[j] if third else FDI14[j])] + summary_labels()
            allh = []
            for c in cases:
                e = c['edges'] if third else [e for e in c['edges'] if e[0] % 10 != 8 and e[1] % 10 != 8]
                if not e:
                    raise ValueError('No non-M3 edge')
                allh.append([hull(e, {int(tag)} if tag.isdigit() else teeth_for(j, tag, third)) for tag in dims])
            a = np.array(allh)
            hulls[j + ('_full16' if third else '_conditional14')] = dict(labels=dims, cohort_hull=np.mean(a, axis=0).tolist(), max_width=float(np.max(np.mean(a[:, :, 1] - a[:, :, 0], axis=0))), level='POPULATION', meaning='Sharp scalar extrema separately, not independent simultaneous box')
    for c in cases:
        e = [e for e in c['edges'] if e[0] % 10 != 8 and e[1] % 10 != 8]
        w = [c['area'][edge] for edge in e]
        a = {j: aggregation(e, w, j) for j in FDI}
        d = {j: aggregation(e, [Q(1)] * len(e), j) for j in FDI}
        for j in FDI:
            degree[j].append(d[j])
        usable = a['upper'] is not None
        if usable:
            keep.append(c['case'])
            for j in FDI:
                proxy[j].append(a[j])
        for j in FDI:
            for (t, v, deg) in zip(FDI14[j], a[j] or [None] * 14, d[j]):
                rows.append(dict(case=c['case'], jaw=j, fdi=t, area_share=v, incidence_share=deg, area_usable=usable, level='PER_TOOTH', force_status='UNMEASURED_PROXY'))
    with open(R / 'raw/GEOMETRY_PROFILES.csv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    delta = 0.0
    for c in cases:
        d = {j: {t: 0.0 for t in FDI14[j]} for j in FDI}
        total = 0.0
        for (edge, weight) in c['area'].items():
            if edge[0] % 10 == 8 or edge[1] % 10 == 8:
                continue
            v = float(weight)
            total += v
            d['upper'][edge[0]] += v
            d['lower'][edge[1]] += v
        if total == 0:
            continue
        k = keep.index(c['case'])
        for j in FDI:
            delta = max(delta, max((abs(d[j][t] / total - proxy[j][k][i]) for (i, t) in enumerate(FDI14[j]))))
    pred = {j: dict(fdi=FDI14[j], area_mean=np.mean(proxy[j], axis=0).tolist(), incidence_mean=np.mean(degree[j], axis=0).tolist(), default_uniform=[1 / 14] * 14) for j in FDI}
    for j in FDI:
        x = np.array(proxy[j])
        boot = np.array([np.mean(x[rng.integers(len(x), size=len(x))], axis=0) for _ in range(2000)])
        pred[j]['area_mean_bootstrap95'] = np.quantile(boot, [0.025, 0.975], axis=0).T.tolist()
    frozen = dict(schema='X94-predictions-v1', prereg_sha256=sha(R / 'PREREG_R1.json'), data_sha256=sha(R / 'raw/GEOMETRY_PROFILES.csv'), frozen_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(), external_search_performed_in_X94=False, external_known_from_predecessors=['Ferrato 2017', 'Hattori 1996'], n_cases=len(cases), n_area_usable=len(keep), profiles=pred, hulls=hulls, source_sha256=sources)
    fp = R / 'FROZEN_PREDICTIONS.json'
    if fp.exists():
        old = json.load(open(fp))
        assert old['profiles'] == pred and old['hulls'] == hulls and (old['source_sha256'] == sources)
    else:
        write('FROZEN_PREDICTIONS.json', frozen)
        (R / 'FROZEN_PREDICTIONS.json.sha256').write_text(sha(fp) + '  FROZEN_PREDICTIONS.json\n')
    partition = manifest['inherited_partition']
    write('raw/PATIENT_PARTITION.json', dict(partition=partition, fit='NONE', geometry_test_reuse=True, patient_force_pairing='NONE', duplicates='Inherited IDs are groups; biological donor identity not independently verified'))
    out = dict(n_cases=len(cases), n_area_usable=len(keep), area_dropout=dict(n=len(cases) - len(keep), fraction=(len(cases) - len(keep)) / len(cases), reason='zero non-M3 absolute-near area'), independent_accumulation_max_error=delta, hull_usefulness='FAIL' if max((x['max_width'] for x in hulls.values())) > 0.2 else 'PASS', labelled_face_fraction_mean={j: float(np.mean([c['stats'][j]['labelled_face_fraction'] for c in cases])) for j in FDI}, upstream_geometry_compute_wall_s=sum((c['cost']['wall_s'] for c in cases)), conditional_geometry=True, physical_force_identified=False, predictions_sha256=sha(fp), wall_s=time.perf_counter() - start, max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    write('raw/GEOMETRY_SUMMARY.json', out)
    write('CURRENT_WORK_STATE.json', dict(lane='X94-force-validation', status='OWN_GEOMETRY_PREDICTIONS_FROZEN', latest_gate=out['hull_usefulness'], next_operation='Run exact identical-geometry support counterexample, then search primary measurement distributions', frozen_predictions_sha256=sha(fp)))
    print(json.dumps({k: v for (k, v) in out.items() if k not in ['labelled_face_fraction_mean']}))
if __name__ == '__main__':
    main()
