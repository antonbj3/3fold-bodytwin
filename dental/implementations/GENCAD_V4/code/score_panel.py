from common import *
import gzip, collections, io, csv
sys.path.insert(0, str(DATA / 'participant_code'))
from task_io import attach
from legacy.checks import all_checks
from legacy.generators import submit
from quality import roof_metrics, nondominated
METRICS = ['anatomy_rmse_mm', 'contact_symdiff_mm2', 'negative_gap_area_mm2']
TOLS = [0.01, 0.05, 0.05]

def verified_npz(path, item):
    blob = Path(path).read_bytes()
    if hashlib.sha256(blob).hexdigest() != item['sha256']:
        raise ValueError('INPUT_HASH_MISMATCH ' + str(path))
    with np.load(io.BytesIO(blob), allow_pickle=False) as a:
        return dict(a)

def boot(v, seed=20261003):
    v = np.asarray(v, float)
    if len(v) == 0:
        return dict(mean=None, CI95=None, patients=0)
    rng = np.random.default_rng(seed)
    b = v[rng.integers(0, len(v), size=(2000, len(v)))].mean(1)
    return dict(mean=float(v.mean()), CI95=np.quantile(b, [0.025, 0.975]).tolist(), patients=len(v), degenerate=bool(np.ptp(b) == 0))

def run():
    from integrity import verify
    verify()
    start = time.perf_counter()
    cfg = read(ROOT / 'PREREG_R2.json')
    gf = read(ROOT / 'FROZEN_GENERATOR.json')
    pf = read(ROOT / 'FROZEN_PREDICTIONS.json')
    rows = []
    cohort = read(DATA / 'private/COHORT.json')
    private = {r['relative'].removeprefix('payload/'): r for r in read(ROOT / 'SOURCE_REUSE.json')['files']}
    public = gf['files']
    names = cfg['participants']
    for (ci, c) in enumerate(cohort):
        key = c['case_key']
        path = DATA / 'public/tasks' / f'{key}.json'
        blob = path.read_bytes()
        if hashlib.sha256(blob).hexdigest() != public['public/tasks/' + path.name]['sha256']:
            raise ValueError('task hash')
        ts = json.loads(blob)
        rr = 'private/references/' + key + '.npz'
        ref = verified_npz(DATA / rr, private[rr])
        preds = {}
        for name in names:
            rel = name + '/' + key + '.npz'
            a = verified_npz(DATA / 'predictions' / rel, pf['files'][rel])
            preds[name] = a
            if a['task_ids'].tolist() != [t['task_id'] for t in ts]:
                raise ValueError('Wrong task order')
            if len(a['status']) != len(ts):
                raise ValueError('Wrong status length')
            allowed = {'task_ids', 'status'}
            for (ii, ss) in enumerate(a['status']):
                if ss not in ['DESIGN', 'ABSTAIN', 'ERROR', 'UNKNOWN_SITE']:
                    raise ValueError('Unknown status')
                if ss == 'UNKNOWN_SITE' and ts[ii]['status'] == 'READY':
                    raise ValueError('Invented unknown site')
                if ss == 'DESIGN':
                    allowed.update(['outer_' + str(ii), 'inner_' + str(ii)])
            if set(a) != allowed:
                raise ValueError('Extra/missing prediction arrays')
        scenes = {}
        cache = {}
        for (i, t0) in enumerate(ts):
            t = t0
            if t['status'] == 'READY':
                rel = 'public/' + t['geometry_file']
                if rel not in scenes:
                    scenes[rel] = verified_npz(DATA / rel, public[rel])
                t = attach(t, scenes[rel])
            taskrows = []
            for (name, a) in preds.items():
                row = dict(case_key=key, patient_group=c['group'], dataset=c['dataset'], split=c['split'], task_id=t['task_id'], family=t['family'], level=t['level'], participant=name)
                s = str(a['status'][i])
                row['status'] = s
                if t['status'] != 'READY':
                    row.update(L1='UNKNOWN_SITE', reason=t['site_error'])
                elif s in ['ABSTAIN', 'ERROR']:
                    row.update(L1=s)
                elif s != 'DESIGN':
                    raise ValueError('Unsupported status ' + s)
                else:
                    z = a['outer_' + str(i)]
                    inn = a['inner_' + str(i)]
                    if z.shape != (len(t['xy']),) or inn.shape != z.shape:
                        raise ValueError('Prediction dimensions')
                    row['L1'] = all_checks(t, submit(t, z, inn))['validity']
                    ck = (t['family'], z.tobytes(), inn.tobytes())
                    if ck not in cache:
                        cache[ck] = roof_metrics(t, z, inn, ref[t['family']])
                    row.update(cache[ck])
                    row.update(cement_physical='UNKNOWN_NO_MATCHED_PROCESS_OR_SPATIAL_MEASUREMENT', fracture_margin='UNKNOWN_SETUP_MISMATCH', proximal='UNKNOWN_ROOF_HAS_NO_AXIAL_CONTACT_SURFACE', removal='UNKNOWN_ROOF_HAS_NO_COMPLETE_UNPREPARED_VOLUME')
                row['pareto_on_identified_metrics'] = None
                rows.append(row)
                taskrows.append(row)
            good = [r for r in taskrows if r['L1'] == 'PASS' and all((r.get(k) is not None for k in METRICS))]
            if good:
                front = nondominated([[r[k] for k in METRICS] for r in good], TOLS)
                for (r, b) in zip(good, front):
                    r['pareto_on_identified_metrics'] = bool(b)
        if ci % 8 == 0:
            print('scored', ci + 1, 'seconds', round(time.perf_counter() - start), flush=True)
    groups = collections.defaultdict(list)
    for r in rows:
        groups[r['dataset'], r['split'], r['family'], r['participant']].append(r)
    summaries = []
    for ((ds, sp, fam, name), rs) in sorted(groups.items()):
        counts = dict(collections.Counter((r['L1'] for r in rs)))
        clusters = collections.defaultdict(list)
        for r in rs:
            clusters[r['patient_group']].append(r)
        ms = {}
        for k in METRICS + ['contact_area_error_mm2']:
            vals = [np.mean([r[k] for r in v if r['L1'] == 'PASS' and r.get(k) is not None]) for v in clusters.values() if any((r['L1'] == 'PASS' and r.get(k) is not None for r in v))]
            ms[k] = boot(vals)
        summaries.append(dict(dataset=ds, split=sp, family=fam, participant=name, requested=len(rs), counts=counts, pass_fraction=boot([np.mean([r['L1'] == 'PASS' for r in v]) for v in clusters.values()]), quality_on_pass=ms, pareto_fraction_of_all_requested=sum((r['pareto_on_identified_metrics'] is True for r in rs)) / len(rs)))
    pairs = []
    bytask = collections.defaultdict(dict)
    for r in rows:
        bytask[r['dataset'], r['split'], r['family'], r['task_id']][r['participant']] = r
    for (ds, sp) in [('Bite2Text', 'test'), ('Bite2Text', 'dev'), ('Bits2Bites', 'auxiliary')]:
        for fam in sorted({r['family'] for r in rows}):
            ts = [p for ((d, s, f, t), p) in bytask.items() if (d, s, f) == (ds, sp, fam)]
            for (ni, a) in enumerate(names):
                for b in names[ni + 1:]:
                    for metric in METRICS:
                        cl = collections.defaultdict(list)
                        n = 0
                        for p in ts:
                            (x, y) = (p[a], p[b])
                            if x['L1'] == y['L1'] == 'PASS' and x.get(metric) is not None and (y.get(metric) is not None):
                                cl[x['patient_group']].append(x[metric] - y[metric])
                                n += 1
                        stat = boot([np.mean(v) for v in cl.values()])
                        pairs.append(dict(dataset=ds, split=sp, family=fam, a=a, b=b, metric=metric, paired_tasks=n, **stat))
    qualified = [p for p in pairs if p['split'] == 'test' and p['CI95'] and (p['patients'] >= 5) and (p['CI95'][0] > 0 or p['CI95'][1] < 0)]
    out = dict(claim_type='capability', decision='SPATIAL_QUALITY_DISCRIMINATES' if qualified else 'NO_IDENTIFIED_QUALITY_CONTRAST', case_count=len(cohort), requested_tasks=len(rows) // len(names), participant_rows=len(rows), summaries=summaries, paired_contrasts=pairs, nonzero_test_contrasts=len(qualified), test_intervals_scope='Descriptive, no multiplicity correction or algorithm superiority inference; patient-case resampling, conditional on both designs passing L1', dropout={'not_selected_from_v3': 1191 - len(cohort), 'selected_unknown_site_tasks': sum((r['L1'] == 'UNKNOWN_SITE' for r in rows)) // len(names), 'participant_ERROR_rows': sum((r['L1'] == 'ERROR' for r in rows)), 'original_v3_source_exclusions': {'empty_STL': 1, 'duplicate_arch_group': 2}}, external_referent=cfg['external_referent'], cost_seconds=time.perf_counter() - start, identified_quality_axes=METRICS, not_identified=['full-crown proximal contacts in roof panel', 'substance removal in roof panel', 'physical cement field', 'generated-crown fracture margin'], reference_scope='Geometric contact against original scan in supplied pose. Registration is measured but physical loaded contact is not. PL discretization and unknown scanner/FDI errors remain.')
    dump(ROOT / 'rounds/R2.json', out)
    keys = list(dict.fromkeys((k for r in rows for k in r)))
    s = io.StringIO()
    w = csv.DictWriter(s, keys)
    w.writeheader()
    w.writerows(rows)
    (ROOT / 'raw/QUALITY_ROWS.csv.gz').write_bytes(gzip.compress(s.getvalue().encode(), mtime=0))
    dump(ROOT / 'raw/QUALITY_ROWS.json', rows)
    state('R2_DECIDED', out['decision'], 'Construct and score axial whole-crown surfaces; retain missing physical endpoints')
    print(out['decision'], out['nonzero_test_contrasts'])
if __name__ == '__main__':
    run()
