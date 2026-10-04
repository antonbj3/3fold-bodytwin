"""All official score and submission entry points construct a verified Bundle."""
import collections, gzip, csv, io, time, resource
from pathlib import Path
import numpy as np
from integrity import Bundle, IntegrityError, predictions, safe_bytes, h, load_npz
from util import ROOT, PAYLOAD, FAMILIES, LEVELS, dump, digest, read
from task_io import attach
from legacy.checks import all_checks
from fast_geometry import feasibility
from case_statistics import cluster_summary, paired_case_difference

def decode(arrays, tasks):
    if arrays.get('task_ids', np.array([])).tolist() != [t['task_id'] for t in tasks]:
        raise IntegrityError('prediction IDs/order/duplicates')
    status = arrays.get('status', np.array([]))
    if status.shape != (len(tasks),) or status.dtype.kind != 'U':
        raise IntegrityError('status dimensions/type')
    allowed = {'task_ids', 'status'}
    for (i, s) in enumerate(status):
        if s not in ['DESIGN', 'ABSTAIN', 'UNKNOWN_SITE']:
            raise IntegrityError('unknown submission status')
        if s == 'DESIGN':
            allowed.update(['outer_' + str(i), 'inner_' + str(i)])
        if s == 'UNKNOWN_SITE' and tasks[i]['status'] == 'READY':
            raise IntegrityError('invented unknown site')
    if set(arrays) != allowed:
        raise IntegrityError('extra/missing prediction arrays')
    return status

def design(a, i, t, status):
    if status != 'DESIGN':
        return dict(task_id=t['task_id'], status='ABSTAIN', reason='frozen participant abstention')
    z = a['outer_' + str(i)]
    inn = a['inner_' + str(i)]
    n = len(t['xy'])
    if z.shape != (n,) or inn.shape != (n,) or z.dtype.kind not in 'fi' or (inn.dtype.kind not in 'fi') or (not np.isfinite(z).all()) or (not np.isfinite(inn).all()):
        raise IntegrityError('invalid coordinate shape/type/nonfinite')
    return dict(task_id=t['task_id'], status='DESIGN', units='mm', frame=t['frame'], outer_vertices=np.c_[t['xy'], z], inner_vertices=np.c_[t['xy'], inn], faces=t['faces'])

def metrics(t, z, ref, band):
    valid = np.isfinite(ref)
    w = t['weights']
    den = float(w[valid].sum())
    if den <= 0:
        return (None, None, 0.0)
    try:
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            rmse = float(np.sqrt(np.sum(w[valid] * (z[valid] - ref[valid]) ** 2) / den))
            ceiling = t['ceiling']
            ok = valid & np.isfinite(ceiling)
            native_gap = ceiling - ref
            pred_gap = ceiling - z
            gold = ok & (native_gap >= 0) & (native_gap <= band)
            pred = ok & (pred_gap >= 0) & (pred_gap <= band)
            area = abs(float(w[pred].sum() - w[gold].sum())) if ok.any() else None
    except FloatingPointError as e:
        raise IntegrityError('unrepresentable metric arithmetic') from e
    if not np.isfinite(rmse) or (area is not None and (not np.isfinite(area))):
        raise IntegrityError('nonfinite reported metric')
    return (rmse, area, float(den / w.sum()))

def evaluate(bundle, predroot, pred_anchor, output):
    start = time.perf_counter()
    from release_anchor import BENCHMARK_SHA256
    destination = Path(output).resolve()
    for prefix in bundle.lock['closed_directories']:
        protected = (bundle.root / prefix).resolve()
        if destination == protected or protected in destination.parents or destination in protected.parents:
            raise IntegrityError('output overlaps immutable input directory')
    for name in bundle.files:
        (base, rel) = bundle.location(name)
        p = (base / rel).resolve()
        if destination == p or destination in p.parents:
            raise IntegrityError('output would overwrite immutable input tree')
    protected = Path(predroot).resolve()
    if destination == protected or destination in protected.parents or protected in destination.parents:
        raise IntegrityError('output overlaps frozen predictions')
    case_files = list(bundle.tasks())
    cases = [Path(p).stem for (p, _) in case_files]
    manifest = predictions(predroot, pred_anchor, cases, BENCHMARK_SHA256)
    names = manifest['participants']
    policy = bundle.json('SCORER_PARAMETERS.json')
    rows = []
    groups = {}
    cohort = bundle.json('payload/private/COHORT.json')['payload']['cases']
    meta = {r['case_key']: r for r in cohort}
    totals = collections.defaultdict(collections.Counter)
    fc = collections.Counter()
    reject = collections.Counter()
    for (path, tasks) in case_files:
        case = Path(path).stem
        ref = bundle.npz('payload/private/references/' + case + '.npz')
        preds = {}
        for name in names:
            rel = name + '/' + case + '.npz'
            blob = safe_bytes(predroot, rel, 2000000)
            if h(blob) != manifest['files'][rel]['sha256']:
                raise IntegrityError('prediction raced after verification')
            a = load_npz(blob, max_expanded=2000000)
            status = decode(a, tasks)
            preds[name] = (a, status)
        cached = {}
        for (i, t0) in enumerate(tasks):
            t = t0
            if t['status'] == 'READY':
                scene = t['geometry_file']
                if scene not in cached:
                    cached[scene] = bundle.npz('payload/public/' + scene)
                t = attach(t, cached[scene])
                f = feasibility(t)
                f.pop('z', None)
            else:
                f = {'status': 'UNKNOWN_SITE'}
                reject[t.get('site_error', 'unknown')] += 1
            fc[f['status']] += 1
            for (name, (a, status)) in preds.items():
                row = dict(case_key=case, patient_group=meta[case]['group'], dataset=meta[case]['dataset'], split=t['split'], family=t['family'], level=t['level'], participant=name, feasibility=f['status'], height_rmse_mm=None, contact_area_error_mm2=None, reference_coverage=None, correct_abstention=False)
                if t['status'] != 'READY':
                    row['verdict'] = 'UNKNOWN_SITE'
                else:
                    d = design(a, i, t, status[i])
                    checks = all_checks(t, d)
                    row['verdict'] = checks['validity']
                    if row['verdict'] == 'PASS' and f['status'] == 'INFEASIBLE':
                        raise IntegrityError('accepted impossible design')
                    if status[i] == 'DESIGN':
                        (row['height_rmse_mm'], row['contact_area_error_mm2'], row['reference_coverage']) = metrics(t, a['outer_' + str(i)], ref[t['family']], policy['contact_band_mm'])
                    row['correct_abstention'] = row['verdict'] == 'ABSTAIN' and f['status'] == 'INFEASIBLE'
                row['pass'] = int(row['verdict'] == 'PASS')
                totals[name][row['verdict']] += 1
                rows.append(row)
                groups.setdefault((row['dataset'], row['split'], row['family'], name), {}).setdefault(case, {})[row['level']] = row
    summaries = []
    matrices = {}
    for (key, rr) in sorted(groups.items()):
        ordered = sorted(rr)
        a = np.asarray([[rr[c][level]['pass'] for level in LEVELS] for c in ordered])
        matrices[key] = a
        stats = cluster_summary(a, replicates=policy['bootstrap_replicates'])
        flat = [rr[c][l] for c in ordered for l in LEVELS]
        rmse = [r['height_rmse_mm'] for r in flat if r['height_rmse_mm'] is not None]
        summaries.append(dict(dataset=key[0], split=key[1], family=key[2], participant=key[3], **stats, counts=dict(collections.Counter((r['verdict'] for r in flat))), height_rmse_median_mm=float(np.median(rmse)) if rmse else None, rmse_scored_rows=len(rmse), rmse_missing_rows=len(flat) - len(rmse), correct_abstentions=sum((r['correct_abstention'] for r in flat))))
    contrasts = []
    for (ds, sp, fam, name) in matrices:
        if name != 'constraint_optimizer' or (ds, sp, fam, 'parametric') not in matrices:
            continue
        contrasts.append(dict(dataset=ds, split=sp, family=fam, **paired_case_difference(matrices[ds, sp, fam, name], matrices[ds, sp, fam, 'parametric'])))
    bundle.verify()
    predictions(predroot, pred_anchor, cases, BENCHMARK_SHA256)
    outcome = dict(claim_type='capability', benchmark_sha256=BENCHMARK_SHA256, prediction_freeze_sha256=pred_anchor, task_count=sum((len(ts) for (_, ts) in case_files)), case_count=len(cases), totals=dict(totals), feasibility_counts=dict(fc), summaries=summaries, paired_contrasts=contrasts, dropout=dict(task_site_reasons=dict(reject), task_site_rejected=sum(reject.values()), task_site_fraction=sum(reject.values()) / sum((len(ts) for (_, ts) in case_files)), source_records=bundle.json('raw/COHORT_SUMMARY.json')), external_referent=bundle.json('PREREG_R1.json')['external_referent'], scope=policy['scope'], seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
    (output / 'SCORED_ROWS.csv.gz').write_bytes(gzip.compress(buf.getvalue().encode(), mtime=0))
    dump(output / 'SCORES.json', outcome)
    semantic = {k: v for (k, v) in outcome.items() if k not in ['seconds', 'peak_rss_MiB']}
    dump(output / 'SCORE_DIGEST.json', dict(sha256=digest(semantic)))
    return outcome

def score(root=None, output=None):
    b = Bundle(root)
    from release_anchor import PREDICTIONS_SHA256
    return evaluate(b, b.payload / 'predictions', PREDICTIONS_SHA256, output or b.payload / 'evaluation')

def score_submission(root, prediction_directory, trusted_freeze_sha256, output):
    """Referee API: the receipt hash comes from trusted submission intake, never the entry."""
    b = Bundle(root)
    if not trusted_freeze_sha256:
        raise IntegrityError('referee freeze receipt required')
    return evaluate(b, Path(prediction_directory), trusted_freeze_sha256, output)
