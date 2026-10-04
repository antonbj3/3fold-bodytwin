"""Read-only, lineage-preserving statistical reanalysis of the pinned dental package.

All endpoint selection lives in adapters(), not in a heuristic numeric-key search.
Bootstrap intervals condition on delivered predictions; no refitting is implied.
"""
from dental_release.paths import expand as _release_expand
import argparse
import collections
import csv
import datetime
import hashlib
import gzip
import itertools
import json
import math
import os
import resource
import time
from pathlib import Path
for _var in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[_var] = '4'
import numpy as np
from scipy import optimize, stats
ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT.parent / 'DEMO48_PACKAGE'
SEED = 6103501
B = 10000
CONSUMED = {}
RAW = []
PIN = {}

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def read(p):
    p = Path(p)
    data = p.read_bytes()
    CONSUMED[str(p)] = dict(path=str(p), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
    if str(p) in PIN and CONSUMED[str(p)]['sha256'] != PIN[str(p)]:
        raise ValueError(f'Source hash drift: {p}')
    if p.name.endswith('.csv.gz'):
        return list(csv.DictReader(gzip.decompress(data).decode().splitlines()))
    if p.suffix == '.csv':
        return list(csv.DictReader(data.decode().splitlines()))
    if p.suffix == '.jsonl':
        return [json.loads(s) for s in data.decode().splitlines() if s.strip()]
    return json.loads(data)

def dump(p, d):

    def convert(x):
        if isinstance(x, np.generic):
            return x.item()
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, Path):
            return str(x)
        raise TypeError(type(x).__name__)
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    Path(p).write_text(json.dumps(d, indent=2, ensure_ascii=False, default=convert, allow_nan=False) + '\n')

def pointer(d, path):
    for k in path.strip('/').split('/'):
        d = d[int(k)] if isinstance(d, list) else d[k]
    return d

def ci(values):
    a = np.asarray(values, float)
    a = a[np.isfinite(a)]
    return [float(v) for v in np.quantile(a, [0.025, 0.975])] if len(a) else None

def mde_d(k):
    """Exact noncentral-t design MDE in independent-cluster SD units."""
    if k < 2:
        return None
    t = stats.t.ppf(0.975, k - 1)

    def power(d):
        nc = d * math.sqrt(k)
        return stats.nct.sf(t, k - 1, nc) + stats.nct.cdf(-t, k - 1, nc)
    upper = 1.0
    while power(upper) < 0.8:
        upper *= 2
    return float(optimize.brentq(lambda d: power(d) - 0.8, 1e-05, upper))

def endpoint(demo, label, clusters, values, unit, kind='mean', weights='rows', note='', primary=False, denominator=None):
    a = np.asarray(values, float)
    ids = np.asarray([str(x) for x in clusters])
    if len(a) != len(ids) or not len(a) or (not np.isfinite(a).all()):
        raise ValueError(f'{demo}/{label}: missing, unequal, or nonfinite raw records')
    if a.ndim == 1:
        a = a[:, None]
    (unique, inv) = np.unique(ids, return_inverse=True)
    (k, n) = (len(unique), len(a))
    count = np.bincount(inv)
    sums = np.array([np.bincount(inv, a[:, j]) for j in range(a.shape[1])]).T
    if kind == 'rmse':
        sums = np.array([np.bincount(inv, a[:, j] ** 2) for j in range(a.shape[1])]).T
    rng = np.random.default_rng(SEED + int(hashlib.sha256((demo + label).encode()).hexdigest()[:8], 16))

    def evaluate(s, m):
        if weights == 'equal_cluster':
            v = (s / m[..., None]).mean(axis=-2)
        else:
            v = s.sum(axis=-2) / m.sum(axis=-1)[..., None]
        if kind == 'rmse':
            return np.sqrt(v[..., 0])
        if kind == 'ratio':
            return np.divide(v[..., 0], v[..., 1], out=np.full_like(v[..., 0], np.nan), where=v[..., 1] != 0)
        return v[..., 0]
    if kind in ['median', 'median_ratio']:

        def med(x):
            return np.median(x[:, 0]) if kind == 'median' else np.median(x[:, 0]) / np.median(x[:, 1])
        point = float(med(a))
        iid = []
        block = []
        groups = [np.flatnonzero(inv == j) for j in range(k)]
        for _ in range(B):
            iid.append(med(a[rng.integers(n, size=n)]))
            chosen = rng.integers(k, size=k)
            block.append(med(a[np.concatenate([groups[j] for j in chosen])]))
        native_mde = None
    else:
        point = float(evaluate(sums, count))
        block = []
        iid = []
        raw_power = a ** 2 if kind == 'rmse' else a
        fixed_row_weights = 1 / count[inv].astype(float)
        for start in range(0, B, 100):
            j = rng.integers(k, size=(min(100, B - start), k))
            block.extend(evaluate(sums[j], count[j]).tolist())
            j = rng.integers(n, size=(min(100, B - start), n))
            if weights == 'equal_cluster':
                mass = fixed_row_weights[j]
                v = (raw_power[j] * mass[..., None]).sum(axis=1) / mass.sum(axis=1)[..., None]
                iid.extend((np.sqrt(v[:, 0]) if kind == 'rmse' else v[:, 0]).tolist())
            else:
                iid.extend(evaluate(raw_power[j], np.ones(j.shape)).tolist())
        native_mde = None
        if kind == 'mean' and k >= 2 and (weights == 'equal_cluster' or np.all(count == count[0])):
            sd = float(np.std(sums[:, 0] / count, ddof=1))
            if sd > 0:
                native_mde = mde_d(k) * sd
    before_point = point
    unweighted_point = float(np.sqrt(np.mean(a[:, 0] ** 2))) if kind == 'rmse' else float(a[:, 0].mean())
    kish = float(n * n / np.sum(count.astype(float) ** 2)) if weights == 'rows' else float(k)
    result = dict(demo=demo, metric=label, primary=primary, point=point, unit=unit, kind=kind, weighting=weights, n_records=n, k_clusters=k, cluster_size_min=int(min(count)), cluster_size_max=int(max(count)), kish_cluster_equivalent=kish, kish_definition='(sum m)^2/sum m^2; rho=1 concentration sensitivity, NOT estimated independent specimen n', naive_row_ci95=ci(iid), naive_row_point=before_point, cluster_ci95=ci(block), mde_80pct_alpha05_cluster_SD=mde_d(k), mde_native_approx=native_mde, mde_assumptions='Independent Gaussian cluster means; planning SD scenario from frozen records. For medians/RMSE/ratios native MDE UNKNOWN.', inference='PILOT_BOOTSTRAP_ONLY_K_LT10' if k < 10 else 'CONDITIONAL_DESCRIPTIVE', uncertainty_scope='Frozen predictions; conditional empirical sampling, excludes fitting, scanner/label bias, selection and external transport', note=note, denominator=denominator, invalid_bootstrap_draws=int(B - np.isfinite(block).sum()), resolution_level='POPULATION', input_resolution='PER_TOOTH')
    if weights == 'equal_cluster':
        result['naive_method'] = 'IID rows retain frozen inverse original within-block counts; normalize sampled weight mass. Same equal-cluster point estimand.'
        result['unweighted_row_point_sensitivity'] = unweighted_point
        result['note'] += ' R4 row-IID retains original study/case weight masses. Prior unweighted sensitivity preserved in RESULTS_R3_COMPLETE.json; it is not the headline estimand.'
    if kind == 'mean' and k >= 2:
        cm = sums[:, 0] / count
        if weights == 'equal_cluster' or np.all(count == count[0]):
            half = stats.t.ppf(0.975, k - 1) * np.std(cm, ddof=1) / math.sqrt(k)
            result['cluster_t95_Gaussian_sensitivity'] = [float(cm.mean() - half), float(cm.mean() + half)]
    RAW.append(dict(demo=demo, metric=label, cluster_ids=ids.tolist(), values=a.tolist(), unit=unit, kind=kind, weighting=weights))
    return result

def categorical(demo, label, rows, pred_key, ref_key, case_key='case_id', unit='κ'):
    classes = sorted({str(r[pred_key]) for r in rows} | {str(r[ref_key]) for r in rows})
    lut = {v: i for (i, v) in enumerate(classes)}
    l = len(classes)
    values = np.zeros((len(rows), l * l))
    for (j, r) in enumerate(rows):
        values[j, lut[str(r[ref_key])] * l + lut[str(r[pred_key])]] = 1
    ids = [str(r[case_key]) for r in rows]
    (u, inv) = np.unique(ids, return_inverse=True)
    k = len(u)
    sums = np.stack([np.bincount(inv, values[:, j], minlength=k) for j in range(l * l)], axis=1)

    def calc(v):
        m = v.reshape(v.shape[:-1] + (l, l))
        n = m.sum(axis=(-1, -2))
        pa = np.trace(m, axis1=-2, axis2=-1) / n
        pe = (m.sum(axis=-1) * m.sum(axis=-2)).sum(axis=-1) / n ** 2
        return np.divide(pa - pe, 1 - pe, out=np.full_like(pa, np.nan), where=np.abs(1 - pe) > 1e-14)
    rng = np.random.default_rng(SEED)
    bs = []
    na = []
    for start in range(0, B, 100):
        size = min(100, B - start)
        j = rng.integers(k, size=(size, k))
        bs.extend(calc(sums[j].sum(axis=1)))
        j = rng.integers(len(rows), size=(size, len(rows)))
        na.extend(calc(values[j].sum(axis=1)))
    count = np.bincount(inv)
    RAW.append(dict(demo=demo, metric=label, cluster_ids=ids, reference=[r[ref_key] for r in rows], prediction=[r[pred_key] for r in rows], kind='kappa'))
    return dict(demo=demo, metric=label, point=float(calc(sums.sum(axis=0))), unit=unit, primary=False, n_records=len(rows), k_clusters=k, kish_cluster_equivalent=float(len(rows) ** 2 / sum(count ** 2)), naive_row_ci95=ci(na), cluster_ci95=ci(bs), mde_80pct_alpha05_cluster_SD=mde_d(k), mde_native_approx=None, inference='CONDITIONAL_DESCRIPTIVE', resolution_level='POPULATION', note='Categorical report agreement, no millimeter or diagnosis accuracy. κ MDE in native units UNKNOWN; design multiplier applies to a scalar cluster outcome, not κ.')

def missing(demo, label, point, unit, reason, n=None, k=None, primary=True):
    return dict(demo=demo, metric=label, point=point, unit=unit, primary=primary, n_records=n, k_clusters=k, kish_cluster_equivalent=None, naive_row_ci95=None, cluster_ci95=None, mde_80pct_alpha05_cluster_SD=mde_d(k) if k and 'UNKNOWN' in reason else None, mde_native_approx=None, inference='UNKNOWN' if 'UNKNOWN' in reason else 'NOT_APPLICABLE', note=reason, resolution_level='PHENOMENOLOGICAL' if 'scenario' in reason.lower() else 'POPULATION')

def auc_endpoint(demo, truth, predictions):
    byid = {r['sample_id']: r['probability'] for r in predictions}
    rr = [r for r in truth if r['split'] == 'test']
    y = np.array([r['class'] != 'H0' for r in rr])
    score = np.array([byid[r['sample_id']] for r in rr])
    ids = [r['tid'].split('_')[0] for r in rr]
    (u, inv) = np.unique(ids, return_inverse=True)
    k = len(u)
    n = len(rr)
    (vals, ties) = np.unique(score, return_inverse=True)
    v = len(vals)
    pos = np.stack([np.bincount(ties[(inv == j) & y], minlength=v) for j in range(k)])
    neg = np.stack([np.bincount(ties[(inv == j) & ~y], minlength=v) for j in range(k)])

    def calc(p, q):
        nc = q.cumsum(axis=-1) - 0.5 * q
        return (p * nc).sum(axis=-1) / (p.sum(axis=-1) * q.sum(axis=-1))
    rng = np.random.default_rng(SEED)
    bs = []
    iid = []
    for start in range(0, B, 50):
        size = min(50, B - start)
        j = rng.integers(k, size=(size, k))
        bs.extend(calc(pos[j].sum(axis=1), neg[j].sum(axis=1)))
        j = rng.integers(n, size=(size, n))
        pp = np.stack([np.bincount(ties[row][y[row]], minlength=v) for row in j])
        qq = np.stack([np.bincount(ties[row][~y[row]], minlength=v) for row in j])
        iid.extend(calc(pp, qq))
    count = np.bincount(inv)
    RAW.append(dict(demo=demo, metric='test AUC', cluster_ids=ids, values=score.tolist(), outcome=y.astype(int).tolist(), kind='AUC'))
    return dict(demo=demo, metric='test AUC', point=float(calc(pos.sum(axis=0), neg.sum(axis=0))), unit='AUC', primary=True, n_records=n, k_clusters=k, kish_cluster_equivalent=float(n * n / sum(count ** 2)), naive_row_ci95=ci(iid), cluster_ci95=ci(bs), mde_80pct_alpha05_cluster_SD=mde_d(k), mde_native_approx=None, inference='PILOT_BOOTSTRAP_ONLY_K_LT10', resolution_level='POPULATION', note='Four held-out STS case prefixes, one tooth per case; repeated synthetic lesions/views/noise are not independent patients. Synthetic truth is not external physical caries validation.')

def adapters(inventory, facit):
    output = []
    demos = []
    controls = []
    for item in inventory:
        id = item['id']
        p = Path(item['path'])
        d = read(p / 'results.json')
        f = facit[id]
        metrics = []
        c = 0
        primary = f.get('primary_pointer')
        primary_val = f.get('primary_value')
        units = f.get('units', '')
        if primary:
            actual = pointer(d, primary)
            tol = max(1e-09, 1e-08 * abs(float(primary_val)))
            ok = abs(float(actual) - float(primary_val)) <= tol
            bad = float(actual) + max(10, abs(float(actual)) * 0.5)
            controls.append(dict(demo=id, kind='Frozen FACIT pointer', pointer=primary, expected=primary_val, observed=actual, tolerance=tol, pass_=ok, injected_value=bad, injected_rejected=abs(bad - float(primary_val)) > tol))
            if not ok:
                raise ValueError(f'{id}: package FACIT drift')
        if id == 'X7':
            predictions = read(p / 'raw/PREDICTIONS_R1.json')
            lab = read(p / 'raw/TEST_LABELS.json')
            for field in ['overbite', 'overjet', 'crossbite', 'midlines', 'spee', 'molar_right', 'molar_left', 'canine_right', 'canine_left']:
                rows = []
                for r in predictions:
                    if r['split'] != 'test' or r.get('error') or r['case_id'] not in lab:
                        continue
                    reports = lab[r['case_id']]
                    y = reports[0]['labels'].get(field) if reports else None
                    if y is None:
                        continue
                    rows.append(dict(case_id=r['case_id'], y=y, pred=r['fields'][field]['point'], rule=r['fields'][field]['threshold']))
                metrics.append(endpoint(id, field + ' accuracy', [r['case_id'] for r in rows], [r['y'] == r['pred'] for r in rows], 'andel', primary=field == 'overbite', note='One field value per case; no extra clustering penalty within this endpoint. Reused internal test set.'))
                metrics.append(categorical(id, field + ' kappa', rows, 'pred', 'y'))
                metrics.append(endpoint(id, field + ' paired accuracy gain vs rule', [r['case_id'] for r in rows], [int(r['y'] == r['pred']) - int(r['y'] == r['rule']) for r in rows], 'andel'))
            disagreement = read(p / 'raw/REPORT_DISAGREEMENT.json')
            for (field, q) in disagreement.items():
                metrics.append(categorical(id, field + ' report-report kappa', q['pairs'], 'second', 'first'))
        elif id == 'X21':
            rows = read(p / 'raw/COMPARISON_ROWS_R2.json')
            for field in sorted(set((r['field'] for r in rows))):
                rr = [r for r in rows if r['field'] == field]
                metrics.append(endpoint(id, field + ' first-report accuracy', [r['case_id'] for r in rr], [r['prediction'] == r['reference_first'] for r in rr], 'andel', primary=field == 'overbite', note='Same case cohort reused from X7; categorical report endpoint.'))
                metrics.append(categorical(id, field + ' kappa', rr, 'prediction', 'reference_first'))
                metrics.append(endpoint(id, field + ' paired accuracy gain vs X7', [r['case_id'] for r in rr], [int(r['prediction'] == r['reference_first']) - int(r['x7_prediction'] == r['reference_first']) for r in rr], 'andel'))
        elif id == 'X3':
            rows = read(p / 'PER_TOOTH_R3.csv')
            rr = [r for r in rows if r['method'] == 'collar_sdf']
            metrics.append(endpoint(id, 'R3 upper-surface p95 mean', [r['case'] for r in rr], [float(r['occlusal_p95_mm']) for r in rr], 'mm', weights='equal_cluster', primary=True, note='32 teeth, four new test cases. Exact original equal-patient estimand; absolute .50 mm gate FAIL.'))
            maps = {(r['case'], r['jaw'], r['fdi'], r['method']): float(r['occlusal_p95_mm']) for r in rows}
            for ctrl in ['mirror', 'collar_mesh']:
                metrics.append(endpoint(id, 'R3 paired p95 difference vs ' + ctrl, [r['case'] for r in rr], [float(r['occlusal_p95_mm']) - maps[r['case'], r['jaw'], r['fdi'], ctrl] for r in rr], 'mm', weights='equal_cluster'))
        elif id == 'X13':
            rows = read(p / 'FROZEN_PREDICTIONS_R3.json')['predictions']
            lookup = {r['row_id']: r for r in read(p / 'measurements.json')}
            er = [r['prediction_um'] - lookup[r['row_id']]['measured_mean_um'] for r in rows]
            metrics.append(endpoint(id, 'R3 study-balanced RMSE', [r['study'] for r in rows], er, 'µm', kind='rmse', weights='equal_cluster', primary=True, note='Three selected held-out aggregate cells in two studies; cannot estimate specimen prediction error or fit uncertainty.'))
        elif id == 'X5':
            metrics.append(missing(id, 'constructed decision counterexamples', 6, 'motprov', 'NOT_APPLICABLE: six selected counterexamples establish existence, not a sampled error/prevalence rate. No meaningful population CI or MDE for this count.', 6, 6))
        elif id == 'X14':
            rows = read(p / 'raw/R1_test.json')
            metrics.append(missing(id, 'R1 angle RMSE', primary_val, units, 'UNKNOWN: six table aggregates, one study, no angle SD/raw specimen values; two manufacturers are fixed conditions, not six independent samples.', 6, 1))
        elif id == 'GENCAD_V2':
            rows = read(p / 'raw/SCORED_ROWS_R1.json')
            for method in sorted(set((r['participant'] for r in rows))):
                rr = [r for r in rows if r['participant'] == method and r['split'] == 'test']
                values = [int(r.get('checks', {}).get('validity') == 'PASS') for r in rr]
                metrics.append(endpoint(id, 'test digital PASS ' + method, [r['case_key'] for r in rr], values, 'andel', primary=method == 'constraint_optimizer', note='All assigned test tasks retained; missing-site, abstention and UNKNOWN are not PASS. 576 all-split tasks != independent test n. Cross-dataset person identity unresolved.'))
        elif id == 'PROOF_LANE':
            tasks = read(p / 'data/tasks.json')
            tm = {r['task_id']: r for r in tasks}
            rows = d['rows']
            metrics.append(missing(id, 'submitted designs', primary_val, units, 'NOT_APPLICABLE: exact software workload count; not 72 independent anatomical observations.', 72, 2))
            for method in d['summary']:
                rr = [r for r in rows if r['generator'] == method]
                if not rr:
                    continue
                metrics.append(endpoint(id, 'full digital PASS ' + method, [tm[r['task_id']]['subject_group'] for r in rr], [r['level1']['status'] == 'PASS' for r in rr], 'andel', note='Only two held-out subjects; antagonist input is absent, so complete design remains UNKNOWN.'))
        elif id == 'GENCAD_V3':
            rows = [r for r in read(p / 'payload/evaluation/SCORED_ROWS.csv.gz') if r['dataset'] == 'Bite2Text' and r['split'] == 'test']
            summaries = {(r['family'], r['participant']): r for r in d['surface_benchmark']['summaries'] if r['dataset'] == 'Bite2Text' and r['split'] == 'test'}
            for (family, method) in sorted(summaries):
                rr = [r for r in rows if r['family'] == family and r['participant'] == method]
                e = endpoint(id, family + ' digital PASS ' + method, [r['patient_group'] for r in rr], [r['verdict'] == 'PASS' for r in rr], 'andel', weights='equal_cluster', primary=method == 'constraint_optimizer', note='Four dependent levels per patient-case. Native scan surface facit, not prepared/milled/seated geometry. Within-dataset588 cases; cross-dataset patient identity UNKNOWN.')
                original = summaries[family, method]
                e['original_reported_ci95'] = original['case_cluster_bootstrap_95']
                e['original_interval_locator'] = '/surface_benchmark/summaries/' + family + '/' + method
                e['original_reported_n_eff_tasks_ICC'] = original['n_eff_tasks_estimated']
                if e['k_clusters'] != 588:
                    raise ValueError('V3 patient-group count is not588')
                expected = original['mean_case_pass_fraction']
                tol = 1e-08
                controls.append(dict(demo=id, kind=family + '/' + method + '/original patient fraction', expected=expected, observed=e['point'], tolerance=tol, pass_=abs(expected - e['point']) < tol, injected_value=e['point'] + 0.1, injected_rejected=abs(expected - e['point'] - 0.1) > tol))
                metrics.append(e)
        elif id == 'X33':
            for name in ['low_peak_C', 'high_peak_C', 'uniform_peak_C']:
                e = missing(id, 'same-heat ' + name, d['R2'][name], '°C', 'NOT_APPLICABLE: deterministic pulsed-model maximum; no new independent specimen, measurement uncertainty or certified continuum error.', 0, 0)
                e['resolution_level'] = 'PER_POINT'
                metrics.append(e)
            metrics.append(missing(id, 'external proxy-transfer MAE', d['external_validation']['proxy_transfer_MAE_C'], '°C', 'UNKNOWN: three held protocol-group means from one source study, missing matched pulp/sensor quantity and specimen lineage.', 3, 1))
            e = missing(id, 'published sensor maximum group mean', 0.84, '°C', 'UNKNOWN: reported SD0.55C is specimen dispersion, not uncertainty for the local all-pulp-field maximum. Independent specimen n/batch and observation mapping unavailable.', None, 1)
            e['reported_group_SD_C'] = 0.55
            metrics.append(e)
        elif id == 'X4':
            rows = [r for r in read(p / 'table_cases.csv') if r['round'] == 'R4' and r['method'] == 'dp' and r['primary_p95_mm']]
            metrics.append(endpoint(id, 'R4 diagnostic median p95', [r['case'] for r in rows], [float(r['primary_p95_mm']) for r in rows], 'mm', kind='median', primary=True, note='Published virtual target; original geometry and fit gates FAIL. Sampling CI cannot rescue target mismatch.'))
        elif id == 'X8':
            rows = [r for r in read(p / 'PER_SITE_FULL_GEOMETRY_RISK.csv') if r['guide'] == 'fully_guided']
            metrics.append(endpoint(id, 'whole-body vs old gap decision change', [r['case'] for r in rows], [int((float(r['old_K3_gap_mm']) >= 2) != (float(r['whole_cylinder_label_gap_lower_mm']) >= 2)) for r in rows], 'andel', primary=True, note='Virtual site decisions, no observed operation or nerve outcome.'))
        elif id == 'X9':
            allrows = read(p / 'artifacts/thickness_by_tooth_side_height.csv')
            rows = [r for r in allrows if float(r['height_fraction']) == 0.5 and r['mean_mm']]
            metrics.append(endpoint(id, 'R1 middle-height median thickness', [r['tid'].split('_')[0] for r in rows], [float(r['mean_mm']) for r in rows], 'mm', kind='median', primary=True, note='Case prefix from STS; mean ray thickness per tooth-side, then median. Not safe IPR.'))
        elif id == 'X12':
            metrics.append(missing(id, 'nontruncated annotation tooth rows', primary_val, units, 'NOT_APPLICABLE: exact annotation census, not an uncertain patient parameter.', 3449, 282))
            rows = [r for r in read(p / 'raw/R8_teeth.csv') if r['domain_truncated'] == 'False']
            metrics.append(endpoint(id, 'hard tissue horn minimum median', [r['case'] for r in rows], [float(r['horn_boundary_min_mm']) for r in rows], 'mm', kind='median', note='Whole hard tissue against annotated pulp; not dentin, scanner trueness or a restoration thickness.'))
        elif id == 'X15':
            rows = [r for r in read(p / 'RAW_R4.jsonl') if r['R3_query']['conditional_geometric_status'] == 'WITHIN_SECTION_CONSTRAINTS']
            if rows:
                metrics.append(endpoint(id, 'R4 previous admission survival', [r['case'] for r in rows], [r['R3_admission_survived'] for r in rows], 'andel', primary=True, note='Repeated motion queries of only six scans; label-surface check, no biological endpoint.'))
        elif id == 'X17':
            metrics.append(missing(id, 'area absolute-error ratio', primary_val, units, 'UNKNOWN: derived from two population group means in one external study; no joint AP/lateral/area/volume covariance or individual prediction errors.', 2, 1))
        elif id == 'X18':
            rows = [r for r in read(p / 'raw/RESULTS_R2_ROWS.json') if r['status'] == 'SCORED']
            a = []
            for r in rows:
                q = r['arms']['0.1']
                a.append([q['informed']['errors']['absolute_area_error_mm2'], q['practice']['errors']['absolute_area_error_mm2']])
            metrics.append(endpoint(id, 'R2 ratio median area errors', [r['case'] for r in rows], a, 'ratio', kind='median_ratio', primary=True, note='Ratio of TWO medians, not median of individual ratios or ratio of means; fixed annotated tooth surface.'))
        elif id == 'X19':
            metrics.append(missing(id, 'formed relative force error', primary_val, units, 'UNKNOWN: three external tooth-type means; same setup, no pairing with simulated case forces. Specimen SD is not SD of case-wise prediction error.', 3, 1))
        elif id == 'X22':
            truth = read(_release_expand('@DENTAL_WORK_ROOT@/X22-caries-synthetic/private/truth_index.json'))
            metrics.append(auc_endpoint(id, truth, read(p / 'raw/R1_reference_logistic_predictions.json')))
        elif id == 'X20':
            metrics.append(missing(id, 'clinical UNKNOWN queries', primary_val, units, 'NOT_APPLICABLE: exact count of declared query answers; selected locations do not identify clinical prevalence.', 28, None))
        elif id == 'X26':
            rows = read(p / 'RAW_R2_CANAL.json')
            metrics.append(endpoint(id, 'conditional class 1 fraction', [r['case'] for r in rows], [r['class'] == 1 for r in rows], 'andel', primary=True, note='Two dependent selection policies per site. Numerical feasibility under assumed Gaussian error floor; physical certificates remain zero.'))
        elif id == 'X18B':
            metrics.append(missing(id, 'maximum source surface export error', 0, 'mm', 'NOT_APPLICABLE: deterministic identity of reused measured input triangles; zero is no accuracy interval for new preparations.', 24, None))
        elif id == 'X30':
            rows = read(p / 'raw/R2_topology.csv')
            metrics.append(endpoint(id, 'one-voxel topology stable fraction', [r['case'] for r in rows], [r['scenario_stable'] == 'True' for r in rows], 'andel', primary=True, note='3449 teeth from 282 CT cases; digital perturbation stability, not an anatomical CI.'))
        elif id == 'X24':
            rows = [r for r in d['rounds']['R2']['rows'] if r['variant'] == 'quadratic3_candidate' and r['device'] != 'Varian']
            metrics.append(missing(id, 'maximum finite held-out calibration error', max((abs(r['error_HU']) for r in rows)), 'HU_ref', 'NOT_APPLICABLE: finite-set maximum, 48 regions on three test scanners plus one reference device; no upper population error guarantee.', 48, 3))
            metrics.append(endpoint(id, 'calibration mean absolute held-out error', [r['device'] for r in rows], [abs(r['error_HU']) for r in rows], 'HU_ref', primary=False, note='Three fixed test scanner devices, separate Varian reference; cluster at device. No randomized scanner population sample.'))
        elif id == 'X31':
            rows = read(p / 'raw/GUIDE_PER_SITE.csv')
            for guide in sorted(set((r['guide'] for r in rows))):
                rr = [r for r in rows if r['guide'] == guide]
                metrics.append(endpoint(id, guide + ' apex-scenario decision change', [r['case'] for r in rr], [r['two_mm_accept'] != r['local_apex_scenario_accept'] for r in rr], 'andel', primary=guide == 'fully_guided', note='Uncalibrated scenario; case resampling does not quantify guide-parameter uncertainty or clinical clearance.'))
            rows = [r for r in read(p / 'raw/REDUCTION_PER_TOOTH.csv') if r['tooth_type'] == 'molar3' and float(r['reduction_mm']) == 2 and (float(r['annotation_sensitivity_per_surface_mm']) == 0)]
            if rows:
                metrics.append(endpoint(id, 'molar3 total hard tissue remaining median', [r['case'] for r in rows], [float(r['hard_tissue_proxy_remaining_mm']) for r in rows], 'mm', kind='median', note='Enamel+dentin proxy; true dentin partially identified, not measured.'))
        elif id == 'X32':
            rows = read(p / 'raw/decision_ledger.csv')
            for did in ['X8', 'X5', 'X13']:
                rr = [r for r in rows if r['demo'] == did]
                metrics.append(endpoint(id, did + ' changed decisions per100', [r['cluster_id'] for r in rr], [100 * (r['changed'] == 'True') for r in rr], 'per100', primary=did == 'X8', note='Same source rows as X8/X5/X13, no additional independent evidence. X13 study families, no patient measurements.'))
                events = {c: any((r['changed'] == 'True' for r in rr if r['cluster_id'] == c)) for c in sorted(set((r['cluster_id'] for r in rr)))}
                metrics.append(endpoint(id, did + ' any-change clusters per100', list(events), [100 * int(v) for v in events.values()], 'per100', note='One case/family event per source cluster, distinct from per-site rate; no additional source cohort.'))
        elif id == 'X23':
            rows = read(p / 'raw/HELDOUT_SCORES.json')
            metrics.append(endpoint(id, 'held-out process q05 loss', [r['group'] for r in rows], [r['losses']['process_q0.05'] for r in rows], 'MPa', primary=True, note='Cluster by process group: shared process/batch uncertainty; only one source study. Individual bar bootstrap alone conditions on process family. CV fits overlap.'))
            metrics.append(endpoint(id, 'q05 process-minus-nominal loss', [r['group'] for r in rows], [r['losses']['process_q0.05'] - r['losses']['nominal_q0.05'] for r in rows], 'MPa', note='Conditional comparison within one study; process groups do not replicate source studies.'))
        elif id == 'X25':
            rows = read(p / 'raw/specimens.json')
            rr = [r for r in rows if r['configuration'] == 'MUSLA 04020' and r['max_load_N'] == 100 and (r['cycles'] >= 5000000)]
            n = len(rr)
            fail = sum((r['failure'] for r in rr))
            upper = 1 - 0.05 ** (1 / n) if n and fail == 0 else None
            metrics.append(missing(id, '100N complete-horizon failures', fail / n if n else None, 'andel', 'UNKNOWN: specimen/binomial inference conditional on independent production; batch/donor lineage absent. Exact one-sided 95% upper risk=' + str(upper), n, None))
            metrics[-1]['independent_specimen_scenario_n'] = n
            metrics[-1]['one_sided95_binomial_upper_if_IID'] = upper
        elif id == 'X16':
            rr = read(p / 'raw/R1_RESULTS.json')['folds']
            metrics.append(endpoint(id, 'LOSO study-family balanced absolute temperature error', [r['study'] for r in rr], [r['MAE_C'] for r in rr], '°C', weights='rows', primary=True, note='Seven study-family folds share only FIVE source studies. Preserve equal-fold original point while bootstrapping whole source studies; future individual tooth error and refit variance UNKNOWN.'))
        elif id == 'X1B':
            q = read(p / 'raw/MATCHED_R2.json')
            rr = [dict(study=f['held_out_study'], **r) for f in q['folds'] for r in f['rows']]
            metrics.append(missing(id, 'prospective M2 double ratio Q', primary_val, units, 'NOT_APPLICABLE: frozen deterministic FE prediction, no manufactured crown response. Rejection window [0.351,0.640] is NOT a confidence interval.', 0, 0))
            metrics.append(endpoint(id, 'R2 study-balanced log ratio RMSE', [r['study'] for r in rr], [r['log_error'] for r in rr], 'log ratio', kind='rmse', weights='equal_cluster', note='Retrospective source-held-out aggregate study ratio; not independent manufactured crown validation.'))
        elif id == 'X29':
            rr = read(p / 'raw/R3_MEASUREMENTS.json')
            r = next((r for r in rr if r['id'] == 'PMC10788321_3Y_angle30'))
            n = int(r['n'])
            v = float(r['mean_N'])
            sd = float(r['sd_N'])
            half = stats.t.ppf(0.975, n - 1) * sd / math.sqrt(n)
            e = missing(id, 'published aged inclined bridge mean fracture load', v, 'N', 'UNKNOWN: published n/mean/SD allow a t interval only under independent specimens. Batch IDs and specimen raw values unavailable.', n, 1)
            e['reported_specimen_t95_if_IID'] = [float(v - half), float(v + half)]
            e['mde_native_if_IID_specimens'] = mde_d(n) * sd
            metrics.append(e)
        elif id == 'X28':
            rr = read(p / 'raw/canine_group_measurements.json')
            metrics.append(missing(id, '500-minus-1000rpm viability', 37.59, 'pp', 'UNKNOWN: four reported group samples; biological donor versus aliquot/technical repeat lineage unavailable, plus cross-group covariance. No valid donor-cluster bootstrap.', 8, None))
        elif id == 'X27':
            metrics.append(missing(id, 'whole-arch neighbor-force groups passing', 2 / 6, 'andel', 'UNKNOWN: six adjacent tooth means share one setup/calibration, not six independent experiments; raw repeat/donor identifiers missing.', 6, 1))
        elif id == 'X10':
            metrics.append(missing(id, 'retrospective positional median reduction', primary_val, units, 'UNKNOWN: ratio of separately published medians, no paired scanner-run values/covariance. 28 scanner×scanbody summary cells are not specimen n.', 28, 1))
        elif id == 'X2':
            metrics.append(missing(id, 'Ferrato mechanics regional MAE', primary_val, units, 'UNKNOWN: regional population force means vs unpaired geometry cohort; published force SD cannot identify prediction-error variance or matched-case uncertainty.', 200, None))
        if not metrics:
            metrics.append(missing(id, 'package headline', primary_val, units, 'UNKNOWN: no audited independent raw-record adapter; exact missing record mapping retained.'))
        if primary and id in ['X7', 'X3', 'X13', 'X4', 'X8', 'X9', 'X15', 'X18', 'X21', 'X22']:
            em = next((x for x in metrics if x.get('primary') and x['inference'] not in ['UNKNOWN', 'NOT_APPLICABLE']), None)
            if em:
                tol = max(1e-09, 1e-08 * abs(float(primary_val)))
                ok = abs(em['point'] - primary_val) <= tol
                controls.append(dict(demo=id, kind='Raw-record headline reconstruction', expected=primary_val, observed=em['point'], tolerance=tol, pass_=ok, injected_rejected=True))
                if not ok:
                    raise ValueError(f"{id}: RAW reconstruction mismatch {em['point']} vs {primary_val}")
        demos.append(dict(id=id, path=str(p), source=item['source'], title=item['title'], metrics=metrics, external_reference=f, claim_type=d.get('claim_type', 'capability'), support='SEE_METRIC_SCOPE'))
        output.extend(metrics)
        dump(ROOT / 'R1_PARTIAL.json', dict(demos=demos, metrics=output, controls=controls, completed_demos=len(demos)))
        dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X35-package-statistics', phase='R1_CLUSTER_BOOTSTRAP_RUNNING', latest_demo=id, completed_demos=len(demos), last_gate='Raw headline and FACIT checks passed for completed demos', next_operation='Continue remaining demo adapters, then small-cluster exact sensitivity', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        print(id, 'complete', len(metrics), flush=True)
    return (demos, output, controls)

def main():
    global PIN
    ap = argparse.ArgumentParser()
    ap.add_argument('--prepare', action='store_true')
    args = ap.parse_args()
    pinpath = ROOT / ('SOURCE_PIN_FINAL.json' if (ROOT / 'SOURCE_PIN_FINAL.json').exists() else 'SOURCE_PIN.json')
    if pinpath.exists():
        PIN = {r['path']: r['sha256'] for r in json.loads(pinpath.read_text())}
    started = time.perf_counter()
    inventory = read(ROOT / 'INPUT_INVENTORY_R3.json')
    facit = read(ROOT / 'FACIT_SNAPSHOT_R3.json')
    (demos, metrics, controls) = adapters(inventory, facit)
    import report_metadata
    metrics = report_metadata.enrich(demos, read, missing, mde_d)
    dump(ROOT / 'raw/ENDPOINT_RECORDS.json', RAW)
    dump(ROOT / 'INPUT_LOCK.json', list(CONSUMED.values()))
    if not pinpath.exists():
        dump(pinpath, list(CONSUMED.values()))
    result = dict(lane='X35-package-statistics', claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), bootstrap_draws=B, seed=SEED, demos=demos, metrics=metrics, n_demos=len(demos), n_metrics=len(metrics), controls=controls, runtime_s=time.perf_counter() - started, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    import round2, checklists, verify, render_report
    result['round2'] = round2.run(RAW)
    result['validation'] = verify.run(result, RAW, result['round2'])
    result['checklists'] = checklists.build(demos)
    render_report.build(result, result['checklists'])
    result['runtime_s'] = time.perf_counter() - started
    result['peak_RSS_MiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    result['raw_data_manifest'] = [dict(path=str(ROOT / 'raw/ENDPOINT_RECORDS.json'), bytes=(ROOT / 'raw/ENDPOINT_RECORDS.json').stat().st_size, sha256=sha(ROOT / 'raw/ENDPOINT_RECORDS.json'))]
    result['protocol_hashes'] = {name: sha(ROOT / name) for name in ['PREREG_R1.json', 'PREREG_R2.json', 'PREREG_R3.json', 'PREREG_R4.json', 'FROZEN_PREDICTIONS.json', 'DECOMPOSITION.json']}
    result['package_snapshot'] = dict(utc=json.loads((ROOT / 'PREREG_R3.json').read_text())['created_utc'], active_demos=34, retained_including_archived=35, archived=['GENCAD_V2'], facit_sha256=sha(ROOT / 'FACIT_SNAPSHOT_R3.json'), dynamic_package_changes='Not silently added; refresh requires a new versioned scope/prereg')
    result['full_cost'] = dict(preparation='Targeted source/lineage/DOI audit; human-equivalent research time not benchmarked', fit=0, discovery='Retrospective endpoint/schema repair; original outcomes known', validation_controls=result['validation']['n_controls'], queries='Official reporting documents browsed;0 new dataset/lab queries', fallback='UNKNOWN when independent specimens/batches or matched quantities absent', compute_wall_s=result['runtime_s'], resource_scope='CPU <=4threads; no GPU; RSS measured; no arrays>50MB')
    dump(ROOT / 'results.json', result)
    if not (ROOT / 'SOURCE_PIN_FINAL.json').exists():
        dump(ROOT / 'SOURCE_PIN_FINAL.json', list(CONSUMED.values()))
    feedback_path = ROOT / 'GRAPH_FEEDBACK.json'
    if feedback_path.exists():
        feedback = json.loads(feedback_path.read_text())
        feedback.update(result_sha256=sha(ROOT / 'results.json'), review_state='PENDING_INDEPENDENT_REVIEW', status='PENDING_INDEPENDENT_REVIEW')
        dump(feedback_path, feedback)
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X35-package-statistics', phase='AUDIT_R1_R4_COMPLETE_PENDING_INDEPENDENT_REVIEW', last_gate=result['validation']['status'], n_demos=len(demos), n_endpoints=len(metrics), package_active=34, package_including_archived=35, latest_result_sha256=sha(ROOT / 'results.json'), next_operation='Independent review and a new case/patient/batch measurement register', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    print(json.dumps(dict(demos=len(demos), endpoints=len(metrics), controls=result['validation']['n_controls'], wall_s=time.perf_counter() - started)))
if __name__ == '__main__':
    main()
