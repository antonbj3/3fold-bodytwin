"""Separate evaluator. Never imports or executes a participant."""
import time, collections
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .common import *
from .tasks import load_task
from .checks import all_checks, validate, feasibility
from .geometry import export_stl, shell

def metric(z, inner, t, ref):
    truth = ref['reference']
    w = np.asarray(t['weights'])
    valid = np.isfinite(truth)
    den = float(w[valid].sum())
    if den <= 0:
        return dict(status='UNKNOWN_NO_REFERENCE')
    error = z[valid] - truth[valid]
    out = dict(status='SCORED', reference_coverage=float(w[valid].sum() / w.sum()), height_rmse_mm=float(np.sqrt(np.sum(w[valid] * error ** 2) / den)), height_p95_mm=float(np.quantile(abs(error), 0.95)), nominal_removed_volume_mm3=float(np.sum(w[valid] * np.maximum(0, truth[valid] - np.asarray(t['preparation_z'])[valid]))), restoration_volume_mm3=float(np.sum(w * (z - inner))), nominal_film_mae_from_80um=float(np.mean(abs(inner - np.asarray(t['preparation_z']) - 0.08)) * 1000))
    ceiling = np.asarray(t['ceiling'])
    contacts = {}
    ok = valid & np.isfinite(ceiling)
    for band in [0.05, 0.1, 0.2]:
        if not ok.any():
            contacts[str(band)] = dict(status='UNKNOWN_NO_ANTAGONIST')
            continue
        native_gap = ceiling - truth
        gap = ceiling - z
        gold = ok & (native_gap >= 0) & (native_gap <= band)
        pred = ok & (gap >= 0) & (gap <= band)
        union = float(w[gold | pred].sum())
        inter = float(w[gold & pred].sum())
        contacts[str(band)] = dict(status='SCORED', reference_area_mm2=float(w[gold].sum()), predicted_area_mm2=float(w[pred].sum()), iou=inter / union if union > 0 else None, empty_both=union == 0, area_abs_error_mm2=abs(float(w[pred].sum() - w[gold].sum())), native_penetration_sample_area_mm2=float(w[ok & (native_gap < 0)].sum()), pose_area_mm2={str(d): float(w[ok & (gap + d >= 0) & (gap + d <= band)].sum()) for d in [-0.05, 0.05]})
    out['contacts'] = contacts
    neighbors = np.asarray(t['neighbor_points'])
    pts = np.c_[t['xy'], z]
    out['approximal_sample_min_mm'] = float(cKDTree(neighbors).query(pts, workers=1)[0].min()) if len(neighbors) else None
    out['approximal_scope'] = 'nearest sampled other-tooth surface point; not signed collision or actual contact certificate'
    return out

def pareto(rows):
    """Within-task partial order; abstentions and invalid tasks cannot dominate designs."""
    vectors = {}
    for r in rows:
        if r.get('checks', {}).get('validity') != 'PASS' or r.get('metrics', {}).get('status') != 'SCORED':
            continue
        m = r['metrics']
        co = m['contacts'].get('0.1', {})
        e = co.get('area_abs_error_mm2')
        if e is None:
            continue
        vectors[r['participant']] = np.array([m['height_rmse_mm'], e, m['nominal_film_mae_from_80um'], m['restoration_volume_mm3']])
    return {name: not any((np.all(b <= a + 1e-09) and np.any(b < a - 1e-09) for (other, b) in vectors.items() if other != name)) for (name, a) in vectors.items()}

def score(round_name='R1', participants=None):
    start = time.perf_counter()
    frozen = read(ROOT / f'FROZEN_PREDICTIONS_{round_name}.json')
    payload = frozen['payload']
    if digest(payload) != frozen['payload_sha256']:
        raise ValueError('prediction manifest altered')
    for (rel, h) in payload['files'].items():
        p = DATA / rel
        if p.is_symlink() or sha(p) != h:
            raise ValueError('prediction drift ' + rel)
    b = read(ROOT / 'data/BENCHMARK_LOCK.json')
    ref_lookup = {r['task_id']: r for r in read(DATA / 'private/references.json')}
    rows = []
    tasks = []
    exported = []
    participants = participants or payload['participants']
    for path in sorted((DATA / 'public').glob('*.json')):
        t = load_task(path)
        tid = t['task_id']
        taskrow = {k: t[k] for k in ['task_id', 'case_key', 'split', 'family', 'level', 'status']}
        group = []
        if t['status'] == 'READY':
            q = ref_lookup[tid]['reference']
            rp = DATA / 'private' / q['file']
            if sha(rp) != q['sha256']:
                raise ValueError('reference hash drift')
            with np.load(rp, allow_pickle=False) as a:
                ref = {k: a[k] for k in a.files}
            f = feasibility(t)
            f.pop('z', None)
            taskrow['feasibility'] = f
        else:
            taskrow['feasibility'] = dict(status='UNKNOWN_SITE')
        tasks.append(taskrow)
        for name in participants:
            r = dict(**taskrow, participant=name)
            p = DATA / payload['prediction_dirs'][name] / (tid + '.json')
            if not p.exists():
                r.update(checks=dict(validity='MISSING', digital_score=0), metrics=dict(status='UNKNOWN_MISSING_SUBMISSION'))
                group.append(r)
                continue
            d = read(p)
            if t['status'] != 'READY':
                r.update(checks=dict(validity='UNKNOWN_SITE', digital_score=None), metrics=dict(status='UNKNOWN_SITE'))
                group.append(r)
                continue
            check = all_checks(t, d)
            r['checks'] = check
            if check['validity'] not in ['INVALID', 'ABSTAIN']:
                (z, inner) = validate(t, d)
                r['metrics'] = metric(z, inner, t, ref)
                key = (t['family'], name)
                if t['level'] == 'normal' and key not in exported:
                    dest = DATA / 'exports' / name
                    dest.mkdir(parents=True, exist_ok=True)
                    (v, faces) = shell(np.asarray(t['xy']), z, inner, np.asarray(t['faces']))
                    export_stl(dest / (t['family'] + '.stl'), v, faces)
                    exported.append(key)
            else:
                r['metrics'] = dict(status='UNKNOWN_' + check['validity'])
            r['abstention_correct'] = taskrow['feasibility']['status'] == 'INFEASIBLE' if check['validity'] == 'ABSTAIN' else None
            if check['validity'] == 'PASS' and taskrow['feasibility']['status'] == 'INFEASIBLE':
                raise ValueError('false acceptance contradicts geometric obstruction ' + tid)
            group.append(r)
        front = pareto(group)
        for r in group:
            r['pareto_nondominated'] = front.get(r['participant'])
            rows.append(r)
    summaries = []
    for split in ['dev', 'test']:
        for family in FAMILIES:
            for level in LEVELS:
                for name in participants:
                    rr = [r for r in rows if r['split'] == split and r['family'] == family and (r['level'] == level) and (r['participant'] == name)]
                    counts = dict(collections.Counter((r['checks']['validity'] for r in rr)))
                    mm = [r['metrics'] for r in rr if r['metrics'].get('status') == 'SCORED']
                    summaries.append(dict(split=split, family=family, level=level, participant=name, n=len(rr), counts=counts, height_rmse_median_mm=float(np.median([m['height_rmse_mm'] for m in mm])) if mm else None, contact_area_error_median_mm2=float(np.median([m['contacts']['0.1']['area_abs_error_mm2'] for m in mm if m['contacts']['0.1']['status'] == 'SCORED'])) if any((m['contacts']['0.1']['status'] == 'SCORED' for m in mm)) else None, correct_abstentions=sum((r.get('abstention_correct') is True for r in rr)), pareto_count=sum((r.get('pareto_nondominated') is True for r in rr))))
    totals = {name: dict(collections.Counter((r['checks']['validity'] for r in rows if r['participant'] == name))) for name in participants}
    out = dict(round=round_name, claim_type='capability', task_count=len(tasks), participant_count=len(participants), totals=totals, feasibility_counts=dict(collections.Counter((t['feasibility']['status'] for t in tasks))), summaries=summaries, complete_restoration_validation='UNKNOWN', uncertainty='sampled contact areas/approximal distances; native registration accuracy and target FDI error uncalibrated; mathematical surface-model checks only', external_referent=read(ROOT / 'PREREG_R1.json')['external_referent'], seconds=time.perf_counter() - start, frozen_predictions_sha256=sha(ROOT / f'FROZEN_PREDICTIONS_{round_name}.json'))
    dump(ROOT / f'raw/SCORED_ROWS_{round_name}.json', rows)
    dump(ROOT / f'raw/TASK_FEASIBILITY_{round_name}.json', tasks)
    dump(ROOT / f'rounds/{round_name}.json', out)
    return out
