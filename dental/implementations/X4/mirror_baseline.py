"""All-case Pre-only mirror baseline, including cases with no graft rail."""
import json, time, resource
from pathlib import Path
import numpy as np
from freeze import ROOT, DATA, write, sha, now
from planner import load, reflection, rigid_fit, transform, dist, stats

def run():
    reg = json.loads((ROOT / 'PREREG_R1.json').read_text())
    frozen = ROOT / 'FROZEN_MIRROR_BASELINE.json'
    records = []
    start = time.perf_counter()
    if not frozen.exists():
        for cid in reg['population']:
            (pre, a, meta) = load(cid, 'Pre')
            (mirror, plane) = reflection(pre, a)
            p = DATA / 'mirror_all' / f'{cid}.npz'
            p.parent.mkdir(exist_ok=True)
            np.savez_compressed(p, points=mirror.astype('f4'))
            records.append(dict(case=cid, path=str(p), sha256=sha(p), input=meta, plane=plane))
        write(frozen, dict(frozen_utc=now(), predictions=records, predictor_code_sha256=sha(__file__), information='target Pre only; coverage supplement, no changed graft prediction or frozen gate', prediction_seconds=time.perf_counter() - start))
    rows = []
    for row in json.loads(frozen.read_text())['predictions']:
        cid = row['case']
        assert sha(row['path']) == row['sha256']
        saved = ROOT / 'mirror_all' / f'{cid}.json'
        if saved.exists():
            rows.append(json.loads(saved.read_text()))
            continue
        (pre, a, _) = load(cid, 'Pre')
        (post, pa, _) = load(cid, 'Post')
        (T, fit) = rigid_fit(a[:4000], pa[:30000])
        q = transform(post, np.linalg.inv(T))
        d = dist(q, pre)
        roi = q[d > 3.0]
        with np.load(row['path']) as z:
            m = z['points'].astype(float)
        candidate = m[dist(m, pre) > 3.0]
        fd = np.minimum(d[d > 3.0], dist(roi, m))
        rd = dist(candidate, q)
        result = dict(case=cid, roi_n=len(roi), primary_p95_mm=float(max(np.quantile(fd, 0.95), np.quantile(rd, 0.95))) if len(fd) and len(rd) else None, forward=stats(fd), reverse=stats(rd), baseline='raw skull-plane Pre mirror proxy, not a fabricated graft', anchor_fit=fit, frozen_prediction=row, external_target='virtual Post', information='Pre-only prediction; Post evaluator')
        rawpath = DATA / 'mirror_all' / f'{cid}_distances.npz'
        np.savez_compressed(rawpath, forward=fd.astype('f4'), reverse=rd.astype('f4'))
        result['raw'] = dict(path=str(rawpath), sha256=sha(rawpath))
        write(saved, result)
        rows.append(result)
        if len(rows) % 20 == 0:
            print('MIRROR_ALL', len(rows), flush=True)
    vals = [r['primary_p95_mm'] for r in rows if r['primary_p95_mm'] is not None]
    out = dict(n=len(rows), n_available=len(vals), diagnostic_median_p95_mm=float(np.median(vals)) if vals else None, complete_cohort_median_p95_mm=float(np.median(vals)) if len(vals) == len(rows) else None, wall_seconds=time.perf_counter() - start, peak_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, cases=rows, meaning='coverage supplement for clinical-standard proxy; does not revise any prereg or graft prediction')
    write(ROOT / 'RESULTS_MIRROR_ALL.json', out)
if __name__ == '__main__':
    run()
