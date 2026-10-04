from common import *
import trimesh, resource
from scipy.spatial import cKDTree

def mesh(t):
    return trimesh.Trimesh(t.reshape(-1, 3), np.arange(len(t) * 3).reshape(-1, 3), process=False)

def component(t):
    m = trimesh.Trimesh(t.reshape(-1, 3), np.arange(len(t) * 3).reshape(-1, 3), process=True)
    parts = m.split(only_watertight=False)
    best = max(parts, key=lambda p: p.area)
    return (best.triangles.copy(), float(1 - best.area / m.area), len(parts))

def sample(t, n):
    a = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2
    i = np.searchsorted(np.cumsum(a), (np.arange(n) + 0.5) * a.sum() / n)
    return t[i].mean(1)

def distance(t, p):
    m = mesh(t)
    return np.concatenate([trimesh.proximity.closest_point(m, q)[1] for q in np.array_split(p, max(1, int(np.ceil(len(p) / 64))))])

def fit(s, t):
    cs = s.mean(0)
    ct = t.mean(0)
    ss = s - cs
    tt = t - ct
    tree = cKDTree(tt)
    out = []
    for angle in [0, np.pi / 2, np.pi, 3 * np.pi / 2]:
        c = np.cos(angle)
        sn = np.sin(angle)
        R = np.array([[c, -sn, 0], [sn, c, 0], [0, 0, 1.0]])
        u = np.zeros(3)
        for _ in range(60):
            a = ss @ R.T + u
            (d, j) = tree.query(a, workers=1)
            (dr, ir) = cKDTree(a).query(tt, workers=1)
            x = np.r_[ss, ss[ir]]
            y = np.r_[tt[j], tt]
            cx = x.mean(0)
            cy = y.mean(0)
            (U, _, V) = np.linalg.svd((x - cx).T @ (y - cy))
            Rn = V.T @ U.T
            if np.linalg.det(Rn) < 0:
                V[-1] *= -1
                Rn = V.T @ U.T
            R = Rn
            u = cy - cx @ R.T
        a = ss @ R.T + u
        err = float(np.sqrt(np.mean(np.r_[tree.query(a)[0], cKDTree(a).query(tt)[0]] ** 2)))
        out.append((err, R, ct + u - cs @ R.T))
    best = min(out, key=lambda x: x[0])
    return (best[1], best[2], [x[0] for x in out])

def metrics(a, b, n):
    da = distance(b, sample(a, n))
    db = distance(a, sample(b, n))
    return dict(p95_mm=float(max(np.quantile(da, 0.95), np.quantile(db, 0.95))), rms_equal_surface_mm=float(np.sqrt(np.mean(np.r_[da, db] ** 2))), sampled_max_mm=float(max(da.max(), db.max())), directed_p95_mm=[float(np.quantile(x, 0.95)) for x in [da, db]])

def run():
    start = time.perf_counter()
    pr = read(ROOT / 'PREREG_LOCAL.json')
    rows = []
    for (p, h) in pr['input_hashes'].items():
        assert sha(p) == h, p
    for r in pr['cohort']:
        key = r['key']
        p = npz(PUBLIC / (key + '.npz'))
        a = npz(V4 / 'payload/whole_private' / key / 'reference.npz')
        st = time.perf_counter()
        row = dict(key=key, family=r['family'], case_key=r['case_key'], dataset=r['dataset'], resolution='PER_TOOTH', status='FAILED')
        try:
            (s, drop_s, n_s) = component(p['homolog_triangles'])
            (t, drop_t, n_t) = component(a['source_triangles'])
            s[:, :, 0] *= -1
            (R, u, starts) = fit(sample(s, 1536), sample(t, 1536))
            aligned = s @ R.T + u
            m = float(p['margin_z'])
            aa = aligned[aligned.mean(1)[:, 2] >= m]
            bb = t[t.mean(1)[:, 2] >= m]
            mm = {str(n): metrics(aa, bb, n) for n in pr['paired_protocol']['evaluate_samples']}
            selfmax = float(distance(bb, sample(bb, 2048)).max())
            row.update(status='MEASURED', metrics=mm, registration_start_rms_mm=starts, rotation=R, translation_mm=u, homolog_dropped_area_fraction=drop_s, target_dropped_area_fraction=drop_t, components=[n_s, n_t], retained_faces=[len(aa), len(bb)], self_max_mm=selfmax, identity_control_pass=selfmax <= 1e-10, clinical_acceptance='UNKNOWN; no same-preparation designs', probe_sensitivity_p95_mm=abs(mm['2048']['p95_mm'] - mm['8192']['p95_mm']))
            np.savez_compressed(DATA / (key + '_aligned.npz'), homolog=aa, target=bb)
        except Exception as e:
            row['reason'] = repr(e)
        row['seconds'] = time.perf_counter() - st
        rows.append(row)
        dump(ROOT / 'raw/LOCAL_PAIRS.json', dict(rows=rows, seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
        print(key, row['status'], row.get('metrics', {}).get('8192'), flush=True)
        state('LOCAL_MEASURING', str(len(rows)) + '/18 pairs processed', 'Complete pairs; then source comparison')
    return rows
if __name__ == '__main__':
    run()
