"""Private research adapter to pinned published Bite2Text geometry; point proximity is a proxy."""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '2'
import sys, io, hashlib, time, zipfile, json, traceback, resource
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / 'sources/upstream/src'))
from bite2text.geom import canonicalize, extract_arch_profile, CaseGeometry, measure_case

def load_stl_bytes(blob, max_points=200000):
    if len(blob) >= 84:
        n = int.from_bytes(blob[80:84], 'little')
        if 84 + 50 * n == len(blob):
            dt = np.dtype([('normal', '<f4', (3,)), ('vertices', '<f4', (3, 3)), ('attr', '<u2')])
            tri = np.frombuffer(blob, dtype=dt, offset=84, count=n)['vertices']
            pts = tri.reshape(-1, 3)
        else:
            pts = None
    else:
        pts = None
    if pts is None:
        import trimesh
        pts = np.asarray(trimesh.load(io.BytesIO(blob), file_type='stl', process=False).vertices)
    raw_n = len(pts)
    if raw_n > max_points:
        pts = pts[np.linspace(0, raw_n - 1, max_points, dtype=np.int64)]
    pts = np.asarray(pts, dtype=np.float64)
    pts = pts[np.isfinite(pts).all(axis=1)]
    if len(pts) < 1000:
        raise ValueError('Unusable mesh: fewer than 1000 finite points')
    return (pts, raw_n)

def geometry(upper, lower, case_id, frame=None):
    frame = frame or canonicalize(upper, lower)
    (up, lw) = (frame.apply(upper), frame.apply(lower))
    lp = extract_arch_profile(lw, 'lower')
    uprof = extract_arch_profile(up, 'upper', pole=lp.pole)
    geom = CaseGeometry(case_id, frame, up, lw, uprof, lp)
    return (measure_case(geom).features(), geom)

def clean(d):
    if isinstance(d, dict):
        return {k: clean(v) for (k, v) in d.items()}
    if isinstance(d, (list, tuple)):
        return [clean(v) for v in d]
    if isinstance(d, np.ndarray):
        return clean(d.tolist())
    if isinstance(d, (np.floating, float)):
        return float(d) if np.isfinite(d) else None
    if isinstance(d, np.integer):
        return int(d)
    return d

def proximity(up, lw):
    si = np.linspace(0, len(up) - 1, min(30000, len(up)), dtype=np.int64)
    (d, ind) = cKDTree(lw).query(up[si], workers=2)
    nearest = np.argsort(d)[:100]
    return dict(metric='sampled upper-vertex to lower-vertex Euclidean distance; surface-distance upper bound', min_mm=float(d.min()), counts_below_0_1mm=int((d < 0.1).sum()), counts_below_0_3mm=int((d < 0.3).sum()), sample_count=len(si), candidates=[dict(upper_xyz_mm=up[si[k]].tolist(), lower_xyz_mm=lw[ind[k]].tolist(), distance_mm=float(d[k])) for k in nearest], load_N=None, physical_contact_validated=False)

def extract():
    manifest = json.loads((P / 'raw/DATA_MANIFEST.json').read_text())
    ids = manifest['train'] + manifest['calibration'] + manifest['test']
    dest = P / 'raw/geometry.jsonl'
    done = {json.loads(l)['case_id'] for l in dest.read_text().splitlines()} if dest.exists() else set()
    start = time.perf_counter()
    print('Extracting', len(ids), 'resume', len(done), flush=True)
    with zipfile.ZipFile(manifest['zip']) as z, dest.open('a') as out:
        for (i, c) in enumerate(ids):
            if c in done:
                continue
            t = time.perf_counter()
            row = {'case_id': c}
            try:
                bup = z.read(c + '/ios/ios_upper.stl')
                blw = z.read(c + '/ios/ios_lower.stl')
                (u, nu) = load_stl_bytes(bup)
                (l, nl) = load_stl_bytes(blw)
                (f, g) = geometry(u, l, c)
                row.update(features=clean(f), frame=clean(g.frame.__dict__), n_raw_vertices=[nu, nl], n_sampled_vertices=[len(u), len(l)], member_sha256=[hashlib.sha256(bup).hexdigest(), hashlib.sha256(blw).hexdigest()], contact=clean(proximity(g.upper_points, g.lower_points)), error=None)
            except Exception as e:
                row['error'] = repr(e)
                row['traceback'] = traceback.format_exc()
            row['wall_seconds'] = time.perf_counter() - t
            out.write(json.dumps(row, allow_nan=False) + '\n')
            out.flush()
            if (i + 1) % 20 == 0:
                print(i + 1, c, row['error'], round(time.perf_counter() - start, 1), 's', flush=True)
                state = dict(lane='X7-bite2text', status='R1_GEOMETRY_EXTRACTION', updated_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(), cases_written=i + 1, latest_gate='PREREG_R1 frozen; geometry not compared with test labels', next_operation='Freeze model predictions then parse test reports')
                (P / 'CURRENT_WORK_STATE.json').write_text(json.dumps(state, indent=2))
    print('finished', round(time.perf_counter() - start, 2), 's peakRSS', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, flush=True)
if __name__ == '__main__':
    extract()
