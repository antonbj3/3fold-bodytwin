"""Vector frame with explicit anatomical-closure status, unlike discrete axis snapping."""
import numpy as np
from geometry import canonicalize, geometry, P
from bite2text.geom.canonical import CanonicalFrame, _axis_closure_score

def vector_frame(u, l):
    s = u.mean(axis=0) - l.mean(axis=0)
    s /= np.linalg.norm(s)
    q = np.eye(3) - np.outer(s, s)
    lp = l - l.mean(axis=0)
    cov = q @ (lp.T @ lp / len(lp)) @ q
    (w, V) = np.linalg.eigh(cov)
    opts = []
    for j in (-1, -2):
        a = V[:, j]
        a -= s * np.dot(a, s)
        a /= np.linalg.norm(a)
        right = np.cross(a, s)
        v = np.column_stack([lp @ right, lp @ a, lp @ s])
        (lo, hi) = _axis_closure_score(v, 1, 0)
        opts.append((abs(hi - lo), a * (1 if hi > lo else -1)))
    (margin, a) = max(opts, key=lambda p: p[0])
    right = np.cross(a, s)
    R = np.vstack([right, a, s])
    return CanonicalFrame(R, -1, -1, 1, 1, float(np.dot(u.mean(axis=0) - l.mean(axis=0), s)), float(margin))

def verify_rotation(frame):
    R = frame.rotation
    return bool(np.max(np.abs(R @ R.T - np.eye(3))) <= 1e-12 and abs(np.linalg.det(R) - 1) <= 1e-12)

def run():
    import json, zipfile, time, hashlib
    from scipy.spatial.transform import Rotation
    from geometry import load_stl_bytes, clean
    from benchmark import write
    r = json.load(open(P / 'PREREG_R3.json'))
    m = json.load(open(P / 'raw/DATA_MANIFEST.json'))
    rows = []
    start = time.perf_counter()
    with zipfile.ZipFile(m['zip']) as z:
        for c in r['cases']:
            (u, _) = load_stl_bytes(z.read(c + '/ios/ios_upper.stl'), 100000)
            (l, _) = load_stl_bytes(z.read(c + '/ios/ios_lower.stl'), 100000)
            origin = l.mean(axis=0)
            u = u - origin
            l = l - origin
            for meth in ['axis_published', 'vector_candidate', 'vector_conventional_control']:
                fr = canonicalize(u, l) if meth == 'axis_published' else vector_frame(u, l)
                assert verify_rotation(fr)
                (f, _) = geometry(u, l, c, fr)
                for e in r['rotations']:
                    rot = Rotation.from_euler('xyz', e, degrees=True).as_matrix()
                    uu = u @ rot.T
                    ll = l @ rot.T
                    fr2 = canonicalize(uu, ll) if meth == 'axis_published' else vector_frame(uu, ll)
                    (f2, _) = geometry(uu, ll, c, fr2)
                    es = {k: abs(f2[k] - v) for (k, v) in f.items() if (k.endswith('_mm') and 'width' not in k and ('wilson' not in k) and ('irregularity' not in k) and (k != 'arch_separation_mm')) and np.isfinite(v) and np.isfinite(f2[k])}
                    lost = [k for (k, v) in f.items() if k.endswith('_mm') and np.isfinite(v) and (not np.isfinite(f2[k]))]
                    rows.append(dict(case_id=c, method=meth, euler_xyz_deg=e, max_error_mm=max(es.values()) if es else None, errors_mm=es, lost_finite_features=lost, reference_features=clean(f), rotated_features=clean(f2), pass_gate=bool(es and max(es.values()) <= 0.05 and (not lost))))
    bad = vector_frame(u, l)
    bad.rotation[0] *= 2
    assert not verify_rotation(bad)
    result = dict(rows=rows, candidate_pass=all((v['pass_gate'] for v in rows if v['method'] == 'vector_candidate')), published_pass=all((v['pass_gate'] for v in rows if v['method'] == 'axis_published')), max_candidate_error_mm=max((v['max_error_mm'] for v in rows if v['method'] == 'vector_candidate')), max_published_error_mm=max((v['max_error_mm'] for v in rows if v['method'] == 'axis_published')), matched_control='TIE; identical conventional vector primitive', injected_bad_frame=dict(injection='scale first basis row by2', accepted=verify_rotation(bad), gate='FAIL', refutes=True), clinical_anatomical_accuracy='UNKNOWN ; no numerical land reference', wall_seconds=time.perf_counter() - start, external_referent=r['external_referent'])
    write('raw/RESULTS_R3.json', result)
    print('R3', result['candidate_pass'], result['max_candidate_error_mm'], 'axis', result['max_published_error_mm'], flush=True)
if __name__ == '__main__':
    run()
