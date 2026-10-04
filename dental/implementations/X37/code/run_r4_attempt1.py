import json, time, resource, datetime, zipfile, numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from geometry import *
start = time.monotonic()
rows = []
roi_clouds = {}
cache_files = []

def orient(patient, t):
    xyz = [90, 0, 0]
    if t == 8:
        xyz = [180, 0, 0] if patient == '3485' else [0, 90, 180]
    if patient == '6457' and t == 4:
        xyz = [-100, 180, 40]
    Rx = Rotation.from_euler('x', xyz[0], degrees=True).as_matrix()
    Ry = Rotation.from_euler('y', xyz[1], degrees=True).as_matrix()
    Rz = Rotation.from_euler('z', xyz[2], degrees=True).as_matrix()
    return mat((Rz @ Ry @ Rx).T)
for p in ['3485', '6457']:
    z = zipfile.ZipFile(DATA / f'{p}_Clinical_Trial.zip')
    for t in [0, 9]:
        name = f'{p}_Clinical_Trial/Sirona/T{t}/{p}_OnyxCeph3_Export_OK.stl'
        cache = CACHE / f'{p}_full_upper_T{t}_100k.npz'
        if cache.exists():
            d = np.load(cache)
            pts = d['points']
            norms = d['normals']
            sha = str(d['source_sha256'])
        else:
            blob = z.read(name)
            H = initial_frame(z, p, 'Sirona', t, 'OK') @ orient(p, t)
            (pts, norms, c, area) = sample_mesh(blob, H, int(digest(name.encode())[:8], 16), n=100000)
            sha = digest(blob)
            np.savez_compressed(cache, points=pts.astype('float32'), normals=norms.astype('float32'), source_sha256=np.array(sha))
        ids = [f for f in range(11, 28) if f % 10 in range(1, 8) and member(p, 'Sirona', t, f) in z.namelist()]
        crowns = {f: load_cloud(z, p, 'Sirona', t, f) for f in ids}
        (cp, cn) = union(crowns, ids, n=4000)
        tree = cKDTree(cp)
        (dist, _) = tree.query(pts, workers=4)
        zi = np.mean([crowns[f]['centroid'][2] for f in ids if f % 10 == 1])
        zm = np.mean([crowns[f]['centroid'][2] for f in ids if f % 10 == 6])
        lo = min(zm + 3, zi - 5)
        hi = max(zm + 3, zi - 5)
        mask = (np.abs(pts[:, 0]) <= 8) & (pts[:, 2] >= lo) & (pts[:, 2] <= hi) & (dist >= 3)
        rp = pts[mask]
        rn = norms[mask]
        ext = np.ptp(rp, axis=0).tolist() if len(rp) else [0.0, 0.0, 0.0]
        coverage = len(rp) >= 500 and ext[0] >= 8 and (ext[2] >= 8)
        rows.append({'patient': p, 'week': t, 'jaw': 'OK', 'resolution': 'PER_SURFACE_REGION', 'original_scan_member': name, 'source_sha256': sha, 'sample_points': len(pts), 'candidate_roi_points': len(rp), 'roi_extent_mm': ext, 'z_bounds_mm': [lo, hi], 'coverage_pass': bool(coverage), 'anatomical_identification': 'UNKNOWN_AUTOMATIC_INNER_SURFACE_NOT_ANNOTATED_RUGAE'})
        roi_clouds[p, t] = (rp, rn)
        cache_files.append(str(cache))
        print(p, t, len(rp), ext, coverage, flush=True)
fits = []
for p in ['3485', '6457']:
    rr = [r for r in rows if r['patient'] == p]
    if not all((r['coverage_pass'] for r in rr)):
        fits.append({'patient': p, 'status': 'REJECTED_MISSING_ROI_COVERAGE', 'absolute_movement': 'UNKNOWN'})
        continue
    (a, an) = roi_clouds[p, 9]
    (b, bn) = roi_clouds[p, 0]
    a = a[:4000]
    (H, e) = icp(a, b, bn, initial=mat(t=b.mean(0) - a.mean(0)))
    center = b.mean(0)
    J = np.c_[np.cross(b - center, bn) / 10, bn]
    sv = np.linalg.svd(J, compute_uv=False)
    condition = float(sv[0] / sv[-1]) if sv[-1] > 1e-12 else None
    left = a[:, 0] <= np.median(a[:, 0])
    Hs = []
    for mask in [left, ~left]:
        (hh, ee) = icp(a[mask], b, bn, initial=H)
        Hs.append(hh)
    splitdiff = float(np.linalg.norm(apply(Hs[0], a.mean(0)) - apply(Hs[1], a.mean(0))))
    (control, ce) = icp(a, b, bn, initial=H, method='point')
    passed = e['trimmed_plane_rms_mm'] <= 0.1 and splitdiff <= 0.1 and (condition is not None) and (condition <= 10000.0)
    fits.append({'patient': p, 'status': 'NUMERICAL_REFERENCE_PASS' if passed else 'REJECTED_REFERENCE_GATE', 'fit': e, 'scaled_6dof_jacobian_condition': condition, 'left_right_centroid_difference_mm': splitdiff, 'point_control_centroid_difference_mm': float(np.linalg.norm(apply(H, a.mean(0)) - apply(control, a.mean(0)))), 'wrong_1mm_residual_rejected': 1.0 > 0.1, 'absolute_movement': 'UNKNOWN_NO_ANNOTATED_AND_INDEPENDENTLY_STABLE_REGION'})
result = {'round': 'R4', 'claim_type': 'information_link', 'external_referent': json.loads((ROOT / 'PREREG_R4.json').read_text())['external_referent'], 'mesh_cases': len(rows), 'coverage_pass': sum((r['coverage_pass'] for r in rows)), 'coverage_rejected': sum((not r['coverage_pass'] for r in rows)), 'rows': rows, 'fits': fits, 'absolute_movement': 'UNKNOWN', 'cache_files': cache_files, 'elapsed_s': time.monotonic() - start, 'peak_rss_mib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024}
(ROOT / 'raw/R4_RESULTS.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
state = json.loads((ROOT / 'CURRENT_WORK_STATE.json').read_text())
state.update(status='R4_MEASURED', latest_gate='Raw inner-surface reference feasibility measured; absolute region provenance/stability unknown', next_operation='Package reproducible demo, figure, dropout, graph proposal and next minimal matched measurement', updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
(ROOT / 'CURRENT_WORK_STATE.json').write_text(json.dumps(state, indent=2) + '\n')
