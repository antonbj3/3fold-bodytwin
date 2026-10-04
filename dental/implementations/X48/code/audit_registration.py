"""Independent repeat fitting to existing measured crowns, not repeat acquisition."""
from common import *
import sys, importlib.util, zipfile, time
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('x37_geometry_readonly', X37 / 'code/geometry.py')
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)

def run():
    start = time.monotonic()
    (rows, ix) = pose_index()
    audits = []
    sources = []
    for p in ['3485', '6457']:
        with zipfile.ZipFile(DATA / f'{p}_Clinical_Trial.zip') as z:
            for f in [17, 27]:
                clouds = []
                for t in [0, 8]:
                    n = g.member(p, 'Sirona', t, f)
                    blob = z.read(n)
                    d = np.load(CACHE / f'{p}_Sirona_T{t}_Z{f}.npz')
                    H = g.initial_frame(z, p, 'Sirona', t, 'OK')
                    seed = int(hashlib.sha256(n.encode()).hexdigest()[:8], 16)
                    (pts, ns, c, area) = g.sample_mesh(blob, H, seed)
                    assert hashlib.sha256(blob).hexdigest() == str(d['source_sha256']), 'raw STL source mismatch'
                    assert np.array_equal(pts.astype('float32'), d['points']) and np.array_equal(ns.astype('float32'), d['normals']), 'cached samples do not match actual STL'
                    sources.append(dict(patient=p, fdi=f, week=t, member=n, bytes=len(blob), sha256=hashlib.sha256(blob).hexdigest(), points_recomputed_from_external_STL=True))
                    clouds.append(dict(points=pts, normals=ns, centroid=c))
                a = clouds[0]['points'][:4000]
                b = clouds[1]['points']
                nb = clouds[1]['normals']
                initial = g.mat(t=b.mean(0) - a.mean(0))
                (F, ef) = g.icp(a, b, nb, initial=initial)
                (C, ec) = g.icp(a, b, nb, initial=initial, method='point')
                (B, eb) = g.icp(b[:4000], clouds[0]['points'], clouds[0]['normals'], initial=np.linalg.inv(F))
                c = clouds[0]['centroid']
                inverse_error = float(np.linalg.norm(g.apply(B @ F, c) - c))
                control_error = float(np.linalg.norm(g.apply(F, c) - g.apply(C, c)))
                partitions = []
                for k in range(4):
                    (h, e) = g.icp(a[k::4], b, nb, initial=F, iterations=25)
                    partitions.append(dict(centroid_mm=float(np.linalg.norm(g.apply(h, c) - g.apply(F, c))), rotation_deg=g.angle(h @ np.linalg.inv(F))))
                audits.append(dict(patient=p, fdi=f, interval='T0_to_T8', resolution='PER_TOOTH', plane_fit_RMS_mm=ef['trimmed_plane_rms_mm'], point_fit_RMS_mm=ec['trimmed_point_rms_mm'], inverse_centroid_error_mm=inverse_error, point_plane_centroid_delta_mm=control_error, point_plane_rotation_delta_deg=g.angle(C @ np.linalg.inv(F)), partition_max_centroid_mm=max((x['centroid_mm'] for x in partitions)), partition_max_rotation_deg=max((x['rotation_deg'] for x in partitions)), fit_pass=ef['trimmed_plane_rms_mm'] <= 0.1, patient_repeat_acquisition_error='UNKNOWN', interpretation='Fit/inverse/partition floor is numerical and rigid-shape sensitivity; not scanner repeatability and not anatomical stability'))
    a = clouds[0]['points'][:4000]
    n = clouds[0]['normals'][:4000]
    known = g.mat(Rotation.from_euler('xyz', [2, -3, 5], degrees=True).as_matrix(), np.array([0.35, -0.2, 0.15]))
    target = g.apply(known, a)
    tn = n @ known[:3, :3].T
    controls = []
    for method in ['point', 'plane']:
        (h, e) = g.icp(a, target, tn, initial=g.mat(t=target.mean(0) - a.mean(0)), method=method)
        err = h @ np.linalg.inv(known)

        def accepts(H):
            return np.linalg.norm(H[:3, 3]) <= 0.01 and g.angle(H) <= 0.1
        badt = err.copy()
        badt[0, 3] += 1.0
        badr = g.mat(Rotation.from_euler('x', 10, degrees=True).as_matrix()) @ err
        controls.append(dict(method=method, translation_error_mm=float(np.linalg.norm(err[:3, 3])), rotation_error_deg=g.angle(err), passed=bool(accepts(err)), wrong_1mm_rejected=bool(not accepts(badt)), wrong_10deg_rejected=bool(not accepts(badr))))
    out = dict(audits=audits, raw_members=sources, controls=controls, code_locator=str(X37 / 'code/geometry.py'), code_sha256=sha(X37 / 'code/geometry.py'), source_scanner_mean_um=50.0, source_scanner_SD_um=25.0, source_locator='10.3390/oral4040039 Table1 p494', physical_scanner_pose_floor='UNKNOWN', wall_s=time.monotonic() - start)
    assert all((x['passed'] and x['wrong_1mm_rejected'] and x['wrong_10deg_rejected'] for x in controls))
    write(ROOT / 'raw/REGISTRATION_AUDIT.json', out)
    print('Independent raw-STL registration audit', len(audits), round(out['wall_s'], 2))
    return out
if __name__ == '__main__':
    run()
