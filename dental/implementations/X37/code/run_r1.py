import csv, json, hashlib, time, resource, zipfile, datetime
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from geometry import *
t0 = time.monotonic()
records = []
poses = {}
measurements = []
rejected = []
controls = []
source_hashes = []

def save(name, x):
    (ROOT / 'raw' / name).write_text(json.dumps(x, indent=2) + '\n')

def refs(ids):
    return {'all_crowns': ids, 'bilateral_molars': [i for i in ids if i % 10 in [6, 7]], 'left_molars': [i for i in ids if i // 10 in [2, 3] and i % 10 in [6, 7]], 'right_molars': [i for i in ids if i // 10 in [1, 4] and i % 10 in [6, 7]]}
z0 = zipfile.ZipFile(DATA / '3485_Clinical_Trial.zip')
probe = load_cloud(z0, '3485', 'Sirona', 0, 11)
a = probe['points'][:4000].astype(float)
n = probe['normals'][:4000]
known = mat(Rotation.from_euler('xyz', [2, -3, 5], degrees=True).as_matrix(), np.array([0.35, -0.2, 0.15]))
b = apply(known, a)
bn = n @ known[:3, :3].T
(out, e) = icp(a, b, bn, initial=mat(t=b.mean(0) - a.mean(0)))
err = out @ np.linalg.inv(known)
num = {'translation_error_mm': float(np.linalg.norm(err[:3, 3])), 'rotation_error_deg': angle(err), 'pass': bool(np.linalg.norm(err[:3, 3]) <= 0.01 and angle(err) <= 0.1), 'external_referent': {'kind': 'closed_form', 'locator': 'SE(3) q=Rp+t; https://doi.org/10.1107/S0567739476001873', 'compared_quantity': 'known rigid transform recovery', 'refutes_us': True}}
num['injected_1mm_rejected'] = bool(np.linalg.norm(err[:3, 3] + np.array([1, 0, 0])) > 0.01)
num['injected_10deg_rejected'] = angle(mat(Rotation.from_euler('x', 10, degrees=True).as_matrix()) @ err) > 0.1
save('R1_NUMERICAL_CONTROLS.json', num)
for patient in ['3485', '6457']:
    z = zipfile.ZipFile(DATA / f'{patient}_Clinical_Trial.zip')
    names = set(z.namelist())
    for (jaw, quad) in [('OK', [1, 2]), ('UK', [3, 4])]:
        ids = [f for f in range(11, 48) if f // 10 in quad and f % 10 in range(1, 8) and (member(patient, 'Sirona', 0, f) in names)]
        base = {f: load_cloud(z, patient, 'Sirona', 0, f) for f in ids}
        rp = refs(ids)
        (x0, n0) = union(base, ids, 1000)
        for src in ['Sirona', 'Bottmedical']:
            for t in range(10):
                avail = [f for f in ids if member(patient, src, t, f) in names]
                for f in set(ids) - set(avail):
                    rejected.append({'patient': patient, 'source': src, 'week': t, 'fdi': f, 'reason': 'missing_segment'})
                target = {f: load_cloud(z, patient, src, t, f) for f in avail}
                (target_all, tn) = union(target, avail, 1000)
                (G_all, gerr) = icp(target_all, x0, n0, initial=mat(t=x0.mean(0) - target_all.mean(0)))
                Gs = {'all_crowns': G_all}
                for (name, rs) in rp.items():
                    if name == 'all_crowns':
                        continue
                    use = [f for f in rs if f in avail]
                    (xr, nr) = union(base, use)
                    (xt, nt) = union(target, use)
                    (Gs[name], eg) = icp(xt, xr, nr, initial=G_all)
                poses[patient, jaw, src, t] = {}
                for f in avail:
                    cloud = target[f]
                    aa = base[f]['points'][:4000]
                    if src == 'Sirona' and t == 0:
                        L = np.eye(4)
                        ee = {'iterations': 0, 'trimmed_plane_rms_mm': 0.0, 'trimmed_point_rms_mm': 0.0}
                    else:
                        (L, ee) = icp(aa, cloud['points'], cloud['normals'], initial=np.linalg.inv(G_all))
                    alternatives = {name: G @ L for (name, G) in Gs.items()}
                    boot = []
                    for s in range(4):
                        (ll, be) = icp(aa[s::4], cloud['points'], cloud['normals'], initial=L, iterations=20)
                        boot.append(Gs['bilateral_molars'] @ ll)
                    (C, ce) = icp(aa, cloud['points'], cloud['normals'], initial=L, method='point', iterations=30)
                    cc = base[f]['centroid']
                    hb = alternatives['bilateral_molars']
                    dr = [np.linalg.norm(apply(h, cc) - apply(hb, cc)) for h in alternatives.values()]
                    da = [angle(h @ np.linalg.inv(hb)) for h in alternatives.values()]
                    u = max((np.linalg.norm(apply(h, cc) - apply(hb, cc)) for h in boot))
                    good = ee['trimmed_plane_rms_mm'] <= 0.1
                    row = {'patient': patient, 'jaw': jaw, 'source': src, 'week': t, 'fdi': f, 'resolution': 'PER_TOOTH', 'surface_fit_pass': good, 'surface_error': ee, 'reference_spread_mm': float(max(dr)), 'reference_rotation_spread_deg': float(max(da)), 'sampling_spread_mm': float(u), 'conditional_reference_numerical_pass': bool(max(dr) <= 0.1 and max(da) <= 1), 'absolute_movement_status': 'UNKNOWN_NO_INDEPENDENT_STABLE_REFERENCE', 'H': {k: h.tolist() for (k, h) in alternatives.items()}, 'point_control_centroid_delta_mm': float(np.linalg.norm(apply(Gs['bilateral_molars'] @ C, cc) - apply(hb, cc))), 'point_control_rotation_delta_deg': angle(C @ np.linalg.inv(L)), 'centroid_base_mm': cc.tolist()}
                    records.append(row)
                    poses[patient, jaw, src, t][f] = row
                    source_hashes.append({'name': member(patient, src, t, f), 'sha256': str(cloud['source_sha256']), 'cache': str(CACHE / f'{patient}_{src}_T{t}_Z{f}.npz')})
                print(patient, jaw, src, t, 'teeth', len(avail), 'elapsed_s', round(time.monotonic() - t0, 1), flush=True)
                save('R1_POSE_CHECKPOINT.json', {'last': [patient, jaw, src, t], 'rows': len(records), 'rejected': rejected, 'elapsed_s': time.monotonic() - t0})
    for (jaw, _) in [('OK', [1, 2]), ('UK', [3, 4])]:
        ids = list(poses[patient, jaw, 'Sirona', 0])
        for t in range(1, 10):
            for f in ids:
                rr = [poses.get((patient, jaw, s, w), {}).get(f) for (s, w) in [('Sirona', t), ('Sirona', t - 1), ('Bottmedical', t), ('Bottmedical', t - 1), ('Bottmedical', 0)]]
                if any((x is None for x in rr)):
                    continue
                (st, sp, pt, pp, p0) = rr
                c = np.array(st['centroid_base_mm'])
                variant = {}
                for name in refs(ids):
                    A = np.array(st['H'][name])
                    B = np.array(sp['H'][name])
                    P = np.array(pt['H'][name]) @ np.linalg.inv(np.array(p0['H'][name]))
                    Q = np.array(pp['H'][name]) @ np.linalg.inv(np.array(p0['H'][name]))
                    variant[name] = {'observed_weekly_delta_mm': (apply(A, c) - apply(B, c)).tolist(), 'planned_weekly_delta_mm': (apply(P, c) - apply(Q, c)).tolist(), 'observed_cumulative_mm': float(np.linalg.norm(apply(A, c) - c)), 'planned_cumulative_mm': float(np.linalg.norm(apply(P, c) - c)), 'tracking_error_mm': float(np.linalg.norm(apply(A, c) - apply(P, c))), 'weekly_rotation_deg': angle(A @ np.linalg.inv(B)), 'planned_weekly_rotation_deg': angle(P @ np.linalg.inv(Q)), 'tracking_rotation_deg': angle(A @ np.linalg.inv(P))}
                vv = variant['bilateral_molars']
                mv = np.array(vv['observed_weekly_delta_mm'])
                pv = np.array(vv['planned_weekly_delta_mm'])
                u = st['reference_spread_mm'] + sp['reference_spread_mm'] + st['sampling_spread_mm'] + sp['sampling_spread_mm'] + 0.2
                ue = st['reference_spread_mm'] + pt['reference_spread_mm'] + p0['reference_spread_mm'] + st['sampling_spread_mm'] + pt['sampling_spread_mm'] + 0.2
                measurements.append({'patient': patient, 'jaw': jaw, 'fdi': f, 'week': t, 'resolution': 'PER_TOOTH', 'observed_weekly_mm': float(np.linalg.norm(mv)), 'planned_weekly_mm': float(np.linalg.norm(pv)), 'observed_weekly_vector_mm': mv.tolist(), 'planned_weekly_vector_mm': pv.tolist(), 'observed_cumulative_mm': vv['observed_cumulative_mm'], 'planned_cumulative_mm': vv['planned_cumulative_mm'], 'tracking_error_mm': vv['tracking_error_mm'], 'weekly_rotation_deg': vv['weekly_rotation_deg'], 'planned_weekly_rotation_deg': vv['planned_weekly_rotation_deg'], 'tracking_rotation_deg': vv['tracking_rotation_deg'], 'weekly_sensitivity_budget_mm': u, 'tracking_sensitivity_budget_mm': ue, 'budget_kind': 'operational worst-sum with declared PHENOMENOLOGICAL floor, not CI', 'surface_fit_pass': all((x['surface_fit_pass'] for x in rr)), 'absolute_movement_status': 'UNKNOWN_NO_INDEPENDENT_STABLE_REFERENCE', 'conditional_tracking_gap_resolved': bool(vv['tracking_error_mm'] > ue and all((x['surface_fit_pass'] for x in rr))), 'reference_variants': variant})
save('R1_POSES.json', records)
save('R1_WEEKLY_MEASUREMENTS.json', measurements)
save('SOURCE_SURFACE_HASHES.json', source_hashes)
flat = [{k: v for (k, v) in m.items() if k not in ['reference_variants', 'observed_weekly_vector_mm', 'planned_weekly_vector_mm']} for m in measurements]
with (ROOT / 'raw' / 'weekly_movement.csv').open('w') as o:
    w = csv.DictWriter(o, fieldnames=list(flat[0]))
    w.writeheader()
    w.writerows(flat)
r = {'claim_type': 'information_link', 'round': 'R1', 'external_referent': json.loads((ROOT / 'PREREG_R1.json').read_text())['external_referent'], 'pose_rows': len(records), 'weekly_tooth_rows': len(measurements), 'patients': 2, 'missing_segments': rejected, 'absolute_reference_gate': 'FAIL_NO_INDEPENDENT_STABLE_REFERENCE', 'surface_fit_rejections': sum((not x['surface_fit_pass'] for x in records)), 'conditional_reference_gate_failures': sum((not x['conditional_reference_numerical_pass'] for x in records)), 'conditional_tracking_gap_resolved': sum((x['conditional_tracking_gap_resolved'] for x in measurements)), 'control': {'max_point_control_centroid_delta_mm': max((x['point_control_centroid_delta_mm'] for x in records)), 'median_point_control_centroid_delta_mm': float(np.median([x['point_control_centroid_delta_mm'] for x in records]))}, 'numerical_controls': num, 'elapsed_s': time.monotonic() - t0, 'peak_rss_mib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'by_patient': {p: {'rows': sum((x['patient'] == p for x in measurements)), 'resolved_tracking_rows': sum((x['patient'] == p and x['conditional_tracking_gap_resolved'] for x in measurements)), 'median_observed_weekly_mm': float(np.median([x['observed_weekly_mm'] for x in measurements if x['patient'] == p])), 'median_planned_weekly_mm': float(np.median([x['planned_weekly_mm'] for x in measurements if x['patient'] == p])), 'max_tracking_mm': max((x['tracking_error_mm'] for x in measurements if x['patient'] == p))} for p in ['3485', '6457']}}
save('R1_RESULTS.json', r)
print(json.dumps(r, indent=2))
state = json.loads((ROOT / 'CURRENT_WORK_STATE.json').read_text())
state.update(status='R1_MEASURED', latest_gate=r['absolute_reference_gate'], next_operation='Freeze R2 gauge-invariant pairwise observable; preserve R1 failure', updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
(ROOT / 'CURRENT_WORK_STATE.json').write_text(json.dumps(state, indent=2) + '\n')
