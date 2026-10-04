import json, csv, time, resource, zipfile, datetime
import numpy as np
from geometry import *
from region import load_full, roi, heldout
start = time.monotonic()
poses = json.loads((ROOT / 'raw/R1_POSES.json').read_text())
ix = {(r['patient'], r['jaw'], r['source'], r['week'], r['fdi']): r for r in poses}
out = []
regs = []
for p in ['3485', '6457']:
    z = zipfile.ZipFile(DATA / f'{p}_Clinical_Trial.zip')
    ids = sorted({r['fdi'] for r in poses if r['patient'] == p and r['jaw'] == 'OK'})
    base = {f: load_cloud(z, p, 'Sirona', 0, f) for f in ids}
    (bp, bn) = union(base, ids, 1000)
    (full0, ns0, sha0) = load_full(z, p, 0)
    (r0, n0, cov0) = roi(full0, ns0, base)
    for t in range(10):
        clouds = {f: load_cloud(z, p, 'Sirona', t, f) for f in ids}
        (pts, ns, sha) = load_full(z, p, t)
        (r, n, cov) = roi(pts, ns, clouds)
        if not cov['coverage_pass'] or not cov0['coverage_pass']:
            regs.append({'patient': p, 'week': t, 'coverage': cov, 'gate_pass': False, 'reason': 'coverage'})
            continue
        a = r[:4000]
        (H, e) = icp(a, r0, n0, initial=mat(t=r0.mean(0) - a.mean(0)))
        left = a[:, 0] <= np.median(a[:, 0])
        Hs = []
        held = []
        for mask in [left, ~left]:
            (hh, ee) = icp(a[mask], r0, n0, initial=H)
            Hs.append(hh)
            held.append(heldout(hh, a[~mask], r0, n0))
        J = np.c_[np.cross(r0 - r0.mean(0), n0) / 10, n0]
        sv = np.linalg.svd(J, compute_uv=False)
        condition = float(sv[0] / sv[-1]) if sv[-1] > 1e-12 else np.inf
        good = e['trimmed_plane_rms_mm'] <= 0.1 and max(held) <= 0.1 and (condition <= 10000.0)
        (cp, cn) = union(clouds, ids, 1000)
        (G, _) = icp(cp, bp, bn, initial=mat(t=bp.mean(0) - cp.mean(0)))
        regs.append({'patient': p, 'week': t, 'coverage': cov, 'fit': e, 'held_out_rms_mm': held, 'jacobian_condition': condition, 'gate_pass': bool(good), 'region_to_base_H': H.tolist(), 'source_sha256': sha, 'wrong1mm_reference_rejected': float(np.linalg.norm(apply(mat(t=np.array([1.0, 0.0, 0.0])) @ H, a.mean(0)) - apply(H, a.mean(0)))) > 0.1, 'anatomy_stability': 'UNKNOWN'})
        for f in ids:
            rr = ix[p, 'OK', 'Sirona', t, f]
            L = np.linalg.inv(G) @ np.array(rr['H']['all_crowns'])
            C = H @ L
            c = base[f]['centroid']
            split = max((float(np.linalg.norm(apply(h @ L, c) - apply(C, c))) for h in Hs))
            out.append({'patient': p, 'jaw': 'OK', 'fdi': f, 'week': t, 'resolution': 'PER_TOOTH', 'conditional_region_frame_pose_H': C.tolist(), 'material_center_mm': apply(C, c).tolist(), 'cumulative_motion_mm': float(np.linalg.norm(apply(C, c) - c)), 'reference_split_error_mm': split, 'sampling_error_mm': rr['sampling_spread_mm'], 'reference_gate_pass': bool(good), 'per_tooth_split_gate_pass': split <= 0.1, 'surface_gate_pass': rr['surface_fit_pass'], 'absolute_movement': 'UNKNOWN_REGION_STABILITY', 'frame': 'Candidate central non-dental surface, anatomy unverified'})
        print(p, t, 'ROI', len(r), 'heldout', round(max(held), 4), 'fitpass', good, flush=True)
oi = {(x['patient'], x['week'], x['fdi']): x for x in out}
weekly = []
for r in out:
    if r['week'] == 0:
        continue
    q = oi.get((r['patient'], r['week'] - 1, r['fdi']))
    if q is None:
        continue
    v = np.array(r['material_center_mm']) - np.array(q['material_center_mm'])
    u = 0.2 + r['reference_split_error_mm'] + q['reference_split_error_mm'] + r['sampling_error_mm'] + q['sampling_error_mm']
    good = all((x['reference_gate_pass'] and x['per_tooth_split_gate_pass'] and x['surface_gate_pass'] for x in [r, q]))
    H = np.array(r['conditional_region_frame_pose_H'])
    P = np.array(q['conditional_region_frame_pose_H'])
    weekly.append({'patient': r['patient'], 'jaw': 'OK', 'fdi': r['fdi'], 'week': r['week'], 'resolution': 'PER_TOOTH', 'weekly_vector_mm': v.tolist(), 'weekly_norm_mm': float(np.linalg.norm(v)), 'weekly_rotation_deg': angle(H @ np.linalg.inv(P)), 'cumulative_motion_mm': r['cumulative_motion_mm'], 'sensitivity_budget_mm': u, 'all_numerical_gates_pass': bool(good), 'conditional_weekly_motion_resolved': bool(good and np.linalg.norm(v) > u), 'absolute_movement': 'UNKNOWN_REGION_STABILITY', 'budget_kind': 'PHENOMENOLOGICAL two-scan floor plus operational worst-sum, notCI', 'plan_absolute_comparison': 'UNKNOWN_PLAN_SKULL_FRAME'})
res = {'round': 'R5', 'claim_type': 'information_link', 'external_referent': json.loads((ROOT / 'PREREG_R5.json').read_text())['external_referent'], 'region_cases': len(regs), 'region_gate_passes': sum((x['gate_pass'] for x in regs)), 'pose_rows': len(out), 'weekly_rows': len(weekly), 'weekly_numerical_passes': sum((x['all_numerical_gates_pass'] for x in weekly)), 'conditional_weekly_resolved': sum((x['conditional_weekly_motion_resolved'] for x in weekly)), 'absolute_anatomical_movement': 'UNKNOWN', 'failed_per_tooth_split': sum((not x['per_tooth_split_gate_pass'] for x in out)), 'elapsed_s': time.monotonic() - start, 'peak_rss_mib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'by_patient': {p: {'weekly_rows': sum((x['patient'] == p for x in weekly)), 'median_weekly_mm': float(np.median([x['weekly_norm_mm'] for x in weekly if x['patient'] == p])), 'max_cumulative_mm': max((x['cumulative_motion_mm'] for x in weekly if x['patient'] == p)), 'max_split_mm': max((x['reference_split_error_mm'] for x in out if x['patient'] == p))} for p in ['3485', '6457']}}
for (n, d) in [('R5_REGIONS.json', regs), ('R5_POSES.json', out), ('R5_WEEKLY.json', weekly), ('R5_RESULTS.json', res)]:
    (ROOT / 'raw' / n).write_text(json.dumps(d, indent=2, allow_nan=False) + '\n')
with (ROOT / 'raw/upper_region_referenced_movement.csv').open('w') as f:
    flat = [{k: v for (k, v) in r.items() if k != 'weekly_vector_mm'} for r in weekly]
    w = csv.DictWriter(f, fieldnames=list(flat[0]))
    w.writeheader()
    w.writerows(flat)
print(json.dumps(res, indent=2))
state = json.loads((ROOT / 'CURRENT_WORK_STATE.json').read_text())
state.update(status='R5_MEASURED', latest_gate='Candidate region weekly poses measured; absolute anatomical status UNKNOWN', next_operation='Package and rerun full demo; next independent repeat-scan plus force acquisition', updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
(ROOT / 'CURRENT_WORK_STATE.json').write_text(json.dumps(state, indent=2) + '\n')
