from common import *
import time, resource, csv

def rel(ix, p, j, src, t, f, a):
    Hi = np.array(ix[p, j, src, t, f]['H']['all_crowns'])
    Ha = np.array(ix[p, j, src, t, a]['H']['all_crowns'])
    return np.linalg.inv(Ha) @ Hi

def rms(v):
    return float(np.sqrt(np.mean(np.sum(v * v, axis=1))))

def point_radius(ix, p, j, src, t, f, a, x, eps=0.1, theta=1.0):
    Hi = np.array(ix[p, j, src, t, f]['H']['all_crowns'])
    Ha = np.array(ix[p, j, src, t, a]['H']['all_crowns'])
    ci = np.array(ix[p, j, src, 0, f]['centroid_base_mm'])
    ca = np.array(ix[p, j, src, 0, a]['centroid_base_mm'])
    chord = 2 * np.sin(np.deg2rad(theta) / 2)
    return 2 * eps + chord * np.linalg.norm(x - ci, axis=1) + chord * np.linalg.norm(apply(Hi, x) - apply(Ha, ca), axis=1)

def fraction_interval(A, P, etaA, etaP):
    (a, p) = (rms(A), rms(P))
    if p <= etaP:
        return None
    n = float(np.mean(np.sum(A * P, axis=1)))
    dn = etaA * p + etaP * a + etaA * etaP
    D = [(p - etaP) ** 2, (p + etaP) ** 2]
    v = [q / d for q in [n - dn, n + dn] for d in D]
    return [float(min(v)), float(max(v))]

def direct_relative_point(ix, p, j, src, t, f, a, x):
    i = np.array(ix[p, j, src, t, f]['H']['all_crowns'])
    h = np.array(ix[p, j, src, t, a]['H']['all_crowns'])
    return (apply(i, x) - h[:3, 3]) @ h[:3, :3]

def run():
    start = time.monotonic()
    (rows, ix) = pose_index()
    reg = read(ROOT / 'PREREG_R2.json')
    out = []
    maxg = 0.0
    maxdirect = 0.0
    maxrot = 0.0
    rng = np.random.default_rng(4802)
    for p in reg['selection']['patients']:
        for (j, anchors) in reg['selection']['anchors'].items():
            ids = sorted({r['fdi'] for r in rows if r['patient'] == p and r['jaw'] == j})
            for f in ids:
                cloud = np.load(CACHE / f'{p}_Sirona_T0_Z{f}.npz')
                x = cloud['points'][:reg['selection']['point_samples']].astype(float)
                c = cloud['centroid']
                ax = axes(ix, p, j, f)
                for a in anchors:
                    if a == f:
                        continue
                    for t in reg['selection']['weeks']:
                        for (interval, s) in [('weekly', t - 1), ('cumulative', 0)]:
                            (ast, asp, pst, psp) = [rel(ix, p, j, src, w, f, a) for (src, w) in [('Sirona', t), ('Sirona', s), ('Bottmedical', t), ('Bottmedical', s)]]
                            ua = apply(ast, x) - apply(asp, x)
                            up = apply(pst, x) - apply(psp, x)
                            gap = rms(ua - up)
                            P = rms(up)
                            A = rms(ua)
                            (ad, at, aw) = delta(ast, asp, c, ax)
                            (pd, pt, pw) = delta(pst, psp, c, ax)
                            Da = ast[:3, :3] @ asp[:3, :3].T
                            Dp = pst[:3, :3] @ psp[:3, :3].T
                            thetaA = float(np.degrees(Rotation.from_matrix(Da).magnitude()))
                            thetaP = float(np.degrees(Rotation.from_matrix(Dp).magnitude()))
                            dg = float(np.degrees(Rotation.from_matrix(Da @ Dp.T).magnitude()))
                            errA = point_radius(ix, p, j, 'Sirona', t, f, a, x) + point_radius(ix, p, j, 'Sirona', s, f, a, x)
                            errP = point_radius(ix, p, j, 'Bottmedical', t, f, a, x) + point_radius(ix, p, j, 'Bottmedical', s, f, a, x)
                            etaA = float(np.sqrt(np.mean(errA ** 2)))
                            etaP = float(np.sqrt(np.mean(errP ** 2)))
                            u = rms(ua - up)
                            fraction = float(np.mean(np.sum(ua * up, axis=1)) / (P * P)) if P > 1e-12 else None
                            rotfrac = float(aw @ pw / (pw @ pw)) if pw @ pw > 1e-12 else None
                            magnitude_interval = [max(0.0, thetaA - 4) / (thetaP + 4), (thetaA + 4) / (thetaP - 4)] if thetaP > 4 else None
                            conditional_lag = thetaA + 4 < thetaP - 4
                            directA = direct_relative_point(ix, p, j, 'Sirona', t, f, a, x) - direct_relative_point(ix, p, j, 'Sirona', s, f, a, x)
                            directP = direct_relative_point(ix, p, j, 'Bottmedical', t, f, a, x) - direct_relative_point(ix, p, j, 'Bottmedical', s, f, a, x)
                            de = float(max(np.max(np.abs(directA - ua)), np.max(np.abs(directP - up))))
                            maxdirect = max(maxdirect, de)
                            mats = [np.array(ix[p, j, src, w, k]['H']['all_crowns'])[:3, :3] for (src, w) in [('Sirona', t), ('Sirona', s), ('Bottmedical', t), ('Bottmedical', s)] for k in [f, a]]
                            (R1, R2, R3, R4) = [mats[k + 1].T @ mats[k] for k in range(0, 8, 2)]
                            directgap = float(np.degrees(Rotation.from_matrix(R1 @ R2.T @ (R3 @ R4.T).T).magnitude()))
                            maxrot = max(maxrot, abs(directgap - dg))
                            qe = 0.0
                            for (src, w) in [('Sirona', t), ('Sirona', s), ('Bottmedical', t), ('Bottmedical', s)]:
                                G = np.eye(4)
                                G[:3, :3] = Rotation.from_rotvec(rng.normal(size=3)).as_matrix()
                                G[:3, 3] = rng.normal(size=3) * 100
                                Hi = G @ np.array(ix[p, j, src, w, f]['H']['all_crowns'])
                                Ha = G @ np.array(ix[p, j, src, w, a]['H']['all_crowns'])
                                qr = np.linalg.inv(Ha) @ Hi
                                qe = max(qe, float(np.max(np.abs(apply(qr, x) - apply(rel(ix, p, j, src, w, f, a), x)))))
                            maxg = max(maxg, qe)
                            signed_error = {k: ad[k] - pd[k] for k in ad}
                            out.append(dict(patient=p, jaw=j, fdi=f, anchor_fdi=a, week=t, interval=interval, resolution='PER_TOOTH', timescale='HANDOVER_WEEK', reference='Named moving central incisor; no fixed anatomy assumption', observed_centroid_vector_mm=at.tolist(), planned_centroid_vector_mm=pt.tolist(), observed_rotation_vector_deg=aw.tolist(), planned_rotation_vector_deg=pw.tolist(), signed_components={'achieved': ad, 'planned': pd, 'error': signed_error}, observed_angular_deg=thetaA, planned_angular_deg=thetaP, angular_tracking_gap_deg=dg, conditional_angular_gap_lower_deg=max(0.0, dg - 8), conditional_angular_tracking_rejected=bool(dg > 8), conditional_angular_lag=bool(conditional_lag), angular_magnitude_ratio=thetaA / thetaP if thetaP > 1e-12 else None, conditional_angular_magnitude_ratio_interval=magnitude_interval, signed_angular_projection_fraction=rotfrac, signed_angular_component_bound='UNKNOWN; not inferred from scalar SO3 magnitude bound', point_tracking_rms_mm=gap, observed_rms_mm=A, planned_rms_mm=P, signed_field_completion_fraction=fraction, conditional_completion_interval=fraction_interval(ua, up, etaA, etaP), conditional_field_tracking_rejected=bool(gap > etaA + etaP), observed_field_radius_mm=etaA, planned_field_radius_mm=etaP, pose_radii=dict(centroid_mm=0.1, rotation_deg=1.0, resolution='PHENOMENOLOGICAL', empirical_bound=False), required_uniform_orientation_radius_deg_for_tracking_rejection=dg / 8, required_uniform_orientation_radius_deg_for_lag=max(0.0, (thetaP - thetaA) / 8), gauge_point_error_mm=qe, direct_point_error_mm=de, direct_rotation_error_deg=abs(directgap - dg), wrong_point_plus1mm_rejected=bool(abs(de + 1) > 1e-09), wrong_rotation_plus10deg_rejected=bool(abs(directgap + 10 - dg) > 1e-09), absolute_motion='UNKNOWN_NO_FIXED_REFERENCE', bodily='UNKNOWN_ROOT_NOT_OBSERVED'))
    targets = [r for r in out if r['jaw'] == 'OK' and r['fdi'] in [17, 27] and (r['week'] == 8) and (r['interval'] == 'cumulative')]
    target_pass = {f'{p}_{f}': all((r['conditional_angular_lag'] for r in targets if r['patient'] == p and r['fdi'] == f)) for p in ['3485', '6457'] for f in [17, 27]}
    condrot = sum((r['conditional_angular_tracking_rejected'] for r in out))
    condlag = sum((r['conditional_angular_lag'] for r in out))
    result = dict(round='R2', claim_type='information_link', external_referent=reg['external_referent'], rows=len(out), patient_count=2, tooth_count=52, pair_count=len(out) // 18, conditional_angular_tracking_rejections=condrot, conditional_angular_lag_rows=condlag, conditional_field_tracking_rejections=sum((r['conditional_field_tracking_rejected'] for r in out)), conditional_completion_intervals=sum((r['conditional_completion_interval'] is not None for r in out)), unidentified_completion_denominators=sum((r['conditional_completion_interval'] is None for r in out)), paper_four_target_conditional_lag=target_pass, paper_target_gate='PASS' if all(target_pass.values()) else 'FAIL', gauge_error_mm=maxg, direct_point_control_error_mm=maxdirect, direct_rotation_control_error_deg=maxrot, physical_pose_error_bound='UNKNOWN_NO_REPEAT_PATIENT_SCANS', bounds='rigorous finite-rotation conditional enclosures under declared .1 mm /1deg per-pose radii; input radii PHENOMENOLOGICAL, not validated', wall_s=time.monotonic() - start, peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    write(ROOT / 'raw/R2_RELATIVE_POSE.json', out)
    write(ROOT / 'raw/R2_PAPER_MOLAR_COMPARISON.json', targets)
    write(ROOT / 'raw/R2_RESULTS.json', result)
    fields = ['patient', 'jaw', 'fdi', 'anchor_fdi', 'week', 'interval', 'resolution', 'timescale', 'observed_angular_deg', 'planned_angular_deg', 'angular_tracking_gap_deg', 'conditional_angular_tracking_rejected', 'conditional_angular_lag', 'signed_angular_projection_fraction', 'point_tracking_rms_mm', 'signed_field_completion_fraction', 'observed_field_radius_mm', 'planned_field_radius_mm', 'required_uniform_orientation_radius_deg_for_tracking_rejection', 'required_uniform_orientation_radius_deg_for_lag']
    with (ROOT / 'raw/relative_motion.csv').open('w') as o:
        w = csv.DictWriter(o, fieldnames=fields)
        w.writeheader()
        w.writerows([{k: r[k] for k in fields} for r in out])
    txt = f"R2: {len(out)} relative pair/interval rows, angular plan-tracking rejected in {condrot} under the conditional 1deg per-crown error model. Angular lag in {condlag}. Fixed anatomy never inferred. Paper four molars: {target_pass}, gate {result['paper_target_gate']}. Absolute patient uncertainty remains UNKNOWN.\nNext: compare frozen X19 same-patient geometric force scenarios to tooth-level relative lag; do not fit force-to-biology.\n"
    (ROOT / 'HANDOFF_R2.md').write_text(txt)
    (ROOT / 'HANDOFF.md').write_text(txt)
    state('R2_ADJUDICATED', result['paper_target_gate'], 'Compare frozen X19 geometric proxy with relative crown lag, report material mismatch and n=2 limits')
    print(result)
    return result
if __name__ == '__main__':
    run()
