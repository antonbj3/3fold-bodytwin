from common import *
from run_r2 import rel, point_radius
import time, csv, resource

def wrap(x):
    return (x + 180.0) % 360.0 - 180.0

def projection(Q, unit, plane_x, plane_y, delta_deg=2.0):
    v = Q[:3, :3] @ unit
    x = float(v @ plane_x)
    y = float(v @ plane_y)
    rho = float(np.hypot(x, y))
    theta = float(np.degrees(np.arctan2(x, y)))
    chord = 2 * np.sin(np.deg2rad(delta_deg) / 2)
    error = float(np.degrees(np.arcsin(chord / rho))) if rho > chord else None
    return (theta, error, rho)

def type_observations(Q1, Q0, ax):
    e = -ax['intrusion']
    m = ax['mesial']
    b = ax['buccal']
    out = {}
    for (name, u, x, y) in [('buccolingual_tip', e, b, e), ('mesiodistal_tip', e, m, e), ('axial_rotation', m, b, m)]:
        (t1, e1, r1) = projection(Q1, u, x, y)
        (t0, e0, r0) = projection(Q0, u, x, y)
        d = wrap(t1 - t0)
        err = e1 + e0 if e1 is not None and e0 is not None else None
        if err is not None and abs(d) + err >= 180.0:
            err = None
        out[name] = dict(value=d, radius=err, min_projection_norm=min(r1, r0), unit='deg')
    return out

def control():
    rng = np.random.default_rng(484)
    checks = []
    maximum = 0.0
    rejected = False
    e = np.array([0.0, 1.0, 0.0])
    x = np.array([1.0, 0.0, 0.0])
    y = e
    for i in range(256):
        R = Rotation.from_rotvec(rng.normal(size=3)).as_matrix()
        Q = np.eye(4)
        Q[:3, :3] = R
        (t, err, rho) = projection(Q, e, x, y)
        v = Rotation.from_matrix(R).apply(e)
        direct = float(np.degrees(np.arctan2(v @ x, v @ y)))
        d = abs(wrap(t - direct))
        maximum = max(maximum, d)
        rejected = rejected or abs(wrap(t + 10 - direct)) > 1e-09
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        theta = rng.uniform(0, 2)
        Rp = Rotation.from_rotvec(direction * np.deg2rad(theta)).as_matrix() @ R
        P = np.eye(4)
        P[:3, :3] = Rp
        (tp, ep, rhop) = projection(P, e, x, y)
        within = True if err is None else abs(wrap(tp - t)) <= err + 1e-10
        checks.append(within)
    return dict(cases=len(checks), finite_chord_disk_enclosure_pass=all(checks), direct_quaternion_axis_error_deg=maximum, wrong_plus10deg_rejected=rejected)

def run():
    start = time.monotonic()
    (rows, ix) = pose_index()
    reg = read(ROOT / 'PREREG_R4.json')
    out = []
    for p in reg['selection']['patients']:
        for (j, anchors) in reg['selection']['anchors'].items():
            ids = sorted({r['fdi'] for r in rows if r['patient'] == p and r['jaw'] == j})
            for f in ids:
                ax = axes(ix, p, j, f)
                c = np.array(ix[p, j, 'Sirona', 0, f]['centroid_base_mm'])
                for a in anchors:
                    if a == f:
                        continue
                    for t in reg['selection']['weeks']:
                        for (interval, s) in [('weekly', t - 1), ('cumulative', 0)]:
                            (A1, A0, P1, P0) = [rel(ix, p, j, src, w, f, a) for (src, w) in [('Sirona', t), ('Sirona', s), ('Bottmedical', t), ('Bottmedical', s)]]
                            ao = type_observations(A1, A0, ax)
                            po = type_observations(P1, P0, ax)
                            e = -ax['intrusion']
                            av = float((apply(A1, c) - apply(A0, c)) @ e)
                            pv = float((apply(P1, c) - apply(P0, c)) @ e)
                            ua = float(point_radius(ix, p, j, 'Sirona', t, f, a, c[None, :])[0] + point_radius(ix, p, j, 'Sirona', s, f, a, c[None, :])[0])
                            up = float(point_radius(ix, p, j, 'Bottmedical', t, f, a, c[None, :])[0] + point_radius(ix, p, j, 'Bottmedical', s, f, a, c[None, :])[0])
                            ao['vertical'] = dict(value=av, radius=ua, unit='mm')
                            po['vertical'] = dict(value=pv, radius=up, unit='mm')
                            for k in ao:
                                (aa, pp) = (ao[k], po[k])
                                (uA, uP) = (aa['radius'], pp['radius'])
                                (a0, p0) = (aa['value'], pp['value'])
                                ratio = ratio_box(a0, p0, uA, uP) if uA is not None and uP is not None else None
                                intended = 1 if p0 >= 0 else -1
                                gap_reject = False
                                lag = False
                                if uA is not None and uP is not None:
                                    gap_reject = bool(abs(a0 - p0) > uA + uP)
                                    lag = bool(intended * p0 - uP > 0 and intended * a0 + uA < intended * p0 - uP)
                                name = k if k != 'vertical' else 'extrusion' if p0 >= 0 else 'intrusion'
                                out.append(dict(patient=p, jaw=j, fdi=f, anchor_fdi=a, week=t, interval=interval, motion=name, observation=k, resolution='PER_TOOTH', timescale='HANDOVER_WEEK', achieved=a0, planned=p0, unit=aa['unit'], achieved_radius=uA, planned_radius=uP, conditional_ratio_interval=ratio, signed_ratio=a0 / p0 if abs(p0) > 1e-12 else None, conditional_tracking_rejected=gap_reject, conditional_lag=lag, projection_condition='PASS' if uA is not None and uP is not None else 'UNKNOWN_PROJECTION_OR_WRAP', radius_resolution='PHENOMENOLOGICAL', input_radii_empirically_validated=False, quantity='transported crown coordinate-axis projected angle / relative crown-centroid vertical translation', fixed_reference=False, body='UNKNOWN_ROOT_NOT_OBSERVED'))
    numerical = control()
    assert numerical['finite_chord_disk_enclosure_pass'] and numerical['direct_quaternion_axis_error_deg'] <= 1e-09 and numerical['wrong_plus10deg_rejected']
    res = dict(round='R4', claim_type='information_link', external_referent=reg['external_referent'], component_rows=len(out), conditional_ratio_intervals=sum((r['conditional_ratio_interval'] is not None for r in out)), unidentified_ratios=sum((r['conditional_ratio_interval'] is None for r in out)), conditional_tracking_rejections=sum((r['conditional_tracking_rejected'] for r in out)), conditional_lag_rows=sum((r['conditional_lag'] for r in out)), by_type={k: dict(rows=sum((r['motion'] == k for r in out)), conditional_ratio_intervals=sum((r['motion'] == k and r['conditional_ratio_interval'] is not None for r in out)), conditional_tracking_rejections=sum((r['motion'] == k and r['conditional_tracking_rejected'] for r in out))) for k in sorted({r['motion'] for r in out})}, controls=numerical, absolute_type_accuracy='UNKNOWN_NO_FIXED_REFERENCE_OR_MATCHED_CLINICAL_AXIS', paper_total_vs_tip_correction='R2 total-angular gate is a cross-quantity challenge, not a direct replication/refutation of source buccolingual crown-tip result; R4 names the projected-axis quantity explicitly', wall_s=time.monotonic() - start, peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    write(ROOT / 'raw/R4_TYPE_ENCLOSURES.json', out)
    write(ROOT / 'raw/R4_RESULTS.json', res)
    with (ROOT / 'raw/type_enclosures.csv').open('w') as o:
        fields = list(out[0])
        w = csv.DictWriter(o, fieldnames=fields)
        w.writeheader()
        w.writerows(out)
    txt = f"R4: {len(out)} signed type rows, {res['conditional_ratio_intervals']} nonlinear conditional intervals; input pose radii still PHENOMENOLOGICAL. Absolute type accuracy UNKNOWN. Transported crown-axis angle is not clinical FACC/root torque. {numerical}. R2 total SO3 lag must not be advertised as a direct replication of absolute source buccolingual tip.\nNext construction: replace assumed pose radii by independent same-week repeat scans against a fixed oriented reference; then matched Naturaligner wrench/time-history and heldout patient outcome.\n"
    (ROOT / 'HANDOFF_R4.md').write_text(txt)
    (ROOT / 'HANDOFF.md').write_text(txt)
    state('R4_ADJUDICATED', 'Conditional nonlinear type enclosures computed; anatomical accuracy UNKNOWN', 'Verify/package replay; next round requires new fixed-reference repeat-scan and material/wrench measurements')
    print(res)
    return res
if __name__ == '__main__':
    run()
