from dental_release.paths import expand as _release_expand
from common import *
from regional import fit_response, predict, inverse
from geometry import attach, basis, roof_checks, ball_access, shell, export_stl
from round01 import full_control
import time, resource, hashlib, copy

def simulated_calibration(t, z, inner, B):
    xy = t['xy']
    c = xy.mean(0)
    span = max(np.ptp(xy[:, 0]), np.ptp(xy[:, 1]))
    (xlo, xhi) = (xy[:, 0].min(), xy[:, 0].max())
    points = np.array([[xlo + 0.25 * (xhi - xlo), c[1]], [xlo + 0.75 * (xhi - xlo), c[1]], c + [-span, -span], c + [-span, span], c + [span, -span], c + [span, span]])
    A = np.c_[np.ones(6), (points - c) / 30]
    C = np.eye(6) / 400
    K = 400 * (np.eye(6) - A @ np.linalg.inv(A.T @ A) @ A.T)
    P = np.zeros((6, 2))
    P[:2] = [[0.75, 0.25], [0.25, 0.75]]
    f0 = np.array([12.0, 8.0, 16.0, 16.0, 16.0, 16.0])
    J = K @ P
    s = 0.025
    ghash = hashlib.sha256(np.c_[xy, z, inner].astype('<f8').tobytes() + t['faces'].astype('<i8').tobytes()).hexdigest()
    bhash = hashlib.sha256(B.astype('<f8').tobytes()).hexdigest()
    probes = dict(schema='regional-actuation-probes-v1', units={'height': 'mm', 'force': 'N'}, case=t['case_key'], frame=t['frame'], geometry_sha256=ghash, basis_sha256=bhash, kind='SIMULATION', observation_operator='calibrated_regional_axial_reactions', linear_support_closure=True, matched_load_rate_history=True, height_error_mm=0.0, step_mm=s, raw_force_channel_bound_N=0.05, baseline_force_N=f0, plus_force_N=np.array([f0 + J[:, j] * s for j in range(2)]), minus_force_N=np.array([f0 - J[:, j] * s for j in range(2)]), support_positions_xy_mm=points, actuation_at_supports=P)
    model = fit_response(clean(probes))
    return (model, probes, A, C, P, f0)

def load_item(r):
    k = r['case_key']
    metas = read(V4 / 'payload/public/tasks' / f'{k}.json')
    i0 = next((i for (i, t) in enumerate(metas) if t['family'] == r['family'] and t['level'] == r['level']))
    m = metas[i0]
    with np.load(V4 / 'payload/public' / m['geometry_file'], allow_pickle=False) as a:
        t = attach(m, dict(a))
    if r['participant'] == 'X1B_morphology':
        path = V4 / 'payload/predictions/X1B_morphology' / f'{k}.npz'
    else:
        path = Path(_release_expand('@DENTAL_WORK_ROOT@/X60/predictions/volume_1.0')) / f'{k}.npz'
    with np.load(path, allow_pickle=False) as p:
        i = list(p['task_ids']).index(m['task_id'])
        if 'outer_' + str(i) not in p.files:
            return (t, None, None)
        z = p['outer_' + str(i)].copy()
        inn = p['inner_' + str(i)].copy()
    return (t, z, inn)

def verify_inputs():
    for (p, item) in read(ROOT / 'INPUT_LOCK_R2.json')['files'].items():
        if sha(p) != item['sha256']:
            raise ValueError('Input drift: ' + p)

def run():
    verify_inputs()
    pr = read(ROOT / 'PREREG_R2.json')
    start = time.perf_counter()
    rows = []
    cache = {}
    predictions = []
    DATA.mkdir(parents=True, exist_ok=True)
    for (idx, r0) in enumerate(pr['cohort']):
        r = dict(r0)
        r['uid'] = r['case_key'] + '_' + r['level'] + '_' + r['participant']
        (t, z, inner) = load_item(r)
        if z is None:
            r.update(status='NO_GENERATED_SURFACE', force_result=None)
            rows.append(r)
            continue
        B = basis(t['xy'])
        (model, probes, A, C, P, f0) = simulated_calibration(t, z, inner, B)
        force_only = inverse(model, B, t['weights'], z - inner - t['requirements']['wall_mm'], additional_geometry=False)
        constrained = inverse(model, B, t['weights'], z - inner - t['requirements']['wall_mm'] + 1e-06)
        r.update(baseline_geometry=roof_checks(t, z, inner), force_only=force_only, force_result=constrained, model=model, source_fdi='NOT_TRANSFERRED_FROM_X82', resolution='PER_SURFACE_REGION')
        a = np.asarray(force_only['coefficients_mm']) if force_only['status'] == 'CONDITIONAL_TARGET_DESIGN' else np.zeros(2)
        final = z + B @ a
        path = DATA / (r['uid'] + '.npz')
        np.savez_compressed(path, xy=t['xy'], faces=t['faces'], unadjusted_z=z, inner=inner, adjusted_z=final, basis=B)
        r.update(artifact_path=str(path), artifact_sha256=sha(path), force_only_geometry=roof_checks(t, final, inner), planned_relief_max_um=float(np.max(-B @ a) * 1000), planned_relief_mean_um=float(t['weights'] @ (-B @ a) / sum(t['weights']) * 1000), milling=ball_access(t['xy'], t['faces'], final))
        if constrained['status'] == 'CONDITIONAL_TARGET_DESIGN':
            aa = np.array(constrained['coefficients_mm'])
            zz = z + B @ aa
            r['adjusted_geometry'] = roof_checks(t, zz, inner)
            r['adjusted_milling'] = ball_access(t['xy'], t['faces'], zz)
            r['geometric_design_status'] = 'CONDITIONAL_ROOF_PASS' if r['adjusted_geometry']['wall_pass'] and r['adjusted_geometry']['antagonist_pass'] and r['adjusted_geometry']['intaglio_milling_pass'] and (r['adjusted_milling']['status'] == 'CONDITIONAL_IDEAL_TOOL_PASS') else 'UNKNOWN_OR_FAILED_GEOMETRY'
            r['design_coefficients_mm'] = aa
            np.savez_compressed(path, xy=t['xy'], faces=t['faces'], unadjusted_z=z, inner=inner, adjusted_z=zz, basis=B)
            r['artifact_sha256'] = sha(path)
        else:
            r['geometric_design_status'] = 'INFEASIBLE_OR_UNKNOWN'
        dump(ROOT / 'raw/probes' / f"{r['uid']}.json", probes)
        predictions.append(dict(uid=r['uid'], force_only=force_only, constrained=constrained, artifact_path=str(path), artifact_sha256=r['artifact_sha256'], geometry_sha256=model['geometry_sha256'], basis_sha256=model['basis_sha256']))
        cache[r['uid']] = (t, z, inner, B, model, A, C, P, f0)
        rows.append(r)
    fp = ROOT / 'FROZEN_PREDICTIONS.json'
    if fp.exists():
        predpath = ROOT / 'raw/FROZEN_PREDICTIONS_REPLAY.json'
        dump(predpath, dict(frozen_utc=now(), predictions=predictions, physical_measurement='NOT_RUN'))
    else:
        freeze(fp, dict(claim_type='capability', prereg_sha256=sha(ROOT / 'PREREG_R2.json'), predictions=predictions, physical_measurement='NOT_RUN', external_source_contact_queries_before_freeze=0, force_fixture='Own declared simulation; not external force facit'))
    from review_integrity import assert_frozen_payload
    assert_frozen_payload(ROOT, 'FROZEN_PREDICTIONS.json', 'predictions', clean(predictions))
    q = v6_contact()
    controls = []
    external = []
    for r in rows:
        if r['uid'] not in cache:
            continue
        (t, z, inner, B, model, A, C, P, f0) = cache[r['uid']]
        w = A.T @ f0
        held = []
        acts = [np.array([-0.017, -0.009]), np.array([-0.006, -0.022]), np.array([0.008, -0.004]), np.array([-0.019, -0.003])]
        if r['force_only']['status'] == 'CONDITIONAL_TARGET_DESIGN':
            acts.append(np.array(r['force_only']['coefficients_mm']))
        for a in acts:
            truth = full_control(C, f0, P @ a, A, w)
            p = predict(model, a)
            lohi = np.array(p['force_interval_N'])
            nom = np.array(p['force_N'])
            held.append(dict(coefficients_mm=a, QP_force_N=truth, maximum_nominal_error_N=float(abs(truth - nom).max()), enclosed=bool(np.all(truth >= lohi[:, 0] - 1e-07) and np.all(truth <= lohi[:, 1] + 1e-07)), injected_plus3_N_rejected=bool(np.any(truth + np.r_[3.0, np.zeros(5)] > lohi[:, 1]))))
        wrong = copy.deepcopy(model)
        wrong['status'] = 'UNKNOWN'
        wrong['reason'] = 'Missing calibrated force observation'
        c = dict(uid=r['uid'], heldouts=held, missing_observation_rejected=inverse(wrong, B, t['weights'], z - inner - 0.5)['status'] == 'UNKNOWN', wall_fault_rejected=not roof_checks(t, np.minimum(z, inner + t['requirements']['wall_mm'] - 0.2), inner)['wall_pass'])
        with np.load(r['artifact_path'], allow_pickle=False) as a:
            final = a['adjusted_z']
        with np.load(V4 / 'payload/private/references' / f"{r['case_key']}.npz", allow_pickle=False) as ref:
            rg = t['ceiling'] - ref[r['family']]
            pg = t['ceiling'] - z
            ag = t['ceiling'] - final
        before = q.compare(t['xy'], t['faces'], pg, rg)
        after = q.compare(t['xy'], t['faces'], ag, rg)
        changed = q.compare(t['xy'], t['faces'], ag, pg)
        parity = q.Q.contact_map(t['xy'], t['faces'], ag, rg)
        c['unchanged_v4_polygon_parity_error_mm2'] = abs(parity['contact_symdiff_mm2'] - after['symdiff_mm2']) if after['status'] == 'SCORED' else None
        if after['status'] == 'SCORED':
            c['contact_fault_rejected'] = abs(after['symdiff_mm2'] + 3 - parity['contact_symdiff_mm2']) > 1e-07
        controls.append(c)
        r.update(status='EVALUATED', contact_before=before, contact_after=after, contact_change=changed, contact_preservation='PASS' if after['status'] == 'SCORED' and after['symdiff_mm2'] <= 1 + 1e-09 and (after['coverage'] >= 0.95) else 'FAIL_OR_UNKNOWN', physical_force='UNKNOWN_NO_MATCHED_MEASUREMENTS', chairside_adjustment_saved_um=None, model_additional_adjustment_um=0.0 if r['force_result']['status'] == 'CONDITIONAL_TARGET_DESIGN' else None)
        if r['force_result']['status'] == 'CONDITIONAL_TARGET_DESIGN':
            (v, f) = shell(t['xy'], t['faces'], final, inner)
            dest = ROOT / 'exports' / r['uid']
            dest.mkdir(parents=True, exist_ok=True)
            export_stl(dest / 'adjusted_roof.stl', v, f)
            (v, f) = shell(t['xy'], t['faces'], z, inner)
            export_stl(dest / 'unadjusted_roof.stl', v, f)
        external.append(dict(uid=r['uid'], native_area_mm2=after.get('reference', {}).get('area_mm2'), before_symdiff_mm2=before.get('symdiff_mm2'), after_symdiff_mm2=after.get('symdiff_mm2'), compared_quantity='Projected contact polygons against registered source IOS, not force'))
    xy = np.array([[-1.0, 0.0], [0.0, 0.0], [1.0, 0.0], [-1.0, 1.0], [0.0, 1.0], [1.0, 1.0]])
    ff = np.array([[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]])
    zz = np.array([1.0, 0.0, 1.0, 1.0, 0.0, 1.0])
    toolfault = ball_access(xy, ff, zz)
    controls.append(dict(uid='INJECTED_CONCAVE_MILLING', false_milling_pass_rejected=toolfault['status'] != 'CONDITIONAL_IDEAL_TOOL_PASS', witness=toolfault))
    dump(ROOT / 'raw/R2_ROWS.json', rows)
    dump(ROOT / 'raw/R2_CONTROLS.json', controls)
    dump(ROOT / 'raw/R2_EXTERNAL_GEOMETRY.json', external)
    summary = dict(requested=len(rows), generated=sum((r.get('status') == 'EVALUATED' for r in rows)), force_only_targets=sum(((r.get('force_only') or {}).get('status') == 'CONDITIONAL_TARGET_DESIGN' for r in rows)), wall_constrained_targets=sum(((r.get('force_result') or {}).get('status') == 'CONDITIONAL_TARGET_DESIGN' for r in rows)), full_conditional_roof_passes=sum((r.get('geometric_design_status') == 'CONDITIONAL_ROOF_PASS' for r in rows)), contact_preservation_passes=sum((r.get('contact_preservation') == 'PASS' for r in rows)), physical_force_validated=0, seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(ROOT / 'raw/R2_SUMMARY.json', summary)
    state('R2_DECIDED', str(summary), 'Change joint geometry design at the active wall/access obstruction; preserve all failures')
    (ROOT / 'history/HANDOFF_R2.md').write_text('R2 decided. ' + json.dumps(summary) + '\nRegional inverse was calibrated on SIMULATED probes, not specimen data. Continuous wall and ideal exterior tool guards remain separate; source contact is independently scored. Next construction: quantify minimum extra uniform preparation that makes the force target compatible with the fixed roof wall. Complete crown CAM/force still UNKNOWN.\n')
    attempts = read(ROOT / 'ATTEMPTS.json')
    attempts.append(dict(round='R2', outcome=summary, evidence='raw/R2_ROWS.json', next_operation='Force-target minimum-preparation co-design, with existing tool failures retained'))
    dump(ROOT / 'ATTEMPTS.json', attempts)
    print(json.dumps(clean(summary)))
if __name__ == '__main__':
    run()
