from common import *
from round2 import load_item, simulated_calibration
from geometry import basis, roof_checks, ball_access, shell, export_stl
from regional import inverse
from plugin import generate
from round01 import full_control
import time

def run():
    start = time.perf_counter()
    pending = []
    cache = {}
    for (i, r0) in enumerate(read(ROOT / 'PREREG_R4.json')['cohort']):
        (t, _, _) = load_item(r0)
        r = dict(r0)
        r['uid'] = 'R4_' + r['case_key'] + '_' + r['level'] + '_' + r['participant']
        if t['status'] != 'READY':
            r['status'] = 'UNKNOWN_SITE'
            pending.append(r)
            continue
        inner = t['preparation_z'] + t['requirements']['film_min_mm']
        z = inner + t['requirements']['wall_mm'] + 0.05
        B = basis(t['xy'])
        (cal, probes, A, C, P, f0) = simulated_calibration(t, z, inner, B)
        tt = dict(t, regional_calibration=cal, actuation_basis=B, unadjusted_z=z, inner_z=inner, geometry_sha256=cal['geometry_sha256'], basis_sha256=cal['basis_sha256'])
        design = generate(tt)
        r['status'] = design['status']
        r['diagnostics'] = design.get('diagnostics')
        if design['status'] == 'DESIGN':
            final = np.array(design['outer_vertices'])[:, 2]
            path = DATA / (r['uid'] + '.npz')
            np.savez_compressed(path, xy=t['xy'], faces=t['faces'], unadjusted_z=z, inner=inner, adjusted_z=final, basis=B)
            r.update(artifact_path=str(path), artifact_sha256=sha(path), planned_relief_mean_um=design['diagnostics']['inverse']['objective_mean_removal_mm'] * 1000, planned_relief_max_um=float(np.max(z - final) * 1000), force_model='SIMULATED_REFITTED_PLANE_PROBES', physical_force='UNKNOWN')
            dest = ROOT / 'exports' / r['uid']
            dest.mkdir(parents=True, exist_ok=True)
            for (tag, top) in [('unadjusted', z), ('adjusted', final)]:
                (v, f) = shell(t['xy'], t['faces'], top, inner)
                export_stl(dest / (tag + '_roof.stl'), v, f)
            dump(dest / 'GenCAD_design.json', design)
            cache[r['uid']] = (t, z, inner, B, A, C, P, f0, final)
        pending.append(r)
    if not (ROOT / 'FROZEN_PREDICTIONS_R4.json').exists():
        freeze(ROOT / 'FROZEN_PREDICTIONS_R4.json', dict(predictions=pending, private_source_contact_queries_before_freeze=0, physical_measurement='NOT_RUN'))
    else:
        dump(ROOT / 'raw/FROZEN_PREDICTIONS_R4_REPLAY.json', dict(predictions=pending, frozen_utc=now(), physical_measurement='NOT_RUN'))
    from review_integrity import assert_frozen_payload
    assert_frozen_payload(ROOT, 'FROZEN_PREDICTIONS_R4.json', 'predictions', clean(pending))
    q = v6_contact()
    controls = []
    for r in pending:
        if r['uid'] not in cache:
            continue
        (t, z, inner, B, A, C, P, f0, final) = cache[r['uid']]
        a = np.array(r['diagnostics']['inverse']['coefficients_mm'])
        truth = full_control(C, f0, P @ a, A, A.T @ f0)
        f = np.array(r['diagnostics']['inverse']['at_selected_height']['force_N'])
        bounds = np.array(r['diagnostics']['inverse']['at_selected_height']['force_interval_N'])
        with np.load(V4 / 'payload/private/references' / f"{r['case_key']}.npz", allow_pickle=False) as ref:
            rg = t['ceiling'] - ref[r['family']]
        before = q.compare(t['xy'], t['faces'], t['ceiling'] - z, rg)
        after = q.compare(t['xy'], t['faces'], t['ceiling'] - final, rg)
        change = q.compare(t['xy'], t['faces'], t['ceiling'] - final, t['ceiling'] - z)
        r.update(contact_before=before, contact_after=after, contact_change=change, contact_retention='PASS' if after['status'] == 'SCORED' and after['symdiff_mm2'] <= 1 + 1e-09 and (after['coverage'] >= 0.95) else 'FAIL_OR_UNKNOWN', full_crown_eligibility='UNKNOWN_AXIAL_ANATOMY_AND_PHYSICAL_SUPPORT', chairside_adjustment_saved_um=None)
        controls.append(dict(uid=r['uid'], full_force_QP_max_error_N=float(abs(f - truth).max()), enclosed=bool(np.all(truth >= bounds[:, 0]) and np.all(truth <= bounds[:, 1])), injected_force_plus3_rejected=bool(truth[0] + 3 > bounds[0, 1]), milling_plane_pass=r['diagnostics']['milling']['status'] == 'CONDITIONAL_IDEAL_TOOL_PASS'))
    full = [r for r in pending if r.get('diagnostics', {}).get('manufacturing_status') == 'CONDITIONAL_ROOF_ONLY']
    summary = dict(requested=len(pending), designs=sum((r['status'] == 'DESIGN' for r in pending)), ideal_tool_certificates=sum((c['milling_plane_pass'] for c in controls)), force_wall_antagonist_ideal_tool_passes=len(full), also_native_contact_passes=sum((r['contact_retention'] == 'PASS' for r in full)), physical_validated_crowns=0, seconds=time.perf_counter() - start)
    dump(ROOT / 'raw/R4_ROWS.json', pending)
    dump(ROOT / 'raw/R4_CONTROLS.json', controls)
    dump(ROOT / 'raw/R4_SUMMARY.json', summary)
    state('R4_DECIDED', str(summary), 'Prepare actual full-crown transfer refusal and10-crown measurement protocol; seal replay')
    print(json.dumps(summary))
if __name__ == '__main__':
    run()
