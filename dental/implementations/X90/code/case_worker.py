"""One case, one process, immutable package functions, source-owned arrays only."""
import sys
import os
sys.dont_write_bytecode = True
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS']:
    os.environ[key] = '4'
os.environ['MPLBACKEND'] = 'Agg'
import time
import datetime
import resource
import argparse
from pathlib import Path
import numpy as np
import common as local_common
from common import read, dump, sha, freeze, decision, stable, validate_manifest, validate_decisions, validate_force_binding, validate_pulp_observation
from operators import load

def contracts(package, edge_ids):
    net = read(package / 'batch23/demos/NET_R4/CONSTRAINT_NET_DENTAL_R4.json')
    cs = read(package / 'batch23/demos/NET_R4/sources/r3_contracts.json')['edges']
    rows = []
    for edge in edge_ids:
        e = next((x for x in net['edges'] if x['id'] == edge))
        c = cs[edge]
        rows.append(dict(edge=edge, resolution=e['resolution_level'], timescale=e['timescale'], measurement=c.get('measurement', c.get('required_operation')), measurement_ids=c.get('measurement_ids', []), would_decide=c.get('would_decide', c.get('source_required_measurement')), locator='NET_R4/CONSTRAINT_NET_DENTAL_R4.json#/edges/' + edge))
    return rows

def run(args):
    started = time.perf_counter()
    package = Path(args.package).resolve()
    out = Path(args.out).resolve()
    m = validate_manifest(read(args.manifest))
    manifest_root = Path(args.manifest).resolve().parent

    def input_path(path):
        q = Path(path)
        return q if q.is_absolute() else manifest_root / q
    out.mkdir(parents=True, exist_ok=True)
    epoch = read(args.epoch)
    (geom, gen, adapt, source, check, force, pulp, nerve, cap, force_bridge) = load(package, args.nerve)
    source.ZIP = Path(args.dataset_root) / 'Bits2Bites/Bits2Bites_v01.zip'
    source.LABEL = Path(args.labels_root) / 'Bits2Bites'
    (arches, sources) = source.pair(int(m['case']))
    meta = read(source.LABEL / m['case'] / 'labels+landmarks.json')
    fr = meta['arches']['lower']['frame']
    rotation = np.column_stack([fr['right_unit'], fr['anterior_unit'], fr['superior_unit']])
    center = np.array(fr['center_mm'])
    if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-10) or np.linalg.det(rotation) <= 0:
        raise ValueError('INVALID_SOURCE_FRAME')
    for jaw in arches.values():
        jaw['v'] = (jaw['v'] - center) @ rotation
        jaw['tri'] = jaw['v'][jaw['f']]
    geometry_hash = __import__('hashlib').sha256(''.join((sources[j]['source_sha256'] for j in ['upper', 'lower'])).encode()).hexdigest()
    src = [{k: v for (k, v) in s.items() if k != 'label_file'} for s in sources.values()]
    source_locators = ['Bits2Bites_v01.zip::' + s['member'] + '#sha256=' + s['source_sha256'] for s in src]
    fdi = m['fdi']
    low = arches['lower']
    (target, drop, components) = geom.component(low['tri'][low['owner'] == fdi])
    (donor, donor_drop, donor_components) = geom.component(low['tri'][low['owner'] == fdi + 10])
    z = float(np.quantile(target.vertices[:, 2], 0.25))
    dz = float(np.quantile(donor.vertices[:, 2], 0.25))
    latest_failure = None
    proxy_control = None
    try:
        (closed_target, target_cap_info) = cap(target.copy())
        (closed_donor, donor_cap_info) = cap(donor.copy())
        target_cut = geom.cut(closed_target, z)
        donor_cut = geom.cut(closed_donor, dz)
        (target_cut, cut_drop, cut_components) = geom.component(target_cut.triangles)
        (donor_cut, donor_cut_drop, donor_cut_components) = geom.component(donor_cut.triangles)
        dump(out / 'PREPARATION_DIAGNOSTICS.json', dict(target_cap=target_cap_info, donor_cap=donor_cap_info, source_cut_components=cut_components, discarded_cut_area_fraction=cut_drop, donor_cut_components=donor_cut_components, discarded_donor_cut_area_fraction=donor_cut_drop, target_loops=len(geom.loops(target_cut)), donor_loops=len(geom.loops(donor_cut)), virtual_margin_z_mm=z, virtual_quantile=0.25))
        ids = geom.boundary(target_cut, z)
        donor_ids = geom.boundary(donor_cut, dz)
        curve = target_cut.vertices[ids]
        collar = target.triangles[(target.triangles[:, :, 2].max(1) <= z + 0.5) & (target.triangles[:, :, 2].min(1) >= z - 1)]
        donor_collar = donor.triangles[(donor.triangles[:, :, 2].max(1) <= dz + 0.5) & (donor.triangles[:, :, 2].min(1) >= dz - 1)]
        if len(collar) < 30 or len(donor_collar) < 30:
            raise ValueError('INSUFFICIENT_CERVICAL_OBSERVATION')
        dc = donor_collar.copy()
        dc[:, :, 0] *= -1
        dt = donor_cut.triangles.copy()
        dt[:, :, 0] *= -1
        dt = dt[:, ::-1]
        fit_start = time.perf_counter()
        (rot, translation, errors) = geom.fit(geom.sample(dc, 1024), geom.sample(collar, 1024))
        ext = geom.mesh(dt @ rot.T + translation, True)
        p = dict(margin_curve=curve, cavity_center=curve[:, :2].mean(0), cavity_scale=0.82, cavity_apex=np.r_[curve[:, :2].mean(0), z + 0.65 * (ext.vertices[:, 2].max() - z)])
        (fitted, trace) = gen.exact_trace(ext, curve)
        (final, roles) = gen.exact_shell(fitted, p)
        (ctrl, ctrl_roles, ctrl_info) = sys.modules['generate'].shell(ext.copy(), p)
        generator_name = 'reviewed FULL_CROWN_R4 exact_margin; virtual cervical adapter'
        generator_scope = 'VIRTUAL_COMPLETE_SHELL; physical crown UNKNOWN'
    except ValueError as exc:
        latest_failure = str(exc)
        generator_name = 'reviewed X18/Patient360 roof generator'
        generator_scope = 'CONDITIONAL_ROOF_PROXY; complete crown UNKNOWN'
        for (name, default) in [('target_cap_info', {}), ('donor_cap_info', {}), ('cut_drop', None), ('cut_components', None), ('donor_cut_drop', None), ('donor_cut_components', None)]:
            if name not in locals():
                pass
        target_cap_info = locals().get('target_cap_info', {})
        donor_cap_info = locals().get('donor_cap_info', {})
        cut_drop = locals().get('cut_drop', None)
        cut_components = locals().get('cut_components', None)
        donor_cut_drop = locals().get('donor_cut_drop', None)
        donor_cut_components = locals().get('donor_cut_components', None)
        site = source.site(low, fdi, 25)
        U = arches['upper']['tri']
        (ceiling, _) = source.query_height(U, site['xy'], True)
        donor_prior_path = package / 'batch9/demos/PATIENT360/dependencies/build/X18_crown_antagonist/donor_template.npz'
        with np.load(donor_prior_path, allow_pickle=False) as prior:
            template = prior[str(fdi % 10)]
        if len(template) != len(site['xy']):
            raise ValueError('FROZEN_DONOR_GRID_MISMATCH')
        fit_start = time.perf_counter()
        z0 = site['z_top'] + template
        roof = np.minimum(z0, ceiling - 0.025)
        (roof, repair) = source.repair(U, site, roof, 0.025)
        control = z0 - np.maximum(z0 - ceiling + 0.025, 0)
        (control, control_repair) = source.repair(U, site, control, 0.025)
        if not np.isfinite(roof).all() or not np.isfinite(control).all():
            raise ValueError('ROOF_NONFINITE')
        (v, f) = source.shell(site['xy'], roof, site['faces'], 1.2)
        (cv, cf) = source.shell(site['xy'], control, site['faces'], 1.2)
        final = geom.trimesh.Trimesh(v, f, process=False)
        ctrl = geom.trimesh.Trimesh(cv, cf, process=False)
        before_winding = bool(final.is_winding_consistent)
        geom.trimesh.repair.fix_normals(final, multibody=True)
        geom.trimesh.repair.fix_normals(ctrl, multibody=True)
        if final.volume < 0:
            final.invert()
        if ctrl.volume < 0:
            ctrl.invert()
        nfaces = len(site['faces'])
        roles = np.r_[np.zeros(nfaces, int), np.ones(nfaces, int), np.full(len(f) - 2 * nfaces, 2, int)]
        ctrl_roles = roles.copy()
        trace = dict(target_knot_identity_error_mm=None, virtual_margin='X18 roof sidewall, not anatomical finish line', repair=stable(repair), original_roof_winding_consistent=before_winding, exported_winding_consistent=bool(final.is_winding_consistent))
        proxy_control = dict(max_roof_difference_mm=float(np.max(np.abs(roof - control))), continuous_control_repair=stable(control_repair))
        dump(out / 'GENERATOR_REFUSAL.json', dict(latest_generator='FULL_CROWN_R4', reason=latest_failure, fallback=generator_name, fallback_scope=generator_scope, missing_measurement='Independently owned whole-tooth surface and actual preparation finish line'))
    fit_seconds = time.perf_counter() - fit_start
    native = final.vertices @ rotation.T + center
    source.write_stl(out / 'crown.stl', native, final.faces)
    np.savez_compressed(out / 'crown_local.npz', vertices=final.vertices, faces=final.faces, roles=roles, rotation=rotation, center=center)
    source.write_stl(out / 'crown_local.stl', final.vertices, final.faces)
    inner = final.triangles[roles == 1]
    source.write_stl(out / 'virtual_prep.stl', inner.reshape(-1, 3), np.arange(3 * len(inner)).reshape(-1, 3)[:, ::-1])
    antagonist = arches['upper']['tri']
    lo = final.bounds[0] - 2
    hi = final.bounds[1] + 2
    sel = np.all(antagonist[:, :, :2].max(1) >= lo[:2], axis=1) & np.all(antagonist[:, :, :2].min(1) <= hi[:2], axis=1)
    antagonist = antagonist[sel]
    source.write_stl(out / 'antagonist_local.stl', antagonist.reshape(-1, 3), np.arange(3 * len(antagonist)).reshape(-1, 3))
    paths = dict(crown=out / 'crown_local.stl', prep=out / 'virtual_prep.stl', antagonist=out / 'antagonist_local.stl')
    contract = dict(units='mm', common_frame_confirmed=True, contact_axis=[0, 0, 1], input_sha256={k: sha(v) for (k, v) in paths.items()}, regions=dict(intaglio=np.flatnonzero(roles == 1).tolist(), exterior=np.flatnonzero(roles == 0).tolist(), occlusal=np.flatnonzero(roles == 0).tolist(), preparation=list(range(len(inner))), antagonist=list(range(len(antagonist)))), ifu_profile='katana-ht', indication='posterior', ball_radius_final_mm=0.5, crown_solid_integrity_asserted=True, physical_geometry_status='VIRTUAL_PREPARATION; digital supplied pose', patient_id=m['patient'], geometry_sha256=geometry_hash, occlusal_gap_limits_mm=[-0.03, 0.1], occlusal_limits_source='Frozen illustrative laboratory band; PHENOMENOLOGICAL')
    dump(out / 'design_contract.json', contract)
    gate = check(paths, '3Y', out / 'design_contract.json', round_version='R2')
    dump(out / 'DESIGN_GATE_RAW.json', gate)
    b = stable(gate['rules']['ball_milling'])
    wall = stable(gate['rules']['material_wall'])
    gap = stable(gate['rules']['occlusal_contact'])
    gap_value = gap.get('measurement', {}).get('minimum_gap_mm')
    links = {'crown': contracts(package, ['D-E-K36']), 'milling': contracts(package, ['D-E-K40', 'D-E-K43']), 'force': contracts(package, ['D-E-DEFORMABLE-FORCES', 'D-E-K21']), 'pulp': contracts(package, ['D-E-K03']), 'nerve': contracts(package, ['D-E-K37'])}
    rows = [decision('crown', 'Crown design and material wall', {'FDI': fdi, 'design_verdict': gate['verdict'], 'wall_rule': wall, 'watertight': bool(final.is_watertight), 'virtual_margin_knot_error_mm': trace['target_knot_identity_error_mm']}, 'mm', 'Virtual preparation ; scanner error , the injectivity of form and physical fit Unknown . Distance enclosure is available in the marked X34 region; floating point enclosure is missing.', 'MODELLED', 'PER_SURFACE_REGION', gate['verdict'], source_locators, links['crown'], generator=generator_name, generator_scope=generator_scope, latest_generator_refusal=latest_failure, unsuitable_newer_operations=['fixed_prep_projection: fixed preparation is virtual', 'measured_neighbours: no signed loaded neighbour calibration']), decision('milling', 'Millability with final radius 0,5 mm', b, 'mm', 'X34 : s bullet can detect obstacles. Absence of obstacles does not prove the accessibility of the entire tool path; shrinkage/tool calibration UNKNOWN .', 'MODELLED', 'PER_SURFACE_REGION', b['status'], ['DEMO48_PACKAGE/X34/designgate'], links['milling']), decision('occlusal_height', 'Occlusal distance in delivered scan pose', gap_value, 'mm', 'Digital geometric measure in delivered pose; loaded pose and scanner repeatability UNKNOWN. No physical contact force can be identified from height alone.', 'MODELLED' if gap_value is not None else 'UNKNOWN', 'PER_POINT', gap['status'], source_locators, links['force'])]
    force_answer = dict(status='UNKNOWN', reason='No same-subject loaded height/force acquisition', force_interval_N=None)
    if m.get('force_state'):
        state = read(input_path(m['force_state']))
        validate_force_binding(state, m, geometry_hash)
        force_answer = force.predict(state, m.get('height_change_mm', 0), geometry_sha256=geometry_hash)
    else:
        force_answer.update(force_bridge.field_request({'case': m['case']}, {}, 'NATIVE_X21_GAPS'))
        force_answer['consumer'] = 'X82.consumer_bridge.field_request: executed native-state refusal'
    observed_fdi = sorted({int(k) for jaw in arches.values() for k in np.unique(jaw['labels']) if k > 0})
    tooth_forces = [dict(fdi=k, force_interval_N=None, evidence='UNKNOWN', resolution='PER_TOOTH') for k in observed_fdi]
    if force_answer.get('force_interval_N') is not None:
        tooth_forces = [dict(fdi=k, force_interval_N=v, evidence='MODELLED', resolution='PER_TOOTH') for (k, v) in zip(force_answer['predicted_fdi'], force_answer['force_interval_N'])]
    rows.append(decision('force', 'Force interval per tooth when changing altitude', force_answer.get('force_interval_N'), 'N', 'UNKNOWN without three loaded height probes with calibrated regional forces, the same geometry, wrench and contact branch. X82 has analytical model inclusion; general floating point/rare memory closure is missing.', 'UNKNOWN' if force_answer.get('force_interval_N') is None else 'MODELLED', 'PER_TOOTH', force_answer['status'], ['X82/crown_corridor.py'], links['force'], operator=force_answer, per_tooth=tooth_forces))
    pulp_value = None
    pulp_result = None
    if m.get('cbct') and m['cbct'].get('pulp_observation'):
        c = m['cbct']
        obs = read(input_path(c['pulp_observation']))
        validate_pulp_observation(obs, m, input_path(c['geometry_file']))
        pulp_result = pulp.uniform_shell_decision(obs['distance_mm'], obs['digital_two_surface_radius_mm'], 0.5, 0.5)
        pulp_value = pulp_result
    rows.append(decision('pulp', 'Pulp distance and preparation budget', pulp_value, 'mm', 'The same individual pulp -/ DEJ - boundary is missing.' if pulp_value is None else 'Digital surface geometry with voxelradie; total hard tissue , not separate dentin . Anatomical boundary accuracy UNKNOWN .', 'UNKNOWN' if pulp_value is None else 'MODELLED', 'PER_POINT', 'UNKNOWN' if pulp_value is None else pulp_result['decision'], ['Dataset modality manifest: no verified same-subject CBCT'] if pulp_value is None else [obs['locator']], links['pulp'], coordinate_frame='CBCT_NATIVE_MM' if pulp_value is not None else 'UNKNOWN', crown_coordinate_transfer='UNKNOWN; native CT observation not registered to the IOS preparation'))
    profiles = read(Path(args.nerve) / 'inputs/guide_profiles.json')['profiles']
    guide_rows = []
    for profile in profiles:
        for target_coverage in [0.9, 0.95]:
            failure = (1 - target_coverage) / 2
            budget = nerve.guide(profile, failure)
            direct_control = nerve.controls(profile, failure, budget)
            guide_rows.append(dict(guide=profile['id'], target=target_coverage, **budget, control=direct_control, resolution='POPULATION', physical_margin_mm=None, physical_coverage='UNKNOWN', source=profile['source']))
    rows.append(decision('nerve', 'Nerve margin per guide type', {'patient_margin_mm': None, 'population_guide_error_budgets': guide_rows}, 'mm', "The patient's anatomic canal wall/implant pose and common guide/wall error model are missing. The population's elements are treated as disclosure; the budgets are not patient margins. Rigorous floating point enclosure missing.", 'MODELLED', 'POPULATION', 'UNKNOWN_PATIENT_MARGIN', [p['source'] for p in profiles], links['nerve'], patient_resolution='PER_POINT', source_review='X87 PENDING_INDEPENDENT_REVIEW'))
    first_freeze_utc = epoch.get('per_case_frozen_utc') or datetime.datetime.now(datetime.timezone.utc).isoformat()
    predictions = dict(schema='dental-case-frozen-v1', frozen_utc=first_freeze_utc, patient=m['patient'], freeze_time_semantics='Actual first freeze of this exact input/code prediction; reruns reuse only its timestamp and must reproduce the complete payload hash', adapter_code_sha256=__import__('hashlib').sha256(''.join((sha(local_common.ROOT / 'code' / name) for name in ['common.py', 'operators.py', 'case_worker.py'])).encode()).hexdigest(), source_geometry_sha256=geometry_hash, source_members=src, input_manifest_sha256=sha(args.manifest), prereg_sha256=sha(local_common.ROOT / 'PREREG_R1.json'), changed_construction_prereg_sha256=sha(local_common.ROOT / 'PREREG_R2.json'), source_domain_prereg_sha256=sha(local_common.ROOT / 'PREREG_R3.json'), fallback_prereg_sha256=sha(local_common.ROOT / 'PREREG_R4.json'), crown_stl_sha256=sha(out / 'crown.stl'), trace=trace, force=force_answer, pulp=pulp_result, decision_payload_sha256=__import__('hashlib').sha256(local_common.encoded(rows)).hexdigest(), physical_measurement='NOT_RUN', prior_exposure='Retrospective digital reconstruction; X11 labels and some cases previously evaluated', prospective_targets=['as-built surface', 'seated regional film', 'loaded region force'], prospective_numeric_predictions='UNKNOWN until matched preparation, material batch and loaded calibration exist')
    freeze(out / 'FROZEN_PREDICTIONS.json', predictions)
    measured_target_cut = geom.cut(target, z)
    metrics = geom.metrics(final.triangles[roles == 0], measured_target_cut.triangles, 8192)
    ctrl_metrics = geom.metrics(ctrl.triangles[ctrl_roles == 0], measured_target_cut.triangles, 8192)
    readback = source_geometry_readback(out / 'crown.stl')
    export_error = float(np.max(np.abs(readback - native[final.faces])))
    if export_error > 1e-05:
        raise ValueError('EXPORT_TRANSPORT_GATE')
    facit = dict(kind='published_dataset', locator=source_locators, compared_quantity='Original measured tooth surface above virtual cervical plane, bidirectional 8192-point point-to-triangle distance; retrospective source reconstruction, not a blinded independent prediction', refutes_us=metrics['p95_mm'] > 0.5, threshold_status='0.5mm descriptive prereg reference only; no clinical acceptance', resolution='PER_TOOTH', uncertainty='Sampled surface metric; no scanner, anatomical or rigorous sampling/roundoff enclosure', measured_original=metrics, equally_informed_rigid_control=ctrl_metrics, proxy_control=proxy_control, data_exposure='Original archive is readable by generator/scorer worker. Numeric scoring occurs after freeze, but no blind target-access separation is claimed.')
    validate_decisions(rows)
    case = dict(schema='dental-case-report-v1', claim_type='capability', patient=m['patient'], dataset='Bits2Bites', licence='CC BY-NC-SA according to dataset brief', source_geometry_sha256=geometry_hash, source_members=src, geometry=dict(frame='lower X11 rigid frame; supplied registered pair', rotation=rotation, origin_mm=center, units='mm', fdi=fdi, target_component_discarded_area_fraction=drop, target_components=components, donor_component_discarded_area_fraction=donor_drop, donor_components=donor_components, target_virtual_closure=target_cap_info, donor_virtual_closure=donor_cap_info, cut_component_discarded_area_fraction=cut_drop, cut_components=cut_components, donor_cut_discarded_area_fraction=donor_cut_drop, donor_cut_components=donor_cut_components, export_float32_max_mm=export_error, virtual_cut_quantile=0.25, virtual_cavity_scale=0.82, virtual_margin_z_mm=z), decisions=rows, external_referent=facit, frozen_predictions_sha256=sha(out / 'FROZEN_PREDICTIONS.json'), clinical_recommendation=False, lab_release='ABSTAIN', physical_measurement='NOT_RUN', review_state='PENDING_INDEPENDENT_REVIEW', source_review='Reviewed package operators; new adapter and X87 pending review')
    dump(out / 'case.json', case)
    figure(out / 'case_figure.png', measured_target_cut.triangles, final.triangles[roles == 0], fdi)
    dump(out / 'COSTS.json', dict(wall_seconds=time.perf_counter() - started, fit_seconds=fit_seconds, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, cpu_seconds=time.process_time(), threads=4, gpu=False, inherited_label_training='UNKNOWN', physical_validation='NOT_RUN', acquisition_and_human_discovery='UNKNOWN'))
    print(m['patient'], gate['verdict'], b['status'], metrics['p95_mm'], flush=True)

def source_geometry_readback(path):
    b = Path(path).read_bytes()
    n = int.from_bytes(b[80:84], 'little')
    dt = np.dtype([('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')])
    return np.frombuffer(b, dt, count=n, offset=84)['v'].astype(float)

def figure(path, original, generated, fdi):
    import matplotlib.pyplot as plt
    (fig, axes) = plt.subplots(1, 2, figsize=(9, 3.8))
    zlo = min(original[:, :, 2].min(), generated[:, :, 2].min())
    zhi = max(original[:, :, 2].max(), generated[:, :, 2].max())
    xlo = min(original[:, :, 0].min(), generated[:, :, 0].min())
    xhi = max(original[:, :, 0].max(), generated[:, :, 0].max())
    ylo = min(original[:, :, 1].min(), generated[:, :, 1].min())
    yhi = max(original[:, :, 1].max(), generated[:, :, 1].max())
    for (ax, tri, title) in zip(axes, [original, generated], ['Measured Original Surface', 'Genererad forskningskrona']):
        pts = tri.mean(1)
        stride = max(1, len(pts) // 5000)
        q = pts[::stride]
        sc = ax.scatter(q[:, 0], q[:, 1], c=q[:, 2], s=2, cmap='viridis', vmin=zlo, vmax=zhi)
        ax.set_xlim(xlo - 0.2, xhi + 0.2)
        ax.set_ylim(ylo - 0.2, yhi + 0.2)
        ax.set_aspect('equal')
        ax.set_title(title)
        ax.set_xlabel('local x [mm ]')
        ax.set_ylabel('local y [ mm ]')
    fig.subplots_adjust(left=0.08, right=0.88, bottom=0.18, top=0.82, wspace=0.45)
    fig.colorbar(sc, ax=axes, label='local z [mm ]', shrink=0.8)
    fig.suptitle('FDI ' + str(fdi) + ' · virtual preparation · physical fit unknown')
    fig.savefig(path, dpi=130, metadata={'Software': 'X90 deterministic case report'})
    plt.close(fig)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    for name in ['manifest', 'out', 'package', 'nerve', 'dataset-root', 'labels-root', 'epoch']:
        ap.add_argument('--' + name, required=True)
    run(ap.parse_args())
