import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '4'
import copy, json, math, resource, time, shutil
import numpy as np
import trimesh
from common import ROOT, dump, sha, stamp, check_frozen, inputs, state
from sources import read_sources
from compare_lab import score

def local_frame(point, normal):
    z = np.array(normal, float)
    z /= np.linalg.norm(z)
    a = np.array([1.0, 0.0, 0.0])
    if abs(a @ z) > 0.8:
        a = np.array([0.0, 1.0, 0.0])
    x = np.cross(a, z)
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    T = np.eye(4)
    T[:3, :3] = np.column_stack([x, y, z])
    T[:3, 3] = point
    return T

def geometry_check(path, frame, expected):
    mesh = trimesh.load_mesh(path, process=True)
    verts = trimesh.transform_points(mesh.vertices, np.linalg.inv(frame))
    bounds = verts.max(axis=0) - verts.min(axis=0)
    error = float(abs(bounds - expected).max())
    return {'dimension_error_mm': error, 'local_dimensions_mm': bounds.tolist(), 'watertight': bool(mesh.is_watertight), 'positive_volume': bool(mesh.volume > 0), 'pass': bool(error <= 1e-05 and mesh.is_watertight and (mesh.volume > 0))}

def main():
    started = time.perf_counter()
    p = check_frozen('PREREG_R4.json')
    inputs()
    outdir = ROOT / 'exports'
    outdir.mkdir(exist_ok=True)
    (bond, optical, target) = read_sources()
    protocol_b = {'zirconia_product': '3M Lava Esthetic', 'air_abrasion_particle_um': 50.0, 'air_abrasion_pressure_bar': 3.5, 'air_abrasion_seconds': 15.0, 'air_abrasion_standoff_mm': 10.0, 'curing_irradiance_mW_cm2': 1514.0, 'curing_wavelength_nm': [430, 480], 'curing_seconds_per_direction': 20.0, 'curing_directions': 5, 'thermocycles': 5000, 'thermal_temperatures_C': [5, 55], 'thermal_dwell_seconds': 30.0, 'thermal_transfer_seconds': 5.0, 'test_crosshead_mm_min': 0.5, 'source_locator': 'doi:10.1016/j.jds.2023.05.011'}
    protocol_o = {'zirconia_product': 'Aidite 3D Pro Zir A1', 'cement_product': 'GC G-CEM LinkForce', 'spectrophotometer': 'OptiShade StyleItaliano', 'reference_Lab': target, 'seating_force_N': 4.9, 'seating_minutes': 5.0, 'curing_exposures': 2, 'curing_seconds_per_exposure': 10.0, 'source_locator': 'doi:10.1016/j.heliyon.2023.e23046'}
    payload = {'bond': bond, 'optical': optical, 'bond_protocol': protocol_b, 'optical_protocol': protocol_o, 'prospective_status': 'Future physical source replication NOT_RUN; values from pre-existing independent measurements. No patient physical forecast.'}
    pred = ROOT / 'FROZEN_LAB_PREDICTIONS.json'
    if pred.exists():
        assert check_frozen(pred.name)['predictions'] == payload
    else:
        dump(pred, {'frozen_utc': stamp(), 'prereg_sha256': sha(ROOT / 'PREREG_R4.json'), 'predictions': payload})
        (ROOT / (pred.name + '.sha256')).write_text(sha(pred) + '\n')
    r2 = json.loads((ROOT / 'rounds/R2.json').read_text())
    point = r2['inverse_measurement_targets'][0]
    rays = np.load(ROOT / 'raw/rays.npz')
    visible = np.flatnonzero(rays['visible'])
    k = int(visible[0])
    frames = {'bond': local_frame(point['coordinate_mm'], point['normal']), 'optical': local_frame(rays['points'][k], rays['normals'][k])}
    specs = [('bond_substrate_disc', trimesh.creation.cylinder(radius=7.0, height=2.5, sections=128), frames['bond'], np.array([14.0, 14.0, 2.5])), ('optical_ceramic_tile', trimesh.creation.box(extents=[8.0, 8.0, 1.0]), frames['optical'], np.array([8.0, 8.0, 1.0])), ('optical_substrate_tile', trimesh.creation.box(extents=[10.0, 10.0, 2.0]), frames['optical'], np.array([10.0, 10.0, 2.0]))]
    checks = []
    for (name, mesh, T, expected) in specs:
        mesh.apply_transform(T)
        path = outdir / (name + '.stl')
        mesh.export(path)
        check = geometry_check(path, T, expected)
        assert check['pass'], check
        bad = trimesh.load_mesh(path, process=True)
        v = trimesh.transform_points(bad.vertices, np.linalg.inv(T))
        v[:, 2] *= 1.1
        bad.vertices = trimesh.transform_points(v, T)
        badpath = ROOT / 'raw' / ('INJECTED_' + name + '.stl')
        bad.export(badpath)
        injected = geometry_check(badpath, T, expected)
        assert not injected['pass']
        checks.append({'name': name, 'file': str(path.relative_to(ROOT)), 'expected_dimensions_mm': expected.tolist(), 'frame_patient_mm': T.tolist(), 'source_locator': p['source_geometry']['bond_locator' if name.startswith('bond') else 'optical_locator'], 'control': check, 'actual_thickness_x1_1_fault': injected, 'sha256': sha(path), 'resolution_level': 'PER_SURFACE_REGION'})
    shutil.copyfile(ROOT / 'raw/crown.stl', outdir / 'STS_D1_patient_reference.stl')
    assert sha(outdir / 'STS_D1_patient_reference.stl') == sha(ROOT / 'raw/crown.stl')
    bfixture = {'kind': 'bond_coupon', 'evidence_type': 'our_own_fixture', 'material_batch_id': 'SYNTHETIC_SCORER_ONLY', 'system': 'RXU', 'with_cover': True, 'adhesive_product': 'Scotchbond Universal', 'measured_cover_thickness_bounds_mm': [1.0, 1.0], 'measured_cover_internal_space_mm': 0.3, 'cover_shape_scan_locator': 'our_own_fixture', **protocol_b, 'specimens': [{'specimen_id': f'SYNTHETIC_B{i:02d}', 'peak_force_N': 8.35 * math.pi * (4.91 / 2) ** 2, 'measured_bonded_area_mm2': math.pi * (4.91 / 2) ** 2, 'failure_mode': 'MIXED', 'failure_origin_record': 'SYNTHETIC', 'raw_observation_locator': 'our_own_fixture'} for i in range(10)]}
    ofixture = {'kind': 'optical_stack', 'evidence_type': 'our_own_fixture', 'material_batch_id': 'SYNTHETIC_SCORER_ONLY', 'substrate': 'CR', 'cement': 'BL', **protocol_o, 'specimens': [{'specimen_id': f'SYNTHETIC_O{i:02d}', 'Lab': [74.52, -0.93, 11.84], 'measured_tile_thickness_mm': 1.0, 'measured_cement_thickness_um': 100.0, 'raw_observation_locator': 'our_own_fixture'} for i in range(10)]}
    faults = []
    for (label, original, mutate) in [('force_x10', bfixture, lambda x: [r.update(peak_force_N=r['peak_force_N'] * 10) for r in x['specimens']]), ('ceramic_failure', bfixture, lambda x: x['specimens'][0].update(failure_mode='CERAMIC_FRACTURE')), ('missing_area', bfixture, lambda x: x['specimens'][0].pop('measured_bonded_area_mm2')), ('wrong_cycles', bfixture, lambda x: x.update(thermocycles=500)), ('wrong_cure', bfixture, lambda x: x.update(curing_directions=1)), ('missing_adhesive_identity', bfixture, lambda x: x.update(system='RXU_ALONE')), ('missing_adhesive_product', bfixture, lambda x: x.pop('adhesive_product')), ('wrong_cover_thickness', bfixture, lambda x: x.update(measured_cover_thickness_bounds_mm=[0.8, 0.8])), ('missing_cover_dome_scan', bfixture, lambda x: x.pop('cover_shape_scan_locator')), ('missing_batch', bfixture, lambda x: x.pop('material_batch_id')), ('Lab_L_plus20', ofixture, lambda x: [r['Lab'].__setitem__(0, r['Lab'][0] + 20) for r in x['specimens']]), ('wrong_tile_thickness', ofixture, lambda x: x['specimens'][0].update(measured_tile_thickness_mm=1.2)), ('missing_film', ofixture, lambda x: x['specimens'][0].pop('measured_cement_thickness_um')), ('wrong_backing', ofixture, lambda x: x.update(substrate='NATURAL_DENTIN')), ('different_zirconia_product', ofixture, lambda x: x.update(zirconia_product='3M Lava Esthetic')), ('unknown_specimen_type', ofixture, lambda x: x.update(kind='PATIENT_CROWN'))]:
        x = copy.deepcopy(original)
        mutate(x)
        result = score(x, payload)
        assert result['gate'] == 'FAIL', label
        faults.append({'fault': label, 'actual_scorer_result': result})
    assert score(bfixture, payload)['gate'] == score(ofixture, payload)['gate'] == 'PASS_SOFTWARE_FIXTURE'
    dump(ROOT / 'raw/SCORER_BOND_FIXTURE.json', bfixture)
    dump(ROOT / 'raw/SCORER_OPTICAL_FIXTURE.json', ofixture)
    dump(ROOT / 'raw/SCORER_FAULTS.json', faults)
    acquisition = {'geometry_checks': checks, 'local_patient_frames': {k: v.tolist() for (k, v) in frames.items()}, 'source_matched_geometry_scope': 'complete dimensions of flat optical tiles and bond substrate disc, STL mm; rigid placement only, not patient anatomy transfer', 'not_exported': ['size-5 capsule dome', '1-mm zirconia cover: requires actual capsule/cast scan; internal space .3 mm'], 'source_patient_products_joinable': False, 'minimum_missing_measurements': ['Procure capsule #5; measure actual bonded internal area; scan its dome/cast and fabricate source-specific 1-mm cover.', 'For curved-crown transfer: matched local stack Lab/spectra and cement-film thickness on the chosen same-product crown.', 'For mechanics: force-displacement and failure origin on the same-product crown; peak SBS alone leaves cohesive tangent/energy/mixed-mode law unknown.'], 'physical_trials': 'NOT_RUN', 'clinical_output': 'NONE'}
    dump(ROOT / 'LAB_REQUESTS.json', acquisition)
    out = {'round': 'R4', 'claim_type': 'capability', 'gate': 'PASS_EXECUTABLE_SOURCE_REPLICATION_PACKET', 'physical_validation': 'NOT_RUN', 'patient_transfer': 'UNKNOWN', 'geometry_checks': checks, 'scorer_fault_rejections': len(faults), 'external_referent': {'kind': 'independent_measurement', 'locator': ['doi:10.1016/j.jds.2023.05.011', 'doi:10.1016/j.heliyon.2023.e23046'], 'compared_quantity': 'complete flat specimen dimensions and source-replication expectation; no new physical observation', 'refutes_us': False}, 'synthetic_fixture_role': 'Scores correctness only; not external reference or measured agreement', 'cost': {'elapsed_s': time.perf_counter() - started, 'max_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'threads_max': 4}}
    dump(ROOT / 'rounds/R4.json', out)
    state(status='R4_DECIDED', latest_gate={'export_and_scoring': 'PASS', 'physical_source_replication': 'NOT_RUN', 'patient_transfer': 'UNKNOWN'}, next_operation='Compile one-command replay, independent controls, selection/dropout and reviewed handoff; physical follow-up needs manufactured matched measurements')
    print(json.dumps({'exported': len(checks), 'scorer_fault_rejections': len(faults), 'physical': 'NOT_RUN', 'patient': 'UNKNOWN'}))
if __name__ == '__main__':
    main()
