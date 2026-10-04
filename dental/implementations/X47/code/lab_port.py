"""Read-only preparation port and disjoint 72-specimen assay design."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, csv, hashlib, time, math
import numpy as np
import trimesh
ROOT = Path(__file__).resolve().parents[1]
X1B = ROOT.parent / _release_expand('X1B')

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def allocation():
    rng = np.random.default_rng(470003)
    rows = []
    for d in ['D1', 'M1', 'M2']:
        for block in [0, 30]:
            assignments = [('fracture', 'unaged')] * 8 + [('retention', 'sham37')] * 2 + [('retention', 'TC10000')] * 2
            rng.shuffle(assignments)
            for (i, (endpoint, age)) in enumerate(assignments, 1):
                rows.append({'specimen_id': f'{d}-A{block}-N{i:02d}', 'design_id': d, 'legacy_angle_block_deg': block, 'ultimate_endpoint': endpoint, 'age_state': age, 'tensile_axis_deg': 0 if endpoint == 'retention' else '', 'thermal_cycles': 10000 if age == 'TC10000' else 0})
    return rows

def verify_allocation(rows):
    if len(rows) != 72 or len({r['specimen_id'] for r in rows}) != 72:
        raise ValueError('72_UNIQUE_SPECIMENS_REQUIRED')
    expected = allocation()
    if rows != expected:
        raise ValueError('FROZEN_ALLOCATION_CHANGED')
    for r in rows:
        if r['ultimate_endpoint'] == 'retention' and float(r['tensile_axis_deg']) != 0:
            raise ValueError('RETENTION_MUST_BE_AXIAL')
    return True

def geometry():
    out = []
    manifest = []
    for d in ['D1', 'M1', 'M2']:
        stl = X1B / f'exports/{d}/preparation.stl'
        gp = X1B / f'inputs/geometry/{d}_grid.json'
        g = json.loads(gp.read_text())
        margin = g['z_m'] + g.get('margin_shift_mm', 0)
        mesh = trimesh.load_mesh(stl, process=True)
        cent = mesh.triangles_center
        normal = mesh.face_normals
        area = mesh.area_faces
        if not mesh.is_watertight or not np.isfinite(cent).all() or np.any(area <= 0):
            raise ValueError('INVALID_PREPARATION_MESH')
        z = cent[:, 2] - margin
        region = np.where(z < 0, 0, np.where(z <= 0.3, 1, np.where(abs(normal[:, 2]) < math.sqrt(0.5), 2, 3)))
        signed = np.degrees(np.arcsin(np.clip(normal[:, 2], -1, 1)))
        rnames = ['BELOW_MARGIN', 'CHAMFER_BAND', 'AXIAL', 'OCCLUSAL']
        face_file = ROOT / f'raw/{d}_interface_faces.npz'
        np.savez_compressed(face_file, centroid_mm=cent, normal=normal, area_mm2=area, height_above_margin_mm=z, region_index=region, signed_normal_tilt_deg=signed)
        regions = {name: float(area[region == k].sum()) for (k, name) in enumerate(rnames)}
        closure = abs(sum(regions.values()) - float(area.sum()))
        ax = region == 2
        quantiles = np.quantile(signed[ax], [0, 0.1, 0.5, 0.9, 1]).tolist() if ax.any() else []
        out.append({'design_id': d, 'face_count': len(area), 'units': {'area': 'mm2', 'height': 'mm', 'normal_tilt': 'deg'}, 'resolution': 'PER_POINT triangles; aggregated PER_SURFACE_REGION descriptive areas', 'region_area_mm2': regions, 'area_partition_identity_error_mm2': closure, 'descriptive_height_above_margin_max_mm': float(mesh.vertices[:, 2].max() - margin), 'axial_signed_normal_tilt_unweighted_quantiles_deg': quantiles, 'summary_decisive': False, 'global_TOC_deg': 'UNKNOWN_NONCONICAL_PREPARATION', 'effective_bonded_area_mm2': 'UNKNOWN_NOT_MEASURED', 'local_face_artifact': str(face_file.relative_to(ROOT)), 'local_face_sha256': sha(face_file), 'absolute_pull_off_prediction_N': 'UNKNOWN_FAILED_TRANSFER_AND_COMPOSITE_SUBSTRATE', 'source_stl_sha256': sha(stl), 'source_grid_sha256': sha(gp), 'source_stl': str(stl), 'source_grid': str(gp), 'margin_z_mm': margin})
        manifest.extend([{'path': str(p), 'sha256': sha(p), 'bytes': p.stat().st_size} for p in [stl, gp]])
    return (out, manifest)

def run():
    start = time.monotonic()
    a = allocation()
    verify_allocation(a)
    with (ROOT / 'LAB_ALLOCATION.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(a[0]))
        w.writeheader()
        w.writerows(a)
    keys = list(a[0]) + ['force_N', 'force_unit', 'failure_mode', 'material_batch', 'cement_product', 'cement_batch', 'surface_treatment', 'die_material_batch', 'actual_TC_cycles', 'TC_low_C', 'TC_high_C', 'dwell_s', 'transfer_s', 'post_storage_hours', 'pull_speed_mm_min', 'measured_pull_axis_deg', 'force_trace_sha256', 'image_sha256', 'notes']
    with (ROOT / 'LAB_MEASUREMENTS_TEMPLATE.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(a)
    (geo, manifest) = geometry()
    legacy = json.loads((X1B / 'FROZEN_PREDICTIONS.json').read_text())
    result = {'round': 'R3', 'claim_type': 'capability', 'allocation_counts': {'fracture_unaged': 48, 'retention_sham': 12, 'retention_TC10000': 12, 'total': 72}, 'fracture_n_per_design_angle': 8, 'retention_n_per_design_aging': 4, 'original_M2_planning_min_n': 7, 'M2_planning_n_gate_passes': True, 'M1_still_underpowered': True, 'retention_power': 'UNKNOWN; no cross-study variance transfer justified', 'geometry': geo, 'source_manifest': manifest, 'timing_wall_s': time.monotonic() - start, 'physical_measurements_performed': False, 'rigorous_physical_error_enclosure': 'MISSING; no affine sensitivity reported', 'external_referent': {'kind': 'independent_measurement', 'locator': 'https://doi.org/10.1155/2023/6664894 Table1, methods; Choi https://doi.org/10.4047/jap.2023.15.2.63 Table1', 'compared_quantity': 'Published disjoint whole-crown retention/fracture endpoints and local-feature force contrast; new 72-crown assay unmeasured', 'refutes_us': True}}
    (ROOT / 'rounds/RESULTS_R3.json').write_text(json.dumps(result, indent=2) + '\n')
    return result
if __name__ == '__main__':
    r = run()
    print(r['allocation_counts'])
    print([(x['design_id'], x['region_area_mm2'], x['descriptive_height_above_margin_max_mm']) for x in r['geometry']])
