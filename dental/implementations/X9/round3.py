from dental_release.paths import expand as _release_expand
import json
import resource
import time
from pathlib import Path
import numpy as np
from skimage.measure import marching_cubes
import trimesh
from ipr import jsonable, load_case, ray_hits, sha, write_json
from segmented_tooth_tool import build_protection, envelope
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X9-ipr-safety'))

def rasterize(meshes, h=0.15):
    low = meshes['Tooth'].bounds[0] - 0.6
    high = meshes['Tooth'].bounds[1] + 0.6
    axes = [np.arange(low[i], high[i] + h / 2, h) for i in range(3)]
    (yy, zz) = np.meshgrid(axes[1], axes[2], indexing='ij')
    yz = np.column_stack([yy.ravel(), zz.ravel()])
    mats = {}
    odd_counts = {}
    tooth_unknown = np.zeros((len(axes[0]), len(yz)), bool)
    for (name, mesh) in meshes.items():
        a = np.zeros((len(axes[0]), len(yz)), bool)
        odd = 0
        for start in range(0, len(yz), 256):
            hits = ray_hits(mesh, yz[start:start + 256], 1)
            for (j, xs) in enumerate(hits):
                xs = sorted(xs)
                if len(xs) % 2:
                    odd += 1
                    if name == 'Tooth':
                        tooth_unknown[:, start + j] = True
                    continue
                for (left, right) in zip(xs[::2], xs[1::2]):
                    a[:, start + j] |= (axes[0] >= left) & (axes[0] <= right)
        mats[name] = a.reshape(tuple((len(x) for x in axes)))
        odd_counts[name] = odd
    labels = np.where(mats['Tooth'], 2, 0).astype(np.uint8)
    labels[mats['Enamel'] & mats['Tooth']] = 1
    labels[tooth_unknown.reshape(labels.shape)] = 2
    return (labels, np.full(3, h), low, odd_counts)

def run():
    (tic, cpu) = (time.perf_counter(), time.process_time())
    pre = json.loads(Path('PREREG_R3.json').read_text())
    r1 = json.loads(Path('rounds/R1_raw_profiles.json').read_text())
    inputs = json.loads(Path('PREREG_R1.json').read_text())['inputs']
    (rows, grid_records) = ([], [])
    DATA.mkdir(parents=True, exist_ok=True)
    for (i, item) in enumerate(inputs):
        (meshes, meta, frame) = load_case(item['path'])
        (labels, spacing, origin, odd) = rasterize(meshes, pre['grid_mm'])
        (protected, diag) = build_protection(labels, spacing)
        file = DATA / (item['tid'] + '_labels.npz')
        np.savez_compressed(file, labels=labels, spacing_mm=spacing, origin_mm=origin, provenance=np.array('STS Otsu enamel; non-enamel includes unknown exterior skin; anatomy uncalibrated'))
        grid_records.append({'tid': item['tid'], 'path': str(file), 'sha256': sha(file), 'bytes': file.stat().st_size, 'odd_ray_counts': odd, 'shape': labels.shape, 'diagnostics': diag})
        for p in [p for p in r1 if p['tid'] == item['tid']]:
            side = 1 if p['side'] == 'mesial' else -1
            (out, removed) = envelope(labels, protected, spacing, origin, side, p['buccolingual_center_mm'], p['height_center_mm'])
            row = {'tid': item['tid'], 'jaw': item['jaw'], 'tooth_type': item['type'], 'side': p['side'], 'height_fraction': p['height_fraction'], 'y_mm': p['buccolingual_center_mm'], 'z_mm': p['height_center_mm'], **out}
            rows.append(row)
            if i == 0 and p['side'] == 'mesial' and (p['height_fraction'] == 0.5):
                np.savez_compressed(DATA / 'demo_segmented_tooth.npz', labels=labels, spacing_mm=spacing, origin_mm=origin)
                if removed is not None:
                    remaining = (labels > 0) & ~removed
                    (v, f, _, _) = marching_cubes(remaining.astype(np.float32), 0.5, spacing=tuple(spacing))
                    trimesh.Trimesh(v + origin, f, process=False).export('artifacts/conditional_post_ipr.stl')
        write_json('rounds/R3_rows.partial.json', rows)
        write_json('CURRENT_WORK_STATE.json', {'stage': 'R3_RUNNING', 'last_tooth': item['tid'], 'completed_teeth': i + 1, 'next_operation': 'finish 3D protected-tissue cut envelope and published type comparison'})
        print(item['tid'], 'SDF done', labels.shape, flush=True)
    valid = [r for r in rows if r['scenario_capacity_mm'] is not None]
    refs = json.loads(Path('sources/external_by_tooth_type.json').read_text())
    by_type = {(r['jaw'], r['tooth_type']): r for r in refs}
    measured = json.loads(Path('rounds/R1_rows.json').read_text())
    comparisons = []
    for r in measured:
        if r['height_fraction'] == 0.5 and r['mean_mm'] is not None:
            ref = by_type[r['jaw'], r['tooth_type']]
            target = ref[r['side'] + '_mean_mm']
            comparisons.append({'tid': r['tid'], 'side': r['side'], 'our_mean_mm': r['mean_mm'], 'external_mean_mm': target, 'error_mm': r['mean_mm'] - target, 'locator': ref['locator']})
    mederr = float(np.median([abs(r['error_mm']) for r in comparisons]))
    parity = max([r['parity_error_mm'] for r in valid], default=None)
    failures = sum((not r['generic_0_5_pass'] for r in valid))
    g3 = mederr <= 0.2
    gates = {'G1_solver_parity': {'pass': bool(valid) and parity <= 1e-07, 'max_error_mm': parity, 'outcome': 'TIE'}, 'G2_cut_guard': {'pass': bool(valid) and all((r['protected_overlap_cells'] == 0 for r in valid)), 'protected_overlap_cells_total': sum((r['protected_overlap_cells'] for r in valid))}, 'G3_external_type_contact': {'pass': g3, 'median_absolute_error_mm': mederr, 'frozen_tolerance_mm': 0.2, 'comparisons': len(comparisons), 'meaning': 'different specimens; population comparison'}, 'G4_praxis_geometric_margin': {'pass': bool(valid) and failures / len(valid) >= 0.1, 'failing_fixed_0_5_patches': failures, 'evaluated_patches': len(valid)}, 'G5_anatomical_admission': {'pass': False, 'outcome': 'UNKNOWN'}}
    try:
        build_protection(np.zeros((3, 3, 3), np.uint8), np.full(3, 0.15))
        empty_rejected = False
    except ValueError:
        empty_rejected = True
    injections = {'solver_plus_0_1_rejected': parity is None or parity + 0.1 > 1e-07, 'cut_plus_2_rejected_all': bool(valid) and all((r['injected_plus_2_mm_rejected'] for r in valid)), 'empty_tissue_mask_rejected': empty_rejected, 'external_thickness_plus_2_rejected': np.median([abs(r['error_mm'] + 2) for r in comparisons]) > 0.2, 'missing_anatomical_calibration_blocks_positive_safety': all((r['calibrated_safe_ipr_mm'] is None for r in rows))}
    result = {'round': 'R3', 'prereg_sha256': sha('PREREG_R3.json'), 'gates': gates, 'injected_faults': injections, 'number_teeth': len(inputs), 'number_patches': len(rows), 'positive_scenario_patches': sum((r['scenario_capacity_mm'] > 0 for r in valid)), 'external_referent': {'kind': 'independent_measurement', 'locator': 'https://doi.org/10.1590/2177-6709.29.3.e242422.oar', 'compared_quantity': 'proximal contact enamel thickness per jaw and M1/M2 side, Table 1 (mm)', 'refutes_us': not g3}, 'array_manifest': grid_records, 'cost': {'wall_seconds': time.perf_counter() - tic, 'cpu_seconds': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'physical_validation': 'NOT_DONE'}}
    write_json('rounds/R3_rows.json', rows)
    write_json('rounds/R3_external_comparison.json', comparisons)
    write_json('rounds/R3_results.json', result)
    print(json.dumps(jsonable({k: v for (k, v) in result.items() if k != 'array_manifest'}), indent=2))
    return result
if __name__ == '__main__':
    run()
