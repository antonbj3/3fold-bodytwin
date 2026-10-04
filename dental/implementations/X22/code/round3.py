from dental_release.paths import expand as _release_expand
import json
import resource
import time
from pathlib import Path
import nibabel as nib
import numpy as np
from scipy import ndimage as ndi
from common import *

def native_coordinates(local, center, R, croplo, spacing):
    return (local @ R + center + croplo * spacing) / spacing

def sample_scene(source, shape, origin, h, center, R, croplo, spacing):
    scene = np.zeros(shape, bool)
    (x, y) = np.meshgrid(np.arange(shape[0]) * h + origin[0], np.arange(shape[1]) * h + origin[1], indexing='ij')
    for k in range(shape[2]):
        z = np.full_like(x, k * h + origin[2])
        local = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
        xyz = native_coordinates(local, center, R, croplo, spacing)
        scene[:, :, k] = ndi.map_coordinates(source, xyz.T, order=0, mode='constant', cval=0, prefilter=False).reshape(x.shape) == 1
    return scene

def run():
    tic = time.perf_counter()
    cpu = time.process_time()
    pre = json.loads((ROOT / 'PREREG_R3_NATIVE_TOOTH_CONTEXT.json').read_text())
    r1 = json.loads((ROOT / 'PREREG_R1_SYNTHETIC_BENCHMARK.json').read_text())
    outdir = DATA / 'native_context'
    outdir.mkdir(exist_ok=True)
    results = []
    manifest = []
    rng = np.random.default_rng(32202)
    for tid in pre['inputs']:
        inp = next((r for r in r1['inputs'] if r['tid'] == tid))
        q = np.load(inp['path'])
        target = q['labels'] > 0
        h = 0.15
        old_origin = q['origin_mm']
        parent = Path(_release_expand('@DENTAL_INPUT_ROOT@/artifacts/CROWN_pop/cases')) / (tid + '.npz')
        p = np.load(parent)
        meta = json.loads(str(p['meta']))
        center = p['molar_centroid']
        m = p['mesial_centroid'] - center
        m[2] = 0
        m /= np.linalg.norm(m)
        z = np.array([0, 0, 1 if meta['jaw'] == 'lower' else -1])
        b = np.cross(z, m)
        R = np.stack([m, b, z])
        spacing = np.array(meta['zooms'])
        croplo = np.array(meta['crop_lo_vox'])
        src = Path(_release_expand('@DENTAL_CORPUS_ROOT@/geometry/STS-3D-Tooth/ROI/Labeled/Mask')) / f"ROI_L_{meta['vol']:03d}.nii.gz"
        hashpass = sha(src) == meta['sha256'][f"Mask/ROI_L_{meta['vol']:03d}.nii.gz"]
        source = np.asanyarray(nib.load(src).dataobj)
        pad = np.array([80, 20, 100])
        shape = np.array(target.shape) + 2 * pad
        origin = old_origin - pad * h
        scene = sample_scene(source, tuple(shape), origin, h, center, R, croplo, spacing)
        sl = tuple((slice(int(i), int(i + n)) for (i, n) in zip(pad, target.shape)))
        occupied = np.mean(scene[sl][target])
        local = np.argwhere(target) * h + old_origin
        native = native_coordinates(local, center, R, croplo, spacing)
        inverse = (native * spacing - croplo * spacing - center) @ R.T
        frame_error = float(np.abs(inverse - local).max())
        wrong_local = local.copy()
        wrong_local[:, 0] += 10
        wrong_native = native_coordinates(wrong_local, center, R, croplo, spacing)
        wrong_occupancy = float(np.mean(ndi.map_coordinates(source, wrong_native.T, order=0, mode='constant', cval=0) == 1))
        wrong_inverse = (wrong_native * spacing - croplo * spacing - center) @ R.T
        wrong_frame_error = float(np.abs(wrong_inverse - local).max())
        del source
        (labels, dt) = tissue(target, h, 1.2)
        target_mu = coefficient(labels)
        mu = np.where(scene, pre['operation']['background_mu_cm_inv'] / 10, 0.0)
        tmp = mu[sl]
        tmp[target] = target_mu[target]
        embedded = np.zeros(tuple(shape), bool)
        embedded[sl] = target
        other = scene & ~embedded
        extra_ratio = float(other.sum() / target.sum())
        zcoords = np.arange(shape[2]) * h + origin[2]
        halves = [int(scene[:, :, zcoords < 0].sum()), int(scene[:, :, zcoords >= 0].sum())]
        survey_path = Path(_release_expand('@DENTAL_INPUT_ROOT@/artifacts/CROWN_pop/survey')) / f"vol_{meta['vol']:03d}.json"
        survey = json.loads(survey_path.read_text())
        visible = []
        for tooth in survey['instances']:
            point = (np.array(tooth['cen_mm']) - center) @ R.T
            grid = (point - origin) / h
            if np.all(grid >= 0) and np.all(grid < shape):
                visible.append({'k': tooth['k'], 'upper_proxy': tooth['upper'], 'lower_proxy': tooth['lower'], 'centroid_local_mm': point.tolist()})
        base = {angle: project(mu, h, angle) + 0.02683 * 10 for angle in pre['operation']['angles_deg']}
        base_file = outdir / f'{tid}_scene.npz'
        np.savez_compressed(base_file, source_tooth_union=scene, target_mask=embedded, mu_mm_inv=mu, spacing_mm=h, origin_mm=origin)
        records = []
        max_linear = 0.0
        same_loss = True
        faults = []
        for site in pre['operation']['sites']:
            for name in pre['operation']['lesions']:
                original = DATA / 'private' / f'{tid}_s1.2_{site}_{name}.npz'
                old = np.load(original)
                basis = old['unit_loss']
                loss = 0.4 * basis
                delta = np.zeros(mu.shape)
                delta[sl] = loss * (target_mu - 0.02683)
                field_path = outdir / f'{tid}_{site}_{name}_truth.npz'
                np.savez_compressed(field_path, local_loss=loss, spacing_mm=h, scene_offset_vox=pad)
                read = np.load(field_path)
                same_loss &= np.array_equal(read['local_loss'], loss)
                views = []
                for angle in pre['operation']['angles_deg']:
                    deficit = project(delta, h, angle)
                    tau = base[angle] - deficit
                    ex = transmission(tau, h)
                    healthy = transmission(base[angle], h)
                    image = noisy_image(ex, 20000, rng)
                    direct = project(mu - delta, h, angle) + 0.02683 * 10
                    error = float(np.max(np.abs(direct - tau)))
                    max_linear = max(max_linear, error)
                    if angle == 0:
                        errsum = float(np.max(np.abs(deficit - delta.sum(axis=1) * h)))
                        max_linear = max(max_linear, errsum)
                        bad_error = float(np.max(np.abs(deficit * 10 - delta.sum(axis=1) * h)))
                        faults.append({'site': site, 'class': name, 'deficit_units_x10_rejected': bad_error > 1e-08, 'error': bad_error})
                    iso_base = project(target_mu, h, angle)
                    iso_delta = project(loss * (target_mu - 0.02683), h, angle)
                    isolated_healthy = transmission(iso_base, h)
                    isolated = transmission(iso_base - iso_delta, h)
                    datafile = outdir / f'{tid}_{site}_{name}_{angle}deg.npz'
                    np.savez_compressed(datafile, image=image.astype(np.float32), healthy_expected=healthy, lesion_expected=ex, deficit_tau=deficit, isolated_healthy=isolated_healthy, isolated_lesion=isolated)
                    signal_counts = (ex - healthy) * 20000
                    variance = (ex + healthy) * 20000 + 18
                    dprime = float(np.sqrt(np.sum(signal_counts ** 2 / variance)))
                    views.append({'angle_deg': angle, 'path': str(datafile), 'superposition_error': error, 'max_scene_transmission_change': float(np.max(ex - healthy)), 'max_isolated_transmission_change': float(np.max(isolated - isolated_healthy)), 'paired_healthy_template_dprime_gaussian_approx': dprime, 'meaning': 'Known healthy template and lesion template ideal-observer bound, not image-only detector result'})
                    manifest.append({'path': str(datafile), 'sha256': sha(datafile), 'bytes': datafile.stat().st_size})
                records.append({'site': site, 'class': name, 'views': views, 'original_loss_file': str(original), 'original_loss_sha256': sha(original), 'embedded_loss_file': str(field_path), 'lesion_voxels': int((loss > 0).sum())})
                manifest.append({'path': str(field_path), 'sha256': sha(field_path), 'bytes': field_path.stat().st_size})
        scene_read = np.load(base_file)
        scene_persistence_error = float(np.max(np.abs(scene_read['mu_mm_inv'] - mu)))
        bad_loss = loss.copy()
        bad_loss.flat[np.flatnonzero(bad_loss)[0]] += 0.01
        empty_context = np.zeros_like(other)
        results.append({'tid': tid, 'source_mask': str(src), 'source_mask_sha256': sha(src), 'parent_case': str(parent), 'parent_case_sha256': sha(parent), 'frame': {'local_rotation_rows': R, 'target_center_parent_mm': center, 'parent_crop_lo_vox': croplo, 'native_spacing_mm': spacing, 'scene_origin_local_mm': origin, 'scene_grid_shape': shape, 'roundtrip_error_mm': frame_error}, 'target_source_containment_fraction': float(occupied), 'other_tooth_volume_ratio': extra_ratio, 'source_foreground_axis_halfspace_voxels': halves, 'visible_centroids_proxies': visible, 'source_geometry_status': 'Actual published mask anatomy; tooth/jaw type at centroids is predecessor geometric proxy', 'scene_export_max_mu_error_mm_inv': scene_persistence_error, 'max_superposition_error': max_linear, 'all_loss_arrays_bitidentical': bool(same_loss), 'gates': {'F_frame': bool(hashpass and frame_error <= 1e-09), 'M_mask': bool(occupied >= 0.8), 'C_context': bool(extra_ratio >= 1 and min(halves) > 0), 'L_linear_projection': bool(max_linear <= 1e-08), 'S_same_lesion': bool(same_loss), 'R_realism': 'UNKNOWN'}, 'faults': {'native_translation10mm_frame_rejected': wrong_frame_error > 1e-09, 'native_translation10mm_mask_rejected': wrong_occupancy < 0.8, 'wrong_mask_containment': wrong_occupancy, 'wrong_frame_error_mm': wrong_frame_error, 'modified_loss_array_rejected': not np.array_equal(bad_loss, loss), 'empty_context_rejected': float(empty_context.sum() / target.sum()) < 1.0, 'deficit_faults': faults}, 'lesions': records})
        manifest.append({'path': str(base_file), 'sha256': sha(base_file), 'bytes': base_file.stat().st_size})
        print(tid, results[-1]['gates'], 'occupancy', occupied, 'context ratio', extra_ratio, flush=True)
        state('R3_RUNNING', f'{tid} native context built', 'finish second source and calibration package')
    report = {'round': 'R3', 'claim_type': 'capability', 'prereg_sha256': sha(ROOT / 'PREREG_R3_NATIVE_TOOTH_CONTEXT.json'), 'cases': results, 'array_manifest': manifest, 'gates': {key: all((r['gates'][key] for r in results)) for key in ['F_frame', 'M_mask', 'C_context', 'L_linear_projection', 'S_same_lesion']}, 'realism': 'UNKNOWN', 'external_referent': {'kind': 'published_dataset', 'locator': 'https://zenodo.org/records/10597292', 'compared_quantity': 'Same-case original tooth-union mask occupancy in reconstructednativeframe', 'refutes_us': not all((r['gates']['M_mask'] and r['gates']['F_frame'] for r in results))}, 'cost': {'wall_s': time.perf_counter() - tic, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'GPU_s': 0}, 'scope': 'Bounded actual multi-tooth anatomy; homogeneous other-tooth material, constant10mmwater background, uncalibrated parallel beam. No bone/soft-tissue morphology or full clinical detector validation.'}
    write(ROOT / 'raw/R3_results.json', report)
    state('R3_COMPLETE', str(report['gates']), 'write reviewable demo, frozen lab predictions and final handoff')
    return report
if __name__ == '__main__':
    run()
