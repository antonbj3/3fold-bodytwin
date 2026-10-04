"""Joint interfaces from raw local CBCT measurements, then sampled tool envelope."""
from dental_release.paths import expand as _release_expand
import json
import resource
import sys
import time
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.optimize import least_squares
from scipy.special import ndtr
from ipr import jsonable, sha, write_json
sys.path.insert(0, _release_expand('@DENTAL_INPUT_ROOT@/workspace/cells/design'))
import crown_case_sts3d as STS
RAY_IDS = [0, 4, 8, 36, 40, 44, 72, 76, 80]

def signal(s, p):
    (u, t, w, a, rise, fall) = p
    return a + rise * ndtr((s - u) / w) - fall * ndtr((s - u - t) / w)

def fit_profile(s, y, seed_t):
    yy = np.asarray(y) / 1000
    a = float(np.clip(np.median(yy[:20]), -1.49, 2.99))
    high = float(np.max(yy))
    tail = float(np.median(yy[-20:]))
    lower = [-1.5, 0.1, 0.05, -1.5, 0, 0]
    upper = [1.5, 3.0, 0.35, 3.0, 6, 6]
    fits = []
    for u in [0.0, -0.4]:
        p0 = [u, np.clip(seed_t, 0.2, 2.5), 0.15, a, np.clip(high - a, 0.1, 5.9), np.clip(high - tail, 0.15, 5.9)]
        fit = least_squares(lambda p: signal(s, p) - yy, p0, bounds=(lower, upper), x_scale='jac', max_nfev=200)
        fits.append(fit)
    fit = min(fits, key=lambda f: np.mean(f.fun ** 2))
    p = fit.x
    norm_rmse = float(np.sqrt(np.mean(fit.fun ** 2)) / max(p[4], p[5], 1e-06))
    bound_ok = all((p[i] - lower[i] > 0.01 and upper[i] - p[i] > 0.01 for i in [0, 1]))
    reasons = []
    if not fit.success:
        reasons.append('solver_not_converged')
    if norm_rmse > 0.1:
        reasons.append('normalized_RMSE')
    if p[4] * 1000 < 500:
        reasons.append('outer_contrast')
    if p[5] * 1000 < 150:
        reasons.append('DEJ_contrast')
    if p[1] / p[2] < 2:
        reasons.append('unresolved_plateau')
    if not bound_ok:
        reasons.append('interface_at_bound')
    if np.max(y) >= 3090:
        reasons.append('saturation')
    return {'parameters': p, 'thickness_mm': float(p[1]), 'outer_shift_mm': float(p[0]), 'PSF_sigma_mm': float(p[2]), 'normalized_RMSE': norm_rmse, 'accepted': not reasons, 'rejections': reasons, 'nfev': sum((f.nfev for f in fits)), 'fit_grey': signal(s, p) * 1000}

def plane_capacity(outer, inner, e_outer=0.15, e_dej=0.15, residual=0.5, tool=0.05):
    (outer, inner) = (np.asarray(outer), np.asarray(inner))
    q = max(np.max(inner + e_dej + residual), np.max((outer + inner + e_outer + e_dej) / 2)) + tool
    return (max(0.0, float(np.max(outer) - q)), float(q))

def plane_feasible(r, outer, inner, e_outer=0.15, e_dej=0.15, residual=0.5, tool=0.05):
    q = np.max(outer) - r - tool
    active = q < np.asarray(outer) + e_outer
    if not np.any(active):
        return True
    return bool(np.all(q - np.asarray(inner)[active] - e_dej >= residual - 1e-12) and np.all(q >= (np.asarray(outer)[active] + np.asarray(inner)[active] + e_outer + e_dej) / 2 - 1e-12))

def plane_bisection(outer, inner):
    (lo, hi) = (0.0, 5.0)
    if not plane_feasible(0, outer, inner):
        return 0.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if plane_feasible(mid, outer, inner):
            lo = mid
        else:
            hi = mid
    return lo

def run():
    (tic, cpu) = (time.perf_counter(), time.process_time())
    pre = json.loads(Path('PREREG_R2.json').read_text())
    r1 = json.loads(Path('rounds/R1_raw_profiles.json').read_text())
    selection = {s['tid']: s for s in json.loads(Path('PREREG_R1.json').read_text())['inputs']}
    (raw, patches) = ([], [])
    (cached_vol, img) = (None, None)
    s = np.arange(-2.0, 5.00001, 0.025)
    source_hashes = {}
    for patch in r1:
        item = selection[patch['tid']]
        with np.load(item['path']) as d:
            meta = json.loads(str(d['meta']))
        vol = meta['vol']
        if cached_vol != vol:
            img = None
            fn = f'{STS.STS_ROI}/Image/ROI_L_{vol:03d}.nii.gz'
            (img, zooms, affine, image_sha) = STS.read_nifti(fn)
            img = img.astype(np.float32)
            source_hashes[fn] = image_sha
            assert image_sha == meta['sha256'][f'Image/ROI_L_{vol:03d}.nii.gz']
            cached_vol = vol
        side = 1 if patch['side'] == 'mesial' else -1
        R = np.array(patch['frame']['rotation'])
        center = np.array(patch['frame']['center_mm'])
        (outer, inner, accepted) = ([], [], [])
        for i in RAY_IDS:
            if patch['thickness_mm'][i] is None:
                raw.append({'tid': patch['tid'], 'side': patch['side'], 'height_fraction': patch['height_fraction'], 'ray_id': i, 'accepted': False, 'rejections': ['R1_missing_enamel']})
                continue
            x0 = patch['enamel_outer_x_mm'][i]
            (y, z) = patch['yz_mm'][i]
            local = np.column_stack([x0 - side * s, np.full(len(s), y), np.full(len(s), z)])
            roi_idx = (local @ R + center) / zooms + np.array(meta['crop_lo_vox'])
            bounds_ok = np.all(roi_idx >= 1) and np.all(roi_idx < np.array(img.shape) - 2)
            grey = map_coordinates(img, roi_idx.T, order=1, mode='nearest', prefilter=False)
            fit = fit_profile(s, grey, patch['thickness_mm'][i])
            if not bounds_ok:
                fit['accepted'] = False
                fit['rejections'].append('outside_image')
            out = side * x0 - fit['outer_shift_mm']
            inn = out - fit['thickness_mm']
            record = {'tid': patch['tid'], 'side': patch['side'], 'height_fraction': patch['height_fraction'], 'ray_id': i, 's_mm': s, 'grey': grey, 'oriented_outer_mm': out, 'oriented_inner_mm': inn, 'original_mesh_thickness_mm': patch['thickness_mm'][i], **fit}
            raw.append(record)
            if fit['accepted']:
                outer.append(out)
                inner.append(inn)
                accepted.append(record)
        row = {'tid': patch['tid'], 'jaw': item['jaw'], 'tooth_type': item['type'], 'side': patch['side'], 'height_fraction': patch['height_fraction'], 'accepted_rays': len(accepted), 'total_rays': 9, 'image_thickness_mean_mm': float(np.mean([r['thickness_mm'] for r in accepted])) if accepted else None, 'positive_plane_output_admitted': len(accepted) == 9, 'calibrated_safe_ipr_mm': None}
        if len(accepted) == 9:
            (cap, q) = plane_capacity(outer, inner)
            ctrl = plane_bisection(outer, inner)
            row.update(scenario_plane_capacity_mm=cap, envelope_q_mm=q, control_capacity_mm=ctrl, parity_error_mm=abs(cap - ctrl), residual_pass=plane_feasible(cap, outer, inner), generic_0_5_pass=plane_feasible(0.5, outer, inner), injected_plus_1_mm_rejected=not plane_feasible(cap + 1, outer, inner))
        else:
            row.update(scenario_plane_capacity_mm=None, reason='incomplete_fit_coverage')
        patches.append(row)
        write_json('rounds/R2_profiles.partial.json', raw)
        write_json('CURRENT_WORK_STATE.json', {'stage': 'R2_RUNNING', 'last_patch': {k: row[k] for k in ['tid', 'side', 'height_fraction', 'accepted_rays']}, 'completed_patches': len(patches), 'next_operation': 'complete joint profile fit then evaluate frozen gates'})
        print(patch['tid'], patch['side'], patch['height_fraction'], 'accepted', len(accepted), flush=True)
    accepted = [r for r in raw if r['accepted']]
    middle = [r['thickness_mm'] for r in accepted if r['height_fraction'] == 0.5]
    median = float(np.median(middle)) if middle else None
    complete = [r for r in patches if r['positive_plane_output_admitted']]
    parity = max([r['parity_error_mm'] for r in complete], default=None)
    bad_counts = {}
    for record in raw:
        for reason in record['rejections']:
            bad_counts[reason] = bad_counts.get(reason, 0) + 1
    g1 = median is not None and 0.96 <= median <= 1.96
    gates = {'G1_external_population': {'pass': g1, 'median_mm': median, 'range_mm': [0.96, 1.96]}, 'G2_equal_information': {'pass': bool(complete) and parity <= 1e-07, 'max_difference_mm': parity}, 'G3_cut_feasibility': {'pass': bool(complete) and all((r['residual_pass'] for r in complete)), 'scope': 'nine sampled directional rays; no whole-surface guarantee'}, 'G4_fit_coverage': {'pass': len(accepted) / len(raw) >= 0.8, 'fraction': len(accepted) / len(raw), 'accepted_rays': len(accepted), 'planned_rays': len(raw), 'complete_patches': len(complete)}, 'G5_anatomical_admission': {'pass': False, 'outcome': 'UNKNOWN'}}
    flat = fit_profile(s, np.full_like(s, 800.0), 1.0)
    injections = {'external_plus_2_mm_rejected': median is None or not 0.96 <= median + 2 <= 1.96, 'control_plus_0_1_rejected': parity is None or abs(parity + 0.1) > 1e-07, 'zero_contrast_profile_rejected': not flat['accepted'], 'cut_plus_1_rejected_all_complete_patches': bool(complete) and all((r['injected_plus_1_mm_rejected'] for r in complete))}
    result = {'round': 'R2', 'prereg_sha256': sha('PREREG_R2.json'), 'gates': gates, 'injected_faults': injections, 'rejection_counts': bad_counts, 'inputs_image_sha256': source_hashes, 'external_referent': {'kind': 'independent_measurement', 'locator': 'https://pubmed.ncbi.nlm.nih.gov/32634888/', 'compared_quantity': 'molar proximal enamel thickness mm, different specimens', 'refutes_us': not g1}, 'cost': {'wall_seconds': time.perf_counter() - tic, 'cpu_seconds': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'fit_function_evaluations': sum((r.get('nfev', 0) for r in raw)), 'physical_validation': 'NOT_DONE'}}
    write_json('rounds/R2_profiles.json', raw)
    write_json('rounds/R2_patches.json', patches)
    write_json('rounds/R2_results.json', result)
    print(json.dumps(jsonable(result), indent=2))
    return result
if __name__ == '__main__':
    run()
