import csv
import json
import resource
import time
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from common import *

def select_threshold(prob, y):
    choices = np.r_[0.0, np.unique(prob), 1.0 + 1e-09]
    feasible = []
    for t in choices:
        sp = float(np.mean(prob[~y] < t))
        se = float(np.mean(prob[y] >= t))
        if sp >= 0.9:
            feasible.append((se, sp, float(t)))
    return max(feasible, key=lambda a: (a[0], a[1], -a[2]))[2]

def score_rows(rows, prob, threshold, split='test'):
    r = [(a, float(p)) for (a, p) in zip(rows, prob) if a['split'] == split]
    negative = [p < threshold for (a, p) in r if a['class'] == 'H0']
    out = {'threshold': float(threshold), 'split': split, 'specificity': float(np.mean(negative)), 'specificity_wilson': wilson(sum(negative), len(negative)), 'healthy_n': len(negative), 'by_depth': {}, 'subgroups': {}, 'uncertainty_warning': 'Wilson intervals treat simulation images independently. Shared geometry and noise replicates are correlated; these are not clinical confidence intervals.'}
    for c in CLASSES[1:]:
        hits = [p >= threshold for (a, p) in r if a['class'] == c]
        out['by_depth'][c] = {'sensitivity': float(np.mean(hits)), 'n': len(hits), 'tp': int(sum(hits)), 'wilson': wilson(sum(hits), len(hits))}
    y = [a['class'] != 'H0' for (a, p) in r]
    p = [p for (a, p) in r]
    out['auc'] = float(roc_auc_score(y, p))
    for group in ['site', 'peak_loss', 'angle_deg', 'photons', 'shell_mm', 'tid']:
        out['subgroups'][group] = []
        for val in sorted(set((a[group] for (a, p) in r))):
            sub = [(a, p) for (a, p) in r if a[group] == val]
            neg = [p < threshold for (a, p) in sub if a['class'] == 'H0']
            dep = {c: float(np.mean([p >= threshold for (a, p) in sub if a['class'] == c])) for c in CLASSES[1:]}
            out['subgroups'][group].append({'value': val, 'n': len(sub), 'specificity': float(np.mean(neg)), **dep})
    return out

def physical_controls(pre):
    h = 0.15
    mu = 0.297
    n = 20
    slab = np.full((80, n, 12), mu)
    errors = []
    for angle in pre['physics']['angles_deg']:
        p = project(slab, h, angle)
        expected = mu * n * h / np.cos(np.deg2rad(angle))
        errors.append(float(np.max(np.abs(p[20:60] - expected))))
    injected = project(slab, h, 0) * 10
    slab_fault = float(np.max(np.abs(injected[20:60] - mu * n * h)))
    rng = np.random.default_rng(24422)
    draws = rng.poisson(500, 20000) + rng.normal(0, 3, 20000)
    meanerr = abs(float(draws.mean()) / 500 - 1)
    varerr = abs(float(draws.var(ddof=1)) / 509 - 1)
    bad = np.full(20000, 500.0) + rng.normal(0, 3, 20000)
    badvar = abs(float(bad.var(ddof=1)) / 509 - 1)
    facit = pre['gates']['G4_external_attenuation']
    model = [2.75 * 0.99, 1.44 * 0.99 + 0.5 * 0.24]
    measured = [facit['enamel_mean_sd_cm_inv'], facit['dentin_mean_sd_cm_inv']]
    standardized = [abs(m - mn) / sd for (m, (mn, sd)) in zip(model, measured)]
    wrong_standardized = [abs(m * 10 - mn) / sd for (m, (mn, sd)) in zip(model, measured)]
    return {'G2_analytic_slab': {'pass': max(errors) <= 1e-08, 'max_error': max(errors), 'quantity': 'optical depth'}, 'G3_noise': {'pass': meanerr <= 0.01 and varerr <= 0.05, 'mean_fraction_error': meanerr, 'variance_fraction_error': varerr}, 'G4_external_attenuation': {'pass': max(standardized) <= 2, 'model_cm_inv': model, 'external_mean_sd_cm_inv': measured, 'absolute_standardized_residual': standardized, 'locator': 'https://pubmed.ncbi.nlm.nih.gov/8819355/', 'limitation': '40keV effective scenario vs70kVp film measurements, broad SD. Compatibility only; not spectrum calibration.'}, 'injections': {'slab_units_x10_rejected': slab_fault > 1e-08, 'slab_units_x10_error': slab_fault, 'non_poisson_noise_rejected': badvar > 0.05, 'bad_noise_variance_fraction_error': badvar, 'attenuation_units_x10_rejected': max(wrong_standardized) > 2}}

def run():
    tic = time.perf_counter()
    cpu = time.process_time()
    pre = json.loads((ROOT / 'PREREG_R1_SYNTHETIC_BENCHMARK.json').read_text())
    for sub in ['public', 'private']:
        (DATA / sub).mkdir(parents=True, exist_ok=True)
    rows = []
    public_rows = []
    patches = []
    truth_checks = []
    rejects = []
    manifest = []
    examples = []
    physical = physical_controls(pre)
    max_parity = 0.0
    rng = np.random.default_rng(pre['physics']['seed'])
    for item in pre['inputs']:
        assert sha(item['path']) == item['sha256'], 'Predecessor hash changed'
        a = np.load(item['path'])
        mask = a['labels'] > 0
        h = float(a['spacing_mm'][0])
        origin = a['origin_mm']
        block = []
        for shell in pre['tissue_model']['enamel_shell_mm']:
            (labels, dt) = tissue(mask, h, shell)
            mu = coefficient(labels)
            gp = DATA / 'private' / f"{item['tid']}_shell{shell:.1f}.npz"
            np.savez_compressed(gp, labels=labels, spacing_mm=h, origin_mm=origin)
            manifest.append({'path': str(gp), 'sha256': sha(gp), 'bytes': gp.stat().st_size, 'role': 'synthetic_tissue_on_real_mask'})
            base = {angle: project(mu, h, angle) for angle in pre['physics']['angles_deg']}
            max_parity = max(max_parity, float(np.max(np.abs(base[0] - mu.sum(axis=1) * h))))
            for site in ['proximal', 'occlusal']:
                try:
                    (axis, center) = choose_site(labels, site, h, origin, shell)
                except ValueError as e:
                    rejects.append({'tid': item['tid'], 'shell_mm': shell, 'site': site, 'reason': str(e)})
                    continue
                for name in CLASSES:
                    unit_loss = np.zeros(labels.shape, np.float32) if name == 'H0' else lesion(labels, h, site, center, name)
                    truth = verify_truth(labels, unit_loss, h, site)
                    passed = truth['actual_class'] == name and (not (truth['pulp_overlap'] or truth['outside_overlap']))
                    truth_checks.append({'tid': item['tid'], 'shell_mm': shell, 'site': site, 'requested_class': name, **truth, 'pass': passed})
                    tp = DATA / 'private' / f"{item['tid']}_s{shell}_{site}_{name}.npz"
                    np.savez_compressed(tp, unit_loss=unit_loss, spacing_mm=h, center_column=center, site=site)
                    manifest.append({'path': str(tp), 'sha256': sha(tp), 'bytes': tp.stat().st_size, 'role': 'voxel_mineral_loss_basis'})
                    saved = np.load(tp)
                    readback = verify_truth(labels, saved['unit_loss'], h, site)
                    if readback != truth:
                        raise ValueError('Saved voxel truth changed during serialization')
                    if not passed:
                        rejects.append(truth_checks[-1])
                        continue
                    delta = unit_loss * (mu - 0.02683)
                    deficit = {angle: project(delta, h, angle) for angle in pre['physics']['angles_deg']}
                    for peak in pre['physics']['mean_peak_fractional_mineral_loss']:
                        for angle in pre['physics']['angles_deg']:
                            tau = base[angle] - peak * deficit[angle]
                            expected = transmission(tau, h)
                            cen = roi_center(center, site, labels, angle)
                            if len(examples) < 6 and item['split'] == 'test' and (shell == 1.2) and (site == 'proximal') and (peak == 0.4) and (angle == 0):
                                ex = DATA / 'private' / f'EXAMPLE_{name}.npz'
                                np.savez_compressed(ex, expected=expected, mu=mu, loss=peak * unit_loss, labels=labels, h=h, center=cen)
                                examples.append({'class': name, 'path': str(ex)})
                            for photons in pre['physics']['incident_photons_per_pixel']:
                                for rep in range(pre['physics']['noise_replicates']):
                                    im = noisy_image(expected, photons, rng).astype(np.float32)
                                    serialized = im.astype(np.float16)
                                    index = len(block)
                                    block.append(serialized)
                                    patches.append(crop(serialized.astype(np.float32), cen))
                                    rid = sha_id = len(rows)
                                    sid = __import__('hashlib').sha256(f'X22-{rid}-20261002'.encode()).hexdigest()[:20]
                                    common = {'sample_id': sid, 'block': f"{item['tid']}.npz", 'array_index': index, 'split': item['split'], 'roi_center_xz_pixel': cen, 'pixel_mm': h, 'photons': photons, 'angle_deg': angle}
                                    public_rows.append(common)
                                    rows.append({**common, 'class': name, 'tid': item['tid'], 'site': site, 'shell_mm': shell, 'peak_loss': peak, 'noise_replicate': rep, 'voxel_truth_path': str(tp), 'actual_depth_mm': truth['max_depth_mm']})
            state('R1_GENERATING', f"{item['tid']}, shell {shell}mm complete; failed placements {len(rejects)}", 'finish benchmark and train frozen reference')
        bp = DATA / 'public' / f"{item['tid']}.npz"
        np.savez_compressed(bp, images=np.stack(block))
        manifest.append({'path': str(bp), 'sha256': sha(bp), 'bytes': bp.stat().st_size, 'role': 'public_images'})
        print(item['tid'], len(block), 'images', flush=True)
    write(DATA / 'public/index.json', public_rows)
    write(DATA / 'private/truth_index.json', rows)
    write(ROOT / 'raw/R1_voxel_truth_checks.json', truth_checks)
    write(ROOT / 'raw/R1_rejected_placements.json', rejects)
    x = np.asarray([features(a) for a in patches])
    y = np.array([r['class'] != 'H0' for r in rows])
    training = np.array([r['split'] == 'train' for r in rows])
    val = np.array([r['split'] == 'validation' for r in rows])
    sc = StandardScaler().fit(x[training])
    xs = sc.transform(x)
    timing = {}
    results = {}
    preds = {}
    for (name, model, xx) in [('reference_logistic', LogisticRegression(C=1, max_iter=1000, random_state=22), xs), ('equal_information_trees', ExtraTreesClassifier(n_estimators=150, min_samples_leaf=4, n_jobs=1, random_state=22), x)]:
        t = time.perf_counter()
        model.fit(xx[training], y[training])
        fit = time.perf_counter() - t
        t = time.perf_counter()
        prob = model.predict_proba(xx)[:, 1]
        query = time.perf_counter() - t
        thresh = select_threshold(prob[val], y[val])
        results[name] = score_rows(rows, prob, thresh)
        timing[name] = {'fit_s': fit, 'query_all_s': query}
        preds[name] = [{'sample_id': r['sample_id'], 'probability': float(p)} for (r, p) in zip(rows, prob)]
        write(ROOT / 'raw' / f'R1_{name}_predictions.json', preds[name])
    reference = results['reference_logistic']
    dep = reference['by_depth']
    proximal = next((r for r in reference['subgroups']['site'] if r['value'] == 'proximal'))
    se_e = np.mean([proximal[c] for c in ['E1', 'E2']])
    se_d = np.mean([proximal[c] for c in ['D1', 'D2', 'D3']])
    clinical = {'locator': 'https://doi.org/10.1016/j.jdsr.2024.02.001', 'enamel_sensitivity': float(se_e), 'dentin_sensitivity': float(se_d), 'specificity': reference['specificity'], 'external_intervals': {'enamel_sensitivity': [0.66, 0.75], 'dentin_sensitivity': [0.8, 0.87], 'specificity': [0.75, 0.96]}, 'pass': bool(0.66 <= se_e <= 0.75 and 0.8 <= se_d <= 0.87 and (0.75 <= proximal['specificity'] <= 0.96)), 'proximal_specificity': proximal['specificity'], 'interpretation': 'Different task, synthetic truth, populations and detector. Gate diagnoses lack of performance resemblance; cannot establish physical realism.'}
    checks = physical.pop('injections')
    tcheck = next((c for c in truth_checks if c['requested_class'] == 'E2' and c['pass']))
    stem = f"{tcheck['tid']}_s{tcheck['shell_mm']}_{tcheck['site']}_E2.npz"
    q = np.load(DATA / 'private' / stem)
    lbl = np.load(DATA / 'private' / f"{tcheck['tid']}_shell{tcheck['shell_mm']:.1f}.npz")['labels']
    replay = verify_truth(lbl, q['unit_loss'], float(q['spacing_mm']), str(q['site']))
    poisoned = q['unit_loss'].copy()
    pulp = np.argwhere(lbl == 3)
    if len(pulp):
        poisoned[tuple(pulp[0])] = 0.4
    badtruth = verify_truth(lbl, poisoned, float(q['spacing_mm']), str(q['site']))
    checks['wrong_E2_label_D3_rejected'] = replay['actual_class'] != 'D3'
    checks['pulp_loss_rejected'] = badtruth['pulp_overlap'] > 0
    checks['axis_projector_x10_rejected'] = float(np.max(np.abs(base[0] * 10 - mu.sum(axis=1) * h))) > 1e-08
    checks['clinical_sensitivity_0_rejected'] = not (0.66 <= 0 <= 0.75 and 0.8 <= se_d <= 0.87)
    bytes_total = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    report = {'round': 'R1', 'claim_type': 'capability', 'prereg_sha256': sha(ROOT / 'PREREG_R1_SYNTHETIC_BENCHMARK.json'), 'n_images': len(rows), 'n_real_geometries': len(pre['inputs']), 'n_lesion_truths': len(truth_checks), 'rejected_placements': rejects, 'gates': {'G1_voxel_truth': {'pass': all((c['pass'] for c in truth_checks)) and (not rejects), 'saved_field_replay': replay}, 'G2_axis_sum': {'pass': max_parity <= 1e-08, 'max_abs_error': max_parity}, **physical, 'G5_realism': {'outcome': 'UNKNOWN', 'pass': False, 'reason': 'Matched specimen mineral loss vs clinical bitewing contrast unavailable; shell, spectrum, PSF and scatter are uncalibrated'}, 'G6_clinical_compatibility': clinical}, 'detectors': results, 'injected_faults': checks, 'all_faults_rejected': all((checks[k] for k in checks if k.endswith('_rejected'))), 'examples': examples, 'array_manifest': manifest, 'disk_bytes': bytes_total, 'cost': {'wall_s': time.perf_counter() - tic, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'detectors': timing, 'preparation_search_s': 'UNKNOWN', 'physical_validation': 'NOT_DONE', 'token_cost': 'UNKNOWN', 'GPU_s': 0}, 'external_referent': {'kind': 'independent_measurement', 'locator': 'https://pubmed.ncbi.nlm.nih.gov/8819355/', 'compared_quantity': 'sound enamel/dentin LAC at70kVp cm^-1 compared with40keV effective model', 'refutes_us': not physical['G4_external_attenuation']['pass']}, 'external_referents_additional': [{'kind': 'independent_measurement', 'locator': clinical['locator'], 'compared_quantity': 'proximal enamel/dentin sensitivity and overall specificity', 'refutes_us': not clinical['pass']}], 'limitations': ['Synthetic DEJ/pulp and lesion classes conditional on closures', 'One-tooth ROI, not full clinical bitewing', 'No sensor, spectrum or paired contrast calibration', 'Not a detector superiority claim', 'No clinical validation or recommendation']}
    write(ROOT / 'raw/R1_results.json', report)
    assert report['all_faults_rejected']
    assert bytes_total < 3000000000
    state('R1_COMPLETE', f"{len(rows)} images; realism UNKNOWN; clinical resemblance {clinical['pass']}", 'Freeze and execute projection-equivalent different-depth lesions')
    print('R1 complete', report['n_images'], 'images; all fault injections rejected:', report['all_faults_rejected'])
    return report
if __name__ == '__main__':
    run()
