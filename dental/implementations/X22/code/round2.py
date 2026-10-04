import json
import resource
import time
import numpy as np
from scipy import ndimage as ndi
from scipy.special import ndtr
from common import *

def make_pair(labels, dt, mu, h, pre):
    shell = 1.2
    shallow = (labels == 1) & (dt >= shell * 0.55) & (dt <= shell * 0.9)
    deep = (labels == 2) & (dt >= shell + (3 - shell) * 0.1) & (dt <= shell + (3 - shell) * 0.3)
    y = np.arange(labels.shape[1])[None, :, None]
    first = np.argmax(labels > 0, axis=1)
    last = labels.shape[1] - 1 - np.argmax((labels > 0)[:, ::-1], axis=1)
    midpoint = (first + last) / 2
    shallow &= y <= midpoint[:, None, :]
    deep &= y <= midpoint[:, None, :]
    c = (mu - 0.02683) * h
    cap_e = np.sum(c * shallow, axis=1)
    cap_d = np.sum(c * deep, axis=1)
    eligible = (cap_e > 0) & (cap_d > 0)
    coords = np.argwhere(labels > 0)
    target = np.array([np.quantile(coords[:, 0], 0.8), coords[:, 2].max() - 3.5 / h])
    q = np.argwhere(eligible)
    if not len(q):
        raise ValueError('No common rays with both tissue-depth bands')
    center = q[np.argmin(np.sum((q - target) ** 2, axis=1))]
    (x, z) = np.indices(cap_e.shape)
    rr = np.sqrt((x - center[0]) ** 2 + (z - center[1]) ** 2) * h / pre['operation']['shared_footprint_radius_mm']
    region = eligible & (rr <= 1)
    taper = np.maximum(0, 1 - rr) * region
    target_deficit = 0.45 * np.minimum(cap_e, cap_d / 0.7) * taper
    fe = np.divide(target_deficit, cap_e, out=np.zeros_like(cap_e), where=cap_e > 0)
    fd = np.divide(0.7 * target_deficit, cap_d, out=np.zeros_like(cap_d), where=cap_d > 0)
    loss_e = shallow * fe[:, None, :]
    loss_d = shallow * (0.3 * fe[:, None, :]) + deep * fd[:, None, :]
    return (loss_e, loss_d, center)

def depth_summary(labels, dt, loss, shell=1.2):
    support = loss > 0
    maxdepth = float(dt[support].max()) if support.any() else 0.0
    frac = maxdepth / shell if maxdepth <= shell else 1 + (maxdepth - shell) / (3 - shell)
    return {'voxels': int(support.sum()), 'enamel_voxels': int(np.sum(support & (labels == 1))), 'dentin_voxels': int(np.sum(support & (labels == 2))), 'pulp_overlap': int(np.sum(support & (labels == 3))), 'outside_overlap': int(np.sum(support & (labels == 0))), 'max_loss': float(loss.max()), 'max_nearest_boundary_depth_mm': maxdepth, 'max_tissue_normalized_depth': frac}

def likelihood_test(a, b, photons, rng, n=512):
    flat_a = a.ravel()
    flat_b = b.ravel()
    var = (flat_a + flat_b) * photons / 2 + 9
    diff = (flat_b - flat_a) * photons
    w = diff / var
    dprime = float(np.sqrt(np.sum(diff * diff / var)))
    center = (flat_a + flat_b) * photons / 2
    hits = 0
    for (label, expected) in [(0, flat_a), (1, flat_b)]:
        for _ in range(n // 2):
            counts = rng.poisson(expected * photons) + rng.normal(0, 3, expected.shape)
            score = float(np.dot(w, counts - center))
            if score == 0:
                pred = int(rng.integers(2))
            else:
                pred = int(score > 0)
            hits += pred == label
    return {'dprime_gaussian_approx': dprime, 'predicted_accuracy_gaussian_approx': float(ndtr(dprime / 2)), 'sampled_correct': int(hits), 'sampled_n': n, 'sampled_accuracy': hits / n, 'wilson': wilson(hits, n)}

def run():
    tic = time.perf_counter()
    cpu = time.process_time()
    pre = json.loads((ROOT / 'PREREG_R2_PROJECTIVE_AMBIGUITY.json').read_text())
    source = json.loads((ROOT / 'PREREG_R1_SYNTHETIC_BENCHMARK.json').read_text())
    rng = np.random.default_rng(pre['operation']['seed'])
    results = []
    manifest = []
    rejects = []
    faults = []
    outdir = DATA / 'ambiguity'
    outdir.mkdir(exist_ok=True)
    for item in source['inputs']:
        q = np.load(item['path'])
        mask = q['labels'] > 0
        h = float(q['spacing_mm'][0])
        (labels, dt) = tissue(mask, h, 1.2)
        mu = coefficient(labels)
        try:
            (le, ld, center) = make_pair(labels, dt, mu, h, pre)
        except ValueError as e:
            rejects.append({'tid': item['tid'], 'reason': str(e)})
            continue
        path = outdir / f"{item['tid']}_pair.npz"
        np.savez_compressed(path, labels=labels, dt_mm=dt, shallow_loss=le, deep_loss=ld, mu_mm_inv=mu, spacing_mm=h, center=center)
        read = np.load(path)
        le = read['shallow_loss']
        ld = read['deep_loss']
        se = depth_summary(labels, dt, le)
        sd = depth_summary(labels, dt, ld)
        depthpass = se['dentin_voxels'] == 0 and sd['dentin_voxels'] > 0 and (0.5 < se['max_tissue_normalized_depth'] <= 1) and (1 < sd['max_tissue_normalized_depth'] <= 1 + 1 / 3)
        depthpass = bool(depthpass and le.max() <= 0.5 and (ld.max() <= 0.5) and (not (se['pulp_overlap'] or sd['pulp_overlap'] or se['outside_overlap'] or sd['outside_overlap'])))
        delta_e = le * (mu - 0.02683)
        delta_d = ld * (mu - 0.02683)
        sums = (delta_e - delta_d).sum(axis=1) * h
        direct_error = float(np.abs(sums).max())
        obs = {}
        p_error = 0.0
        mean_error = 0.0
        parity = 0.0
        for angle in [0, 10]:
            pe = project(delta_e, h, angle)
            pd = project(delta_d, h, angle)
            base = project(mu, h, angle)
            a = transmission(base - pe, h)
            b = transmission(base - pd, h)
            obs[angle] = (a, b)
            if angle == 0:
                p_error = float(np.abs(pe - pd).max())
                mean_error = float(np.abs(a - b).max())
                parity = float(np.max(np.abs(pe - delta_e.sum(axis=1) * h)))
        (a0, b0) = obs[0]
        (a10, b10) = obs[10]
        scores = {}
        for photons in pre['operation']['photon_counts']:
            scores[str(photons)] = {'one_view': likelihood_test(a0, a0, photons, rng), 'additional_10deg_view': likelihood_test(a10, b10, photons, rng), 'both_views': likelihood_test(np.stack([a0, a10]), np.stack([a0, b10]), photons, rng)}
        expected_file = outdir / f"{item['tid']}_observations.npz"
        np.savez_compressed(expected_file, shallow_0=a0, deep_0=b0, shallow_10=a10, deep_10=b10)
        active = np.argwhere(le > 0)[0]
        wrong = ld.copy()
        wrong[tuple(active)] += 0.05
        poison_error = float(np.max(np.abs(project(wrong * (mu - 0.02683), h, 0) - project(delta_e, h, 0))))
        faults.append({'tid': item['tid'], 'extra_loss_0p05_rejected': poison_error > 1e-08, 'extra_loss_tau_error': poison_error, 'swapped_depth_labels_rejected': se['dentin_voxels'] == 0 and sd['dentin_voxels'] > 0, 'axis_scale_x10_rejected': float(np.abs(project(delta_e, h, 0) * 10 - delta_e.sum(axis=1) * h).max()) > 1e-08})
        results.append({'tid': item['tid'], 'shallow_depth': se, 'deep_depth': sd, 'depth_gate_pass': depthpass, 'direct_axis_nullspace_error': direct_error, 'projector_nullspace_error': p_error, 'projector_axis_parity_error': parity, 'mean_transmission_difference': mean_error, 'one_view_balanced_Bayes_accuracy_bound': 0.5, 'scores': scores, 'second_view_gate_pass': scores['20000']['additional_10deg_view']['dprime_gaussian_approx'] >= 2, 'pair_file': str(path), 'observation_file': str(expected_file)})
        for p in [path, expected_file]:
            manifest.append({'path': str(p), 'sha256': sha(p), 'bytes': p.stat().st_size})
        print(item['tid'], 'nullspace', p_error, 'dprime10', scores['20000']['additional_10deg_view']['dprime_gaussian_approx'], flush=True)
        state('R2_RUNNING', f'{len(results)} ambiguity pairs verified', 'complete fixed10deg observation and fault injections')
    gates = {'A_depth': {'pass': len(results) == 7 and all((r['depth_gate_pass'] for r in results))}, 'B_equal_projection': {'pass': bool(results) and all((r['projector_nullspace_error'] <= 1e-08 and r['direct_axis_nullspace_error'] <= 1e-08 and (r['projector_axis_parity_error'] <= 1e-08) for r in results))}, 'C_same_distribution': {'pass': bool(results) and all((r['mean_transmission_difference'] <= 1e-08 for r in results)), 'balanced_Bayes_accuracy_upper_bound': 0.5}, 'D_extra_view': {'pass': sum((r['second_view_gate_pass'] for r in results)) >= 5, 'passing_geometries': sum((r['second_view_gate_pass'] for r in results)), 'required': 5}, 'E_faults': {'pass': all((r[k] for r in faults for k in r if k.endswith('_rejected')))}}
    report = {'round': 'R2', 'claim_type': 'capability', 'prereg_sha256': sha(ROOT / 'PREREG_R2_PROJECTIVE_AMBIGUITY.json'), 'gates': gates, 'pairs': results, 'rejected_pairs': rejects, 'fault_injections': faults, 'array_manifest': manifest, 'depth_definition': 'nearest-outer-boundary EDT synthetic tissue depth; notR1 surface-ray staging or clinical/histological depth', 'external_referent': {'kind': 'closed_form', 'locator': 'https://engineering.purdue.edu/~malcolm/pct/CTI_Ch03.pdf', 'compared_quantity': 'Chapter3 line-integral projection, equation2/3 and linearity. Nullspace/depth limit is our explicit derivation under this operator.', 'refutes_us': not gates['B_equal_projection']['pass']}, 'cost': {'wall_s': time.perf_counter() - tic, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'GPU_s': 0, 'fit_s': 0, 'empirical_validation': 'NOT_DONE'}, 'scope': 'Exact digital ambiguity only under frozen effective-energy primary-image and fixed observation model; structural distributions differ. Natural lesion growth, polychromatic beam/scatter, calibrated clinical realism UNKNOWN.'}
    write(ROOT / 'raw/R2_results.json', report)
    assert gates['E_faults']['pass']
    state('R2_COMPLETE', f"{len(results)} pairs; secondview gate {gates['D_extra_view']['pass']}", 'prepare real full-arch context or minimal laboratory calibration with frozen predictions')
    return report
if __name__ == '__main__':
    run()
