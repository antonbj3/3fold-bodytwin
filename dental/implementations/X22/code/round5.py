import json
import resource
import time
import numpy as np
from common import *
from round4 import phases

def run():
    tic = time.perf_counter()
    cpu = time.process_time()
    r2 = json.loads((ROOT / 'raw/R2_results.json').read_text())
    rows = []
    manifest = []
    out = DATA / 'mineral_ambiguity'
    out.mkdir(exist_ok=True)
    for pair in r2['pairs']:
        q = np.load(pair['pair_file'])
        labels = q['labels']
        h = float(q['spacing_mm'])
        le = q['shallow_loss']
        old_ld = q['deep_loss']
        C = np.choose(labels, [0.0, 2.75, 1.44, 0.0])
        mass_e = (C * le).sum(axis=1)
        shape_d = old_ld * (labels == 2)
        mass_d = (C * shape_d).sum(axis=1)
        scale = np.divide(0.7 * mass_e, mass_d, out=np.zeros_like(mass_d), where=mass_d > 0)
        ld = 0.3 * le + shape_d * scale[:, None, :]
        mass_error = float(np.max(np.abs((C * (le - ld)).sum(axis=1) * h / 10)))
        endpoint = []
        old_errors = []
        for rho in [2.99, 3.15]:
            (*_, mu_e, damaged_e) = phases(labels, le, rho)
            (*_, mu_d, damaged_d) = phases(labels, ld, rho)
            (*_, mu_old, damaged_old) = phases(labels, old_ld, rho)
            base = project(mu_e, h, 0)
            old_error = float(np.max(np.abs(project(damaged_e, h, 0) - project(damaged_old, h, 0))))
            old_errors.append(old_error)
            tau_e = project(damaged_e, h, 0)
            tau_d = project(damaged_d, h, 0)
            equal_error = float(np.max(np.abs(tau_e - tau_d)))
            mean_error = float(np.max(np.abs(transmission(tau_e, h) - transmission(tau_d, h))))
            a = transmission(project(damaged_e, h, 10), h)
            b = transmission(project(damaged_d, h, 10), h)
            dprime = float(np.sqrt(np.sum(((b - a) * 20000) ** 2 / ((a + b) * 10000 + 9))))
            bad = ld.copy()
            idx = tuple(np.argwhere(le > 0)[0])
            bad[idx] += 0.05
            (*_, bad_mu, bad_damaged) = phases(labels, bad, rho)
            fault_mass = float(np.max(np.abs((C * (le - bad)).sum(axis=1) * h / 10)))
            fault_tau = float(np.max(np.abs(project(damaged_e, h, 0) - project(bad_damaged, h, 0))))
            endpoint.append({'rho_mineral_g_cm3': rho, 'optical_depth_difference': equal_error, 'mean_transmission_difference': mean_error, 'extra_10deg_dprime_20000': dprime, 'extra_view_pass': dprime >= 2, 'old_R2_pair_new_material_optical_depth_difference': old_error, 'fault_mass_rejected': fault_mass > 1e-12, 'fault_projection_rejected': fault_tau > 1e-08})
        path = out / f"{pair['tid']}_mass_pair.npz"
        np.savez_compressed(path, labels=labels, shallow_loss=le, deep_loss=ld, C_mineral_g_cm3=C, spacing_mm=h)
        manifest.append({'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size})
        bounded = bool(le.max() <= 0.5 and ld.max() <= 0.5 and (not np.any((ld > 0) & ((labels == 0) | (labels == 3)))))
        rows.append({'tid': pair['tid'], 'mass_projection_error_g_cm2': mass_error, 'bounded_fields_pass': bounded, 'endpoints': endpoint, 'old_R2_exact_equality_transfers_to_new_material': max(old_errors) <= 1e-08})
        print(pair['tid'], 'mass error', mass_error, 'old transfer', max(old_errors), 'new dprime', [e['extra_10deg_dprime_20000'] for e in endpoint], flush=True)
    gates = {'mass_projection': all((r['mass_projection_error_g_cm2'] <= 1e-12 for r in rows)), 'bounded_fields': all((r['bounded_fields_pass'] for r in rows)), 'rho_family': all((e['optical_depth_difference'] <= 1e-08 and e['mean_transmission_difference'] <= 1e-08 for r in rows for e in r['endpoints'])), 'added_view': all((sum((r['endpoints'][i]['extra_view_pass'] for r in rows)) >= 5 for i in [0, 1])), 'faults': all((e['fault_mass_rejected'] and e['fault_projection_rejected'] for r in rows for e in r['endpoints']))}
    result = {'round': 'R5', 'claim_type': 'capability', 'prereg_sha256': sha(ROOT / 'PREREG_R5_MINERAL_AREAL_AMBIGUITY.json'), 'gates': gates, 'rows': rows, 'array_manifest': manifest, 'old_R2_transfer_failure_count': sum((not r['old_R2_exact_equality_transfers_to_new_material'] for r in rows)), 'external_referent': {'kind': 'closed_form', 'locator': 'https://physics.nist.gov/PhysRefData/Xcom/Text/chap3.html', 'compared_quantity': 'NIST mass-mixture attenuation andderivedmineralarealmass invariance across commoncompositionwaterreplacementfamily', 'refutes_us': not gates['rho_family']}, 'balanced_identical_image_Bayes_accuracy_bound': 0.5, 'empirical_realism': 'UNKNOWN', 'cost': {'wall_s': time.perf_counter() - tic, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'GPU_s': 0}}
    write(ROOT / 'raw/R5_results.json', result)
    assert gates['faults']
    state('R5_COMPLETE', str(gates), 'final stable one-command replay and labhandoff')
    return result
if __name__ == '__main__':
    run()
