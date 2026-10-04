import json
import resource
import time
import numpy as np
from common import *

def phases(labels, loss, rho_m):
    C = np.choose(labels, [0.0, 2.75, 1.44, 0.0])
    Co = np.choose(labels, [0.0, 0.02 * 1.4, 0.5, 0.0])
    hard = (labels == 1) | (labels == 2)
    pm = C / rho_m
    po = Co / 1.4
    pw = np.where(hard, 1 - pm - po, np.where(labels == 3, 1.0, 0.0))
    if np.any(pm < 0) or np.any(po < 0) or np.any(pw < 0):
        raise ValueError('Impossible phase volume')
    removed = C * loss
    pm_after = (C - removed) / rho_m
    pw_after = pw + removed / rho_m
    mu = (0.99 * C + 0.24 * Co + 0.2683 * pw) / 10
    damaged = (0.99 * (C - removed) + 0.24 * Co + 0.2683 * pw_after) / 10
    return (C, Co, pm, po, pw, removed, pm_after, pw_after, mu, damaged)

def mass_fraction_mixture(C, Co, pw):
    total = C + Co + pw
    masses = np.stack([C, Co, pw])
    weights = np.divide(masses, total[None], out=np.zeros_like(masses), where=total[None] > 0)
    mac = np.einsum('i,ijkl->jkl', np.array([0.99, 0.24, 0.2683]), weights)
    return total * mac / 10

def run():
    tic = time.perf_counter()
    cpu = time.process_time()
    pre = json.loads((ROOT / 'PREREG_R4_MINERAL_INVENTORY.json').read_text())
    parent = json.loads((ROOT / 'PREREG_R1_SYNTHETIC_BENCHMARK.json').read_text())
    out = DATA / 'mineral_inventory'
    out.mkdir(exist_ok=True)
    rows = []
    manifest = []
    faults = []
    for item in parent['inputs']:
        q = np.load(DATA / 'private' / f"{item['tid']}_shell1.2.npz")
        labels = q['labels']
        h = float(q['spacing_mm'])
        loss = 0.4 * np.load(DATA / 'private' / f"{item['tid']}_s1.2_proximal_D2.npz")['unit_loss']
        volume_cm3 = h ** 3 / 1000
        for rho in [2.99, 3.15]:
            (C, Co, pm, po, pw, removed, pma, pwa, mu, damaged) = phases(labels, loss, rho)
            hard = (labels == 1) | (labels == 2)
            volume_error = float(np.abs((pma + po + pwa)[hard] - 1).max())
            mixture_error = float(np.max(np.abs(damaged - mass_fraction_mixture(C - removed, Co, pwa))))
            direct_delta = removed * (0.99 - 0.2683 / rho) / 10
            delta_error = float(np.max(np.abs(mu - damaged - direct_delta)))
            mineral_out = float(removed.sum() * volume_cm3)
            water_in = float((removed / rho).sum() * volume_cm3)
            mineral_drop = float(np.sum(C - (C - removed)) * volume_cm3)
            water_gain = float(np.sum(pwa - pw) * volume_cm3)
            ledger_error = max(abs(mineral_out - mineral_drop), abs(water_in - water_gain))
            healthy_tau = project(mu, h, 0)
            new_tau = project(damaged, h, 0)
            projection_error = float(np.max(np.abs(new_tau - (healthy_tau - project(direct_delta, h, 0)))))
            old_mu = coefficient(labels)
            old_damaged = old_mu - loss * (old_mu - 0.02683)
            old_delta = project(old_mu - old_damaged, h, 0)
            new_delta = healthy_tau - new_tau
            bad_pwa = pw.copy()
            omitted_water_error = float(np.abs((pma + po + bad_pwa)[hard] - 1).max())
            bad_mu = (0.99 * (C - removed) + 0.24 * Co * 10 + 0.2683 * pwa) / 10
            matrix_fault = float(np.max(np.abs(bad_mu - mass_fraction_mixture(C - removed, Co, pwa))))
            impossible_C = C.copy()
            impossible_C[hard] = 4
            impossible_volume = (impossible_C / rho + po)[hard]
            inventory_fault = abs(mineral_out - mineral_drop * 2)
            faults.append({'tid': item['tid'], 'rho_mineral_g_cm3': rho, 'omitted_water_rejected': omitted_water_error > 1e-12, 'omitted_water_phase_error': omitted_water_error, 'matrix_x10_rejected': matrix_fault > 1e-12, 'matrix_mu_error_mm_inv': matrix_fault, 'impossible_C4_rejected': bool(np.any(impossible_volume > 1)), 'mineral_ledger_x2_rejected': inventory_fault > 1e-12})
            fp = out / f"{item['tid']}_rho{rho:.2f}.npz"
            np.savez_compressed(fp, labels=labels, C_mineral_g_cm3=C, C_removed_g_cm3=removed, phi_water_before=pw, phi_water_after=pwa, mu_healthy_mm_inv=mu, mu_lesion_mm_inv=damaged, primary_healthy=np.exp(-healthy_tau), primary_lesion=np.exp(-new_tau), spacing_mm=h)
            manifest.append({'path': str(fp), 'sha256': sha(fp), 'bytes': fp.stat().st_size})
            saved = np.load(fp)
            persisted_error = float(np.max(np.abs(saved['C_removed_g_cm3'] - removed)))
            sound = [float(mu[labels == 1][0] * 10), float(mu[labels == 2][0] * 10)]
            external_z = [abs(sound[0] - 2.97) / 0.71, abs(sound[1] - 2.12) / 0.92]
            rows.append({'tid': item['tid'], 'rho_mineral_g_cm3': rho, 'removed_mineral_g': mineral_out, 'water_added_g': water_in, 'volume_fraction_error': volume_error, 'mass_ledger_error_g': ledger_error, 'mixture_error_mm_inv': mixture_error, 'delta_law_error_mm_inv': delta_error, 'projection_error': projection_error, 'persisted_mineral_error_g_cm3': persisted_error, 'sound_mu_cm_inv': sound, 'external_abs_SD': external_z, 'max_old_attenuation_deficit': float(old_delta.max()), 'max_new_attenuation_deficit': float(new_delta.max()), 'relative_change_max_deficit': float(new_delta.max() / old_delta.max() - 1), 'gates': {'I_inventory': bool(max(volume_error, ledger_error, persisted_error) <= 1e-12), 'A_attenuation': bool(max(mixture_error, delta_error) <= 1e-12), 'P_projection': bool(projection_error <= 1e-08), 'E_external_sound': bool(max(external_z) <= 2)}})
        print(item['tid'], 'mineral/water inventory checked', flush=True)
    gates = {k: all((r['gates'][k] for r in rows)) for k in ['I_inventory', 'A_attenuation', 'P_projection', 'E_external_sound']}
    gates['F_faults'] = all((v for f in faults for (k, v) in f.items() if k.endswith('_rejected')))
    result = {'round': 'R4', 'claim_type': 'capability', 'prereg_sha256': sha(ROOT / 'PREREG_R4_MINERAL_INVENTORY.json'), 'gates': gates, 'rows': rows, 'fault_injections': faults, 'array_manifest': manifest, 'external_referent': {'kind': 'closed_form', 'locator': 'https://physics.nist.gov/PhysRefData/Xcom/Text/chap3.html', 'compared_quantity': 'NIST mass-weighted compound/mixture attenuation; independentmass-fraction check against phase-density sum', 'refutes_us': not gates['A_attenuation']}, 'empirical_realism': 'UNKNOWN', 'sound_measurement_reference': 'https://pubmed.ncbi.nlm.nih.gov/8819355/', 'cost': {'wall_s': time.perf_counter() - tic, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'GPU_s': 0}, 'scope': 'Mineral/phase closure internally conservative; physicaldensities, composition andtissueboundaries remain scenarios. R1detector scores remain R1closure scores; notsilently updated toR4.'}
    write(ROOT / 'raw/R4_results.json', result)
    assert gates['F_faults']
    state('R4_COMPLETE', str(gates), 'integrate conservative material port and prospective lab predictions into demo')
    return result
if __name__ == '__main__':
    run()
