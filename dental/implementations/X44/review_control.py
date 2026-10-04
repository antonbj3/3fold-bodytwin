"""Corrected X44 fault probes; independent raw calculations, existing observations."""
import sys, json, math, copy
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parent
sys.path.insert(0, str(R / 'code'))
from common import ROOT
from calibration_port import alarm
from phase_observation import observation_query
from material_port import isotropic

def main():
    rows = json.loads((ROOT / 'results/F2_pdl_nonlinear/jepsen2023_fig2b.json').read_text())['points']
    pr = json.loads((R / 'PREREG_K50_R1.json').read_text())
    cal = [rows[i] for i in pr['calibration_points_k']]
    expected = sum((r['deflection_mm'] * r['F_mean_N'] for r in cal)) / sum((r['deflection_mm'] ** 2 for r in cal))
    actual = json.loads((R / 'raw/K50_R1.json').read_text())['posterior']['k_mean_N_per_mm']

    def coefficient_ok(v):
        return math.isfinite(v) and abs(v - expected) <= 1e-10
    p = np.load(R / 'raw/K09_PROBES.npz')
    n = np.load(R / 'local_data/region_field.npz')
    import trimesh
    mesh = trimesh.Trimesh(n['V'], n['faces_2'], process=False)
    pts = p['points_mm'][300:303]
    truth = []
    for q in pts:
        truth.append(trimesh.proximity.closest_point_naive(mesh, q[None])[1][0])
    truth = np.array(truth)
    actual_d = abs(p['sdf_mm'][300:303, 1])

    def distances_ok(v):
        return bool(np.all(np.isfinite(v)) and np.max(abs(v - truth)) <= 1e-08)
    work = sum(((a['F_mean_N'] + b['F_mean_N']) * (b['deflection_mm'] - a['deflection_mm']) / 2 for (a, b) in zip(rows, rows[1:])))
    reported = json.loads((R / 'raw/K50_R2.json').read_text())['protocol_guard']['signed_work_if_plotted_magnitude_is_conjugate_Nmm']

    def work_ok(v):
        return math.isfinite(v) and abs(v - work) <= 1e-12
    model = json.loads((R / 'OBSERVATION_PORT.json').read_text())['model']
    a = observation_query(model, 0.1, 'loading')
    b = observation_query(model, 0.1, 'unloading')
    tests = [dict(control='normal-equation coefficient from original calibration points', baseline=coefficient_ok(actual), wrong_rejected=not coefficient_ok(actual + 10)), dict(control='direct all-triangle distance versus actual probe array', baseline=distances_ok(actual_d), wrong_rejected=not distances_ok(actual_d + 0.01)), dict(control='independent raw closed-polygon quadrature', baseline=work_ok(reported), wrong_rejected=not work_ok(-reported))]
    E = json.loads((R / 'raw/K13.json').read_text())['conditional_E_mid_MPa']
    nu = 0.45
    strain = np.array([0.1, 0.2, -0.3, 0.02, 0.04, 0.05])
    eps = np.array([[strain[0], strain[3] / 2, strain[5] / 2], [strain[3] / 2, strain[1], strain[4] / 2], [strain[5] / 2, strain[4] / 2, strain[2]]])
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    truth_energy = 0.5 * lam * np.trace(eps) ** 2 + mu * np.sum(eps * eps)
    C = isotropic(E, nu)

    def energy_ok(matrix):
        return bool(abs(float(strain @ matrix @ strain / 2) - truth_energy) <= 1e-12 * truth_energy)
    tests.append(dict(control='full-tensor elastic-energy reference', baseline=energy_ok(C), wrong_rejected=not energy_ok(2 * C)))

    def fault_probe(predicate):
        return predicate(actual) and (not predicate(actual + 10))
    always_pass_mutation_detected = not fault_probe(lambda submitted: True)
    out = dict(controls=tests, always_pass_mutation_detected=always_pass_mutation_detected, expected_coefficient_N_per_mm=expected, independent_work_Nmm=work, sufficiency=dict(summary='deflection_mm', states=['loading', 'unloading'], identical_summary=[0.1, 0.1], identity_error=0.0, downstream_force_N=[a['F_mean_N'], b['F_mean_N']], downstream_difference_N=b['F_mean_N'] - a['F_mean_N'], minimum_extension='phase at minimum for this observed two-branch curve; history/rate/sign/paired geometry still required for a physical material law'), conclusion_changed=False, physical_material_update='REFUSED')
    assert all((t['baseline'] and t['wrong_rejected'] for t in tests)) and a['F_mean_N'] != b['F_mean_N']
    (R / 'raw/REVIEW_CONTROL.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out))
if __name__ == '__main__':
    main()
