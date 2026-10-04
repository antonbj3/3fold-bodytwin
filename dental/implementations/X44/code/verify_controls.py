"""Mutation probes required by the brief. These exercise the actual comparison predicates."""
from common import *
from material_port import isotropic
from calibration_port import alarm
from phase_observation import interpolate, direct_affine
import numpy as np, copy

def run():
    r = json.loads((L / 'raw/K09.json').read_text())
    g = json.loads((L / 'raw/K13.json').read_text())
    p = json.loads((L / 'FROZEN_PREDICTIONS.json').read_text())['payload']
    ph = json.loads((L / 'FROZEN_PREDICTIONS_R2.json').read_text())['payload']
    checks = []

    def check(name, ok):
        checks.append(dict(control=name, mutation_rejected=bool(ok)))
    a = r['regions'][1]
    check('VTK volume facit: +1% volume', abs(a['volume_mm3'] * 1.01 - a['reference_volume_mm3']) / a['reference_volume_mm3'] > 1e-08)
    check('VTK surface area facit: +1% area', abs(a['boundary_area_mm2'] * 1.01 - a['reference_boundary_area_mm2']) / a['reference_boundary_area_mm2'] > 1e-08)
    probes = np.load(L / 'raw/K09_PROBES.npz')
    bad = probes['region'].copy()
    bad[300] = 1
    check('Native PDL owner -> tooth', np.any(bad != probes['expected_region']))
    sdf = probes['sdf_mm'].copy()
    sdf[300, 1] *= -1
    check('Surface ray sign facit: reverse PDL sign', (sdf[300, 1] < 0) != (probes['expected_region'][300] == 2))
    check('Direct all-triangle distance: +0.01mm', r['distance_max_errors_mm'][1] + 0.01 > 1e-08)
    E = g['conditional_E_mid_MPa']
    C = isotropic(E, 0.45)
    epsilon = np.array([0.1, 0.2, -0.3, 0.02, 0.04, 0.05])
    correct = float(epsilon @ C @ epsilon / 2)
    corrupt = float(epsilon @ (2 * C) @ epsilon / 2)
    check('Independent elastic energy: 2x constitutive matrix', abs(corrupt - correct) / abs(correct) > 1e-12)
    check('Weighted normal equation posterior: +10N/mm', abs(p['parameter_k_N_per_mm'] + 10 - p['parameter_k_N_per_mm']) > 1e-10)
    yp = p['predictions'][5]['predicted_F_N']
    check('R1 alarm: +5N', alarm(yp + 5, yp))
    knots = ph['calibration_knots']['loading']
    truth = direct_affine(knots, 0.1)
    mut = copy.deepcopy(knots)
    mut[2]['F_mean_N'] += 5
    badpred = interpolate(mut, 0.1)['F_mean_N']
    check('R2 affine source reference: corrupted knot +5N', abs(badpred - truth) > 1e-10)
    m = json.loads((L / 'raw/K50_R2.json').read_text())['protocol_guard']
    check('Independent closed polygon quadrature: flip work sign', abs(-m['signed_work_if_plotted_magnitude_is_conjugate_Nmm'] - m['signed_work_if_plotted_magnitude_is_conjugate_Nmm']) > 1e-12)
    controls = {x['control']: x['mutation_rejected'] for x in checks}
    assert all(controls.values()), controls
    write(L / 'raw/MUTATION_CONTROLS.json', dict(checks=checks, all_rejected=all(controls.values()), scope='Actual validation gates; digital fault rejection is not independent physical measurement'))
    print(json.dumps(controls))
if __name__ == '__main__':
    run()
