from pathlib import Path
import json, hashlib, datetime, itertools, math
from fractions import Fraction
import numpy as np
from motion_port import stiffness, identify, project_intervals, maximum_abs_interval, local_threshold
LANE = Path(__file__).resolve().parents[1]

def main():
    cases = {'spread_support': {'z': [-4.0, 4.0], 'k': [32.0, 32.0]}, 'central_support': {'z': [-1.0, 1.0], 'k': [32.0, 32.0]}}
    marker_z = np.array([-4.0, 4.0])
    W = np.array([[1.0, 0.0], [0.0, 4.0]])
    query = np.array([1.0, 4.0])
    regional_z = np.array([-4.0, 0.0, 4.0])
    frozen = {'round': 'R2', 'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'claim_type': 'capability', 'status': 'PREDICTION_BEFORE_VIRTUAL_MEASUREMENT_ONLY_NO_REAL_NEW_MEASUREMENT', 'fixture_cases': cases, 'wrenches_N_Nmm_columns': W.tolist(), 'held_query': query.tolist(), 'predictions_um': {'spread_support': [-0.0, 15.625, 31.25], 'central_support': [-234.375, 15.625, 265.625]}, 'external_primary_source': 'PMC10792929#Tab1, already measured/published; cannot prospectively freeze that historical measurement'}
    payload = json.dumps(frozen, sort_keys=True, indent=2) + '\n'
    frozen_file = LANE / 'FROZEN_PREDICTIONS_R2.json'
    if frozen_file.exists():
        old = json.loads(frozen_file.read_text())
        assert old['predictions_um'] == frozen['predictions_um'] and old['fixture_cases'] == cases
    else:
        frozen_file.write_text(payload)
    (LANE / 'FROZEN_PREDICTIONS_R2.sha256').write_text(hashlib.sha256(frozen_file.read_bytes()).hexdigest() + '  ' + frozen_file.name + '\n')
    result_cases = {}
    summary = []
    for (name, d) in cases.items():
        K = stiffness(d['z'], d['k'])
        C = np.linalg.inv(K)
        q = C @ query
        H = np.column_stack([np.ones(2), marker_z])
        bone = np.array([[0.125, 0.0625], [0.125, 0.0625]])
        Y = H @ C @ W + bone
        ident = identify(marker_z, W, Y, bone)
        cnew = np.array(ident['C'])
        err = np.linalg.norm(cnew - C) / np.linalg.norm(C)
        marker_noise_mm = 1 / 1024
        qerr = np.array([marker_noise_mm, marker_noise_mm / 4])
        intervals = project_intervals(q - qerr, q + qerr, regional_z)
        peak_interval = (1000 * np.array(maximum_abs_interval(intervals))).tolist()
        pred = 1000 * (np.column_stack([np.ones(3), regional_z]) @ q)
        analytic = frozen['predictions_um'][name]
        assert np.array_equal(pred, np.array(analytic))
        k0 = sum((Fraction(x) for x in d['k']))
        k2 = sum((Fraction(k) * Fraction(z) ** 2 for (k, z) in zip(d['k'], d['z'])))
        rational = [Fraction(1000) * (Fraction(1) / k0 + Fraction(z) * Fraction(4) / k2) for z in regional_z]
        rational_error = max((abs(float(a) - float(b)) for (a, b) in zip(pred, rational)))
        verts = []
        for signs in itertools.product([-1, 1], repeat=2):
            qcorner = q + qerr * np.array(signs)
            vals = np.column_stack([np.ones(3), regional_z]) @ qcorner
            verts.append(vals)
            assert np.all(vals >= intervals[:, 0]) and np.all(vals <= intervals[:, 1])
        bad = Y.copy()
        bad[0, 1] += 1
        neg = identify(marker_z, W, bad, bone)
        summary.append([float(K[0, 0]), float(K[0, 1]), float(C[0, 0])])
        result_cases[name] = {'K': K.tolist(), 'C': C.tolist(), 'scalar_summary': summary[-1], 'scalar_summary_units': ['N/mm', 'N', 'mm/N'], 'K_matrix_units': [['N/mm', 'N/rad'], ['N', 'Nmm/rad']], 'C_matrix_units': [['mm/N', 'mm/Nmm'], ['rad/N', 'rad/Nmm']], 'regional_motion_um': pred.tolist(), 'regional_abs_max_um': float(max(abs(pred))), 'identified': ident, 'relative_recovery_error': float(err), 'rational_solution_identity_error_um': rational_error, 'sampled_peak_interval_um': peak_interval, 'interval_claim': 'outward-enclosed independent q-box for known C; not dental FE enclosure, not uncertain fitted-C enclosure', 'threshold': local_threshold(peak_interval, True), 'biological_outcome': local_threshold(peak_interval, True, True), 'injected_1mm_calibration': neg, 'missing_bone_reference': identify(marker_z, W, Y), 'resolution': 'PER_SURFACE_REGION', 'stiffness_status': 'PHENOMENOLOGICAL; replace with actual independent local bone-reference load/displacement observations', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'code/round2.py dyadic inputs + Fraction exact equilibrium', 'compared_quantity': 'regional displacement and reconstructed matrix', 'refutes_us': False}}
    identity = float(np.max(np.abs(np.array(summary[0]) - np.array(summary[1]))))
    assert identity == 0
    data = json.loads((LANE / 'raw/measurements.json').read_text())
    pick = lambda c: next((d for d in data if d['study'] == 'PMC10792929' and d['condition'] == c))
    (axial_a, axial_b) = (pick('B11_Z'), pick('MT_Z'))
    (total_a, total_b) = (pick('B11_total'), pick('MT_total'))
    witness = {'axial_mean_um': [axial_a['motion_um'], axial_b['motion_um']], 'reported_mean_identity_error_um': abs(axial_a['motion_um'] - axial_b['motion_um']), 'total_mean_um': [total_a['motion_um'], total_b['motion_um']], 'total_mean_difference_um': total_b['motion_um'] - total_a['motion_um'], 'total_mean_ratio': total_b['motion_um'] / total_a['motion_um'], 'units': 'um', 'resolution': 'POPULATION; equal published rounded means, hidden specimen states not identical', 'minimum_extension': 'retain transverse displacement/rotation under specified wrench, and reference both to local bone; axial mean alone loses total motion', 'causal_limit': 'different shapes and resin fixation and spoon moment confound mechanism; no causal geometry coefficient established', 'external_referent': {'kind': 'independent_measurement', 'locator': 'doi:10.1186/s12903-024-03854-1; PMC10792929#Tab1 B11,MT rows (Z and Total)', 'compared_quantity': 'same reported axial displacement mean versus total displacement mean at 50 N', 'refutes_us': True}}
    out = {'round': 'R2', 'claim_type': 'capability', 'verdict': 'SCALAR_SUMMARY_INSUFFICIENT; LINEAR_MEASUREMENT_PORT_DEMONSTRATED_UNDER_ASSUMPTIONS', 'summary_identity_error': identity, 'downstream_peak_difference_um': result_cases['central_support']['regional_abs_max_um'] - result_cases['spread_support']['regional_abs_max_um'], 'downstream_peak_ratio': result_cases['central_support']['regional_abs_max_um'] / result_cases['spread_support']['regional_abs_max_um'], 'minimum_extension': 'K2 for symmetric K1=0 with known K0; general asymmetric support needs K1 and K2 (2x2 compliance) plus actual local bone-reference observation', 'cases': result_cases, 'published_sufficiency_witness': witness, 'external_referent': witness['external_referent'], 'strongest_equal_information_check': 'direct complete spring-matrix solve agrees exactly; no algorithm superiority claimed', 'failed_or_unknown': ['No real local bone-reference matrix measured; source group means only', 'Cannot convert source optical displacement to local interface maximum', 'No rigorous dental FE sensitivity enclosure', 'identification under uncertain measurement intervals needs nonlinear inverse-box enclosure; q projection enclosure only'], 'full_cost': {'preparation': '4 local bone-referenced marker readings (2 markers x 2 independent loads); bone coordinates measured separately', 'fit': '2x2 inversion, 3 reciprocal coefficients', 'validation': 'third held wrench in linear fixture; symmetry/PD and bad-reading rejection', 'fallback': 'missing local reference/rank returns UNKNOWN'}, 'time_scale': 'SIMULTANEOUS for mechanics; later healing HANDOVER not validated'}
    (LANE / 'rounds/R2').mkdir(parents=True, exist_ok=True)
    (LANE / 'rounds/R2/results.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: out[k] for k in ['verdict', 'summary_identity_error', 'downstream_peak_difference_um', 'downstream_peak_ratio', 'published_sufficiency_witness']}, indent=2))
if __name__ == '__main__':
    main()
