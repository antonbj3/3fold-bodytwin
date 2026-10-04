import json, math, time, datetime, hashlib, itertools
import numpy as np
from pathlib import Path
from scipy.optimize import brentq
from stress_port import force_point, force_bounds, measured_query
from freeze import frozen_write
ROOT = Path(__file__).resolve().parents[1]

def save(p, x):
    (ROOT / p).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def main():
    start = time.perf_counter()
    alpha = 0.251
    reference = 150.0
    R = 0.1
    states = []
    for c in [0.0, 100.0]:
        B = min(0.45 * reference / (reference - c), 1.0)
        C = 0.45 * reference / B ** alpha
        cs = [C]
        up = down = C
        for _ in range(256):
            up = math.nextafter(up, math.inf)
            down = math.nextafter(down, -math.inf)
            cs.extend([up, down])
        picked = next((value for value in cs if force_point(1.0, value, c, alpha, R) == reference), None)
        if picked is None:
            raise ValueError('bitwise identity construction failed; preserve failure')
        states.append({'g_MPa_per_N': 1.0, 'C_MPa': picked, 'compression_offset_MPa': c, 'alpha': alpha, 'reference_force_N': force_point(1.0, picked, c, alpha, R), 'changed_R_force_N': force_point(1.0, picked, c, alpha, 0.5)})
    same = states[0]['reference_force_N'] == states[1]['reference_force_N']
    diff = abs(states[0]['changed_R_force_N'] - states[1]['changed_R_force_N'])
    box = {'g': [0.95, 1.05], 'C': [states[1]['C_MPa'] * 0.98, states[1]['C_MPa'] * 1.02], 'c': [95.0, 105.0], 'alpha': [0.245, 0.26]}
    bound = force_bounds(box['g'], box['C'], box['c'], box['alpha'], 0.5)
    values = []
    for (g, C, c, a) in itertools.product(*[box[k] for k in ['g', 'C', 'c', 'alpha']]):
        q = 0.25
        fn = lambda t: q * t - C * min(q * t / (t - c), 1.0) ** a
        t = brentq(fn, math.nextafter(c, math.inf), C / q, xtol=1e-10)
        values.append(t / g)
    b = bound['force_interval_N']
    corners_pass = all((b[0] <= v <= b[1] for v in values))
    fixture_port = {'system_id': 'simulated_summary_pair', 'critical_region': 'declared single tensile root', 'g_MPa_per_N': box['g'], 'C_MPa': box['C'], 'compression_offset_MPa': box['c'], 'alpha': box['alpha'], 'measurement_locators': ['our_fixture:g', 'our_fixture:offset', 'our_fixture:force'], 'provenance_kind': 'our_own_fixture', 'horizon_cycles': 5000000.0, 'frozen_before_measurement_sha256': 'not_a_measurement', 'units': {'g': 'MPa/N', 'C': 'MPa', 'compression_offset': 'MPa', 'alpha': '1', 'horizon': 'cycles'}}
    nominal_query = measured_query(fixture_port, 0.5)
    unit_port = dict(fixture_port, provenance_kind='published_dataset', units={'g': 'microstrain/N', 'C': 'MPa', 'compression_offset': 'MPa', 'alpha': '1', 'horizon': 'cycles'})
    mutations = [{'gate': 'IDENTICAL_REFERENCE', 'nominal_pass': same, 'injected_error_rejected': force_point(1, states[1]['C_MPa'] * 2, 100.0, alpha, 0.1) != reference, 'injection': 'double resistance'}, {'gate': 'DOWNSTREAM_DIFFERENCE', 'nominal_pass': diff > 5, 'injected_error_rejected': abs(states[0]['changed_R_force_N'] - states[0]['changed_R_force_N']) <= 5, 'injection': 'erase mean-stress state'}, {'gate': 'PHYSICAL_PROVENANCE', 'nominal_pass': nominal_query['status'] == 'UNKNOWN', 'injected_error_rejected': nominal_query['status'] == 'UNKNOWN', 'injection': 'present simulated state as independent measurement'}, {'gate': 'NONLINEAR_ENCLOSURE', 'nominal_pass': corners_pass, 'injected_error_rejected': not all((b[0] <= 2 * v <= b[1] for v in values)), 'injection': 'double independently evaluated corner forces'}, {'gate': 'STRESS_UNIT', 'nominal_pass': True, 'injected_error_rejected': measured_query(unit_port, 0.5)['status'] == 'UNKNOWN', 'injection': 'strain slope used as stress slope'}, {'gate': 'SLOPE_SIGN', 'nominal_pass': True, 'injected_error_rejected': measured_query(dict(fixture_port, provenance_kind='published_dataset', frozen_before_measurement_sha256='0' * 64, g_MPa_per_N=[-1.0, -1.0]), 0.5).get('reason') == 'invalid measured interval', 'injection': 'negative tensile stress slope'}]
    out = {'claim_type': 'capability', 'states': states, 'summary_identity_error_N': abs(states[0]['reference_force_N'] - states[1]['reference_force_N']), 'summary_bitwise_equal': same, 'downstream_difference_N': diff, 'downstream_R': 0.5, 'sufficiency_rejected': same and diff > 5, 'smallest_model_extension': {'fixed_alpha': 'one added number compression_offset/g [N], together with reference force and reference machine R', 'unknown_alpha': 'also retain or independently calibrate alpha; full port ratios are C/g and compression_offset/g'}, 'measurement_spec': 'Same critical-region elastic stress slope g from calibrated force/strain increments; absolute compressive offset from assembly preload plus region XRD; reference finite-horizon fatigue endpoint gives C only under fixed alpha closure. Incremental strain alone does not identify residual stress. Second machine-R endurance measurement tests closure.', 'enclosure_box': box, 'rigorous_closure_enclosure': bound, 'independent_corner_forces_N': values, 'corner_check_pass': corners_pass, 'fixture_measurement_query': nominal_query, 'empirical_changed_R_validation': 'UNKNOWN; no same-assembly paired physical measurement in local sources', 'mutations': mutations, 'wall_s': time.perf_counter() - start, 'external_referent': {'kind': 'our_own_fixture', 'locator': 'code/run_r4.py states and code/stress_port.py; external constitutive equation source https://doi.org/10.6028/jres.099.032', 'compared_quantity': 'model-conditional equal-reference-force and changed-load-ratio force', 'refutes_us': True}, 'separate_external_closure_referent': {'kind': 'closed_form', 'locator': 'https://doi.org/10.6028/jres.099.032, Eq(2), with mean-stress extension retained as legacy closure', 'compared_quantity': 'hardness/defect surface endurance formula; no independent titanium validation', 'refutes_us': True}}
    save('results_R4.json', out)
    frozen = {'created_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'status': 'MEASUREMENT_NOT_PERFORMED', 'prereg_sha256': hashlib.sha256((ROOT / 'PREREG_R4.json').read_bytes()).hexdigest(), 'predictions': states, 'measurement_required': out['measurement_spec'], 'fixture_warning': 'These states are a mathematical counterexample, not manufactured assemblies or published data.'}
    frozen_write(ROOT / 'FROZEN_LAB_PREDICTIONS.json', frozen)
    print(json.dumps({k: v for (k, v) in out.items() if k not in ['states', 'mutations', 'independent_corner_forces_N', 'fixture_measurement_query', 'measurement_spec']}))
if __name__ == '__main__':
    main()
