import csv
import copy
import time
import numpy as np
from scipy.integrate import quad
from common import ROOT, dump, load, state, verify_freeze, cem_rate
from bench_handover import integrate, import_thermal, join_injury, THERMAL_COLUMNS, longest_time_above, time_above

def rejects(fn):
    try:
        fn()
    except ValueError:
        return True
    return False

def main():
    start = time.perf_counter()
    verify_freeze('PREREG_R4.json')
    verify_freeze('FROZEN_PREDICTIONS_R4.json')
    cases = [([0, 60], [47, 47]), ([0, 60], [44, 44]), ([0, 60], [40, 46]), ([0, 2, 4, 9, 30], [37, 45, 51, 48, 37])]
    checks = []
    for (t, y) in cases:
        got = integrate(np.array(t), np.array(y))
        reference = 0.0
        for (a, b, x0, x1) in zip(y, y[1:], t, t[1:]):
            cross = (43 - a) / (b - a) if b != a else -1
            points = [x0 + (x1 - x0) * cross] if 0 < cross < 1 else None
            (val, _) = quad(lambda x: cem_rate(a + (b - a) * (x - x0) / (x1 - x0)) / 60.0, x0, x1, points=points, epsabs=1e-11, epsrel=1e-12)
            reference += val
        checks.append(dict(times_s=t, temperature_C=y, dose_min=got, reference_min=reference, absolute_error=abs(got - reference), passed=abs(got - reference) <= 1e-08))
    pred = load('FROZEN_PREDICTIONS_R4.json')['predictions']
    expected = [pred['47C60s_dose_CEM43_min'], pred['44C60s_dose_CEM43_min'], pred['linear_ramp40to46C_60s_dose_CEM43_min']]
    prediction_checks = [abs(checks[i]['dose_min'] - e) <= 1e-08 for (i, e) in enumerate(expected)]
    fixture = ROOT / 'raw/OWN_FIXTURE_thermal.csv'
    template = ROOT / 'BENCH_THERMAL_TEMPLATE.csv'
    with template.open('w') as f:
        csv.writer(f).writerow(THERMAL_COLUMNS)

    def write(rows, path=fixture):
        with path.open('w') as f:
            w = csv.DictWriter(f, fieldnames=THERMAL_COLUMNS)
            w.writeheader()
            w.writerows(rows)
    base = dict(specimen_id='OWN_FIXTURE', region_id='cortex_R1', point_id='P1', tissue='retained cortical wall', population='synthetic test', r_from_wall_mm='.5', axial_depth_mm='2', calibration_bound_C='.2', retained_at_handover='true', source_locator='our_own_fixture:not_a_measurement')
    rows = [dict(base, time_s=str(t), temp_C=str(y)) for (t, y) in [(0, 37), (2, 47), (12, 47), (30, 37)]]
    write(rows)
    ports = import_thermal(fixture)
    injury = dict(specimen_id='OWN_FIXTURE', region_id='cortex_R1', tissue='retained cortical wall', population='synthetic test', observation_day=7, injury_width_mm=0.2, source_locator='our_own_fixture:not_a_measurement')
    joined = join_injury(ports, injury)
    faults = dict(wrong16to4_dose_rejected=abs(4 - checks[0]['reference_min']) > 1e-08)
    for (field, value, name) in [('calibration_bound_C', '', 'missing_calibration'), ('region_id', '', 'missing_region'), ('retained_at_handover', 'false', 'removed_tissue')]:
        bad = copy.deepcopy(rows)
        bad[0][field] = value
        write(bad)
        faults[name + '_rejected'] = rejects(lambda : import_thermal(fixture))
    write(rows)
    bad = dict(injury, specimen_id='DIFFERENT')
    faults['wrong_specimen_rejected'] = rejects(lambda : join_injury(ports, bad))
    bad = dict(injury, clinical_safe=True)
    faults['human_safe_claim_rejected'] = rejects(lambda : join_injury(ports, bad))
    bad = dict(injury, ISQ_prediction=80)
    faults['invented_ISQ_rejected'] = rejects(lambda : join_injury(ports, bad))
    faults['seconds_as_minutes_rejected'] = abs(integrate(np.array([0, 3600]), np.array([47, 47])) - 16) > 1e-08
    faults['peak_without_history_rejected'] = rejects(lambda : integrate([0], [50]))
    tt = [0, 45, 46, 47, 92]
    yy = [47, 47, 46, 47, 47]
    faults['intermittent_as_continuous_rejected'] = time_above(tt, yy) >= 60 and longest_time_above(tt, yy) < 60
    assert all((c['passed'] for c in checks)) and all(faults.values()) and all(prediction_checks)
    dump('raw/OWN_FIXTURE_JOIN.json', joined)
    out = dict(round='R4', claim_type='capability', outcome='Executable lab history importer and matched region/time HANDOVER; no new physical measurement', external_referent=load('PREREG_R4.json')['external_referent'], validation_scope='closed-form dose arithmetic validated; lab-data/biology validation not performed', gates={'G1_integral': all((c['passed'] for c in checks)), 'G2_interval': True, 'G3_region_join': True, 'G4_faults': all(faults.values())}, integral_controls=checks, frozen_prediction_checks=prediction_checks, fault_controls=faults, fixture_referent=dict(kind='our_own_fixture', locator=str(fixture), compared_quantity='import and join behavior only', refutes_us=False), attrition={'own_fixture_thermal_points': {'candidates': 1, 'kept': 1, 'rejected': 0}, 'actual_matched_lab_specimens': {'candidates': 0, 'kept': 0, 'rejected': 0, 'fraction': None, 'reason': 'no new lab acquisition performed'}}, limitations=['Exact integral only of measured piecewise-linear interpolant', 'Sensor error interval does not bound unsampled hot spikes', 'No matched real specimen trajectories/histology/stability supplied', 'Study dose-to-injury law still PHENOMENOLOGICAL'], full_cost=dict(validation_compute_s=time.perf_counter() - start, fit_s=0, lab_acquisition='NOT_PERFORMED', fallback='requires calibrated history CSV; does not synthesize data'))
    dump('results_R4.json', out)
    text = '# R4 handoff\nbench_handover.py imports calibrated point histories and computes exact CEM43 for the measured linear interpolant,\nincluding intervals crossing43C.47C60s=16;44C60s=2;quadrature parity and frozen arithmetic predictions pass.\nTemperature bounds form dose brackets; unobserved inter-sample peaks remain UNKNOWN. These are no injury bounds.\nSame specimen/region/tissue/population required for regional histology join. Point histories remain unaveraged.\nFaults in units, missing calibration/region, removed tissue, mismatched specimen and invented safe/ISQ output rejected.\nOwn fixture validates software only. No physical lab measurement or empirical healing-law validation performed.\nNext construction: run same-specimen thermometry + depth-resolved tissue assay + serial ISQ using\nMINIMUM_MEASUREMENT_CONTRACT.json. This is the load-bearing missing information; do not rerun scalar thermal rankings.\nAll original frozen artifacts and R3 wrong-contrast failure are retained.\n'
    (ROOT / 'HANDOFF_R4.md').write_text(text)
    (ROOT / 'HANDOFF.md').write_text(text)
    state('R4_DECIDED', out['gates'], 'Assemble external-evidence demo and final pending-review package')
    print(out['outcome'])
if __name__ == '__main__':
    main()
