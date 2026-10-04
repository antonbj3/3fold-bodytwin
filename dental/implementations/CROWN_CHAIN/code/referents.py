from common_chain import *
import numpy as np
from scipy.stats import weibull_min
import time

def run():
    verify_freeze(ROOT / 'PREREG_REFERENTS.json')
    started = time.perf_counter()
    (co, matched) = scoped_modules(BASE / 'LANE_X1B_CROWN_LOOP/code', ['common', 'matched_contrasts'])
    groups = matched.derived_groups()
    rows = []
    for g in groups:
        if g['study'] != 'PMC8558575':
            continue
        (F0, m) = (g['characteristic_load_N'], g['Weibull_m'])
        q = F0 * (-np.log1p(-0.05)) ** (1 / m)
        independent = weibull_min.ppf(0.05, m, scale=F0)
        roundtrip = weibull_min.cdf(q, m, scale=F0)
        (flo, fhi) = g['F0_90CI_N']
        (mlo, mhi) = g['m_90CI']
        band = [flo * (-np.log1p(-0.05)) ** (1 / mlo), fhi * (-np.log1p(-0.05)) ** (1 / mhi)]
        wrong = F0 * (-np.log1p(-0.95)) ** (1 / m)
        assert abs(roundtrip - 0.05) < 1e-12 and abs(q - independent) < 1e-08
        rows.append(dict(**g, resolution='POPULATION', q05_N=q, parameter_rectangle_q05_N=band, rectangle_scope='Monotone image of two marginal90% CIs; simultaneous parameter coverage at least80% by union bound if each reported CI is valid; not prediction interval for our crowns', CDF_error=abs(roundtrip - 0.05), wrong_q95_rejected=abs(weibull_min.cdf(wrong, m, scale=F0) - 0.05) > 1e-12, patient_transfer='UNKNOWN_UNMATCHED_CROWN_MATERIAL_CONTACT_DIE_FILM', physical_patient_prediction_N=None))
    port = module('x23_port', BASE / 'LANE_X23_CROWN_WEIBULL_FLOOR/code/port_checks.py')
    x23 = read(BASE / 'LANE_X23_CROWN_WEIBULL_FLOOR/raw/R1_RESULTS.json')['per_group'][0]
    checks = [dict(name='X23_correct_lower_tail', good=port.quantile_ok(x23['q5_MPa'], x23['m'], x23['eta_MPa'], 0.05), fault_rejected=not port.quantile_ok(x23['q95_MPa'], x23['m'], x23['eta_MPa'], 0.05)), dict(name='X23_material_identity', good=port.material_transfer_ok(port.CUSTOM, port.CUSTOM), fault_rejected=not port.material_transfer_ok(port.CUSTOM, '3Y')), dict(name='X23_quantile_is_not_hard_floor', good=port.quantile_ok(x23['q5_MPa'], x23['m'], x23['eta_MPa'], 0.05), fault_rejected=not port.quantile_ok(x23['q5_MPa'], x23['m'], x23['eta_MPa'], 0.05, hard_floor=True))]
    cement = module('chain_cement', BASE / 'LANE_X13_CEMENT_GAP/cement_port.py')
    measurements = read(BASE / 'LANE_X13_CEMENT_GAP/measurements.json')
    cr = [cement.hydraulic_input_from_group_mean(r) for r in measurements if r['region'] in ('marginal', 'axial', 'occlusal')]
    alarm = module('chain_lab_alarm', BASE / 'LANE_X61_LAB_ALARM/code/lab_alarm.py')
    blank = ROOT / 'raw/EMPTY_MEASUREMENTS.csv'
    blank.write_text('specimen_id,design_id,load_angle_deg,fracture_force_N,MG_distance_mean_um\n')
    try:
        alarm.read_csv(blank)
        blank_rejected = False
    except ValueError:
        blank_rejected = True
    import csv
    template = next(csv.DictReader((BASE / 'LANE_X61_LAB_ALARM/inputs/LAB_MEASUREMENTS_TEMPLATE.csv').open()))
    template.update(fracture_force_N='1000', MG_distance_mean_um='50')
    valid = ROOT / 'raw/AUDIT_SYNTHETIC_VALID.csv'

    def fixture(path, row):
        with path.open('w') as f:
            writer = csv.DictWriter(f, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
    fixture(valid, template)
    valid_accepted = len(alarm.read_csv(valid)) == 1
    invalid = ROOT / 'raw/AUDIT_SYNTHETIC_ZERO_FORCE.csv'
    fixture(invalid, dict(template, fracture_force_N='0'))
    try:
        alarm.read_csv(invalid)
        invalid_rejected = False
    except ValueError:
        invalid_rejected = True
    checks.append(dict(name='X61_blank_protocol_is_not_measurement', good=valid_accepted, fault_rejected=blank_rejected and invalid_rejected))
    x82 = module('chain_x82_height', BASE / 'LANE_X82_HEIGHT_TO_FORCE/code/height_force.py') if False else None
    sys.path.insert(0, str(BASE / 'LANE_X82_HEIGHT_TO_FORCE/code'))
    x82 = module('chain_x82_height', BASE / 'LANE_X82_HEIGHT_TO_FORCE/code/height_force.py')
    cal = read(BASE / 'LANE_X82_HEIGHT_TO_FORCE/inputs/F5367_R4_state.json')
    edit = x82.make_edit(cal, cal['predicted_fdi'].index(16), 0)
    good = x82.predict(cal, edit)
    wrong = dict(edit, case=read(ROOT / 'FROZEN_INPUTS.json')['cohort'][0]['case_key'])
    try:
        x82.predict(cal, wrong)
        reject = False
    except ValueError:
        reject = True
    checks.append(dict(name='X82_real_function_rejects_other_patient', good=good.get('force_interval_N') is not None, fault_rejected=reject))
    bad = dict(edit, unit='um')
    try:
        x82.predict(cal, bad)
        reject = False
    except ValueError:
        reject = True
    checks.append(dict(name='X82_height_unit', good=True, fault_rejected=reject))
    assert all((r['good'] and r['fault_rejected'] for r in checks))
    out = dict(claim_type='capability', published_fracture_rows=rows, cement_group_mean_ports=cr, controls=checks, source_geometry_force_calibration=dict(case=cal['case'], source_geometry_sha256=cal['geometry_sha256'], status=good['status'], scope='Inherited simulated instrument ONLY; not any crown in this cohort'), population_to_patient_admitted=0, cement_to_seated_spatial_film_admitted=0, physical_lab_observations=0, cost=dict(seconds=time.perf_counter() - started, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    dump(ROOT / 'raw/REFERENTS.json', out)
    return out
if __name__ == '__main__':
    run()
