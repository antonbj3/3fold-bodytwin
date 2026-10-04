from patient_geometry import P, write, sha
from load_port import bind, fixed_load_decision
import numpy as np, json, copy, time, resource, datetime
from scipy.optimize import linprog

def main():
    tick = time.perf_counter()
    g = json.load(open(P / 'rounds/R1/results.json'))
    external = json.load(open(P / 'inputs/hattori1994_closed_bite.json'))
    fdi = external['fdi']
    ix = fdi.index(36)
    rows = []
    for r in external['rows']:
        f = r['reported_force_N'][ix]
        interval = [f - 0.5, f + 0.5]
        rows.append(dict(trial=r['trial'], fdi=36, reported_force_N=f, display_rounding_interval_N=interval, display_force_load_decision=fixed_load_decision(interval), total_printed_N=r['reported_total_N'], sum_printed_N=sum(r['reported_force_N']), balance_error_N=sum(r['reported_force_N']) - r['reported_total_N'], resolution='PER_TOOTH', measurement_error='UNKNOWN; displayed rounding interval not sensor error bound'))
    m = json.load(open(P / 'raw/SOURCE_MANIFEST.json'))
    identity = sha(P / 'raw/SOURCE_MANIFEST.json')
    all_fdi = sorted((v['fdi'] for v in m['stls']))
    blank = dict(patient_id='LiuHao2023Demo1', geometry_manifest_sha256=identity, measurement_kind='independent_measurement', acquisition_id=None, source_fdi=all_fdi, calibration_locator=None, independent_check_locator=None, loaded_pose_locator=None, directed_registration_error_locator=None, loaded_upper_to_cbct_4x4=None, directed_registration_error_mm=None, sensor_thickness_mm=None, sensor_thickness_error_mm=None, force_unit='N', resolution='PER_POINT', timescale='SIMULTANEOUS', sensor_state='declared_inserted_sensor_state', force_kind='calibrated_normal_contact_force', total_force_interval_N=None, contacts=None)
    write(P / 'exports/DEMO1_MEASUREMENT_TEMPLATE.json', blank)
    demo = bind(blank, blank['patient_id'], identity, expected_source_fdi=all_fdi)
    fx = copy.deepcopy(blank)
    fx.update(measurement_kind='our_own_fixture', acquisition_id='FIXTURE_NOT_PHYSICALLY_MEASURED', calibration_locator='our_own_fixture', independent_check_locator='our_own_fixture', loaded_pose_locator='our_own_fixture', directed_registration_error_locator='our_own_fixture', loaded_upper_to_cbct_4x4=np.eye(4).tolist(), directed_registration_error_mm=0.01, sensor_thickness_mm=0.1, sensor_thickness_error_mm=0.01, contacts=[dict(upper_fdi=r['upper_fdi'], lower_fdi=r['lower_fdi'], force_interval_N=[0.9, 1.1], point_cbct_mm=r['lower_point_cbct_mm'], source_face_locator='our_own_fixture:' + str(r['lower_source_face']), force_error_locator='our_own_fixture') for r in g['pairs']])
    covered = {t for c in fx['contacts'] for t in [c['upper_fdi'], c['lower_fdi']]}
    fx.update(verified_zero_fdi=sorted(set(all_fdi) - covered), zero_detection_bound_locator='our_own_fixture', total_force_interval_N=[26.0, 28.0])
    valid = bind(fx, fx['patient_id'], identity, allow_fixture=True, expected_source_fdi=all_fdi)
    faults = []

    def run(name, mut, allow=True):
        q = copy.deepcopy(fx)
        mut(q)
        try:
            a = bind(q, fx['patient_id'], identity, allow_fixture=allow, expected_source_fdi=all_fdi)
            rejected = a['status'] == 'UNKNOWN'
            reason = a.get('reason')
        except (ValueError, TypeError) as e:
            rejected = True
            reason = str(e)
        faults.append(dict(name=name, rejected=rejected, reason=reason))
    run('wrong_patient', lambda q: q.update(patient_id='AnotherPatient'))
    run('wrong_geometry', lambda q: q.update(geometry_manifest_sha256='00' * 32))
    run('force_unit_wrong', lambda q: q.update(force_unit='kN'))
    run('non_point_force_resolution', lambda q: q.update(resolution='POPULATION'))
    run('handover_mixed_with_simultaneous', lambda q: q.update(timescale='HANDOVER'))
    run('missing_calibration', lambda q: q.update(calibration_locator=None))
    run('missing_independent_check', lambda q: q.update(independent_check_locator=None))
    run('unknown_bite_pose', lambda q: q.update(loaded_upper_to_cbct_4x4=None))
    run('improper_pose', lambda q: q.update(loaded_upper_to_cbct_4x4=np.diag([-1, 1, 1, 1]).tolist()))
    run('negative_directed_error', lambda q: q.update(directed_registration_error_mm=-0.1))
    run('NaN_directed_error', lambda q: q.update(directed_registration_error_mm=float('nan')))
    run('missing_sensor_thickness', lambda q: q.update(sensor_thickness_mm=None))
    run('unknown_force_operator', lambda q: q.update(force_kind='relative_TScan_percentage'))
    run('unclear_sensor_state', lambda q: q.update(sensor_state='sensor_free_assumed'))
    run('injected_total100N', lambda q: q.update(total_force_interval_N=[100, 101]))
    run('reversed_force_interval', lambda q: q['contacts'][0].update(force_interval_N=[2, 1]))
    run('negative_force', lambda q: q['contacts'][0].update(force_interval_N=[-1, 1]))
    run('unknown_point_error', lambda q: q['contacts'][0].update(force_error_locator=None))
    run('absent_source_tooth', lambda q: q['contacts'][0].update(lower_fdi=38))
    run('missing_last_contact_and_zero_proof', lambda q: (q.update(contacts=q['contacts'][:-1]), q.update(verified_zero_fdi=[], zero_detection_bound_locator=None)))
    run('fixture_masquerading_as_patient', lambda q: None, allow=False)
    run('gain_unknown', lambda q: q.update(measurement_kind=None))
    try:
        bind(dict(patient_id=external['patient_id']), blank['patient_id'], identity, expected_source_fdi=all_fdi)
        mismatch = False
    except ValueError:
        mismatch = True
    p = [np.array(r['lower_point_cbct_mm']) for r in g['pairs'] if r['lower_fdi'] == 36]
    v1 = np.array([0.0, 0.0, 100.0])
    momdiff = np.cross(p[1] - p[0], v1)
    suff = dict(summary_total_N=[100.0, 100.0], summary_tooth36_N=[100.0, 100.0], identity_error_N=0.0, contact_point_A_cbct_mm=p[0].tolist(), contact_point_B_cbct_mm=p[1].tolist(), downstream_moment_difference_vector_Nmm=momdiff.tolist(), downstream_moment_difference_norm_Nmm=float(np.linalg.norm(momdiff)), minimal_extension='For resultant wrench: force-weighted contact centroid (two transverse coordinates) + force direction. For local crown stress: per-contact pressure field still needed.', scope='Own force-state identifiability fixture on externally released geometry; not independently measured force', resolution='PER_POINT')
    lows = np.array([c['force_interval_N'][0] for c in fx['contacts']])
    highs = np.array([c['force_interval_N'][1] for c in fx['contacts']])
    A = np.ones((1, len(lows)))
    controls = []
    for row in valid['forces']:
        c = np.array([row['fdi'] in [v['upper_fdi'], v['lower_fdi']] for v in fx['contacts']], float)
        lo = linprog(c, A_ub=np.r_[A, -A], b_ub=[28, -26], bounds=list(zip(lows, highs)), method='highs', options={'threads': 4})
        hi = linprog(-c, A_ub=np.r_[A, -A], b_ub=[28, -26], bounds=list(zip(lows, highs)), method='highs', options={'threads': 4})
        err = max(abs(row['force_interval_N'][0] - lo.fun), abs(row['force_interval_N'][1] + hi.fun))
        controls.append(dict(fdi=row['fdi'], error_N=float(err), pass_gate=err <= 1e-07, injected_endpoint1N_rejected=abs(row['force_interval_N'][1] + 1 + hi.fun) > 1e-07))
    out = dict(round='R3', claim_type='capability', patient_id=blank['patient_id'], actual_patient_port=demo, physical_K16='STILL_MISSING; no matched measured closed bite + force', external_reference_subject=external['patient_id'], external_rows=rows, external_force36_observed_range_N=[min((r['reported_force_N'] for r in rows)), max((r['reported_force_N'] for r in rows))], fixed100N_display_exceedances=sum((r['display_force_load_decision']['status'] == 'EXCEEDS' for r in rows)), published_force_transfer_rejected=mismatch, external_model_validation='Benchmark tests published force consumer and calibration only; no validation of Demo1 patient forces', controls=controls, faults=faults, all_control_gates_pass=all((c['pass_gate'] and c['injected_endpoint1N_rejected'] for c in controls)) and all((f['rejected'] for f in faults)), fixture_valid_result=valid, sufficiency=suff, external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.2186/jjps.38.835; Table2 p838 six clench maps', compared_quantity='Reported per-tooth static force and fixed100N load enclosure on reference subject; not Demo1', refutes_us=True), attrition=dict(reference_printed_force_rows=84, retained_for_descriptive_table=84, accepted_for_Demo1=0, rejected_for_Demo1=84, rejection_fraction=1.0, reason='Different subject, no Demo1 pose/identity or calibrated error; zero signals additionally censored'), cost=dict(wall_s=time.perf_counter() - tick, maxrss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, fit='NONE', upstream_discovery='UNKNOWN', source_transcription='UNKNOWN', physical_acquisition='NOT_RUN'))
    write(P / 'rounds/R3/results.json', out)
    write(P / 'CURRENT_WORK_STATE.json', dict(status='R3_COMPLETE', latest_gate='Executable static port; Demo1 UNKNOWN; external fixed100N refuted for5/6 reported force states', next_operation='Acquire Demo1 loaded bite and independently calibrated registered force map; separately perturb FDI36 height', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    print(json.dumps({k: out[k] for k in ['physical_K16', 'external_force36_observed_range_N', 'fixed100N_display_exceedances', 'published_force_transfer_rejected', 'all_control_gates_pass']}, indent=2))
    print('faults', len(faults), 'moment difference', suff['downstream_moment_difference_norm_Nmm'])
if __name__ == '__main__':
    main()
