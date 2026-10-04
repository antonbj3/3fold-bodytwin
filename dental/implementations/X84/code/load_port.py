"""Measurement-to-static-load port. No physical force invented from geometry."""
import numpy as np

def finite_nonnegative(v, name):
    if v is None:
        raise ValueError('UNKNOWN:' + name)
    a = np.asarray(v, dtype=float)
    if not np.isfinite(a).all() or (a < 0).any():
        raise ValueError('Invalid ' + name)
    return a

def bind(packet, patient_id, source_sha256, allow_fixture=False, expected_source_fdi=None):
    """Strict sufficient contract. A source hash is identity, never measurement truth."""
    if packet.get('patient_id') != patient_id:
        raise ValueError('Patient identity mismatch')
    if packet.get('geometry_manifest_sha256') != source_sha256:
        raise ValueError('Geometry identity mismatch')
    fixture = packet.get('measurement_kind') == 'our_own_fixture'
    if expected_source_fdi is None:
        return dict(status='UNKNOWN', reason='Source inventory not supplied independently by consumer')
    if sorted(packet.get('source_fdi', [])) != sorted(expected_source_fdi):
        return dict(status='UNKNOWN', reason='Source inventory does not match bound manifest')
    if fixture and (not allow_fixture):
        return dict(status='UNKNOWN', reason='A fixture is not a patient measurement')
    if not fixture and packet.get('measurement_kind') != 'independent_measurement':
        return dict(status='UNKNOWN', reason='No independent measurement provenance')
    provenance_keys = ['acquisition_id', 'calibration_locator', 'independent_check_locator', 'loaded_pose_locator', 'directed_registration_error_locator']
    if not fixture and (any(('fixture' in str(packet.get(k, '')).lower() for k in provenance_keys)) or any(('fixture' in str(c.get(k, '')).lower() for c in packet.get('contacts') or [] for k in ['source_face_locator', 'force_error_locator']))):
        return dict(status='UNKNOWN', reason='Synthetic provenance cannot be relabelled independent')
    req = ['acquisition_id', 'calibration_locator', 'independent_check_locator', 'loaded_pose_locator', 'directed_registration_error_locator', 'force_unit', 'resolution', 'timescale', 'sensor_state', 'force_kind']
    for k in req:
        if not packet.get(k):
            return dict(status='UNKNOWN', reason='Missing ' + k)
    if packet['force_unit'] != 'N':
        raise ValueError('Force unit must be N')
    if packet['resolution'] != 'PER_POINT' or packet['timescale'] != 'SIMULTANEOUS':
        raise ValueError('Wrong force resolution or time scale')
    if packet['force_kind'] != 'calibrated_normal_contact_force':
        return dict(status='UNKNOWN', reason='No calibrated normal-force observation operator')
    if packet['sensor_state'] != 'declared_inserted_sensor_state':
        return dict(status='UNKNOWN', reason='Sensor-present versus sensor-free state unidentified')
    for k in ['directed_registration_error_mm', 'sensor_thickness_mm', 'sensor_thickness_error_mm']:
        try:
            a = finite_nonnegative(packet.get(k), k)
            if a.shape != ():
                raise ValueError('Expected scalar ' + k)
        except ValueError as e:
            return dict(status='UNKNOWN', reason=str(e))
    pose = packet.get('loaded_upper_to_cbct_4x4')
    if pose is None:
        return dict(status='UNKNOWN', reason='Loaded bite transform UNKNOWN')
    T = np.array(pose, float)
    if T.shape != (4, 4) or not np.isfinite(T).all() or (not np.allclose(T[3], [0, 0, 0, 1], atol=1e-10, rtol=0)) or (not np.allclose(T[:3, :3].T @ T[:3, :3], np.eye(3), atol=1e-10, rtol=0)) or (abs(np.linalg.det(T[:3, :3]) - 1) > 1e-10):
        raise ValueError('Loaded pose is not a proper rigid transform')
    contacts = packet.get('contacts')
    if not contacts:
        return dict(status='UNKNOWN', reason='No measured contacts')
    fdi_expected = sorted(packet.get('source_fdi', []))
    if len(fdi_expected) != len(set(fdi_expected)) or len(fdi_expected) == 0:
        raise ValueError('Invalid tooth inventory')
    if packet.get('verified_zero_fdi'):
        return dict(status='UNKNOWN', reason='Censored zeros need numerical contact-force bounds; locator alone is insufficient')
    covered = set()
    for c in contacts:
        try:
            if c.get('force_unit', 'N') != 'N':
                raise ValueError('Contact force unit contradicts packet unit')
            if c.get('resolution', 'PER_POINT') != 'PER_POINT':
                raise ValueError('Contact resolution contradicts packet resolution')
            if c['upper_fdi'] // 10 not in [1, 2] or c['lower_fdi'] // 10 not in [3, 4]:
                raise ValueError('Wrong FDI arch')
            for key in ['upper_fdi', 'lower_fdi']:
                if c[key] not in fdi_expected:
                    raise ValueError('Contact FDI absent from source')
            x = finite_nonnegative(c['force_interval_N'], 'force interval')
            if x.shape != (2,) or x[0] > x[1]:
                raise ValueError('Reversed force interval')
            point = np.array(c['point_cbct_mm'], float)
            if point.shape != (3,) or not np.isfinite(point).all():
                raise ValueError('Invalid contact point')
            if not c.get('source_face_locator') or not c.get('force_error_locator'):
                return dict(status='UNKNOWN', reason='Missing point/error provenance')
            covered.update([c['upper_fdi'], c['lower_fdi']])
        except KeyError as e:
            return dict(status='UNKNOWN', reason='Missing contact field ' + str(e))
    if set(fdi_expected) != covered:
        return dict(status='UNKNOWN', reason='Unmeasured teeth cannot be assigned zero', uncovered_fdi=sorted(set(fdi_expected) - covered))
    if packet.get('verified_zero_fdi') and (not packet.get('zero_detection_bound_locator')):
        return dict(status='UNKNOWN', reason='Zeros need subthreshold force bound')
    total = finite_nonnegative(packet.get('total_force_interval_N'), 'total force')
    if total.shape != (2,) or total[0] > total[1]:
        raise ValueError('Invalid total')
    lows = np.array([c['force_interval_N'][0] for c in contacts])
    highs = np.array([c['force_interval_N'][1] for c in contacts])
    tl = max(total[0], lows.sum())
    th = min(total[1], highs.sum())
    if tl > th:
        raise ValueError('Measured total inconsistent with contact intervals')
    forces = []
    for fdi in fdi_expected:
        inc = np.array([fdi in [c['upper_fdi'], c['lower_fdi']] for c in contacts])
        other = ~inc
        lo = max(lows[inc].sum(), tl - highs[other].sum())
        hi = min(highs[inc].sum(), th - lows[other].sum())
        forces.append(dict(fdi=fdi, force_interval_N=[float(lo), float(hi)], resolution='PER_TOOTH'))
    return dict(status='CONDITIONAL_FIXTURE_STATIC_PORT' if fixture else 'CONDITIONAL_MEASURED_STATIC_PORT', forces=forces, total_force_interval_N=[float(tl), float(th)], joint_constraint='Both arches share one contact vector; totals never independently summed', physical_admission='PENDING_INDEPENDENT_REVIEW', height_response='UNKNOWN_NEEDS_PERTURBATION_MEASUREMENTS', force_direction_scope='Normal scalar only; tangential force UNKNOWN', enclosure='Analytic box+sum marginal; no formal float enclosure', acquisition_id=packet['acquisition_id'])

def fixed_load_decision(interval, limit=100):
    if interval is None:
        return dict(status='UNKNOWN', reason='Missing same-patient absolute force')
    (lo, hi) = finite_nonnegative(interval, 'force interval')
    if lo > hi:
        raise ValueError('Reversed interval')
    return dict(status='ENCLOSED' if hi <= limit else 'EXCEEDS' if lo > limit else 'UNDECIDED', fixed_research_load_N=limit, force_interval_N=[float(lo), float(hi)], resolution='PER_TOOTH')
