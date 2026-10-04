"""Prospective observation binding. Metadata alone never establishes physical truth."""
from pathlib import Path
import datetime, hashlib, json, math

def file_sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def iso(s):
    return datetime.datetime.fromisoformat(s.replace('Z', '+00:00'))

def validate(record, freeze, freeze_sha, test_only=False):
    """Same rejecting predicate used for real intake and live fixture mutations."""
    errors = []
    if not test_only and (not freeze.get('physical_geometry_ready', False)):
        return dict(status='UNKNOWN_CONTRACT_NOT_READY', errors=['NEW_FREEZE_WITH_ACTUAL_CAD_HASHES_REQUIRED'], physical_validation=False)
    expected = dict(scan=('um', 'AS_BUILT', 'PER_POINT'), film=('um', 'CEMENTED', 'PER_POINT'), force=('N', 'LOADED', 'PER_SURFACE_REGION'), fracture=('N', 'FRACTURED', 'PER_TOOTH'))
    if record.get('quantity') not in expected:
        return dict(status='REJECT', errors=['UNKNOWN_QUANTITY'])
    (unit, state, level) = expected[record['quantity']]
    for (key, want) in [('unit', unit), ('measurement_state', state), ('resolution_level', level), ('freeze_sha256', freeze_sha)]:
        if record.get(key) != want:
            errors.append('MISMATCH_' + key)
    if record.get('specimen_id') not in freeze['specimen_ids']:
        errors.append('UNREGISTERED_SPECIMEN')
    value = record.get('value')
    bound = record.get('absolute_error_bound')
    if value is None or bound is None:
        return dict(status='UNKNOWN_NOT_MEASURED', errors=errors + ['MISSING_VALUE_OR_METROLOGY_BOUND'])
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        return dict(status='REJECT', errors=errors + ['NONFINITE_VALUE'], physical_validation=False)
    if not isinstance(bound, (int, float)) or not math.isfinite(bound) or bound < 0:
        return dict(status='REJECT', errors=errors + ['BAD_ERROR_BOUND'], physical_validation=False)
    if record['quantity'] in ['film', 'force', 'fracture'] and value < 0:
        errors.append('NEGATIVE_VALUE')
    try:
        if iso(freeze['timestamp_utc']) >= iso(record['acquired_at_utc']):
            errors.append('FREEZE_NOT_BEFORE_OBSERVATION')
    except (ValueError, KeyError, TypeError):
        errors.append('MISSING_OR_BAD_TIMESTAMP')
    for key in ['registration_locator', 'calibration_locator', 'operator_id', 'geometry_sha256']:
        if not record.get(key):
            errors.append('MISSING_' + key)
    if record.get('geometry_sha256') != freeze['geometry_source_sha256']:
        errors.append('GEOMETRY_HASH_MISMATCH')
    origin = record.get('fracture_origin')
    if record['quantity'] == 'fracture' and origin not in ['intaglio_tensile_zone', 'contact_zone', 'cement_interface', 'die', 'unknown']:
        errors.append('INVALID_FRACTURE_ORIGIN')
    declared = record.get('provenance_kind')
    locator = record.get('raw_path', '')
    fixture = declared == 'our_own_fixture' or any((s in locator.lower() for s in ['fixture', 'virtual', 'synthetic']))
    if fixture and (not test_only):
        errors.append('FIXTURE_IS_NOT_PHYSICAL_MEASUREMENT')
    if declared not in ['independent_measurement', 'our_own_fixture']:
        errors.append('BAD_PROVENANCE_KIND')
    try:
        raw = json.loads(Path(locator).read_text())
        if file_sha(locator) != record.get('raw_sha256'):
            errors.append('RAW_HASH_MISMATCH')
        for k in ['value', 'unit', 'specimen_id', 'measurement_state', 'quantity', 'resolution_level']:
            if raw.get(k) != record.get(k):
                errors.append('RAW_RECORD_MISMATCH_' + k)
    except (OSError, ValueError, TypeError):
        errors.append('RAW_ARTIFACT_UNAVAILABLE')
    if errors:
        return dict(status='REJECT', errors=errors, physical_validation=False)
    return dict(status='TEST_ONLY_ACCEPTED' if test_only else 'EVIDENCE_BOUND_PENDING_INDEPENDENT_REVIEW', value_interval=[value - bound, value + bound], mechanism_status='UNKNOWN' if origin == 'unknown' else 'OBSERVED_ORIGIN_NOT_MODEL_VALIDATION', physical_validation=False, errors=[], provenance_scope='Raw hash/record consistency; calibration/registration truth still needs independent review')
