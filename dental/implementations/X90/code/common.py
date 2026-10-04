"""Stable serialization and source-bound contracts for the X90 consumer."""
from pathlib import Path
import hashlib
import json
import math
ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda : stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def clean(value):
    if hasattr(value, 'tolist'):
        return clean(value.tolist())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): clean(v) for (k, v) in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, float) and (not math.isfinite(value)):
        return None
    return value

def encoded(value):
    return (json.dumps(clean(value), ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()

def dump(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_bytes(encoded(value))

def read(path):
    return json.loads(Path(path).read_text())

def stable(value):
    """Runtime telemetry belongs in COSTS, never in immutable numerical predictions."""
    value = clean(value)
    if isinstance(value, dict):
        return {k: stable(v) for (k, v) in value.items() if k not in ['created_utc', 'query_s', 'query_wall_s', 'cost']}
    if isinstance(value, list):
        return [stable(v) for v in value]
    return value

def freeze(path, value):
    path = Path(path)
    blob = encoded(value)
    if path.exists() and path.read_bytes() != blob:
        raise ValueError('FROZEN_CONTENT_CHANGED: ' + path.name)
    path.write_bytes(blob)
    path.with_suffix(path.suffix + '.sha256').write_text(hashlib.sha256(blob).hexdigest() + '\n')

def validate_manifest(m):
    if m.get('units') != 'mm':
        raise ValueError('UNIT_NOT_MM')
    if m.get('source_pose') != 'REGISTERED_DIGITAL_PAIR':
        raise ValueError('FRAME_NOT_BOUND')
    if m.get('dataset') != 'Bits2Bites' or not str(m['case']).isdigit():
        raise ValueError('UNSUPPORTED_SUBJECT_CONTRACT')
    if m.get('cbct') is not None:
        c = m['cbct']
        if c.get('subject') != m['patient']:
            raise ValueError('CROSS_SUBJECT_CBCT')
        if c.get('frame') != 'CBCT_NATIVE_MM' or c.get('units') != 'mm':
            raise ValueError('CBCT_FRAME_OR_UNIT')
    return m

def validate_decisions(rows):
    required = {'id', 'title', 'value', 'unit', 'uncertainty', 'evidence', 'resolution', 'timescale', 'sources', 'would_change', 'status'}
    ids = set()
    for row in rows:
        if required - set(row):
            raise ValueError('DECISION_CONTRACT_MISSING: ' + str(required - set(row)))
        if row['id'] in ids:
            raise ValueError('DUPLICATE_DECISION_ID')
        ids.add(row['id'])
        if row['evidence'] not in ['PROVEN', 'CALIBRATED', 'MODELLED', 'UNKNOWN']:
            raise ValueError('EVIDENCE_LABEL')
        if not row['uncertainty'] or not row['sources'] or (not row['would_change']):
            raise ValueError('EMPTY_DECISION_CONTRACT')
        if row['value'] is None and row['evidence'] != 'UNKNOWN':
            raise ValueError('UNKNOWN_VALUE_PROMOTED')
    return True

def validate_force_binding(state, m, geometry_hash):
    if state.get('case') != m['case']:
        raise ValueError('CROSS_SUBJECT_CALIBRATION')
    if state.get('geometry_sha256') != geometry_hash:
        raise ValueError('FORCE_GEOMETRY_MISMATCH')
    if state.get('gap_convention') != 'MEASURED_LEVELED_REFERENCE':
        raise ValueError('FORCE_LOADED_STATE_MISSING')
    if state.get('acquisition_kind') not in ['PHYSICAL_MEASUREMENT', 'SIMULATED_ACQUISITION']:
        raise ValueError('FORCE_ACQUISITION_UNKNOWN')

def validate_pulp_observation(obs, m, geometry_file):
    if obs.get('patient') != m['patient'] or obs.get('unit') != 'mm' or obs.get('quantity') != 'total_hard_tissue_to_pulp' or (obs.get('fdi') != m['fdi']) or (obs.get('resolution') != 'PER_POINT'):
        raise ValueError('PULP_OBSERVATION_BINDING')
    if not obs.get('locator') or not obs.get('input_sha256'):
        raise ValueError('PULP_SOURCE_NOT_BOUND')
    if sha(geometry_file) != obs['input_sha256']:
        raise ValueError('PULP_GEOMETRY_HASH')
    for key in ['distance_mm', 'digital_two_surface_radius_mm']:
        value = obs.get(key)
        if type(value) not in [int, float] or not math.isfinite(value) or value < 0:
            raise ValueError('PULP_DISTANCE_OR_RADIUS')

def decision(ident, title, value, unit, uncertainty, evidence, resolution, status, sources, would_change, **extra):
    return dict(id=ident, title=title, value=value, unit=unit, uncertainty=uncertainty, evidence=evidence, resolution=resolution, timescale='SIMULTANEOUS', status=status, sources=sources, would_change=would_change, **extra)
