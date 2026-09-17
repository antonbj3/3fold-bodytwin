"""Structural validator for the frozen I2 -> I1 hand-off payload
(``schema == 'i2_region_example_v1'``). It checks units, provenance, ``missing``-vs-``value``
rules, mass bases and region-id uniqueness so that a consumer (the result envelope in
``bodytwin.framework``) can refuse a malformed region/mass/heat record *before* sealing it.

Location: src/bodytwin/geometry/i2_payload_validate_v1.py
No existing symbol is modified and nothing is written.

Design constraints:
  * pure function: no I/O, no imports beyond ``math`` (verified by tests/test_i2_payload_validate_v1.py);
  * structure only -- the example's numeric values are never hardcoded or recomputed;
  * a bare number is never accepted where a quantity object is required;
  * mass quantities (``regions[].mass``, ``totals.mass``) must be in ``kg`` and heat capacities
    in ``J/K``; any other unit label is an error (no silent conversion);
  * ``missing != default``: ``value is None`` requires a nonempty ``missing`` list;
  * optional fields (``area.area_to_volume``, ``regions[].frame``,
    ``regions[].uncertainty``, ``regions[].thickness_distribution``) are validated only when
    the key is present; an absent optional field never produces an error. Frame validation is
    structural only -- axis orthonormality / handedness is ``frame_v1``'s job and is not
    re-checked here.

Cross-reference: the area kinds, mass bases and provenance labels mirror
``region_mass_v1`` (``AREA_KINDS``, ``MASS_BASIS``, ``PROVENANCE``).
"""
from __future__ import annotations

import math

SCHEMA_ID = 'i2_region_example_v1'

# Identical to region_mass_v1.AREA_KINDS (the area kinds the region payload emits).
I2_AREA_KINDS = ('exterior_area', 'material_interface_area', 'functional_exchange_area')

# Identical to region_mass_v1.MASS_BASIS.
I2_MASS_BASIS = ('region_volume', 'active_muscle', 'total_muscle', 'total_body')

# Payload-level provenance set: region_mass_v1.PROVENANCE plus 'computed' (e.g. mass_fraction).
I2_PROVENANCE = ('measured', 'calibrated', 'literal_cited', 'assumption', 'computed', 'UNKNOWN')

QUANTITY_KEYS = ('name', 'value', 'units', 'provenance', 'missing')

# Required unit tags of the mass-moments block (see i1_region_example.py).
MASS_MOMENT_UNITS = {'mass': 'kg', 'centroid': 'mm', 'inertia': 'kg*mm2'}
# Declared units of the mass and heat-capacity quantity objects. A value under any other unit
# label (e.g. grams) is refused, never accepted unchanged and never converted silently.
MASS_UNITS = 'kg'
HEAT_CAPACITY_UNITS = 'J/K'

_UNSET = object()


def _is_number(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _is_finite_number(value) -> bool:
    return _is_number(value) and math.isfinite(float(value))


def _check_quantity(errors, path, obj, *, expected_units=None):
    """Validate one quantity object {name, value, units, provenance, missing}."""
    if not isinstance(obj, dict):
        errors.append(
            f'{path}: quantity must be an object with keys {list(QUANTITY_KEYS)}; got '
            f'{type(obj).__name__} -- a bare numeric mass/heat value is not allowed')
        return
    for key in QUANTITY_KEYS:
        if key not in obj:
            errors.append(f'{path}: missing required key {key!r}')
    name = obj.get('name', _UNSET)
    if name is not _UNSET and (not isinstance(name, str) or not name):
        errors.append(f'{path}.name: nonempty string required')
    units = obj.get('units', _UNSET)
    if units is not _UNSET:
        if not isinstance(units, str) or not units:
            errors.append(f'{path}.units: nonempty string required')
        elif expected_units is not None and units != expected_units:
            errors.append(f'{path}.units: expected {expected_units!r}, got {units!r}')
    provenance = obj.get('provenance', _UNSET)
    if provenance is not _UNSET and provenance not in I2_PROVENANCE:
        errors.append(f'{path}.provenance: {provenance!r} not in {list(I2_PROVENANCE)}')
    missing = obj.get('missing', _UNSET)
    if missing is not _UNSET and not isinstance(missing, list):
        errors.append(f'{path}.missing: list required, got {type(missing).__name__}')
    value = obj.get('value', _UNSET)
    if value is not _UNSET:
        if value is None:
            if not isinstance(missing, list) or not missing:
                errors.append(
                    f'{path}: value None requires a nonempty missing list (missing != default)')
        elif not _is_finite_number(value):
            errors.append(f'{path}.value: finite number or null required, got {value!r}')


def _check_material(errors, path, obj):
    """A bound material value must carry units and provenance, or be an explicit null."""
    if obj is None:
        return  # caller checks that the region declares a matching missing entry
    if not isinstance(obj, dict):
        errors.append(f'{path}: material object (or explicit null) required, '
                      f'got {type(obj).__name__}')
        return
    for key in ('name', 'units', 'provenance'):
        if key not in obj:
            errors.append(f'{path}: missing required key {key!r}')
    name = obj.get('name', _UNSET)
    if name is not _UNSET and (not isinstance(name, str) or not name):
        errors.append(f'{path}.name: nonempty string required')
    units = obj.get('units', _UNSET)
    if units is not _UNSET and (not isinstance(units, str) or not units):
        errors.append(f'{path}.units: nonempty string required')
    provenance = obj.get('provenance', _UNSET)
    if provenance is not _UNSET and provenance not in I2_PROVENANCE:
        errors.append(f'{path}.provenance: {provenance!r} not in {list(I2_PROVENANCE)}')


def _check_area_to_volume(errors, path, obj):
    """Optional ``area.area_to_volume``: ``{value: finite > 0, units: '1/mm'}`` exactly."""
    if not isinstance(obj, dict):
        errors.append(f"{path}: object with keys ['value', 'units'] required, "
                      f'got {type(obj).__name__}')
        return
    if 'value' not in obj:
        errors.append(f'{path}.value: required finite positive number')
    elif not _is_finite_number(obj['value']) or float(obj['value']) <= 0:
        errors.append(f'{path}.value: finite positive number required, got {obj["value"]!r}')
    units = obj.get('units', _UNSET)
    if units is _UNSET:
        errors.append(f"{path}.units: required '1/mm'")
    elif units != '1/mm':
        errors.append(f"{path}.units: expected '1/mm', got {units!r}")


def _check_frame(errors, path, obj):
    """Optional ``regions[].frame``; structural only (no orthonormality test)."""
    if not isinstance(obj, dict):
        errors.append(f'{path}: frame object required, got {type(obj).__name__}')
        return
    for key in ('origin_mm', 'axes', 'pitch_mm', 'units', 'declared'):
        if key not in obj:
            errors.append(f'{path}.{key}: required key missing')

    origin = obj.get('origin_mm', _UNSET)
    if origin is not _UNSET:
        if (not isinstance(origin, list) or len(origin) != 3
                or not all(_is_finite_number(component) for component in origin)):
            errors.append(f'{path}.origin_mm: list of exactly 3 finite numbers required, '
                          f'got {origin!r}')

    axes = obj.get('axes', _UNSET)
    if axes is not _UNSET:
        ok = (isinstance(axes, list) and len(axes) == 3
              and all(isinstance(row, list) and len(row) == 3
                      and all(_is_finite_number(component) for component in row)
                      for row in axes))
        if not ok:
            errors.append(f'{path}.axes: 3x3 nested list of finite numbers required, '
                          f'got {axes!r}')

    pitch = obj.get('pitch_mm', _UNSET)
    if pitch is not _UNSET:
        if not _is_finite_number(pitch) or float(pitch) <= 0:
            errors.append(f'{path}.pitch_mm: finite positive number required, got {pitch!r}')

    units = obj.get('units', _UNSET)
    if units is not _UNSET and units != 'mm':
        errors.append(f"{path}.units: expected 'mm', got {units!r}")

    declared = obj.get('declared', _UNSET)
    if declared is not _UNSET and not isinstance(declared, bool):
        errors.append(f'{path}.declared: bool required, got {type(declared).__name__}')


def _check_uncertainty(errors, path, obj):
    """Optional ``regions[].uncertainty``: nonempty list of uncertainty entries."""
    if not isinstance(obj, list):
        errors.append(f'{path}: nonempty list required, got {type(obj).__name__}')
        return
    if not obj:
        errors.append(f'{path}: nonempty list required (at least one uncertainty entry)')
        return
    for index, entry in enumerate(obj):
        entry_path = f'{path}[{index}]'
        if not isinstance(entry, dict):
            errors.append(f'{entry_path}: uncertainty object required, '
                          f'got {type(entry).__name__}')
            continue
        for key in ('name', 'value', 'units', 'kind', 'provenance'):
            if key not in entry:
                errors.append(f'{entry_path}.{key}: required key missing')

        name = entry.get('name', _UNSET)
        if name is not _UNSET and (not isinstance(name, str) or not name):
            errors.append(f'{entry_path}.name: nonempty string required')

        value = entry.get('value', _UNSET)
        if value is not _UNSET:
            if not _is_finite_number(value) or float(value) < 0:
                errors.append(f'{entry_path}.value: finite number >= 0 required, got {value!r}')

        units = entry.get('units', _UNSET)
        if units is not _UNSET and (not isinstance(units, str) or not units):
            errors.append(f'{entry_path}.units: nonempty string required')

        kind = entry.get('kind', _UNSET)
        if kind is not _UNSET and kind not in ('absolute', 'relative'):
            errors.append(f"{entry_path}.kind: {kind!r} not in ['absolute', 'relative']")

        provenance = entry.get('provenance', _UNSET)
        if provenance is not _UNSET and provenance not in I2_PROVENANCE:
            errors.append(f'{entry_path}.provenance: {provenance!r} not in '
                          f'{list(I2_PROVENANCE)}')

        if kind == 'relative' and units != '1':
            errors.append(f"{entry_path}.units: relative uncertainty must declare units '1', "
                          f'got {units!r}')


def _check_thickness_distribution(errors, path, obj):
    """Optional ``regions[].thickness_distribution``: sample of positive mm thicknesses."""
    if not isinstance(obj, dict):
        errors.append(f'{path}: thickness_distribution object required, '
                      f'got {type(obj).__name__}')
        return
    for key in ('samples_mm', 'units', 'method', 'provenance'):
        if key not in obj:
            errors.append(f'{path}.{key}: required key missing')

    samples = obj.get('samples_mm', _UNSET)
    if samples is not _UNSET:
        ok = (isinstance(samples, list) and bool(samples)
              and all(_is_finite_number(sample) and float(sample) > 0 for sample in samples))
        if not ok:
            errors.append(f'{path}.samples_mm: nonempty list of finite positive numbers '
                          f'required, got {samples!r}')

    units = obj.get('units', _UNSET)
    if units is not _UNSET and units != 'mm':
        errors.append(f"{path}.units: expected 'mm', got {units!r}")

    method = obj.get('method', _UNSET)
    if method is not _UNSET and (not isinstance(method, str) or not method):
        errors.append(f'{path}.method: nonempty string required')

    provenance = obj.get('provenance', _UNSET)
    if provenance is not _UNSET and provenance not in I2_PROVENANCE:
        errors.append(f'{path}.provenance: {provenance!r} not in {list(I2_PROVENANCE)}')


def _check_region(errors, path, region):
    """Validate one region record; return its region_id when it is a usable string."""
    if not isinstance(region, dict):
        errors.append(f'{path}: region must be an object, got {type(region).__name__}')
        return None
    region_id = region.get('region_id', _UNSET)
    if region_id is _UNSET or not isinstance(region_id, str) or not region_id:
        errors.append(f'{path}.region_id: nonempty string required')
        region_id = None
    volume = region.get('volume_mm3', _UNSET)
    if volume is _UNSET:
        errors.append(f'{path}.volume_mm3: required (mm3)')
    elif not _is_finite_number(volume) or float(volume) <= 0:
        errors.append(f'{path}.volume_mm3: finite positive number required, got {volume!r}')
    method = region.get('volume_method', _UNSET)
    if method is _UNSET or not isinstance(method, str) or not method:
        errors.append(f'{path}.volume_method: nonempty string required')
    basis = region.get('mass_basis', _UNSET)
    if basis is _UNSET:
        errors.append(f'{path}.mass_basis: required, one of {list(I2_MASS_BASIS)}')
    elif basis not in I2_MASS_BASIS:
        errors.append(f'{path}.mass_basis: {basis!r} not in {list(I2_MASS_BASIS)}')

    region_missing = region.get('missing', _UNSET)
    if region_missing is _UNSET:
        errors.append(f'{path}.missing: list required (use [] when nothing is missing)')
        region_missing = []
    elif not isinstance(region_missing, list):
        errors.append(f'{path}.missing: list required, got {type(region_missing).__name__}')
        region_missing = []

    for material_key in ('density', 'specific_heat'):
        if material_key not in region:
            errors.append(f'{path}.{material_key}: required (or explicit null with a '
                          f'matching missing entry)')
            continue
        material = region[material_key]
        if material is None:
            if not any(isinstance(entry, str) and material_key in entry
                       for entry in region_missing):
                errors.append(f'{path}.{material_key}: explicit null requires a matching entry '
                              f'in {path}.missing')
        else:
            _check_material(errors, f'{path}.{material_key}', material)

    if 'mass' in region:
        _check_quantity(errors, f'{path}.mass', region['mass'], expected_units=MASS_UNITS)
    else:
        errors.append(f'{path}.mass: required quantity object')
    if 'heat_capacity' in region:
        _check_quantity(errors, f'{path}.heat_capacity', region['heat_capacity'],
                        expected_units='J/K')
    else:
        errors.append(f"{path}.heat_capacity: required quantity object with units 'J/K'")

    if 'frame' in region:
        _check_frame(errors, f'{path}.frame', region['frame'])
    if 'uncertainty' in region:
        _check_uncertainty(errors, f'{path}.uncertainty', region['uncertainty'])
    if 'thickness_distribution' in region:
        _check_thickness_distribution(errors, f'{path}.thickness_distribution',
                                      region['thickness_distribution'])
    return region_id


def _check_mass_moments(errors, path, moments):
    if not isinstance(moments, dict):
        errors.append(f'{path}: mass_moments must be an object, got {type(moments).__name__}')
        return
    units = moments.get('units')
    if not isinstance(units, dict):
        errors.append(f'{path}.units: object required with {MASS_MOMENT_UNITS}')
    else:
        for key, required in MASS_MOMENT_UNITS.items():
            if units.get(key) != required:
                errors.append(f'{path}.units.{key}: expected {required!r}, '
                              f'got {units.get(key)!r}')
    si_units = moments.get('inertia_SI_units')
    if not isinstance(si_units, str) or not si_units:
        errors.append(f'{path}.inertia_SI_units: nonempty string required')


def validate_payload(payload) -> list[str]:
    """Return a list of human-readable errors; an empty list means the payload is valid.

    Structure only: no example value is hardcoded and nothing is recomputed.
    """
    errors: list[str] = []
    if not isinstance(payload, dict):
        return [f'payload: JSON object required, got {type(payload).__name__}']

    if payload.get('schema') != SCHEMA_ID:
        errors.append(f"schema: expected {SCHEMA_ID!r}, got {payload.get('schema')!r}")
    if not isinstance(payload.get('synthetic'), bool):
        errors.append(f'synthetic: bool required, got {type(payload.get("synthetic")).__name__}')

    area = payload.get('area')
    if not isinstance(area, dict):
        errors.append(f'area: object required, got {type(area).__name__}')
    else:
        kind = area.get('kind')
        if kind not in I2_AREA_KINDS:
            errors.append(f'area.kind: {kind!r} not in {list(I2_AREA_KINDS)}')
        if 'value' not in area:
            errors.append('area.value: required')
        elif not _is_finite_number(area['value']):
            errors.append(f'area.value: finite number required, got {area["value"]!r}')
        units = area.get('units')
        if not isinstance(units, str) or not units:
            errors.append('area.units: nonempty string required')
        if 'area_to_volume' in area:
            _check_area_to_volume(errors, 'area.area_to_volume', area['area_to_volume'])

    regions = payload.get('regions')
    if not isinstance(regions, list) or not regions:
        errors.append('regions: nonempty list required')
        regions = []

    region_ids: list[str] = []
    bases: set[str] = set()
    for index, region in enumerate(regions):
        region_id = _check_region(errors, f'regions[{index}]', region)
        if region_id is not None:
            region_ids.append(region_id)
        if isinstance(region, dict) and region.get('mass_basis') in I2_MASS_BASIS:
            bases.add(region['mass_basis'])

    seen: set[str] = set()
    duplicates: list[str] = []
    for region_id in region_ids:
        if region_id in seen and region_id not in duplicates:
            duplicates.append(region_id)
        seen.add(region_id)
    if duplicates:
        errors.append(f'regions: duplicate region_id {duplicates!r}; region ids must be unique')

    totals = payload.get('totals')
    if not isinstance(totals, dict):
        errors.append(f'totals: object required, got {type(totals).__name__}')
    else:
        for key, obj in totals.items():
            _check_quantity(errors, f'totals.{key}', obj,
                            expected_units={'mass': MASS_UNITS,
                                            'heat_capacity': HEAT_CAPACITY_UNITS}.get(key))
        for required in ('mass', 'heat_capacity'):
            if required not in totals:
                errors.append(f'totals.{required}: required quantity object')
        if len(bases) > 1:
            errors.append(f'totals.mass: mixed mass_basis across regions {sorted(bases)}; '
                          f'a total must not mix bases')

    for key, obj in payload.items():
        if isinstance(key, str) and key.startswith('mass_moments'):
            _check_mass_moments(errors, key, obj)

    return errors


def assert_valid_payload(payload) -> dict:
    """Return ``payload`` unchanged when valid; raise ``ValueError`` with all errors otherwise."""
    errors = validate_payload(payload)
    if errors:
        raise ValueError('i2 payload invalid:\n- ' + '\n- '.join(errors))
    return payload
