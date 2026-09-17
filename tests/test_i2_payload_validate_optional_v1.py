#!/usr/bin/env python3
"""Tests for the OPTIONAL payload fields added to the I2 -> I1 structural validator.

Run:
    OMP_NUM_THREADS=1 PYTHONPATH=src python \
        tests/test_i2_payload_validate_optional_v1.py
    OMP_NUM_THREADS=1 PYTHONPATH=src python \
        -m pytest -q tests/test_i2_payload_validate_optional_v1.py

Scope: the four optional fields ``area.area_to_volume``, ``regions[].frame``,
``regions[].uncertainty`` and ``regions[].thickness_distribution``. Absent fields never
produce an error; present-but-malformed fields always do. The required-field rules are not
changed and are re-checked here so a regression would be caught.

Structural exactness / tolerances (documented, fixed for this file)
-------------------------------------------------------------------
* All checks are structural and exact. The validator never recomputes a value, so this file
  applies **no numeric tolerance or epsilon**: values are only classified as finite /
  positive / non-negative with ``math.isfinite`` and ``isinstance``.
* Unit and enum comparisons are exact string equality (e.g. ``'1/mm'``, ``'mm'``,
  ``'absolute' | 'relative'``). No normalisation, no case folding.
* ``origin_mm`` must be a list of exactly 3 finite numbers; ``axes`` a list of exactly 3 rows
  of exactly 3 finite numbers. Orthonormality / handedness is deliberately NOT tested here
  (that is ``frame_v1``'s job).
* ``uncertainty`` must be a nonempty list; ``samples_mm`` a nonempty list of finite strictly
  positive numbers.
* The tests load the real ``examples/geometry/region_example.json``. The extra values used to populate
  optional fields are arbitrary synthetic numbers; no physiological claim is made and no
  number from the example is asserted.
"""
import copy
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / 'examples' / 'geometry' / 'region_example.json'

from bodytwin.geometry.i2_payload_validate_v1 import assert_valid_payload, validate_payload  # noqa: E402


def load_example() -> dict:
    with open(EXAMPLE) as handle:
        return json.load(handle)


def _has_error(errors, needle: str) -> bool:
    return any(needle.lower() in str(error).lower() for error in errors)


# ---------------------------------------------------------------- optional-field fixtures
def _area_to_volume(**overrides) -> dict:
    obj = {'value': 0.1, 'units': '1/mm'}
    obj.update(overrides)
    return obj


def _frame(**overrides) -> dict:
    obj = {
        'origin_mm': [0.0, 0.0, 0.0],
        'axes': [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        'pitch_mm': 1.0,
        'units': 'mm',
        'declared': True,
    }
    obj.update(overrides)
    return obj


def _uncertainty(**overrides) -> list:
    entry = {
        'name': 'density',
        'value': 0.01,
        'units': 'kg/L',
        'kind': 'absolute',
        'provenance': 'assumption',
    }
    entry.update(overrides)
    return [entry]


def _thickness_distribution(**overrides) -> dict:
    obj = {
        'samples_mm': [3.0, 3.5, 4.0],
        'units': 'mm',
        'method': 'synthetic sampling',
        'provenance': 'assumption',
    }
    obj.update(overrides)
    return obj


def _with_optionals(**fields):
    """Return a deep copy of the example with the named optional fields populated."""
    payload = copy.deepcopy(load_example())
    builders = {
        'area_to_volume': ('area', _area_to_volume()),
        'frame': ('region', _frame()),
        'uncertainty': ('region', _uncertainty()),
        'thickness_distribution': ('region', _thickness_distribution()),
    }
    for name, value in fields.items():
        scope, _default = builders[name]
        if scope == 'area':
            payload['area'][name] = value
        else:
            payload['regions'][0][name] = value
    return payload


def _with_all_optionals():
    return _with_optionals(
        area_to_volume=_area_to_volume(),
        frame=_frame(),
        uncertainty=_uncertainty(),
        thickness_distribution=_thickness_distribution(),
    )


# ---------------------------------------------------------------- positive controls
def test_real_example_validates_empty():
    payload = load_example()
    errors = validate_payload(payload)
    assert errors == [], f'unexpected errors on the real example: {errors}'
    assert assert_valid_payload(payload) is payload


def test_all_optional_fields_valid_payload_passes():
    payload = _with_all_optionals()
    errors = validate_payload(payload)
    assert errors == [], f'a valid payload with all optional fields must pass: {errors}'
    assert assert_valid_payload(payload) is payload


def test_each_optional_field_alone_passes():
    single = {
        'area_to_volume': _area_to_volume(),
        'frame': _frame(),
        'uncertainty': _uncertainty(),
        'thickness_distribution': _thickness_distribution(),
    }
    for name, value in single.items():
        payload = _with_optionals(**{name: value})
        errors = validate_payload(payload)
        assert errors == [], f'optional {name!r} alone must pass: {errors}'


def test_uncertainty_relative_units_one_is_valid():
    payload = _with_optionals(
        uncertainty=_uncertainty(kind='relative', units='1', provenance='computed'))
    errors = validate_payload(payload)
    assert errors == [], f"relative uncertainty with units '1' must pass: {errors}"


# ---------------------------------------------------------------- required-field regression
def test_required_rules_unaffected():
    """The pre-existing required-field rules must still fire after this extension."""
    invalid_synthetic = copy.deepcopy(load_example())
    invalid_synthetic['synthetic'] = 'yes'
    assert _has_error(validate_payload(invalid_synthetic), 'synthetic')

    bare_mass = copy.deepcopy(load_example())
    bare_mass['regions'][0]['mass'] = 0.03816
    assert _has_error(validate_payload(bare_mass), 'object')

    duplicate = copy.deepcopy(load_example())
    duplicate['regions'][1]['region_id'] = duplicate['regions'][0]['region_id']
    assert _has_error(validate_payload(duplicate), 'duplicate')


# ---------------------------------------------------------------- area_to_volume negatives
def test_negative_area_to_volume_wrong_units():
    payload = _with_optionals(area_to_volume=_area_to_volume(units='1/m'))
    errors = validate_payload(payload)
    assert errors, "area_to_volume units '1/m' must be rejected"
    assert _has_error(errors, 'area.area_to_volume.units'), errors


def test_negative_area_to_volume_nonpositive_value():
    payload = _with_optionals(area_to_volume=_area_to_volume(value=0.0))
    errors = validate_payload(payload)
    assert errors, 'a non-positive area_to_volume value must be rejected'
    assert _has_error(errors, 'area.area_to_volume.value'), errors


def test_negative_area_to_volume_not_object():
    payload = _with_optionals(area_to_volume=0.1)
    errors = validate_payload(payload)
    assert errors, 'a bare number area_to_volume must be rejected'
    assert _has_error(errors, 'area.area_to_volume'), errors


# ---------------------------------------------------------------- frame negatives
def test_negative_frame_bad_axes_shape():
    payload = _with_optionals(frame=_frame(axes=[[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]))
    errors = validate_payload(payload)
    assert errors, 'a non-3x3 axes matrix must be rejected'
    assert _has_error(errors, 'regions[0].frame.axes'), errors


def test_negative_frame_units_metres():
    payload = _with_optionals(frame=_frame(units='m'))
    errors = validate_payload(payload)
    assert errors, "frame units 'm' must be rejected (mm guard)"
    assert _has_error(errors, 'regions[0].frame.units'), errors


def test_negative_frame_bad_origin_length():
    payload = _with_optionals(frame=_frame(origin_mm=[0.0, 0.0]))
    errors = validate_payload(payload)
    assert errors, 'an origin_mm of length 2 must be rejected'
    assert _has_error(errors, 'regions[0].frame.origin_mm'), errors


def test_negative_frame_nonpositive_pitch():
    payload = _with_optionals(frame=_frame(pitch_mm=0.0))
    errors = validate_payload(payload)
    assert errors, 'a non-positive frame pitch must be rejected'
    assert _has_error(errors, 'regions[0].frame.pitch_mm'), errors


def test_negative_frame_declared_not_bool():
    payload = _with_optionals(frame=_frame(declared='yes'))
    errors = validate_payload(payload)
    assert errors, 'a non-bool declared flag must be rejected'
    assert _has_error(errors, 'regions[0].frame.declared'), errors


# ---------------------------------------------------------------- uncertainty negatives
def test_negative_uncertainty_empty_list():
    payload = _with_optionals(uncertainty=[])
    errors = validate_payload(payload)
    assert errors, 'an empty uncertainty list must be rejected'
    assert _has_error(errors, 'regions[0].uncertainty'), errors


def test_negative_uncertainty_bad_kind():
    payload = _with_optionals(uncertainty=_uncertainty(kind='systematic'))
    errors = validate_payload(payload)
    assert errors, "uncertainty kind 'systematic' must be rejected"
    assert _has_error(errors, 'regions[0].uncertainty[0].kind'), errors


def test_negative_uncertainty_relative_with_mm_units():
    payload = _with_optionals(uncertainty=_uncertainty(kind='relative', units='mm'))
    errors = validate_payload(payload)
    assert errors, "a relative uncertainty with units 'mm' must be rejected"
    assert _has_error(errors, 'regions[0].uncertainty[0].units'), errors


def test_negative_uncertainty_negative_value():
    payload = _with_optionals(uncertainty=_uncertainty(value=-0.01))
    errors = validate_payload(payload)
    assert errors, 'a negative uncertainty value must be rejected'
    assert _has_error(errors, 'regions[0].uncertainty[0].value'), errors


def test_negative_uncertainty_bad_provenance():
    payload = _with_optionals(uncertainty=_uncertainty(provenance='guessed'))
    errors = validate_payload(payload)
    assert errors, 'an unknown uncertainty provenance must be rejected'
    assert _has_error(errors, 'regions[0].uncertainty[0].provenance'), errors


# ---------------------------------------------------------------- thickness negatives
def test_negative_thickness_samples_empty():
    payload = _with_optionals(thickness_distribution=_thickness_distribution(samples_mm=[]))
    errors = validate_payload(payload)
    assert errors, 'empty thickness samples must be rejected'
    assert _has_error(errors, 'regions[0].thickness_distribution.samples_mm'), errors


def test_negative_thickness_zero_sample():
    payload = _with_optionals(
        thickness_distribution=_thickness_distribution(samples_mm=[3.0, 0.0]))
    errors = validate_payload(payload)
    assert errors, 'a zero thickness sample must be rejected'
    assert _has_error(errors, 'regions[0].thickness_distribution.samples_mm'), errors


def test_negative_thickness_units_metres():
    payload = _with_optionals(thickness_distribution=_thickness_distribution(units='m'))
    errors = validate_payload(payload)
    assert errors, "thickness units 'm' must be rejected (mm guard)"
    assert _has_error(errors, 'regions[0].thickness_distribution.units'), errors


def test_negative_thickness_bad_provenance():
    payload = _with_optionals(
        thickness_distribution=_thickness_distribution(provenance='guessed'))
    errors = validate_payload(payload)
    assert errors, 'an unknown thickness provenance must be rejected'
    assert _has_error(errors, 'regions[0].thickness_distribution.provenance'), errors


# ---------------------------------------------------------------- non-object optionals
def test_negative_optionals_wrong_container_type():
    payload = _with_optionals(
        frame='frame',
        uncertainty={'name': 'x'},
        thickness_distribution=[1.0, 2.0],
    )
    errors = validate_payload(payload)
    assert errors, 'wrong optional container types must be rejected'
    assert _has_error(errors, 'regions[0].frame'), errors
    assert _has_error(errors, 'regions[0].uncertainty'), errors
    assert _has_error(errors, 'regions[0].thickness_distribution'), errors


def main() -> int:
    tests = [(name, obj) for name, obj in sorted(globals().items())
             if name.startswith('test_') and callable(obj)]
    failures = []
    for name, test in tests:
        try:
            test()
        except Exception as exc:  # noqa: BLE001 - report every failure, do not stop early
            failures.append((name, exc))
            print(f'FAIL {name}: {exc!r}')
        else:
            print(f'PASS {name}')
    print(f'{len(tests) - len(failures)}/{len(tests)} passed')
    if failures:
        for name, exc in failures:
            print(f'  {name}: {exc!r}')
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
