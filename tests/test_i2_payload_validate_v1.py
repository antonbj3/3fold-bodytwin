#!/usr/bin/env python3
"""Tests for the region-payload validator.

Run:
    OMP_NUM_THREADS=1 PYTHONPATH=src python \
        tests/test_i2_payload_validate_v1.py
    OMP_NUM_THREADS=1 PYTHONPATH=src python \
        -m pytest -q tests/test_i2_payload_validate_v1.py

The tests are structural: they load the real ``examples/geometry/region_example.json`` and apply
negative controls that each remove or corrupt one structural fact. No physiological claim is
made and no numeric value from the example is asserted.
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


# ---------------------------------------------------------------- positive control
def test_real_example_validates_empty():
    payload = load_example()
    errors = validate_payload(payload)
    assert errors == [], f'unexpected errors on the real example: {errors}'
    assert assert_valid_payload(payload) is payload


# ---------------------------------------------------------------- negative controls
def test_negative_missing_units_on_mass_quantity():
    payload = copy.deepcopy(load_example())
    del payload['regions'][0]['mass']['units']
    errors = validate_payload(payload)
    assert errors, 'removing a mass quantity unit must produce an error'
    assert _has_error(errors, 'units'), errors


def test_negative_bare_number_mass():
    payload = copy.deepcopy(load_example())
    payload['regions'][0]['mass'] = 0.03816
    errors = validate_payload(payload)
    assert errors, 'a bare numeric mass must produce an error'
    assert _has_error(errors, 'object') or _has_error(errors, 'bare'), errors


def test_negative_empty_provenance():
    payload = copy.deepcopy(load_example())
    payload['regions'][0]['mass']['provenance'] = ''
    errors = validate_payload(payload)
    assert errors, "an empty provenance must produce an error"
    assert _has_error(errors, 'provenance'), errors


def test_negative_duplicate_region_id():
    payload = copy.deepcopy(load_example())
    payload['regions'][1]['region_id'] = payload['regions'][0]['region_id']
    errors = validate_payload(payload)
    assert errors, 'two regions with the same id must produce an error'
    assert _has_error(errors, 'duplicate'), errors


def test_negative_invalid_mass_basis():
    payload = copy.deepcopy(load_example())
    payload['regions'][0]['mass_basis'] = 'muscle'
    errors = validate_payload(payload)
    assert errors, "mass_basis 'muscle' is not in the I2 enum and must error"
    assert _has_error(errors, 'mass_basis'), errors


def test_negative_invalid_area_kind():
    payload = copy.deepcopy(load_example())
    payload['area']['kind'] = 'exterior_wall'
    errors = validate_payload(payload)
    assert errors, "'exterior_wall' is not an I2 area kind and must error"
    assert _has_error(errors, 'area.kind'), errors


def test_negative_synthetic_missing():
    payload = copy.deepcopy(load_example())
    del payload['synthetic']
    errors = validate_payload(payload)
    assert errors, "'synthetic' is required and must be a bool"
    assert _has_error(errors, 'synthetic'), errors


def test_negative_value_present_with_missing_entry_is_allowed():
    """A present value with a caveat list is structurally allowed (missing != default is
    about value None, not about a present value carrying an annotation)."""
    payload = copy.deepcopy(load_example())
    payload['regions'][0]['mass']['missing'] = ['density']
    errors = validate_payload(payload)
    assert errors == [], f'a present value with a nonempty missing list must stay valid: {errors}'


def test_negative_value_none_with_empty_missing_errors():
    payload = copy.deepcopy(load_example())
    payload['regions'][0]['mass']['value'] = None
    payload['regions'][0]['mass']['missing'] = []
    errors = validate_payload(payload)
    assert errors, 'value None with an empty missing list must error'
    assert _has_error(errors, 'value None') or _has_error(errors, 'missing'), errors


def test_negative_assert_valid_payload_raises():
    payload = copy.deepcopy(load_example())
    payload['synthetic'] = 'yes'
    try:
        assert_valid_payload(payload)
    except ValueError as exc:
        assert 'synthetic' in str(exc), str(exc)
    else:
        raise AssertionError('assert_valid_payload must raise ValueError on an invalid payload')


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


def test_mass_units_must_be_kg():
    """A mass labelled in grams is refused, never accepted unchanged or converted silently."""
    payload = load_example()
    payload['regions'][0]['mass']['units'] = 'g'
    payload['totals']['mass']['units'] = 'g'
    errors = validate_payload(payload)
    assert "regions[0].mass.units: expected 'kg', got 'g'" in errors
    assert "totals.mass.units: expected 'kg', got 'g'" in errors
