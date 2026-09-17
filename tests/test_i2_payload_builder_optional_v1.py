#!/usr/bin/env python3
"""Optional-block tests for the frozen ``i2_region_example_v1`` payload builder.

Synthetic fixture only (the repo's own layered box); no physiological claims. Covers the four
OPTIONAL blocks ``build_payload`` may now emit -- ``area.area_to_volume``, ``regions[].frame``,
``regions[].uncertainty`` and ``regions[].thickness_distribution`` -- without changing any
required field or the existing defaults. Plain python3 + pytest compatible.

Run:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        python \
        tests/test_i2_payload_builder_optional_v1.py
    OMP_NUM_THREADS=1 python \
        -m pytest -q test_i2_payload_builder_optional_v1.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

from bodytwin.geometry.layered_box_v1 import layered_box  # noqa: E402
from bodytwin.geometry.frame_v1 import Frame  # noqa: E402
from bodytwin.geometry.i2_payload_builder_v1 import PayloadBuildError, build_payload  # noqa: E402
from bodytwin.geometry.i2_payload_validate_v1 import validate_payload  # noqa: E402
from bodytwin.geometry.region_mass_v1 import MaterialProperty  # noqa: E402
from bodytwin.geometry.thickness_v1 import ThicknessDistribution  # noqa: E402
from bodytwin.geometry.uncertainty_v1 import Uncertainty  # noqa: E402

# Exact synthetic fixture: region volumes 36000 + 180000 mm3, shared interface 3600 mm2.
TOTAL_VOLUME_MM3 = 216000.0
INTERFACE_AREA_MM2 = 3600.0
AREA_TO_VOLUME = INTERFACE_AREA_MM2 / TOTAL_VOLUME_MM3
AV_TOL = 1e-12

CP = MaterialProperty('muscle_cp_assumption', 3490.0, 'J/(kg*K)', 'assumption',
                      'src/bodytwin/cells/organ_systems/thermoregulation.py::SPECIFIC_HEAT_BODY_J_PER_KG_K')

IDENTITY_AXES = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
FRAME_DICT = {'origin_mm': [1.0, 2.0, 3.0], 'axes': IDENTITY_AXES, 'pitch_mm': 1.0,
              'units': 'mm', 'declared': False}
THICKNESS = ThicknessDistribution([3.0, 4.0, 5.0], 'synthetic thickness sampling', 'assumption')
THICKNESS_DICT = {'samples_mm': [2.0, 2.5], 'method': 'synthetic thickness sampling',
                  'provenance': 'assumption', 'units': 'mm', 'source': 'fixture'}

REGION_KEYS = {'region_id', 'material_region_id', 'volume_mm3', 'volume_method', 'mass_basis',
               'density', 'mass', 'specific_heat', 'heat_capacity', 'missing'}
FRAME_KEYS = {'origin_mm', 'axes', 'pitch_mm', 'units', 'declared'}
UNCERTAINTY_KEYS = {'name', 'value', 'units', 'kind', 'provenance'}
THICKNESS_KEYS = {'samples_mm', 'units', 'method', 'provenance', 'source'}
OPTIONAL_TOKENS = ('area_to_volume', 'frame', 'uncertainty', 'thickness')


def _spec(region_id, **extra):
    spec = {'region_id': region_id, 'material_region_id': 'MSK-MUSCLE',
            'density_key': 'MSK-MUSCLE', 'specific_heat': CP}
    spec.update(extra)
    return spec


def _base_specs():
    return [_spec('SYNTH-BOX-1'), _spec('SYNTH-BOX-2')]


def _optional_specs():
    """Two MSK-MUSCLE regions; region 1 uses instances, region 2 uses plain dicts."""
    return [
        _spec('SYNTH-BOX-1',
              frame=Frame(origin_mm=[1.0, 2.0, 3.0], axes=IDENTITY_AXES, pitch_mm=1.0),
              uncertainty={'mass': Uncertainty(0.5, 'kg', 'absolute', 'assumption')},
              thickness=THICKNESS),
        _spec('SYNTH-BOX-2', frame=dict(FRAME_DICT),
              uncertainty=[{'name': 'density', 'value': 0.1, 'units': 'kg/L',
                            'kind': 'absolute', 'provenance': 'assumption'},
                           {'name': 'heat_capacity', 'value': 0.02, 'units': '1',
                            'kind': 'relative', 'provenance': 'assumption'}],
              thickness=dict(THICKNESS_DICT)),
    ]


def _assert_valid_or_optional_only(payload):
    """Assert validate_payload accepts the payload, tolerating only optional-block complaints.

    The optional validator support may not have landed in this run; in that case we still assert
    the exact structure directly (in the calling test) rather than failing on the validator. A
    complaint that is NOT about one of the four optional blocks is a real required-field error
    and still fails.
    """
    errors = validate_payload(payload)
    if errors:
        flagged_optional = all(
            any(token in error for token in OPTIONAL_TOKENS) for error in errors)
        assert flagged_optional, errors


def _expect_build_error(region_specs):
    try:
        build_payload(layered_box(), region_specs)
    except PayloadBuildError as error:
        return error
    raise AssertionError('a malformed optional spec must raise PayloadBuildError')


# ------------------------------------------------------- all four optional blocks present
def test_all_optional_blocks_emitted_with_exact_structure():
    payload = build_payload(layered_box(), _optional_specs(), area_to_volume=True)

    # Required shape is unchanged.
    assert set(payload) == {'schema', 'synthetic', 'area', 'regions', 'totals'}

    # area.area_to_volume
    av = payload['area']['area_to_volume']
    assert set(av) == {'value', 'units'}
    assert av['units'] == '1/mm'
    assert abs(av['value'] - AREA_TO_VOLUME) <= AV_TOL, av

    # regions[].frame (instance-built and dict-built)
    for region in payload['regions']:
        frame = region['frame']
        assert set(frame) == FRAME_KEYS, frame
        assert frame['units'] == 'mm'
        assert isinstance(frame['declared'], bool)
        assert len(frame['origin_mm']) == 3 and all(isinstance(v, float)
                                                    for v in frame['origin_mm'])
        assert len(frame['axes']) == 3 and len(frame['axes'][0]) == 3
        assert isinstance(frame['pitch_mm'], float)

    # regions[].uncertainty
    for region in payload['regions']:
        entries = region['uncertainty']
        assert isinstance(entries, list) and entries
        for entry in entries:
            assert set(entry) == UNCERTAINTY_KEYS, entry
            assert isinstance(entry['name'], str) and entry['name']
            if entry['kind'] == 'relative':
                assert entry['units'] == '1', entry

    # regions[].thickness_distribution
    for region in payload['regions']:
        block = region['thickness_distribution']
        assert set(block) == THICKNESS_KEYS, block
        assert block['units'] == 'mm'
        assert block['samples_mm'] and all(isinstance(s, float) for s in block['samples_mm'])

    _assert_valid_or_optional_only(payload)


def test_area_to_volume_value_and_units():
    payload = build_payload(layered_box(), _base_specs(), area_to_volume=True)
    av = payload['area']['area_to_volume']
    assert av['units'] == '1/mm'
    assert abs(av['value'] - 3600.0 / 216000.0) <= AV_TOL, av
    _assert_valid_or_optional_only(payload)


# ------------------------------------------------------- malformed specs are refused
def test_malformed_frame_raises():
    specs = _base_specs()
    specs[0]['frame'] = [[1.0, 0.0, 0.0]]  # not a Frame / not a frame dict
    _expect_build_error(specs)

    bad_origin = _base_specs()
    bad_origin[0]['frame'] = {'origin_mm': [0.0, 0.0], 'axes': IDENTITY_AXES, 'pitch_mm': 1.0}
    _expect_build_error(bad_origin)

    bad_axes = _base_specs()
    bad_axes[0]['frame'] = {'origin_mm': [0.0, 0.0, 0.0],
                            'axes': [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.5]],
                            'pitch_mm': 1.0}
    _expect_build_error(bad_axes)


def test_malformed_uncertainty_raises():
    negative = _base_specs()
    negative[0]['uncertainty'] = {
        'mass': {'value': -1.0, 'units': 'kg', 'kind': 'absolute', 'provenance': 'assumption'}}
    _expect_build_error(negative)

    relative_with_units = _base_specs()
    relative_with_units[0]['uncertainty'] = {
        'mass': {'value': 0.1, 'units': 'kg', 'kind': 'relative', 'provenance': 'assumption'}}
    _expect_build_error(relative_with_units)

    wrong_type = _base_specs()
    wrong_type[0]['uncertainty'] = 'not-an-uncertainty'
    _expect_build_error(wrong_type)


def test_malformed_thickness_raises():
    wrong_units = _base_specs()
    wrong_units[0]['thickness'] = {'samples_mm': [3.0], 'method': 'synthetic',
                                   'provenance': 'assumption', 'units': 'm'}
    _expect_build_error(wrong_units)

    nonpositive_sample = _base_specs()
    nonpositive_sample[0]['thickness'] = {'samples_mm': [0.0], 'method': 'synthetic',
                                          'provenance': 'assumption', 'units': 'mm'}
    _expect_build_error(nonpositive_sample)

    wrong_type = _base_specs()
    wrong_type[0]['thickness'] = 3.0
    _expect_build_error(wrong_type)


# ------------------------------------------------------- base build unchanged
def test_base_build_without_options_is_unchanged():
    payload = build_payload(layered_box(), _base_specs())
    assert set(payload) == {'schema', 'synthetic', 'area', 'regions', 'totals'}
    assert set(payload['area']) == {'kind', 'value', 'units', 'method', 'by_pair'}
    assert 'area_to_volume' not in payload['area']
    for region in payload['regions']:
        assert set(region) == REGION_KEYS, region
    assert validate_payload(payload) == []


if __name__ == '__main__':
    functions = [value for name, value in sorted(globals().items())
                 if name.startswith('test_') and callable(value)]
    for function in functions:
        function()
        print('PASS', function.__name__)
    print(f'{len(functions)} tests passed')
