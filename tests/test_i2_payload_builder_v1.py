#!/usr/bin/env python3
"""Tests for the reusable frozen ``i2_region_example_v1`` payload builder.

Synthetic fixture only (the repo's own layered box); no physiological claims. Plain python3
runnable + pytest compatible.

Run:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        PYTHONPATH=src python tests/test_i2_payload_builder_v1.py
    OMP_NUM_THREADS=1 python \
        -m pytest -q test_i2_payload_builder_v1.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

from bodytwin.geometry.layered_box_v1 import layered_box  # noqa: E402
from bodytwin.geometry.i2_payload_builder_v1 import (  # noqa: E402
    PayloadBuildError, build_and_check, build_payload)
from bodytwin.geometry.i2_payload_validate_v1 import validate_payload  # noqa: E402
from bodytwin.geometry.i2_region_registration_v1 import registration_issues  # noqa: E402
from bodytwin.geometry.region_mass_v1 import MaterialProperty  # noqa: E402

# Fixed tolerances: the layered box is exact (analytic volumes/areas), so these only absorb
# floating-point reduction noise.
VOLUME_TOL = 1e-9
AREA_TOL = 1e-9

CP = MaterialProperty('muscle_cp_assumption', 3490.0, 'J/(kg*K)', 'assumption',
                      'src/bodytwin/cells/organ_systems/thermoregulation.py::SPECIFIC_HEAT_BODY_J_PER_KG_K')
SOURCE_SHA256 = {'layered_box_v1.py': 'test-fixture-sha256'}


def _spec(region_id, material_region_id, density_key, **extra):
    spec = {'region_id': region_id, 'material_region_id': material_region_id,
            'density_key': density_key, 'specific_heat': CP}
    spec.update(extra)
    return spec


def _specs(first_density_key='MSK-MUSCLE', first_material_region_id='MSK-MUSCLE',
           with_moments=False):
    return [
        _spec('SYNTH-BOX-1', first_material_region_id, first_density_key),
        _spec('SYNTH-BOX-2', 'MSK-MUSCLE', 'MSK-MUSCLE', with_moments=with_moments),
    ]


# ------------------------------------------------------------------ positive control
def test_build_validates_and_registers():
    payload = build_and_check(layered_box(), _specs(), source_sha256=SOURCE_SHA256)
    assert payload['schema'] == 'i2_region_example_v1'
    assert payload['synthetic'] is True
    assert payload['source_sha256'] == SOURCE_SHA256
    assert validate_payload(payload) == []
    assert registration_issues(payload) == []


def test_region_volumes_and_interface_pairs():
    payload = build_payload(layered_box(), _specs())
    volumes = [region['volume_mm3'] for region in payload['regions']]
    assert abs(volumes[0] - 36000.0) <= VOLUME_TOL * 36000.0, volumes
    assert abs(volumes[1] - 180000.0) <= VOLUME_TOL * 180000.0, volumes
    assert payload['area']['by_pair'] == {'1-2': 3600.0}, payload['area']['by_pair']
    assert abs(payload['area']['value'] - 3600.0) <= AREA_TOL * 3600.0


# ------------------------------------------------------------------ missing material
def test_missing_density_is_null_and_reported_by_roster():
    payload = build_payload(
        layered_box(),
        _specs(first_density_key='MSK-TENDON', first_material_region_id='MSK-TENDON'))
    region = payload['regions'][0]
    assert region['density'] is None, region['density']
    assert 'density' in region['missing'], region['missing']
    assert validate_payload(payload) == []
    issues = registration_issues(payload)
    assert any('MSK-TENDON' in issue and 'MISSING' in issue for issue in issues), issues


def test_missing_material_declaration_is_refused():
    spec = _spec('SYNTH-BOX-1', 'MSK-MUSCLE', 'MSK-MUSCLE')
    del spec['density_key']
    try:
        build_payload(layered_box(), [spec, _spec('SYNTH-BOX-2', 'MSK-MUSCLE', 'MSK-MUSCLE')])
    except PayloadBuildError as error:
        assert 'density' in str(error), error
    else:
        raise AssertionError('an undeclared density input must be refused')


# ------------------------------------------------------------------ duplicates
def test_duplicate_region_ids_are_refused():
    specs = _specs()
    specs[1]['region_id'] = specs[0]['region_id']
    try:
        build_payload(layered_box(), specs)
    except PayloadBuildError as error:
        assert 'duplicate' in str(error).lower(), error
    else:
        raise AssertionError('duplicate region ids must be refused via the validator')


# ------------------------------------------------------------------ moments
def test_with_moments_block_validates():
    payload = build_payload(layered_box(), _specs(with_moments=True))
    keys = [key for key in payload if key.startswith('mass_moments')]
    assert keys, 'with_moments=True must emit a mass_moments block'
    block = payload[keys[0]]
    assert block['units'] == {'mass': 'kg', 'centroid': 'mm', 'inertia': 'kg*mm2'}
    assert block['mass_kg'] is not None
    assert len(block['centroid_mm']) == 3
    assert isinstance(block['inertia_kg_mm2'], list)
    assert block['inertia_SI_units']
    assert validate_payload(payload) == []


if __name__ == '__main__':
    functions = [value for name, value in sorted(globals().items())
                 if name.startswith('test_') and callable(value)]
    for function in functions:
        function()
        print('PASS', function.__name__)
    print(f'{len(functions)} tests passed')
