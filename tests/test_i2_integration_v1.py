#!/usr/bin/env python3
"""Integration test for region geometry: geometry certificate -> material roster ->
region mass/heat -> frozen I2 payload -> structural validator.

Synthetic fixture only (the repo's own layered box). Plain python3 runnable + pytest compatible.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]

_EX_SPEC = importlib.util.spec_from_file_location(
    'i1_region_example', REPO / 'examples' / 'geometry' / 'i1_region_example.py')
_EX_MOD = importlib.util.module_from_spec(_EX_SPEC)
_EX_SPEC.loader.exec_module(_EX_MOD)
build_example = _EX_MOD.build

from bodytwin.geometry.layered_box_v1 import layered_box  # noqa: E402
from bodytwin.geometry.geometry_validity_v1 import (  # noqa: E402
    GeometryContractError, geometry_certificate, require_geometry_certificate)
from bodytwin.geometry.i2_payload_validate_v1 import assert_valid_payload, validate_payload  # noqa: E402
from bodytwin.geometry.material_roster_v1 import (  # noqa: E402
    ABSENT, KNOWN, is_absent, is_known, lookup, roster_table)
from bodytwin.geometry.region_mass_v1 import (  # noqa: E402
    MASS_BASIS, area_mm2, interface_area_by_pair_mm2, mass_moments_surface, lookup_density)


def outward_faces(box, label):
    return np.concatenate([box['faces'][box['back'] == label],
                           box['faces'][box['front'] == label][:, ::-1]])


def test_chain_certificate_roster_moments_payload():
    box = layered_box()
    faces = outward_faces(box, 2)
    cert = geometry_certificate(box['vertices'], faces, self_intersections=False)
    assert cert['accepted'] and cert['conforming'] and cert['watertight']
    assert abs(cert['volume_mm3'] - 180000.0) < 1e-9 * 180000.0
    assert {p['property_name'] for p in roster_table()} >= {
        'density', 'specific_heat', 'thermal_conductivity', 'perfusion_heat_transport'}
    density = lookup('MSK-MUSCLE', 'density')
    assert is_known(density) and not is_absent(density)
    moments = mass_moments_surface(box['vertices'], faces, density=density,
                                   region_id='SYNTH-BOX-2', basis='region_volume')
    assert moments.known and moments.basis in MASS_BASIS
    assert abs(moments.mass_kg - density.as_density_kg_per_mm3().value * 180000.0) < 1e-18
    assert np.linalg.norm(moments.centroid_mm - np.array([30.0, 30.0, 35.0])) < 1e-9
    assert is_absent(lookup('MSK-MUSCLE', 'thermal_conductivity'))  # absence, not zero
    payload = build_example()
    assert validate_payload(payload) == []
    assert assert_valid_payload(payload) is payload
    assert payload['mass_moments_region_2']['centroid_mm'] == [30.0, 30.0, 35.0]
    assert payload['area']['by_pair'] == {'1-2': 3600.0}
    assert lookup_density('MSK-TENDON') is None


def test_chain_refuses_unevaluated_certificate():
    box = layered_box()
    faces = outward_faces(box, 1)
    unknown = geometry_certificate(box['vertices'], faces)  # self_intersections=None
    assert unknown['accepted'] is False  # abstain is not acceptance
    try:
        require_geometry_certificate(box['vertices'], faces)
        raise AssertionError('unevaluated certificate accepted')
    except GeometryContractError:
        pass
    ok = geometry_certificate(box['vertices'], faces, self_intersections=False)
    assert ok['accepted'] and abs(ok['volume_mm3'] - 36000.0) < 1e-9 * 36000.0


def test_chain_payload_negatives_are_caught():
    payload = build_example()
    payload['regions'][1]['region_id'] = payload['regions'][0]['region_id']
    assert any('duplicate' in e for e in validate_payload(payload))
    payload = build_example()
    payload['regions'][0]['mass'] = 0.03816  # bare number, no units
    assert any('bare numeric' in e or 'quantity must be an object' in e
               for e in validate_payload(payload))
    payload = build_example()
    payload['regions'][0]['mass_basis'] = 'muscle'
    assert any('mass_basis' in e for e in validate_payload(payload))
    assert abs(area_mm2(layered_box(), 'material_interface_area')
               - sum(interface_area_by_pair_mm2(layered_box()).values())) < 1e-9


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for fn in fns:
        fn()
        print('PASS', fn.__name__)
    print(f'{len(fns)} tests passed')
