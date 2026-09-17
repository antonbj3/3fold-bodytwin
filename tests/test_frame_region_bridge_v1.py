#!/usr/bin/env python3
"""Acceptance tests for the frame<->region-mass-moments bridge.

Run:
    OMP_NUM_THREADS=1 PYTHONPATH=src python tests/test_frame_region_bridge_v1.py
    OMP_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q tests/test_frame_region_bridge_v1.py

Tolerances were fixed BEFORE the run: exact synthetic box geometry (trimesh) is compared
at relative 1e-9. All data here is synthetic; no physiological claims are made.
"""
import sys
from pathlib import Path

import numpy as np
import trimesh

REPO = Path(__file__).resolve().parents[1]

from bodytwin.geometry.frame_v1 import Frame, transfer_points, transfer_tensors  # noqa: E402
from bodytwin.geometry.region_mass_v1 import MaterialProperty, mass_moments_surface  # noqa: E402
from bodytwin.geometry.frame_region_bridge_v1 import (  # noqa: E402
    FrameBindingError,
    bind_region_moments,
    moments_in_frame,
)

RTOL = 1e-9
ATOL = 1e-9

REGION_ID = 'region_box'
DENSITY = MaterialProperty('muscle_density', 1.06, 'kg/L', 'literal_cited')
CENTROID_SHIFT = np.array([7.0, -3.0, 5.0])
ROT_Z90 = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
EXPECTED_KEYS = {'region_id', 'mass_kg', 'centroid_mm', 'inertia_kg_mm2', 'frame',
                 'units', 'missing'}
EXPECTED_UNITS = {'mass': 'kg', 'centroid': 'mm', 'inertia': 'kg*mm2'}


def assert_close(actual, expected):
    np.testing.assert_allclose(
        np.asarray(actual, dtype=np.float64),
        np.asarray(expected, dtype=np.float64),
        rtol=RTOL,
        atol=ATOL,
    )


def box_vertices_faces():
    """A synthetic 10x20x30 mm box, translated so its centroid is not the origin."""
    mesh = trimesh.creation.box(extents=[10.0, 20.0, 30.0])
    mesh.apply_translation(CENTROID_SHIFT)
    return np.asarray(mesh.vertices, dtype=np.float64), np.asarray(mesh.faces)


def global_moments(density=DENSITY):
    vertices, faces = box_vertices_faces()
    return mass_moments_surface(vertices, faces, density=density, region_id=REGION_ID,
                                basis='region_volume')


def expect_binding_error(call):
    try:
        call()
    except FrameBindingError:
        return
    raise AssertionError('expected FrameBindingError was not raised')


def test_identity_frame_reproduces_global_moments():
    moments = global_moments()
    frame = Frame()
    result = moments_in_frame(moments, frame)

    assert set(result) == EXPECTED_KEYS
    assert result['region_id'] == REGION_ID
    assert result['frame'] is frame
    assert result['units'] == EXPECTED_UNITS
    assert result['missing'] == []
    assert_close(result['mass_kg'], moments.mass_kg)
    assert_close(result['centroid_mm'], moments.centroid_mm)
    assert_close(result['inertia_kg_mm2'], moments.inertia_kg_mm2)
    # default frame is identity: local coordinates are the shared mm coordinates.
    assert_close(result['centroid_mm'], CENTROID_SHIFT)


def test_pure_translation_shifts_centroid_inverse_offset_inertia_unchanged():
    moments = global_moments()
    offset = np.array([4.0, -1.5, 2.0])
    frame = Frame(origin_mm=tuple(offset))
    result = moments_in_frame(moments, frame)

    # same axes => from_global subtracts the origin.
    assert_close(result['centroid_mm'], np.asarray(moments.centroid_mm) - offset)
    # the moments inertia is about the centroid, so translation leaves it unchanged.
    assert_close(result['inertia_kg_mm2'], moments.inertia_kg_mm2)


def test_rotation_only_inertia_is_axes_transpose_I_axes():
    moments = global_moments()
    frame = Frame(axes=ROT_Z90)
    result = moments_in_frame(moments, frame)

    inertia_global = np.asarray(moments.inertia_kg_mm2)
    expected = ROT_Z90.T @ inertia_global @ ROT_Z90
    assert_close(result['inertia_kg_mm2'], expected)
    # congruence by an orthogonal matrix preserves the eigenvalues.
    assert_close(np.linalg.eigvalsh(result['inertia_kg_mm2']),
                 np.linalg.eigvalsh(inertia_global))
    # the centroid is rotated into the frame, not translated.
    assert_close(result['centroid_mm'], CENTROID_SHIFT @ ROT_Z90)


def test_equivalent_to_transfer_tensors_between_identity_and_frame():
    moments = global_moments()
    identity = Frame()
    frame = Frame(origin_mm=(3.0, -2.0, 5.0), axes=ROT_Z90)
    result = moments_in_frame(moments, frame)

    inertia_global = np.asarray(moments.inertia_kg_mm2)
    centroid_global = np.asarray(moments.centroid_mm)
    expected_inertia = transfer_tensors(identity, frame, inertia_global[None, :, :])[0]
    expected_centroid = transfer_points(identity, frame, centroid_global[None, :])[0]
    assert_close(result['inertia_kg_mm2'], expected_inertia)
    assert_close(result['centroid_mm'], expected_centroid)


def test_missing_density_returns_none_moments_and_missing_list_without_raising():
    vertices, faces = box_vertices_faces()
    moments = mass_moments_surface(vertices, faces, density=None, region_id=REGION_ID,
                                   basis='region_volume')
    assert moments.known is False

    result = moments_in_frame(moments, Frame())
    assert result['mass_kg'] is None
    assert result['centroid_mm'] is None
    assert result['inertia_kg_mm2'] is None
    assert isinstance(result['missing'], list)
    assert result['missing']
    assert 'density' in result['missing']

    # the binding helper abstains the same way instead of raising.
    bound = bind_region_moments(REGION_ID, vertices, faces, None, Frame())
    assert bound['mass_kg'] is None
    assert bound['centroid_mm'] is None
    assert bound['inertia_kg_mm2'] is None
    assert bound['missing']


def test_bind_region_moments_matches_moments_in_frame():
    frame = Frame(origin_mm=(1.0, 2.0, -4.0), axes=ROT_Z90)
    vertices, faces = box_vertices_faces()
    bound = bind_region_moments(REGION_ID, vertices, faces, DENSITY, frame)
    expected = moments_in_frame(global_moments(), frame)

    assert set(bound) == EXPECTED_KEYS
    assert bound['region_id'] == REGION_ID
    assert bound['frame'] is frame
    assert bound['units'] == EXPECTED_UNITS
    assert_close(bound['mass_kg'], expected['mass_kg'])
    assert_close(bound['centroid_mm'], expected['centroid_mm'])
    assert_close(bound['inertia_kg_mm2'], expected['inertia_kg_mm2'])


def test_non_frame_argument_raises_frame_binding_error():
    assert issubclass(FrameBindingError, ValueError)

    moments = global_moments()
    for bad in (None, object(), 'frame', ROT_Z90):
        expect_binding_error(lambda bad=bad: moments_in_frame(moments, bad))

    vertices, faces = box_vertices_faces()
    expect_binding_error(
        lambda: bind_region_moments(REGION_ID, vertices, faces, DENSITY, 'not-a-frame'))


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for fn in fns:
        fn()
        print('PASS', fn.__name__)
    print(f'{len(fns)} tests passed')
