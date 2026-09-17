#!/usr/bin/env python3
"""Acceptance tests for the frame / coordinate descriptor.

Run:
    OMP_NUM_THREADS=1 PYTHONPATH=src python tests/test_frame_v1.py
    OMP_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q tests/test_frame_v1.py

Tolerances were fixed BEFORE the run: exact synthetic rotations/translations are
compared at relative AND absolute 1e-9 (float64 analytic). All data here is synthetic;
no physiological claims are made.
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]

from bodytwin.geometry.frame_v1 import (  # noqa: E402
    Frame,
    FrameError,
    from_global,
    to_global,
    transfer_points,
    transfer_tensors,
    transfer_vectors,
)

RTOL = 1e-9
ATOL = 1e-9


def assert_close(actual, expected):
    np.testing.assert_allclose(
        np.asarray(actual, dtype=np.float64),
        np.asarray(expected, dtype=np.float64),
        rtol=RTOL,
        atol=ATOL,
    )


def rot_z(deg):
    t = np.deg2rad(float(deg))
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def expect_frame_error(call):
    try:
        call()
    except FrameError:
        return
    raise AssertionError('expected FrameError was not raised')


def test_identity_round_trip():
    f = Frame()
    p = np.array([[1.0, 2.0, 3.0], [-4.0, 5.0, -6.0]])
    assert_close(to_global(f, p), p)
    assert_close(from_global(f, p), p)
    assert_close(from_global(f, to_global(f, p)), p)
    v = np.array([[0.0, 0.0, 1.0], [2.0, -1.0, 0.5]])
    assert_close(transfer_vectors(f, f, v), v)
    T = np.array([np.diag([1.0, 2.0, 3.0])])
    assert_close(transfer_tensors(f, f, T), T)


def test_pure_translation_points_shift_vectors_unchanged():
    a = Frame()
    b = Frame(origin_mm=(10.0, -2.0, 0.5))
    p = np.array([[1.0, 1.0, 1.0], [0.0, 0.0, 0.0]])
    # a and b share axes, so points shift by exactly -origin_mm (in b coordinates).
    assert_close(transfer_points(a, b, p), p - np.array([10.0, -2.0, 0.5]))
    # vectors carry no translation.
    v = np.array([[1.0, 2.0, 3.0], [-4.0, 5.0, -6.0]])
    assert_close(transfer_vectors(a, b, v), v)


def test_rotation_90_about_z_known_coordinates():
    f = Frame(axes=rot_z(90.0))
    local = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    expected = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert_close(to_global(f, local), expected)
    assert_close(from_global(f, expected), local)


def test_to_from_global_inverse():
    f = Frame(origin_mm=(3.0, -1.0, 2.0), axes=rot_z(37.0))
    p = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, -3.0], [-5.0, 4.0, 0.5]])
    assert_close(from_global(f, to_global(f, p)), p)
    g = np.array([[3.0, -1.0, 2.0], [10.0, 0.0, -7.0]])
    assert_close(to_global(f, from_global(f, g)), g)


def test_a_to_b_to_a_round_trip():
    a = Frame(origin_mm=(1.0, -2.0, 3.0), axes=rot_z(30.0))
    b = Frame(origin_mm=(-4.0, 0.5, 2.0), axes=rot_z(115.0))
    p = np.array([[0.0, 0.0, 0.0], [1.5, -2.0, 0.25], [-3.0, 4.0, 5.0]])
    assert_close(transfer_points(b, a, transfer_points(a, b, p)), p)
    v = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 2.0, 3.0]])
    assert_close(transfer_vectors(b, a, transfer_vectors(a, b, v)), v)
    T = np.array([
        np.diag([1.0, 2.0, 3.0]),
        np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]),
    ])
    assert_close(transfer_tensors(b, a, transfer_tensors(a, b, T)), T)


def test_vector_transfer_translation_invariant_point_not():
    a = Frame()
    b1 = Frame(origin_mm=(10.0, 0.0, 0.0))
    b2 = Frame(origin_mm=(0.0, 10.0, 0.0))
    v = np.array([[1.0, 2.0, 3.0], [-4.0, 5.0, -6.0]])
    assert_close(transfer_vectors(a, b1, v), transfer_vectors(a, b2, v))
    p = np.array([[1.0, 1.0, 1.0]])
    assert not np.allclose(
        transfer_points(a, b1, p), transfer_points(a, b2, p), rtol=RTOL, atol=ATOL)


def test_tensor_covariance_matches_explicit_rotation():
    a = Frame()
    b = Frame(axes=rot_z(90.0))
    # R = frame_b.axes.T @ frame_a.axes; for axes_a = I this is rot_z(-90).
    R_explicit = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert_close(b.axes.T @ a.axes, R_explicit)
    T = np.array([[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]])
    assert_close(transfer_tensors(a, b, T), R_explicit @ T @ R_explicit.T)


def test_axes_stored_read_only():
    f = Frame()
    assert f.axes.flags.writeable is False
    try:
        f.axes[0, 0] = 5.0
    except ValueError:
        return
    raise AssertionError('axes must be read-only')


def test_reject_units_not_mm():
    for bad in ('m', 'cm', ''):
        expect_frame_error(lambda bad=bad: Frame(units=bad))


def test_reject_non_orthonormal_axes():
    expect_frame_error(lambda: Frame(axes=np.diag([1.0, 2.0, 1.0])))
    expect_frame_error(lambda: Frame(axes=1e-3 * np.eye(3)))


def test_reject_left_handed_axes():
    bad = np.eye(3)
    bad[2, 2] = -1.0
    expect_frame_error(lambda: Frame(axes=bad))


def test_reject_negative_or_zero_pitch():
    for bad in (-1.0, 0.0, -0.5):
        expect_frame_error(lambda bad=bad: Frame(pitch_mm=bad))
    expect_frame_error(lambda: Frame(pitch_mm=float('nan')))
    expect_frame_error(lambda: Frame(pitch_mm=float('inf')))


def test_reject_nonunit_pitch_without_declaration():
    expect_frame_error(lambda: Frame(pitch_mm=2.0))
    expect_frame_error(lambda: Frame(pitch_mm=0.5))
    f = Frame(pitch_mm=2.0, pitch_declared=True)
    assert f.pitch_mm == 2.0


def test_reject_wrong_tensor_shape():
    f = Frame()
    for bad in (np.zeros((3, 3)), np.zeros((2, 3, 3, 1)), np.zeros((1, 3))):
        expect_frame_error(lambda bad=bad: transfer_tensors(f, f, bad))


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for fn in fns:
        fn()
        print('PASS', fn.__name__)
    print(f'{len(fns)} tests passed')
