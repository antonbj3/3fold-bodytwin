#!/usr/bin/env python3
"""Acceptance tests for the geometry validity certificate.

Run:
    OMP_NUM_THREADS=1 PYTHONPATH=src python tests/test_geometry_validity_v1.py
    OMP_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q tests/test_geometry_validity_v1.py

Tolerances were fixed BEFORE the run: exact synthetic polytopes are compared at 1e-9
relative (analytic). All geometry here is synthetic; no physiological claims are made.
"""
import sys
from pathlib import Path

import numpy as np
import trimesh

REPO = Path(__file__).resolve().parents[1]

from bodytwin.geometry.geometry_validity_v1 import (  # noqa: E402
    GeometryContractError,
    geometry_certificate,
    require_geometry_certificate,
)

RTOL = 1e-9  # fixed before the run for exact polytopes (float64 / analytic)


def box_vertices_faces(side):
    m = trimesh.creation.box(extents=[float(side)] * 3)
    return np.asarray(m.vertices, dtype=np.float64), np.asarray(m.faces)


def open_box():
    v, f = box_vertices_faces(1.0)
    return v, np.delete(f, 0, axis=0)


def two_cubes_touching():
    """Two separate 1 mm cubes touching at z=1; shared face has duplicate coordinates."""
    a = trimesh.creation.box(extents=[1.0, 1.0, 1.0])
    a.apply_translation([0.5, 0.5, 0.5])
    b = trimesh.creation.box(extents=[1.0, 1.0, 1.0])
    b.apply_translation([0.5, 0.5, 1.5])
    va = np.asarray(a.vertices, dtype=np.float64)
    fa = np.asarray(a.faces)
    vb = np.asarray(b.vertices, dtype=np.float64)
    fb = np.asarray(b.faces)
    return np.vstack([va, vb]), np.vstack([fa, fb + len(va)])


def test_analytic_box_accepted_and_volume():
    v, f = box_vertices_faces(10.0)
    cert = geometry_certificate(v, f, self_intersections=False)
    assert cert['accepted'] is True
    assert cert['units'] == 'mm'
    assert cert['watertight'] is True
    assert cert['winding_consistent'] is True
    assert cert['positive_volume'] is True
    assert cert['conforming'] is True
    assert cert['self_intersections'] is False
    assert abs(cert['volume_mm3'] - 1000.0) <= RTOL * 1000.0
    assert cert['tolerance_mm'] == 1e-6
    assert isinstance(cert['reason'], str) and cert['reason']


def test_open_mesh_rejected():
    v, f = open_box()
    cert = geometry_certificate(v, f, self_intersections=False)
    assert cert['accepted'] is False
    assert cert['watertight'] is False
    assert 'closed' in cert['reason'].lower() or 'watertight' in cert['reason'].lower()


def test_nonconforming_assembly_rejected():
    v, f = two_cubes_touching()
    cert = geometry_certificate(v, f, self_intersections=False)
    # surface_mesh itself accepts the closed assembly (two closed boxes) ...
    assert cert['watertight'] is True
    assert abs(cert['volume_mm3'] - 2.0) <= RTOL * 2.0
    # ... but the coincident split plane makes it non-conforming and therefore rejected.
    assert cert['conforming'] is False
    assert cert['accepted'] is False
    assert 'conform' in cert['reason'].lower()


def test_unevaluated_self_intersection_abstains():
    v, f = box_vertices_faces(2.0)
    cert = geometry_certificate(v, f)  # self_intersections defaults to None
    assert cert['self_intersections'] is None
    assert cert['accepted'] is False
    assert 'self-intersection' in cert['reason'].lower()


def test_self_intersections_true_rejected():
    v, f = box_vertices_faces(1.0)
    cert = geometry_certificate(v, f, self_intersections=True)
    assert cert['accepted'] is False
    assert 'self-intersection' in cert['reason'].lower()


def test_require_geometry_certificate_returns_for_healthy():
    v, f = box_vertices_faces(3.0)
    cert = require_geometry_certificate(v, f, self_intersections=False)
    assert cert['accepted'] is True
    assert abs(cert['volume_mm3'] - 27.0) <= RTOL * 27.0


def test_require_geometry_certificate_raises_for_rejected():
    cases = [
        (open_box(), {'self_intersections': False}),
        (two_cubes_touching(), {'self_intersections': False}),
        (box_vertices_faces(1.0), {}),  # unevaluated self-intersection abstains
    ]
    for (v, f), kwargs in cases:
        try:
            require_geometry_certificate(v, f, **kwargs)
            assert False, 'rejected geometry must raise GeometryContractError'
        except GeometryContractError:
            pass


def test_declared_units_m_not_accepted():
    v, f = box_vertices_faces(1.0)
    cert = geometry_certificate(v, f, self_intersections=False, declared_units='m')
    assert cert['accepted'] is False
    assert 'mm' in cert['reason'].lower()


def test_bad_tolerance_raises():
    v, f = box_vertices_faces(1.0)
    for bad in (0.0, -1.0, float('nan'), float('inf')):
        try:
            geometry_certificate(v, f, tolerance_mm=bad, self_intersections=False)
            assert False, f'tolerance_mm={bad!r} must raise'
        except ValueError:
            pass


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for fn in fns:
        fn()
        print('PASS', fn.__name__)
    print(f'{len(fns)} tests passed')
