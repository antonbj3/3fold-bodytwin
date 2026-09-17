#!/usr/bin/env python3
"""Acceptance tests for *evaluated* self-intersection in the I2 geometry certificate.

Run:
    OMP_NUM_THREADS=1 PYTHONPATH=src python tests/test_geometry_validity_selfintersect_v1.py
    OMP_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q tests/test_geometry_validity_selfintersect_v1.py

All geometry is synthetic; no physiological claims are made. The tri-state contract under test:
only an evaluated ``False`` may accept, while ``None`` (unavailable/unevaluated) abstains.
"""
import sys
from pathlib import Path

import numpy as np
import pytest
import trimesh

REPO = Path(__file__).resolve().parents[1]

from bodytwin.geometry import geometry_validity_v1 as gv  # noqa: E402
from bodytwin.geometry.geometry_validity_v1 import (  # noqa: E402
    GeometryContractError,
    evaluate_self_intersections,
    geometry_certificate,
    require_geometry_certificate,
)


def box_vertices_faces(side):
    m = trimesh.creation.box(extents=[float(side)] * 3)
    return np.asarray(m.vertices, dtype=np.float64), np.asarray(m.faces)


def two_overlapping_boxes():
    """Two individually watertight 2 mm cubes that interpenetrate (no shared vertices)."""
    a = trimesh.creation.box(extents=[2.0, 2.0, 2.0])
    a.apply_translation([1.0, 1.0, 1.0])
    b = trimesh.creation.box(extents=[2.0, 2.0, 2.0])
    b.apply_translation([2.0, 2.0, 2.0])
    va = np.asarray(a.vertices, dtype=np.float64)
    fa = np.asarray(a.faces)
    vb = np.asarray(b.vertices, dtype=np.float64)
    fb = np.asarray(b.faces)
    return np.vstack([va, vb]), np.vstack([fa, fb + len(va)])


def test_direct_evaluator_clean_box_is_false():
    v, f = box_vertices_faces(2.0)
    assert evaluate_self_intersections(v, f) is False


def test_clean_box_evaluated_and_accepted():
    v, f = box_vertices_faces(2.0)
    cert = geometry_certificate(v, f, evaluate=True)
    assert cert['self_intersections'] is False
    assert cert['accepted'] is True
    assert cert['reason'] and 'accepted' in cert['reason'].lower()


def test_evaluate_false_and_no_declaration_abstains():
    v, f = box_vertices_faces(2.0)
    for cert in (geometry_certificate(v, f), geometry_certificate(v, f, evaluate=False)):
        assert cert['self_intersections'] is None
        assert cert['accepted'] is False
        assert 'self-intersection' in cert['reason'].lower()


def test_unavailable_evaluation_abstains_not_accepts(monkeypatch):
    v, f = box_vertices_faces(2.0)
    monkeypatch.setattr(gv, 'evaluate_self_intersections', lambda *a, **k: None)
    cert = geometry_certificate(v, f, evaluate=True)
    assert cert['self_intersections'] is None
    assert cert['accepted'] is False
    with pytest.raises(GeometryContractError):
        require_geometry_certificate(v, f, evaluate=True)


def test_declaration_wins_over_evaluation():
    v, f = box_vertices_faces(2.0)
    # An explicit declaration is never overwritten by evaluation.
    cert = geometry_certificate(v, f, self_intersections=True, evaluate=True)
    assert cert['self_intersections'] is True
    assert cert['accepted'] is False
    clean = geometry_certificate(v, f, self_intersections=False, evaluate=True)
    assert clean['self_intersections'] is False
    assert clean['accepted'] is True


def test_self_intersecting_fixture_rejected_when_evaluated():
    v, f = two_overlapping_boxes()
    # The assembled surface is still watertight/positive-volume, so only the
    # self-intersection evaluation can reject it.
    assert evaluate_self_intersections(v, f) is True
    cert = geometry_certificate(v, f, evaluate=True)
    assert cert['self_intersections'] is True
    assert cert['accepted'] is False
    assert 'self-intersection' in cert['reason'].lower()
    with pytest.raises(GeometryContractError):
        require_geometry_certificate(v, f, evaluate=True)


def test_clean_mesh_require_helper_accepts_when_evaluated():
    v, f = box_vertices_faces(3.0)
    cert = require_geometry_certificate(v, f, evaluate=True)
    assert cert['accepted'] is True
    assert cert['self_intersections'] is False


if __name__ == '__main__':
    import inspect

    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    ran = 0
    for fn in fns:
        if inspect.signature(fn).parameters:
            print('SKIP', fn.__name__, '(requires a pytest runner)')
            continue
        fn()
        ran += 1
        print('PASS', fn.__name__)
    print(f'{ran} tests passed')
