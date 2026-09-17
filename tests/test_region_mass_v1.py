#!/usr/bin/env python3
"""Acceptance tests for region/material/mass/heat-capacity integration.

Run: OMP_NUM_THREADS=1 PYTHONPATH=src python tests/test_region_mass_v1.py
Tolerances were fixed before the run: exact polytopes 1e-9 relative, and every value is
compared against an independent analytic reference or the repo's own fixture.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import trimesh

REPO = Path(__file__).resolve().parents[1]

from bodytwin.geometry.region_mass_v1 import (  # noqa: E402
    DENSITY_REGISTRY, MaterialProperty, Region, UnitError, UndeclaredError,
    anatomical_cross_section_mm2, area_mm2, closed_surface_volume_mm3,
    functional_exchange_area_mm2, integrate_region, lookup_density,
    physiological_cross_section_mm2, region_boundary_volume_mm3, tetra_volume_mm3,
    to_kg_per_mm3, to_mm3, total_quantity, voxel_region_volume_mm3)
from bodytwin.geometry.region_mass_v1 import (  # noqa: E402
    NonConformingError, assert_conforming_interfaces, assert_unique_regions,
    coincident_vertex_groups, interface_area_by_pair_mm2, mass_fraction,
    mass_moments_surface, mass_moments_tetra, material_interface_area_checked)
from bodytwin.geometry.label_interfaces_v1 import label_interfaces  # noqa: E402
from bodytwin.geometry.layered_box_v1 import layered_box  # noqa: E402
from bodytwin.geometry.mesh_ingest_v1 import surface_mesh, tetrahedral_mesh  # noqa: E402
from bodytwin.geometry.tetra_interfaces_v1 import tetra_interfaces  # noqa: E402

RTOL = 1e-9  # fixed before the run for exact polytopes (float64 / analytic)



def _load_c101():
    """Earlier area-descriptor prototype for the area parity check (BODYTWIN_REF_AREA_DESCRIPTOR;
    skipped when unset, failed when set but missing)."""
    from research_data import optional_path
    path = optional_path('BODYTWIN_REF_AREA_DESCRIPTOR')
    spec = importlib.util.spec_from_file_location('geometry_descriptor_v1', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def box(side):
    m = trimesh.creation.box(extents=[side, side, side])
    return np.asarray(m.vertices), np.asarray(m.faces)


def unit_labels():
    labels = np.zeros((1, 1, 2), np.uint8)
    labels[0, 0, 0] = 1
    labels[0, 0, 1] = 2
    return labels


def cube6():
    v = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    t = np.array([[0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6],
                  [0, 7, 4, 6], [0, 4, 5, 6], [0, 5, 1, 6]], np.int64)
    return v, t


def conforming_two_cubes():
    v = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    t = np.array([[0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6],
                  [0, 7, 4, 6], [0, 4, 5, 6], [0, 5, 1, 6]], np.int64)
    V = np.vstack([v, np.array([[2, 0, 0], [2, 1, 0], [2, 0, 1], [2, 1, 1]], float)])
    remap = {0: 1, 1: 8, 2: 9, 3: 2, 4: 5, 5: 10, 6: 11, 7: 6}
    T = np.vstack([t, np.vectorize(remap.get)(t)])
    L = np.array([1] * 6 + [2] * 6, np.int32)
    return V, T, L


def test_analytic_box_volume_and_mass():
    v, f = box(10.0)
    assert abs(closed_surface_volume_mm3(v, f) - 1000.0) < RTOL * 1000.0
    region = Region('MSK-MUSCLE', 1000.0, 'mm3', 'analytic box 10 mm cube')
    state = integrate_region(region, density_key='MSK-MUSCLE')
    expect_mass = 1.06e-6 * 1000.0
    assert state.mass_kg.known and state.mass_kg.units == 'kg'
    assert abs(state.mass_kg.value - expect_mass) < RTOL * expect_mass


def test_analytic_tetra_volume_and_mass():
    v = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    t = np.array([[0, 1, 2, 3]], np.int64)
    assert abs(tetra_volume_mm3(v, t) - 1.0 / 6.0) < RTOL / 6.0
    surface = tetrahedral_mesh(v, t, units='mm')
    assert abs(surface.volume_mm3 - 1.0 / 6.0) < RTOL / 6.0
    region = Region('MSK-MUSCLE', surface.volume_mm3, 'mm3', 'analytic unit tetra')
    state = integrate_region(region, density_key='MSK-MUSCLE')
    expect = 1.06e-6 / 6.0
    assert abs(state.mass_kg.value - expect) < RTOL * expect


def test_mass_scales_as_length_cubed():
    ratios = []
    for side in (0.5, 2.0, 10.0):
        region = Region('MSK-MUSCLE', side ** 3, 'mm3', 'analytic cube')
        mass = integrate_region(region, density_key='MSK-MUSCLE').mass_kg.value
        ratios.append(mass / side ** 3)
    spread = max(ratios) - min(ratios)
    assert spread < RTOL * max(ratios)
    assert abs(ratios[0] - 1.06e-6) < RTOL * 1.06e-6


def test_mm_m_errors_rejected():
    v, f = box(1.0)
    try:
        surface_mesh(v, f, units='m')
        raise AssertionError('metre-tagged mesh accepted')
    except ValueError:
        pass
    for bad in ('metre', 'Mm', 'inch'):
        try:
            to_mm3(1.0, bad)
            raise AssertionError(f'volume unit {bad!r} accepted')
        except UnitError:
            pass
    try:
        to_kg_per_mm3(1.06, 'lb/ft3')
        raise AssertionError('density unit accepted')
    except UnitError:
        pass
    assert abs(to_mm3(1.0, 'm3') - 1.0e9) < 1e-3
    same_number_mm3 = integrate_region(Region('MSK-MUSCLE', 1.0, 'mm3', 'x'),
                                       density_key='MSK-MUSCLE').mass_kg.value
    same_number_m3 = integrate_region(Region('MSK-MUSCLE', 1.0, 'm3', 'x'),
                                      density_key='MSK-MUSCLE').mass_kg.value
    assert abs(same_number_m3 / same_number_mm3 - 1.0e9) < 1.0
    try:
        Region('MSK-MUSCLE', 1.0, 'mm', 'x')
        raise AssertionError('length unit accepted as volume unit')
    except UnitError:
        pass


def test_interfaces_counted_once():
    r = label_interfaces(unit_labels())
    material = area_mm2(r, 'material_interface_area')
    exterior = area_mm2(r, 'exterior_area')
    assert abs(material - 1.0) < RTOL
    assert abs(exterior - 10.0) < RTOL
    assert abs(material + exterior - 11.0) < RTOL

    b = layered_box()
    interface = area_mm2(b, 'material_interface_area')
    assert abs(interface - 3600.0) < RTOL
    v1 = region_boundary_volume_mm3(b, 1)
    v2 = region_boundary_volume_mm3(b, 2)
    assert abs(v1 - 36000.0) < RTOL * 36000.0
    assert abs(v2 - 180000.0) < RTOL * 180000.0
    assert abs((v1 + v2) - 60.0 ** 3) < RTOL * 60.0 ** 3
    exterior = area_mm2(b, 'exterior_area')
    assert abs(exterior - 6 * 3600.0) < RTOL * 6 * 3600.0
    assert abs(exterior + interface - 7 * 3600.0) < RTOL * 7 * 3600.0


def test_missing_density_is_unknown_not_default():
    region = Region('MSK-TENDON', 1000.0, 'mm3', 'analytic box')
    assert lookup_density('MSK-TENDON') is None
    unbound = integrate_region(region)
    assert not unbound.mass_kg.known and 'density' in unbound.mass_kg.missing
    assert 'GENERIC' not in DENSITY_REGISTRY['MSK-MUSCLE'].name
    generic = integrate_region(region, density_key='GENERIC')
    assert generic.mass_kg.known and generic.mass_kg.provenance == 'assumption'
    assert generic.density.name == 'density_default_generic'


def test_heat_capacity_reported_in_j_per_k():
    cp = MaterialProperty('muscle_cp_assumption', 3490.0, 'J/(kg*K)', 'assumption',
                          'thermoregulation lumped value, regime mismatch flagged')
    region = Region('MSK-MUSCLE', 1000.0, 'mm3', 'analytic box')
    state = integrate_region(region, density_key='MSK-MUSCLE', specific_heat=cp)
    expect = 1.06e-3 * 3490.0
    assert state.heat_capacity_j_per_k.units == 'J/K'
    assert abs(state.heat_capacity_j_per_k.value - expect) < RTOL * expect
    assert state.heat_capacity_j_per_k.provenance == 'assumption'
    without_cp = integrate_region(region, density_key='MSK-MUSCLE')
    assert not without_cp.heat_capacity_j_per_k.known
    assert 'specific_heat' in without_cp.heat_capacity_j_per_k.missing
    cp_grams = MaterialProperty('cp_g', 3.49, 'J/(g*K)', 'assumption', 'equivalent')
    assert abs(cp.as_cp_J_per_kgK().value - cp_grams.as_cp_J_per_kgK().value) < RTOL * 3490.0


def test_cross_section_kinds_are_distinct():
    anatomical = anatomical_cross_section_mm2(1000.0, 10.0)
    physiological = physiological_cross_section_mm2(1000.0, 10.0, 30.0, angle_units='deg')
    assert abs(anatomical - 100.0) < RTOL * 100.0
    assert abs(physiological - 100.0 * np.cos(np.deg2rad(30.0))) < RTOL * 100.0
    assert abs(anatomical - physiological) > 1.0
    try:
        physiological_cross_section_mm2(1000.0, 10.0, 30.0, angle_units='turns')
        raise AssertionError('undeclared angle unit accepted')
    except UnitError:
        pass


def test_functional_exchange_area_requires_declaration():
    b = layered_box()
    try:
        area_mm2(b, 'functional_exchange_area')
        raise AssertionError('wall area silently returned as exchange area')
    except UndeclaredError:
        pass
    try:
        functional_exchange_area_mm2(b)
        raise AssertionError('undeclared exchange pairs accepted')
    except UndeclaredError:
        pass
    declared = functional_exchange_area_mm2(b, declared_pairs=[{1, 2}])
    assert abs(declared - 3600.0) < RTOL


def test_parity_with_lane_c_prototype():
    c101 = _load_c101()
    r = label_interfaces(unit_labels())
    for kind in ('exterior_area', 'material_interface_area'):
        assert abs(area_mm2(r, kind) - c101.interface_area_mm2(r, kind=kind)) < 1e-12
        assert abs(area_mm2(r, kind) - c101.interface_area_mm2(r, kind=kind)) < 1e-12
    b = layered_box()
    assert abs(area_mm2(b, 'material_interface_area')
               - c101.interface_area_mm2(b, kind='material_interface_area')) < 1e-9
    V, T, L = conforming_two_cubes()
    t = tetra_interfaces(V, T, L)
    assert abs(area_mm2(t, 'material_interface_area') - 1.0) < 1e-12
    assert abs(area_mm2(t, 'material_interface_area')
               - c101.interface_area_mm2(t, kind='material_interface_area')) < 1e-12


def test_total_partition_conservation():
    cp = MaterialProperty('muscle_cp_assumption', 3490.0, 'J/(kg*K)', 'assumption', 'declared')
    b = layered_box()
    states = [integrate_region(Region('MSK-MUSCLE', region_boundary_volume_mm3(b, region),
                                      'mm3', 'outward boundary', area_kind='material_interface_area',
                                      area_mm2=area_mm2(b, 'material_interface_area')),
                               density_key='MSK-MUSCLE', specific_heat=cp)
              for region in (1, 2)]
    total_mass = total_quantity(states, 'mass_kg')
    total_heat = total_quantity(states, 'heat_capacity_j_per_k')
    assert total_mass.known and abs(total_mass.value - 1.06e-6 * 216000.0) < RTOL * 1.06e-3 * 216
    assert total_heat.units == 'J/K'
    assert abs(total_heat.value - total_mass.value * 3490.0) < RTOL * total_heat.value
    unknown = [integrate_region(Region('MSK-TENDON', 10.0, 'mm3', 'x'))]
    assert not total_quantity(unknown, 'mass_kg').known
    assert total_quantity(unknown, 'mass_kg').missing == ('MSK-TENDON',)


def test_voxel_region_volume_and_pitch_contract():
    labels = unit_labels()
    assert abs(voxel_region_volume_mm3(labels, 1) - 1.0) < RTOL
    assert abs(voxel_region_volume_mm3(labels, 2) - 1.0) < RTOL
    try:
        voxel_region_volume_mm3(labels, 1, pitch_mm=2.0)
        raise AssertionError('undeclared pitch accepted')
    except UnitError:
        pass
    assert abs(voxel_region_volume_mm3(labels, 1, pitch_mm=2.0, pitch_declared=True) - 8.0) < RTOL


def nonconforming_two_cubes():
    """Two 1 mm cubes touching at z=1 with coincident vertices of distinct indices."""
    low = trimesh.creation.box(extents=[1, 1, 1])
    low.apply_translation([0, 0, 0.5])
    high = trimesh.creation.box(extents=[1, 1, 1])
    high.apply_translation([0, 0, 1.5])
    vertices = np.vstack([low.vertices, high.vertices])
    faces = np.vstack([low.faces, high.faces + len(low.vertices)]).astype(np.int32)
    front = np.zeros(len(faces), np.int32)
    back = np.concatenate([np.ones(len(low.faces), np.int32),
                           np.full(len(high.faces), 2, np.int32)])
    return dict(vertices=vertices, faces=faces, front=front, back=back)


def test_conforming_guard_rejects_silent_exterior():
    bad = nonconforming_two_cubes()
    assert area_mm2(bad, 'material_interface_area') == 0.0  # the non-conforming silent failure
    groups, vertices = coincident_vertex_groups(bad['vertices'])
    assert groups == 4 and vertices == 8
    try:
        material_interface_area_checked(bad)
        raise AssertionError('non-conforming assembly accepted')
    except NonConformingError:
        pass
    for good in (layered_box(), label_interfaces(unit_labels())):
        assert assert_conforming_interfaces(good) is True
        assert material_interface_area_checked(good) >= 0.0
    V, T, L = conforming_two_cubes()
    assert assert_conforming_interfaces(tetra_interfaces(V, T, L)) is True


def test_mass_basis_separation_and_fraction():
    region = Region('MSK-MUSCLE', 1000.0, 'mm3', 'analytic box')
    active = integrate_region(region, density_key='MSK-MUSCLE', basis='active_muscle')
    body = integrate_region(region, density_key='MSK-MUSCLE', basis='total_body')
    assert active.basis == 'active_muscle' and body.basis == 'total_body'
    try:
        integrate_region(region, density_key='MSK-MUSCLE', basis='muscle')
        raise AssertionError('undeclared basis accepted')
    except UnitError:
        pass
    fraction = mass_fraction(3.0, 12.0, numerator_basis='active_muscle',
                             denominator_basis='total_body')
    assert fraction.units == '1' and abs(fraction.value - 0.25) < 1e-12
    for args in ({'numerator_basis': 'active_muscle', 'denominator_basis': 'active_muscle'},
                 {'numerator_basis': 'total_muscle', 'denominator_basis': 'total_muscle'}):
        try:
            mass_fraction(3.0, 12.0, **args)
            raise AssertionError('identical bases accepted as a fraction')
        except UnitError:
            pass
    try:
        mass_fraction(13.0, 12.0, numerator_basis='active_muscle', denominator_basis='total_body')
        raise AssertionError('inconsistent bases accepted')
    except UnitError:
        pass
    try:
        integrate = [active, body]
        total_quantity(integrate, 'mass_kg')
        raise AssertionError('mixed bases summed')
    except UnitError:
        pass


def test_mass_invariant_under_rigid_transform():
    v, f = box(10.0)
    base = closed_surface_volume_mm3(v, f)
    angle = 0.7
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0.0],
                         [np.sin(angle), np.cos(angle), 0.0],
                         [0.0, 0.0, 1.0]])
    moved = v @ rotation.T + np.array([1000.0, -250.0, 37.0])
    assert abs(closed_surface_volume_mm3(moved, f) - base) < RTOL * base
    m_base = integrate_region(Region('MSK-MUSCLE', base, 'mm3', 'x'),
                              density_key='MSK-MUSCLE').mass_kg.value
    m_moved = integrate_region(Region('MSK-MUSCLE', closed_surface_volume_mm3(moved, f), 'mm3', 'x'),
                               density_key='MSK-MUSCLE').mass_kg.value
    assert abs(m_base - m_moved) < RTOL * m_base


def test_mass_moments_box_analytic():
    rho = MaterialProperty('muscle_density', 1.06, 'kg/L', 'literal_cited')
    extents = (10.0, 20.0, 30.0)
    m = trimesh.creation.box(extents=list(extents))
    moments = mass_moments_surface(m.vertices, m.faces, density=rho)
    mass = 1.06e-6 * 6000.0
    assert abs(moments.mass_kg - mass) < RTOL * mass
    assert np.linalg.norm(moments.centroid_mm) < 1e-9
    a, b, c = extents
    expected = mass / 12.0 * np.array([b ** 2 + c ** 2, a ** 2 + c ** 2, a ** 2 + b ** 2])
    assert np.allclose(np.diag(moments.inertia_kg_mm2), expected, rtol=RTOL, atol=0)
    off = moments.inertia_kg_mm2 - np.diag(np.diag(moments.inertia_kg_mm2))
    assert np.abs(off).max() < 1e-12 * mass * c ** 2


def test_mass_moments_surface_and_tetra_agree():
    rho = MaterialProperty('muscle_density', 1.06, 'kg/L', 'literal_cited')
    box = trimesh.creation.box(extents=[1, 1, 1])
    box.apply_translation([0.5, 0.5, 0.5])
    surface = mass_moments_surface(box.vertices, box.faces, density=rho)
    v, t = cube6()
    tetra = mass_moments_tetra(v, t, density=rho)
    assert abs(surface.mass_kg - tetra.mass_kg) < RTOL * surface.mass_kg
    assert np.allclose(surface.centroid_mm, tetra.centroid_mm, rtol=RTOL, atol=1e-12)
    assert np.allclose(surface.inertia_kg_mm2, tetra.inertia_kg_mm2, rtol=1e-9, atol=1e-15)
    assert np.allclose(np.diag(tetra.inertia_kg_mm2),
                       1.06e-6 / 6.0 * np.ones(3), rtol=RTOL, atol=0)


def test_mass_moments_missing_density_and_units():
    box = trimesh.creation.box(extents=[1, 1, 1])
    unknown = mass_moments_surface(box.vertices, box.faces)
    assert not unknown.known and unknown.missing == ('density',)
    assert 'kg*m2' in unknown.inertia_SI_units
    rho = MaterialProperty('m', 1.06, 'kg/L', 'literal_cited')
    try:
        mass_moments_surface(box.vertices, box.faces, density=rho, basis='muscle')
        raise AssertionError('undeclared basis accepted')
    except UnitError:
        pass
    try:
        mass_moments_tetra(np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 0]]),
                           np.array([[0, 1, 2, 3]]), density=rho)
        raise AssertionError('degenerate tetra accepted')
    except ValueError:
        pass


def test_mass_moments_rotation_covariance():
    rho = MaterialProperty('m', 1.06, 'kg/L', 'literal_cited')
    m = trimesh.creation.box(extents=[10.0, 20.0, 30.0])
    base = mass_moments_surface(m.vertices, m.faces, density=rho)
    angle = 0.6
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0.0],
                         [np.sin(angle), np.cos(angle), 0.0], [0.0, 0.0, 1.0]])
    moved = np.asarray(m.vertices) @ rotation.T
    rotated = mass_moments_surface(moved, np.asarray(m.faces), density=rho)
    assert np.allclose(rotated.inertia_kg_mm2, rotation @ base.inertia_kg_mm2 @ rotation.T,
                       rtol=1e-9, atol=1e-12)
    assert np.allclose(np.sort(np.linalg.eigvalsh(rotated.inertia_kg_mm2)),
                       np.sort(np.linalg.eigvalsh(base.inertia_kg_mm2)), rtol=1e-9, atol=1e-12)
    assert np.linalg.norm(rotated.centroid_mm) < 1e-9


def test_interface_pairs_and_unique_regions():
    b = layered_box()
    pairs = interface_area_by_pair_mm2(b)
    assert set(pairs) == {(1, 2)}
    assert abs(pairs[(1, 2)] - 3600.0) < RTOL * 3600.0
    assert abs(sum(pairs.values()) - area_mm2(b, 'material_interface_area')) < 1e-9
    grid = label_interfaces(unit_labels())
    assert interface_area_by_pair_mm2(grid) == {(1, 2): 1.0}
    states = [integrate_region(Region('MSK-MUSCLE', 1.0, 'mm3', 'x'), density_key='MSK-MUSCLE')]
    assert assert_unique_regions(states) is True
    try:
        assert_unique_regions(states + states)
        raise AssertionError('duplicate region ids accepted')
    except ValueError:
        pass


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for fn in fns:
        fn()
        print('PASS', fn.__name__)
    print(f'{len(fns)} tests passed')
