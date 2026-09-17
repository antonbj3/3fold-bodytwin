#!/usr/bin/env python3
"""Independent cross-checks for the geometry modules (region_mass_v1 / material_roster_v1 /
i2_payload_validate_v1).

The point of this file is that the module values are compared against *independent*
sources so a module bug cannot hide behind self-consistency.  Independent sources:

  * the audited evidence CSV, parsed here with the stdlib ``csv`` module (roster values);
  * a closed-form analytic box inertia and a hand-written numpy signed-tetra volume
    integral (moments / volume);
  * trimesh's own ``moment_inertia`` (a separate implementation) scaled by the density,
    with a fixed-seed Monte-Carlo volumetric fallback if that public API is absent;
  * direct JSON key inspection plus mutation probes (the payload validator must be
    falsifiable, not vacuously true);
  * hand-derived face counts on a 1x1x2 voxel grid (area parity, not prototype internals).

Fixed tolerances (chosen BEFORE running):
  * REL_TOL_VOLUME   = 1e-9   -- module volume vs numpy signed-tetra vs analytic 10*20*30;
  * REL_TOL_INERTIA  = 1e-9   -- module inertia vs closed-form analytic;
  * REL_TOL_TRIMESH  = 1e-9   -- module inertia vs trimesh ``moment_inertia`` * density;
  * ABS_TOL_INERTIA  = 1e-15  -- kg*mm2, only to absorb float noise on the ~0 off-diagonals;
  * REL_TOL_ROSTER   = 1e-12  -- audited CSV literal vs roster literal;
  * MC_SAMPLES = 200000, MC_SEED = 20260916, MC_REL_TOL = 5e-2 -- fallback only; in this
    environment trimesh 5.1.0 *does* expose ``moment_inertia`` so the MC path is not run.

All geometry/density values are synthetic or audited literals; no physiological claims.
"""
from __future__ import annotations

import copy
import csv
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import trimesh

REPO = Path(__file__).resolve().parents[1]
# The audited property registry is not part of this repository: it is located by
# BODYTWIN_REF_PROPERTY_REGISTRY_DIR (see tests/research_data.py).
PROPERTY_REGISTRY = 'property_registry.csv'
REGION_RESOLUTION = 'registry_regions_resolved.csv'
EXAMPLE_PAYLOAD = REPO / 'examples' / 'geometry' / 'region_example.json'


from bodytwin.geometry.label_interfaces_v1 import label_interfaces  # noqa: E402
from bodytwin.geometry.i2_payload_validate_v1 import validate_payload  # noqa: E402
from bodytwin.geometry.material_roster_v1 import MISSING, REGION_IDS, is_known, lookup  # noqa: E402
from bodytwin.geometry.region_mass_v1 import (  # noqa: E402
    MaterialProperty, area_mm2, closed_surface_volume_mm3, mass_moments_surface)

REL_TOL_VOLUME = 1e-9
REL_TOL_INERTIA = 1e-9
REL_TOL_TRIMESH = 1e-9
ABS_TOL_INERTIA = 1e-15
REL_TOL_ROSTER = 1e-12
MC_SAMPLES = 200000
MC_SEED = 20260916
MC_REL_TOL = 5e-2

# 1 L = 1e6 mm^3, so 1 kg/L = 1e-6 kg/mm^3 (computed here, not taken from the module).
LITRE_MM3 = 1_000_000.0
KG_PER_L_TO_KG_PER_MM3 = 1.0 / LITRE_MM3
BOX_DENSITY_KG_PER_L = 1.06  # audited muscle density literal (same row the roster cites)
BOX_EXTENTS_MM = (10.0, 20.0, 30.0)

# A material-density row starts with ``density_`` (density_muscle, density_default_generic,
# density_brain, density_heart, density_muscle_krogh) or ``apparent_density_`` (the audited
# ``apparent_density_bone`` row that the task explicitly requires for MSK-BONE).  This
# deliberately excludes ``energy_density_ratio`` (a dimensionless exercise ratio) and
# ``air_density_generic`` (air, not a tissue material).  ``specific_heat_body`` is the only
# audit row whose property kind is specific heat.
_DENSITY_PATTERN = re.compile(r'^(density_|apparent_density_)')
_MASS_BASES = ('region_volume', 'active_muscle', 'total_muscle', 'total_body')


def _canonical_units(units):
    """Fold cosmetic unit spellings ('g/cm^3' == 'g/cm3') for comparison only."""
    return str(units).strip().lower().replace('^', '').replace('**', '').replace(' ', '')


def _load_csv(name):
    from research_data import optional_path
    path = optional_path('BODYTWIN_REF_PROPERTY_REGISTRY_DIR') / name
    if not path.is_file():
        import pytest
        pytest.fail(f'{path} missing under BODYTWIN_REF_PROPERTY_REGISTRY_DIR')
    with open(path, newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def _signed_tetra_volume_mm3(vertices, faces):
    """Volume by the divergence/signed-tetrahedron sum, written from scratch with numpy."""
    tri = np.asarray(vertices, dtype=np.float64)[np.asarray(faces, dtype=np.int64)]
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    return float(np.einsum('ij,ij->i', a, np.cross(b, c)).sum() / 6.0)


def _box_10_20_30():
    mesh = trimesh.creation.box(extents=BOX_EXTENTS_MM)
    return (np.asarray(mesh.vertices, dtype=np.float64),
            np.asarray(mesh.faces, dtype=np.int64), mesh)


def _analytic_box_inertia_kg_mm2(mass_kg, a, b, c):
    """Inertia of a solid a x b x c box about its centroid."""
    return np.diag((
        mass_kg / 12.0 * (b * b + c * c),
        mass_kg / 12.0 * (a * a + c * c),
        mass_kg / 12.0 * (a * a + b * b),
    ))


def _monte_carlo_box_inertia(vertices, faces, rho, *, samples=MC_SAMPLES, seed=MC_SEED):
    """Fixed-seed uniform-volume Monte-Carlo inertia (fallback only; not used here)."""
    mesh = trimesh.Trimesh(vertices, faces, process=False)
    lo, hi = mesh.bounds
    rng = np.random.default_rng(seed)
    pts = rng.uniform(lo, hi, size=(samples, 3))
    inside = pts[mesh.contains(pts)]
    assert len(inside) > 0, 'Monte-Carlo found no interior samples'
    q = inside - inside.mean(axis=0)
    cov = (q.T @ q) / len(q)
    unit_inertia = np.trace(cov) * np.eye(3) - cov
    return (rho * _signed_tetra_volume_mm3(vertices, faces)) * unit_inertia


def test_roster_against_audited_csv():
    """material_roster_v1 values/units must equal the audited property_registry.csv rows."""
    rows = _load_csv(PROPERTY_REGISTRY)
    selected = [r for r in rows
                if r['property'] == 'specific_heat_body'
                or _DENSITY_PATTERN.search(r['property'])]
    assert selected, 'audited CSV holds no density_*/specific_heat_body rows'

    # UNMAPPED rows are resolved through the audited region-resolution table (independent data).
    resolution = {r['property']: r['region_id_new'] for r in _load_csv(REGION_RESOLUTION)}
    kind_for = lambda prop: 'specific_heat' if prop == 'specific_heat_body' else 'density'

    resolved_regions = set()
    for row in selected:
        prop = row['property']
        kind = kind_for(prop)
        csv_region = row['region_id']
        region = csv_region if csv_region in REGION_IDS else resolution.get(prop)
        csv_value = float(row['value_num'])
        csv_units = row['units']

        # No silent invention: the property name must not be attached to a roster region
        # other than the one the CSV/resolution attributes it to.
        name_regions = {rid for rid in REGION_IDS
                        if is_known(lookup(rid, kind)) and lookup(rid, kind).name == prop}
        assert not name_regions or name_regions == {region}, (
            f'{prop!r} is attached to roster regions {sorted(name_regions)} '
            f'but the audit resolves it to {region!r}')

        if region is None:
            for rid in REGION_IDS:
                node = lookup(rid, kind)
                if is_known(node):
                    assert not (node.value == csv_value and
                                _canonical_units(node.units) == _canonical_units(csv_units)), (
                        f'unmapped CSV row {prop!r} silently reappears under {rid!r}')
            continue

        node = lookup(region, kind)
        assert is_known(node), f'{region}/{kind} is not KNOWN although the CSV has a row'
        assert math.isclose(node.value, csv_value, rel_tol=REL_TOL_ROSTER,
                            abs_tol=REL_TOL_ROSTER), (
            f'{region}/{kind}: roster {node.value!r} != audited CSV {csv_value!r}')
        assert _canonical_units(node.units) == _canonical_units(csv_units), (
            f'{region}/{kind}: roster units {node.units!r} != audited CSV {csv_units!r}')
        resolved_regions.add(region)

    required = {'MSK-MUSCLE', 'NEU-TISSUE', 'CARD-VESSEL', 'MSK-BONE', 'WHOLE-BODY'}
    missing = required - resolved_regions
    assert not missing, f'required audited mappings never exercised: {sorted(missing)}'


def test_roster_missing_is_never_silently_defaulted():
    """A documented gap stays MISSING (no number) and an unknown region is refused."""
    tendon = lookup('MSK-TENDON', 'density')
    assert isinstance(tendon, dict) and tendon.get('status') == MISSING
    assert not is_known(tendon)
    assert 'value' not in tendon, 'a MISSING record must carry no numeric value'
    try:
        lookup('NOT-A-REGION', 'density')
    except KeyError:
        pass
    else:
        raise AssertionError('an unregistered region was silently accepted')


def test_inertia_against_independent_implementations():
    """Module inertia vs (a) closed form and (b) trimesh's own moment_inertia."""
    vertices, faces, mesh = _box_10_20_30()
    rho = BOX_DENSITY_KG_PER_L * KG_PER_L_TO_KG_PER_MM3  # conversion computed in this file
    density = MaterialProperty('density_muscle', BOX_DENSITY_KG_PER_L, 'kg/L', 'literal_cited')
    moments = mass_moments_surface(vertices, faces, density=density,
                                   region_id='SYNTH-BOX', basis='region_volume')
    assert moments.known and moments.inertia_kg_mm2 is not None

    volume = _signed_tetra_volume_mm3(vertices, faces)  # independent volume
    expected_mass = rho * volume

    assert math.isclose(moments.mass_kg, expected_mass, rel_tol=REL_TOL_INERTIA, abs_tol=0.0), (
        f'mass {moments.mass_kg!r} != independent rho*V {expected_mass!r}')

    analytic = _analytic_box_inertia_kg_mm2(expected_mass, *BOX_EXTENTS_MM)
    assert np.allclose(moments.inertia_kg_mm2, analytic, rtol=REL_TOL_INERTIA,
                       atol=ABS_TOL_INERTIA), (
        f'candidate inertia\n{moments.inertia_kg_mm2}\n!= analytic\n{analytic}')

    trimesh_mi = getattr(mesh, 'moment_inertia', None)
    if trimesh_mi is not None:
        # trimesh's moment_inertia is for unit density, so scale by rho (homogeneous solid).
        scaled = np.asarray(trimesh_mi, dtype=np.float64) * rho
        assert np.allclose(moments.inertia_kg_mm2, scaled, rtol=REL_TOL_TRIMESH,
                           atol=ABS_TOL_INERTIA), (
            f'candidate inertia\n{moments.inertia_kg_mm2}\n'
            f'!= trimesh.moment_inertia * rho\n{scaled}')
    else:  # pragma: no cover - trimesh 5.1.0 exposes the API, so this is not run here
        mc = _monte_carlo_box_inertia(vertices, faces, rho)
        rel = np.linalg.norm(moments.inertia_kg_mm2 - mc) / np.linalg.norm(mc)
        assert rel < MC_REL_TOL, f'Monte-Carlo inertia mismatch rel={rel:.3g}'


def test_volume_without_candidate():
    """Hand-written numpy volume must match the module and the analytic 6000 mm^3."""
    vertices, faces, _ = _box_10_20_30()
    numpy_volume = _signed_tetra_volume_mm3(vertices, faces)
    candidate_volume = closed_surface_volume_mm3(vertices, faces)
    analytic = BOX_EXTENTS_MM[0] * BOX_EXTENTS_MM[1] * BOX_EXTENTS_MM[2]
    assert math.isclose(numpy_volume, analytic, rel_tol=REL_TOL_VOLUME, abs_tol=0.0)
    assert math.isclose(candidate_volume, analytic, rel_tol=REL_TOL_VOLUME, abs_tol=0.0)
    assert math.isclose(candidate_volume, numpy_volume, rel_tol=REL_TOL_VOLUME, abs_tol=0.0)


def test_payload_validator_is_not_self_referential():
    """The example validates, really carries units/provenance/basis, and fails when mutated."""
    payload = json.loads(EXAMPLE_PAYLOAD.read_text(encoding='utf-8'))
    assert validate_payload(payload) == []

    regions = payload.get('regions')
    assert isinstance(regions, list) and regions
    for index, region in enumerate(regions):
        assert isinstance(region.get('mass_basis'), str) and region['mass_basis'], (
            f'regions[{index}].mass_basis missing/empty')
        assert region['mass_basis'] in _MASS_BASES, (
            f'regions[{index}].mass_basis {region["mass_basis"]!r} outside the schema set')
        for key, units in (('mass', 'kg'), ('heat_capacity', 'J/K')):
            quantity = region.get(key)
            assert isinstance(quantity, dict), f'regions[{index}].{key} is not an object'
            assert quantity.get('units') == units, (
                f'regions[{index}].{key}.units != {units!r}')
            assert isinstance(quantity.get('provenance'), str) and quantity['provenance'], (
                f'regions[{index}].{key}.provenance missing/empty')

    # Mutation probes: a validator that always returned [] would fail these.
    def errors_after(mutate):
        broken = copy.deepcopy(payload)
        mutate(broken)
        return validate_payload(broken)

    assert errors_after(lambda p: p['regions'][0]['mass'].pop('provenance'))
    assert errors_after(lambda p: p['regions'][0].pop('mass_basis'))
    assert errors_after(lambda p: p['regions'][0]['heat_capacity'].__setitem__('units', 'J'))
    assert errors_after(lambda p: p['regions'][0]['mass'].__setitem__('value', None))


def test_area_parity_on_voxel_grid():
    """1x1x2 grid of two unit voxels: exterior 10 mm^2, one shared material face 1 mm^2.

    By hand: voxel A (label 1) spans z in [0,1], voxel B (label 2) spans z in [1,2].  Each
    unit cube has surface 6 mm^2 -> 12 mm^2 total.  They share exactly one 1x1 mm face, which
    is removed from both cubes' exterior surface and counted once as material interface:
    exterior = 12 - 2*1 = 10, material_interface = 1.
    """
    labels = np.array([[[1, 2]]], dtype=np.uint8)
    interfaces = label_interfaces(labels)
    assert area_mm2(interfaces, 'material_interface_area') == 1.0
    assert area_mm2(interfaces, 'exterior_area') == 10.0


if __name__ == '__main__':
    functions = [name for name in sorted(globals())
                 if name.startswith('test_') and callable(globals()[name])]
    for name in functions:
        globals()[name]()
        print('PASS', name)
    print(f'{len(functions)} tests passed')
