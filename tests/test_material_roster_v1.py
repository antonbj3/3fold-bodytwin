#!/usr/bin/env python3
"""Acceptance tests for the material roster.

Run:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        PYTHONPATH=src python tests/test_material_roster_v1.py
    PYTHONPATH=src python -m pytest -q tests/test_material_roster_v1.py

Tolerances fixed before the run: the audited values are literals copied from
property_registry.csv, not computed quantities, so KNOWN value/units/provenance/sha256 are
compared exactly (==). No integration or floating-point tolerance is involved.

The audited data (property registry, unit and region resolution tables, material audit) is
not part of this repository.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

from bodytwin.geometry.material_roster_v1 import (  # noqa: E402
    ABSENT, KNOWN, MISSING, PROPERTY_NAMES, REGION_IDS,
    assert_registered, coverage_report, is_absent, is_known, is_missing, lookup, roster_table)
from bodytwin.geometry.region_mass_v1 import MaterialProperty  # noqa: E402

MUSCLE_DENSITY_SHA = 'd08a611c6549d3de86304e6b4fa74525471d7f44650af5ec8a82c0f574e0ed0c'


def test_known_density_muscle():
    p = lookup('MSK-MUSCLE', 'density')
    assert is_known(p) and isinstance(p, MaterialProperty)
    assert p.value == 1.06
    assert p.units == 'kg/L'
    assert p.provenance == 'literal_cited'
    assert p.source == 'src/bodytwin/cells/energy/muscle_energetics.py::MUSCLE_DENSITY_KG_PER_L'
    assert p.source_sha256 == MUSCLE_DENSITY_SHA


def test_known_specific_heat_whole_body():
    p = lookup('WHOLE-BODY', 'specific_heat')
    assert is_known(p) and isinstance(p, MaterialProperty)
    assert p.value == 3490.0
    assert p.units == 'J/(kg*K)'
    assert p.provenance == 'literal_cited'
    assert p.source == 'src/bodytwin/cells/organ_systems/thermoregulation.py::SPECIFIC_HEAT_BODY_J_PER_KG_K'


def test_absent_is_not_zero_and_not_missing():
    for region_id in REGION_IDS:
        for prop in ('thermal_conductivity', 'perfusion_heat_transport'):
            rec = lookup(region_id, prop)
            assert not is_known(rec)
            assert not isinstance(rec, MaterialProperty)
            assert rec is not None
            assert rec != 0          # never a silent zero
            assert rec != 0.0
            assert rec != MISSING    # distinct from MISSING
            assert rec != ABSENT     # it is a record, not the bare sentinel
            assert rec['status'] == ABSENT
            assert is_absent(rec) and not is_missing(rec)
            assert isinstance(rec['reason'], str) and rec['reason'].strip()
            assert 'audit' in rec['reason'].lower()
            assert rec.get('source')
            assert rec.get('source_sha256')


def test_missing_tendon_density_is_unknown_not_default():
    rec = lookup('MSK-TENDON', 'density')
    assert not is_known(rec) and not isinstance(rec, MaterialProperty)
    assert rec is not None
    assert rec != 0          # never a silent zero
    assert rec != 0.0
    assert rec['status'] == MISSING
    assert rec['status'] != ABSENT
    assert is_missing(rec) and not is_absent(rec)
    assert rec['field'] == 'density'
    assert isinstance(rec['reason'], str) and rec['reason'].strip()


def test_missing_is_distinct_from_absent():
    # Per-region gap (MISSING) vs codebase-wide absence (ABSENT) must not be conflated.
    assert lookup('WHOLE-BODY', 'density')['status'] == MISSING
    assert lookup('MSK-MUSCLE', 'specific_heat')['status'] == MISSING
    assert lookup('MSK-MUSCLE', 'thermal_conductivity')['status'] == ABSENT
    assert lookup('MSK-MUSCLE', 'perfusion_heat_transport')['status'] == ABSENT


def test_generic_density_only_under_generic_and_keeps_assumption():
    generic = lookup('GENERIC', 'density')
    assert is_known(generic) and isinstance(generic, MaterialProperty)
    assert generic.name == 'density_default_generic'
    assert generic.value == 1.05
    assert generic.units == 'kg/L'
    assert generic.provenance == 'assumption'  # audited default, not a measurement
    assert lookup('MSK-MUSCLE', 'density').provenance == 'literal_cited'
    for region_id in REGION_IDS:
        value = lookup(region_id, 'density')
        if region_id == 'GENERIC':
            assert isinstance(value, MaterialProperty)
            assert value.name == 'density_default_generic'
        elif isinstance(value, MaterialProperty):
            assert value.name != 'density_default_generic'  # never leak the generic row


def test_assert_registered():
    assert assert_registered('MSK-MUSCLE') is True
    assert assert_registered('MSK-TENDON') is True
    try:
        assert_registered('NOT-A-REGION')
        raise AssertionError('an unregistered region was accepted')
    except KeyError as exc:
        assert 'NOT-A-REGION' in str(exc)


def test_roster_table_and_coverage_report_consistent():
    rows = roster_table()
    assert len(rows) == len(REGION_IDS) * len(PROPERTY_NAMES)
    keys = [(r['region_id'], r['property_name']) for r in rows]
    assert len(set(keys)) == len(keys)  # exactly one row per (region, property)

    known_rows = 0
    for row in rows:
        value = lookup(row['region_id'], row['property_name'])
        if isinstance(value, MaterialProperty):
            assert row['status'] == KNOWN
            assert row['value'] == value.value
            assert row['units'] == value.units
            assert row['provenance'] == value.provenance
            assert row['source'] == value.source
            assert row['source_sha256'] == value.source_sha256
            assert row['reason'] is None
            known_rows += 1
        else:
            assert row['status'] == value['status']
            assert row['value'] is None and row['units'] is None
            assert row['provenance'] is None
            assert isinstance(row['reason'], str) and row['reason'].strip()  # nonempty reason
        if row['status'] == ABSENT:
            assert row['source'] and row['source_sha256']

    cov = coverage_report()
    assert set(cov) == set(PROPERTY_NAMES)
    total_known = 0
    for prop in PROPERTY_NAMES:
        counts = cov[prop]
        assert set(counts) == {KNOWN, ABSENT, MISSING}
        assert sum(counts.values()) == len(REGION_IDS)
        rows_for = [r for r in rows if r['property_name'] == prop]
        for status in (KNOWN, ABSENT, MISSING):
            assert counts[status] == sum(1 for r in rows_for if r['status'] == status)
        total_known += counts[KNOWN]
    assert total_known == known_rows  # KNOWN count equals rows with a MaterialProperty


def test_absent_and_missing_records_are_copies():
    first = lookup('MSK-MUSCLE', 'thermal_conductivity')
    first['reason'] = 'tampered'
    assert lookup('MSK-MUSCLE', 'thermal_conductivity')['reason'] != 'tampered'
    assert lookup('MSK-TENDON', 'density')['reason'] != 'tampered'


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for fn in fns:
        fn()
        print('PASS', fn.__name__)
    print(f'{len(fns)} tests passed')
