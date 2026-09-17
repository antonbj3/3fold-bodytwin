"""Auditable material roster: it distinguishes KNOWN / ABSENT / MISSING per (region, property) so
a consumer can never substitute a silent zero or default for an unmeasured property.

Additive only. It reuses the sibling ``region_mass_v1.MaterialProperty`` for unit handling, and
reuses that module's audited registry literals; it defines no unit-conversion table of its own. It changes no
existing symbol, performs no I/O at import and writes nothing.

Location: src/bodytwin/geometry/material_roster_v1.py

Status semantics -- mutually exclusive, and never conflated with a numeric zero or None:
  * KNOWN   -- an audited registry row supplies a numerical value with units, provenance, source
               and sha256. KNOWN does NOT imply measured: always read ``provenance`` (it may be
               the documented ``assumption`` for the GENERIC row, which is not a measurement).
  * ABSENT  -- a repository audit established the property does not exist in the codebase at all
               (no spatial thermal conduction and no perfusion heat transport), so no value can be
               read. The record carries a reason and source. Distinct from MISSING.
  * MISSING -- no audited value is registered for this particular (region, property); the gap is
               per-region, not a global absence. The record carries the missing field name. A
               consumer must measure or declare it explicitly; never default.

Source data (read, not invented): an audited property registry of the repository (values,
units and region resolution) and a repository material audit (no conduction or perfusion
term; no tendon or ligament density). Neither document is part of this repository.
"""
from __future__ import annotations

from .region_mass_v1 import DENSITY_REGISTRY, SPECIFIC_HEAT_REGISTRY, MaterialProperty

__all__ = [
    'ABSENT', 'KNOWN', 'MISSING', 'PROPERTY_NAMES', 'REGION_IDS', 'STATUSES',
    'assert_registered', 'coverage_report', 'is_absent', 'is_known', 'is_missing',
    'lookup', 'roster_table',
]

KNOWN = 'KNOWN'
ABSENT = 'ABSENT'
MISSING = 'MISSING'
STATUSES = (KNOWN, ABSENT, MISSING)

PROPERTY_NAMES = ('density', 'specific_heat', 'thermal_conductivity', 'perfusion_heat_transport')

# Reuse the audited literals already carried by region_mass_v1 (single source of truth: values,
# units, provenance, source and sha256 all come from property_registry.csv there). The roster
# adds only the status layer and coverage on top.
_DENSITY = dict(DENSITY_REGISTRY)
_SPECIFIC_HEAT = dict(SPECIFIC_HEAT_REGISTRY)

# The only region with no density value anywhere in either repo; kept MISSING rather than ABSENT
# because the density concept exists for other regions and the generic row must never leak here.
_NO_TENDON_DENSITY = 'MSK-TENDON'

REGION_IDS = tuple(sorted(set(_DENSITY) | set(_SPECIFIC_HEAT) | {_NO_TENDON_DENSITY}))

_ABSENT_SOURCE_SHA256 = '4474bfbf7b2bec7209d2283b658ad573a88ef00d94843e7c47870f7d1387dbc8'

_THERMAL_CONDUCTIVITY_ABSENCE = {
    'status': ABSENT,
    'field': 'thermal_conductivity',
    'reason': (
        "Repository material audit (heat from muscle/metabolism): thermoregulation.py is lumped heat "
        "accounting with no spatial conduction; the audited property registry row "
        "'thermal_conductivity_k_ABSENT' records no thermal-conductivity constant in source. Absent in "
        "the audited codebase, not a measured zero."
    ),
    'source': (
        'src/bodytwin/cells/organ_systems/thermoregulation.py (module docstring: lumped heat '
        'accounting, no spatial conduction)'
    ),
    'source_sha256': _ABSENT_SOURCE_SHA256,
}

_PERFUSION_ABSENCE = {
    'status': ABSENT,
    'field': 'perfusion_heat_transport',
    'reason': (
        "Repository material audit (heat from muscle/metabolism): thermoregulation.py has no "
        "perfusion/vascular heat transport; the audited property registry row "
        "'perfusion_heat_transport_ABSENT' records no perfusion heat transport term. Absent in the "
        "audited codebase, not a measured zero."
    ),
    'source': (
        'src/bodytwin/cells/organ_systems/thermoregulation.py (module docstring: no '
        'vascular/perfusion heat transport)'
    ),
    'source_sha256': _ABSENT_SOURCE_SHA256,
}

_GENERIC_MISSING_REASON = (
    'no audited value registered for this (region, property) in the audited property '
    'registry of the repository; measure or declare it explicitly, never '
    'substitute a default or zero.'
)

_TENDON_DENSITY_MISSING_REASON = (
    "no tendon or ligament density with provenance exists in the repository (material audit); "
    "density for MSK-TENDON is unknown, never a generic fallback."
)


def _absent_record(property_name):
    base = (_THERMAL_CONDUCTIVITY_ABSENCE if property_name == 'thermal_conductivity'
            else _PERFUSION_ABSENCE)
    return dict(base)


def _missing_record(region_id, property_name):
    if region_id == _NO_TENDON_DENSITY and property_name == 'density':
        reason = _TENDON_DENSITY_MISSING_REASON
    else:
        reason = _GENERIC_MISSING_REASON
    return {
        'status': MISSING,
        'region_id': region_id,
        'field': property_name,
        'reason': reason,
        'source': None,
        'source_sha256': None,
    }


def _resolve(region_id, property_name):
    """Internal single source of truth; returns a MaterialProperty or a status record."""
    if property_name == 'density':
        value = _DENSITY.get(region_id)
        return value if value is not None else _missing_record(region_id, property_name)
    if property_name == 'specific_heat':
        value = _SPECIFIC_HEAT.get(region_id)
        return value if value is not None else _missing_record(region_id, property_name)
    if property_name in ('thermal_conductivity', 'perfusion_heat_transport'):
        return _absent_record(property_name)
    raise KeyError(
        f'property {property_name!r} is not registered; registered properties are '
        f'{list(PROPERTY_NAMES)}'
    )


def assert_registered(region_id):
    """Return True for a registered region; otherwise raise KeyError with guidance.

    Registration is the only gate that lets a region carry a status at all. An unregistered
    region must not be silently bound to a generic default; add an audited registry row first.
    """
    if region_id not in REGION_IDS:
        raise KeyError(
            f'region {region_id!r} is not registered in the material roster. Registered region '
            f'ids: {list(REGION_IDS)}. Add an audited registry row (with source and sha256) before '
            f'binding material data; do not substitute a generic default.'
        )
    return True


def is_known(value):
    return isinstance(value, MaterialProperty)


def is_absent(value):
    return (not is_known(value)) and isinstance(value, dict) and value.get('status') == ABSENT


def is_missing(value):
    return (not is_known(value)) and isinstance(value, dict) and value.get('status') == MISSING


def lookup(region_id, property_name):
    """Return one of: MaterialProperty (KNOWN), an ABSENT record, or a MISSING record.

    ABSENT records carry ``status``/``reason``/``source``/``source_sha256``; MISSING records carry
    ``status``/``field``/``reason``. Neither is the number 0 nor None. Fresh dicts are returned so
    callers cannot mutate the roster.
    """
    assert_registered(region_id)
    value = _resolve(region_id, property_name)
    if is_known(value):
        return value
    return dict(value)


def roster_table():
    """One row per (region_id, property_name), ordered by region then property.

    Keys: region_id, property_name, status, value, units, provenance, source, source_sha256,
    reason. Non-KNOWN rows have value/units/provenance None and a nonempty reason.
    """
    rows = []
    for region_id in REGION_IDS:
        for property_name in PROPERTY_NAMES:
            value = lookup(region_id, property_name)
            if is_known(value):
                rows.append({
                    'region_id': region_id,
                    'property_name': property_name,
                    'status': KNOWN,
                    'value': value.value,
                    'units': value.units,
                    'provenance': value.provenance,
                    'source': value.source,
                    'source_sha256': value.source_sha256,
                    'reason': None,
                })
            else:
                rows.append({
                    'region_id': region_id,
                    'property_name': property_name,
                    'status': value['status'],
                    'value': None,
                    'units': None,
                    'provenance': None,
                    'source': value.get('source'),
                    'source_sha256': value.get('source_sha256'),
                    'reason': value['reason'],
                })
    return rows


def coverage_report():
    """Count KNOWN / ABSENT / MISSING per property; totals equal len(REGION_IDS)."""
    report = {name: {KNOWN: 0, ABSENT: 0, MISSING: 0} for name in PROPERTY_NAMES}
    for row in roster_table():
        report[row['property_name']][row['status']] += 1
    return report
