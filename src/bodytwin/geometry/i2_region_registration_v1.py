"""Opt-in cross-check that binds the frozen I2 -> I1 payload's region ids to the audited material
roster (``material_roster_v1``), the region-mapping prerequisite for region mass. It does NOT
touch the pure structural validator (``i2_payload_validate_v1``); it only *calls* it and wraps
its errors with a ``'structural: '`` prefix.

Location: src/bodytwin/geometry/i2_region_registration_v1.py
Pure: no I/O, reads nothing from disk, imports no ``json`` and writes nothing. It reuses
``material_roster_v1`` only (never ``region_mass_v1`` directly) and adds a separate,
opt-in gate on top of structural validation.

Design note (region -> material binding)
----------------------------------------
A payload region carries a synthetic geometry label in ``region_id`` (e.g. ``SYNTH-BOX-2``)
and an explicit roster binding in ``material_region_id`` (e.g. ``MSK-MUSCLE``). Registration
resolves ``material_region_id`` first; if it is absent, ``region_id`` itself must be a roster
id. No prefix heuristic and no implicit fallback: an unbound or unknown id is an issue, and no
generic/default row may be substituted.
"""
from __future__ import annotations

import math

from . import i2_payload_validate_v1, material_roster_v1

__all__ = ['assert_regions_registered', 'registration_issues']

_REL_TOL = 1e-12
_ABS_TOL = 0.0


def _is_finite_number(value):
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(float(value)))


def _values_agree(declared, roster):
    return math.isclose(float(declared), float(roster), rel_tol=_REL_TOL, abs_tol=_ABS_TOL)


def _resolve_roster_region_id(region, region_id):
    """Return the roster key for a region, or ``None`` when it cannot be bound.

    Resolution order: an explicit ``material_region_id`` (which must itself be registered),
    otherwise ``region_id`` when it is a registered roster id. No prefix heuristic and no
    density-material guessing; an unbound region returns ``None`` and is reported.
    """
    explicit = region.get('material_region_id')
    if explicit is not None:
        if not isinstance(explicit, str) or not explicit:
            return None
        try:
            material_roster_v1.assert_registered(explicit)
            return explicit
        except KeyError:
            return None
    try:
        material_roster_v1.assert_registered(region_id)
        return region_id
    except KeyError:
        return None


def registration_issues(payload, *, require_known_density=True) -> list[str]:
    """Return registration/cross-check issues for the frozen I2 -> I1 payload.

    Order of checks per region:
      1. structural errors from ``i2_payload_validate_v1.validate_payload`` (prefixed
         ``'structural: '``), always first;
      2. the region id must bind to a material-roster region;
      3. when ``require_known_density`` is true, the roster density status must be KNOWN
         (ABSENT/MISSING are reported and no default may be substituted);
      4. when the region declares a density ``value``/``units`` and the roster row is KNOWN,
         the declared value (relative tolerance 1e-12) and units (exact) must agree.

    Never raises for a malformed payload; an empty list means no issue was found.
    """
    issues = ['structural: ' + error
              for error in i2_payload_validate_v1.validate_payload(payload)]
    if not isinstance(payload, dict):
        return issues
    regions = payload.get('regions')
    if not isinstance(regions, list):
        return issues

    registered = list(material_roster_v1.REGION_IDS)
    for index, region in enumerate(regions):
        path = f'regions[{index}]'
        if not isinstance(region, dict):
            continue
        region_id = region.get('region_id')
        if not isinstance(region_id, str) or not region_id:
            continue  # structural validation already reported the malformed id

        roster_id = _resolve_roster_region_id(region, region_id)
        if roster_id is None:
            binding = region.get('material_region_id')
            label = (f'material_region_id {binding!r}' if binding is not None
                     else f'region_id {region_id!r}')
            issues.append(
                f'{path}: {label} is not registered in the material roster '
                f'(registered region ids: {registered}). Add an audited registry row (with '
                f'source and sha256) before binding material data; do not substitute a generic '
                f'default.')
            continue

        density = material_roster_v1.lookup(roster_id, 'density')
        if not material_roster_v1.is_known(density):
            if require_known_density:
                issues.append(
                    f'{path}.density: material roster has status {density["status"]!r} for '
                    f'{roster_id!r}/density ({density.get("reason")}); no default may be '
                    f'substituted.')
            continue

        declared = region.get('density')
        if not isinstance(declared, dict):
            continue  # explicit null is a structural/binding concern, not a value mismatch
        value = declared.get('value')
        units = declared.get('units')
        if value is not None and _is_finite_number(value) and not _values_agree(value, density.value):
            issues.append(
                f'{path}.density.value: declared {value!r} disagrees with the roster value '
                f'{density.value!r} for {roster_id!r} (rel tol {_REL_TOL:g}).')
        if isinstance(units, str) and units != density.units:
            issues.append(
                f'{path}.density.units: declared {units!r} disagrees with the roster units '
                f'{density.units!r} for {roster_id!r}.')

    return issues


def assert_regions_registered(payload, **kwargs) -> dict:
    """Return ``payload`` unchanged when it has no registration issue; raise otherwise.

    ``**kwargs`` are forwarded to ``registration_issues`` (e.g. ``require_known_density=False``).
    """
    issues = registration_issues(payload, **kwargs)
    if issues:
        raise ValueError('i2 region registration invalid:\n- ' + '\n- '.join(issues))
    return payload
