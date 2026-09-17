"""Lumped heat-capacity boundary for region geometry (algebra only).

Location: src/bodytwin/geometry/thermal_capacity_v1.py
It is additive. It changes no existing symbol, adds no state, performs no I/O and
writes nothing.

Purpose
-------
``region_mass_v1`` already *emits* a lumped heat capacity in J/K (units are declared, and a missing
specific heat stays UNKNOWN -- see ``region_mass_v1``). This module is the strictly algebraic
boundary that consumes such a capacity:

    Q = C * dT  ->  dT = Q / C

Nothing else. In particular this file deliberately does NOT model conduction, perfusion, blood
flow, convection, radiation, spatial gradients or any transport: those are out of scope here. Every
helper here is a pure scalar identity.

ABSTAIN contract
----------------
If a heat capacity is unknown (``None``) or zero, ``temperature_rise_K`` raises
``ThermalCapacityError`` -- it never invents an infinite rise and never falls back to a zero
rise. ``ABSTAIN_note()`` documents that an unknown capacity yields no temperature estimate.

Adiabatic lumped bound
----------------------
``temperature_rise_K`` is the adiabatic, lumped-capacity bound: all heat is stored in one
capacity, no heat leaves the node. It is an upper bound on the real temperature change whenever
any loss path exists.

C / AUDIT boundary and no double counting
-----------------------------------------
``metabolic_heat_J`` is the single place where mechanical work is split into heat. The
complement ``(1 - efficiency)`` is the heat fraction; it is the C/AUDIT boundary. Tendon /
series-elastic dissipation is NOT added again here, because the efficiency already bundles the
metabolic heat fraction by construction. Adding a separate tendon-dissipation term would double
count.

Constants and honesty
---------------------
The only physical constants embedded are the audited lumped literals in the read-only repo
("audited" means the value is present in the repo with a file:line pointer; its primary
literature source is not verified here):
``SPECIFIC_HEAT_BODY_J_PER_KG_K = 3490.0`` (thermoregulation.py) is not re-declared here --
callers pass a ``region_mass_v1.MaterialProperty`` -- and
``LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G = 2426.0`` (thermoregulation.py) is re-declared as a
named constant. The latent heat is a whole-body lumped constant, NOT a per-tissue value. No
physiological claim is made; all test numbers are synthetic.
"""
from __future__ import annotations

import math
import numbers

from . import region_mass_v1

__all__ = [
    'ThermalCapacityError',
    'LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G',
    'LATENT_HEAT_SOURCE',
    'SPECIFIC_HEAT_SOURCE',
    'heat_capacity_J_per_K',
    'temperature_rise_K',
    'metabolic_heat_J',
    'sweat_heat_removal_J',
    'ABSTAIN_note',
]

# Audited lumped constant reused verbatim from the read-only repo, with its source string.
# thermoregulation.py LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G -- "at ~skin temperature";
# Gagge two-node / ISO 7933 Predicted Heat Strain. This is a whole-body lumped value, NOT a per-tissue latent heat.
LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G = 2426.0
LATENT_HEAT_SOURCE = (
    'src/bodytwin/cells/organ_systems/thermoregulation.py::LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G '
    'LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G = 2426.0 J/g '
    '(Gagge two-node / ISO 7933; whole-body lumped, not per-tissue)'
)

# Source of the whole-body specific heat that region_mass_v1.MaterialProperty carries.
SPECIFIC_HEAT_SOURCE = (
    'src/bodytwin/cells/organ_systems/thermoregulation.py::SPECIFIC_HEAT_BODY_J_PER_KG_K '
    'SPECIFIC_HEAT_BODY_J_PER_KG_K = 3490.0 J/(kg*K) '
    '(whole-body lumped; no per-tissue specific heat exists in the repository)'
)

_ABSTAIN_NOTE = (
    'Unknown heat capacity ABSTAINS: with no declared capacity (None) there is no temperature '
    'change estimate. A missing capacity is NOT zero capacity and NOT an infinite rise; '
    'temperature_rise_K raises ThermalCapacityError instead of fabricating a value.'
)


class ThermalCapacityError(ValueError):
    """A heat-capacity, heat-flow, efficiency or latent-heat input violates the I2 contract."""


def _real(value, what):
    """Return ``float(value)`` if it is a finite real number, excluding bool; else raise."""
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ThermalCapacityError(f'{what} must be a finite real number, got {value!r}')
    v = float(value)
    if not math.isfinite(v):
        raise ThermalCapacityError(f'{what} must be finite, got {value!r}')
    return v


def _positive_real(value, what):
    v = _real(value, what)
    if v <= 0.0:
        raise ThermalCapacityError(f'{what} must be finite and > 0, got {value!r}')
    return v


def _nonnegative_real(value, what):
    v = _real(value, what)
    if v < 0.0:
        raise ThermalCapacityError(f'{what} must be finite and >= 0, got {value!r}')
    return v


def heat_capacity_J_per_K(mass_kg, specific_heat):
    """Lumped heat capacity ``C = m * cp`` in J/K.

    ``mass_kg`` must be a finite real number > 0. ``specific_heat`` must be a
    ``region_mass_v1.MaterialProperty`` carrying a specific-heat unit in the audited contract
    (``J/(kg*K)``, ``J/(g*K)`` or ``kJ/(kg*K)``); its value is converted through
    ``as_cp_J_per_kgK()`` so the unit handling is shared, never re-derived here.

    A missing material (``None``) or a raw number instead of a ``MaterialProperty`` is refused
    with ``ThermalCapacityError`` -- an unknown capacity is never defaulted.
    """
    m = _positive_real(mass_kg, 'mass_kg')
    if specific_heat is None:
        raise ThermalCapacityError(
            'specific_heat is None: heat capacity is unknown, refusing to default')
    if not isinstance(specific_heat, region_mass_v1.MaterialProperty):
        raise ThermalCapacityError(
            'specific_heat must be a region_mass_v1.MaterialProperty with explicit units, '
            f'got {type(specific_heat).__name__}')
    try:
        cp = specific_heat.as_cp_J_per_kgK().value
    except region_mass_v1.UnitError as exc:
        raise ThermalCapacityError(f'specific_heat units rejected: {exc}') from exc
    if cp is None:
        raise ThermalCapacityError('specific_heat value is unknown (None)')
    return float(m * cp)


def temperature_rise_K(heat_J, heat_capacity_J_per_K):
    """Adiabatic lumped temperature change ``dT = Q / C`` in K.

    This is the pure algebraic identity only: all heat is stored in one lumped capacity and no
    heat leaves the node (an adiabatic upper bound). It intentionally models no conduction,
    perfusion or other transport -- those are out of scope here.

    ``heat_J`` may be negative (cooling), zero or positive and must be finite.
    ``heat_capacity_J_per_K`` must be finite and strictly > 0. A capacity of ``None`` (unknown)
    or ``0`` raises ``ThermalCapacityError``: capacity unknown is NOT zero capacity and NOT an
    infinite rise, so no temperature is fabricated.
    """
    Q = _real(heat_J, 'heat_J')
    if heat_capacity_J_per_K is None:
        raise ThermalCapacityError(
            'heat_capacity_J_per_K is None: capacity unknown ABSTAINS, no temperature estimate')
    C = _positive_real(heat_capacity_J_per_K, 'heat_capacity_J_per_K')
    return float(Q / C)


def metabolic_heat_J(mechanical_work_J, efficiency):
    """Heat released by metabolism for a mechanical work output: ``Q = (1 - e) * W`` in J.

    ``efficiency`` is the dimensionless fraction of metabolic energy delivered as mechanical
    work and must lie in ``[0, 1)``. A value of exactly ``1.0`` (no heat) or outside the range
    is refused. ``mechanical_work_J`` must be finite and >= 0.

    This function is the C/AUDIT boundary. The heat fraction ``(1 - efficiency)`` already
    includes all metabolic losses attributable to the work, so tendon / series-elastic
    dissipation must NOT be added again here -- doing so would double count.
    """
    work = _nonnegative_real(mechanical_work_J, 'mechanical_work_J')
    e = _real(efficiency, 'efficiency')
    if not (0.0 <= e < 1.0):
        raise ThermalCapacityError(
            f'efficiency must be in [0, 1), got {efficiency!r}')
    return float((1.0 - e) * work)


def sweat_heat_removal_J(evaporated_water_g, latent_heat_J_per_g=LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G):
    """Evaporative heat removed by a mass of sweat: ``Q = m_evap * L_vap`` in J.

    ``evaporated_water_g`` is the *evaporated* mass in grams (finite, >= 0). ``latent_heat_J_per_g``
    defaults to the audited whole-body lumped constant ``2426.0``
    (``LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G``, sourced in ``LATENT_HEAT_SOURCE``). It is a
    single lumped value and NOT a per-tissue latent heat; callers may override it only with an
    equally explicit, sourced constant.
    """
    mass = _nonnegative_real(evaporated_water_g, 'evaporated_water_g')
    latent = _positive_real(latent_heat_J_per_g, 'latent_heat_J_per_g')
    return float(mass * latent)


def ABSTAIN_note():
    """Constant note: an unknown capacity yields no temperature change estimate."""
    return _ABSTAIN_NOTE
