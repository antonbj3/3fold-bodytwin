#!/usr/bin/env python3
"""Acceptance tests for the heat-capacity boundary (thermal_capacity_v1).

Run:
    OMP_NUM_THREADS=1 PYTHONPATH=src python tests/test_thermal_capacity_v1.py
    OMP_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q tests/test_thermal_capacity_v1.py

Tolerances were fixed BEFORE the run: the identities here are exact algebra, so they are
compared at RTOL = 1e-12 relative. Every number is a small synthetic value; no physiological
claim is made or tested. The only physical constants exercised are the audited lumped literals
(whole-body specific heat via region_mass_v1, latent heat 2426 J/g).
"""
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

from bodytwin.geometry.thermal_capacity_v1 import (  # noqa: E402
    ABSTAIN_note,
    LATENT_HEAT_SOURCE,
    LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G,
    ThermalCapacityError,
    heat_capacity_J_per_K,
    metabolic_heat_J,
    sweat_heat_removal_J,
    temperature_rise_K,
)
from bodytwin.geometry import region_mass_v1  # noqa: E402  (reused unit handling / MaterialProperty)

RTOL = 1e-12  # fixed before the run for exact float64 algebraic identities


def approx(a, b, rtol=RTOL):
    return math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=0.0)


def raises(fn):
    try:
        fn()
    except ThermalCapacityError:
        return True
    return False


def _cp(value=3490.0, units='J/(kg*K)'):
    """A synthetic specific-heat MaterialProperty in the audited unit contract."""
    return region_mass_v1.MaterialProperty(
        'specific_heat_body', value, units, 'literal_cited', source='synthetic-test')


# ------------------------------------------------------------------- heat_capacity_J_per_K

def test_heat_capacity_is_mass_times_cp():
    # 2.0 kg * 3490 J/(kg*K) == 6980 J/K
    C = heat_capacity_J_per_K(2.0, _cp(3490.0, 'J/(kg*K)'))
    assert approx(C, 6980.0), C
    assert type(C) is float


def test_heat_capacity_converts_specific_heat_units():
    # 2.0 kg * 3.49 kJ/(kg*K) == 2.0 kg * 3490 J/(kg*K) == 6980 J/K
    C_kj = heat_capacity_J_per_K(2.0, _cp(3.49, 'kJ/(kg*K)'))
    C_g = heat_capacity_J_per_K(2.0, _cp(3.49, 'J/(g*K)'))
    assert approx(C_kj, 6980.0), C_kj
    assert approx(C_g, 6980.0), C_g


def test_heat_capacity_rejects_missing_material():
    assert raises(lambda: heat_capacity_J_per_K(2.0, None))


def test_heat_capacity_rejects_raw_number_instead_of_material_property():
    # a bare cp number carries no units -> refused, not guessed
    assert raises(lambda: heat_capacity_J_per_K(2.0, 3490.0))
    assert raises(lambda: heat_capacity_J_per_K(2.0, '3490 J/(kg*K)'))


def test_heat_capacity_rejects_nonpositive_or_nonfinite_mass():
    for bad in (0.0, -1.0, float('nan'), float('inf'), -float('inf')):
        assert raises(lambda b=bad: heat_capacity_J_per_K(b, _cp())), bad


# -------------------------------------------------------------------- temperature_rise_K

def test_temperature_rise_one_kelvin_from_matching_joules():
    # Q = 6980 J into C = 6980 J/K -> dT = 1.0 K
    dT = temperature_rise_K(6980.0, 6980.0)
    assert approx(dT, 1.0), dT
    assert type(dT) is float


def test_temperature_rise_zero_heat_is_zero():
    assert approx(temperature_rise_K(0.0, 6980.0), 0.0)


def test_temperature_rise_negative_heat_cools():
    assert approx(temperature_rise_K(-6980.0, 6980.0), -1.0)


def test_temperature_rise_unknown_capacity_abstains():
    # None capacity -> refuse; no infinite rise and no zero fallback
    assert raises(lambda: temperature_rise_K(100.0, None))


def test_temperature_rise_zero_capacity_raises():
    assert raises(lambda: temperature_rise_K(100.0, 0.0))


def test_temperature_rise_rejects_negative_or_nonfinite_capacity():
    for bad in (-1.0, float('nan'), float('inf'), -float('inf')):
        assert raises(lambda b=bad: temperature_rise_K(100.0, b)), bad


def test_temperature_rise_rejects_nonfinite_heat():
    for bad in (float('nan'), float('inf'), -float('inf')):
        assert raises(lambda b=bad: temperature_rise_K(b, 6980.0)), bad


def test_temperature_rise_is_algebraic_inverse_of_capacity():
    # dT * C recovers Q exactly (same operand), the pure Q = C dT identity
    Q, C = 1234.5678, 4321.0
    assert approx(temperature_rise_K(Q, C) * C, Q)


# ----------------------------------------------------------------------- metabolic_heat_J

def test_metabolic_heat_quarter_efficiency():
    # (1 - 0.25) * 1000 J == 750 J
    Q = metabolic_heat_J(1000.0, 0.25)
    assert approx(Q, 750.0), Q
    assert type(Q) is float


def test_metabolic_heat_zero_work_is_zero():
    assert approx(metabolic_heat_J(0.0, 0.25), 0.0)


def test_metabolic_heat_rejects_efficiency_at_or_above_one():
    for bad in (1.0, 1.0000001, 2.0):
        assert raises(lambda b=bad: metabolic_heat_J(1000.0, b)), bad


def test_metabolic_heat_rejects_negative_efficiency():
    for bad in (-0.1, -1.0):
        assert raises(lambda b=bad: metabolic_heat_J(1000.0, b)), bad


def test_metabolic_heat_rejects_nonfinite_inputs():
    for bad in (float('nan'), float('inf'), -float('inf')):
        assert raises(lambda b=bad: metabolic_heat_J(1000.0, b)), bad
        assert raises(lambda b=bad: metabolic_heat_J(b, 0.25)), bad


# -------------------------------------------------------------------- sweat_heat_removal_J

def test_sweat_heat_removal_one_gram_uses_audited_latent_heat():
    assert approx(sweat_heat_removal_J(1.0), 2426.0)
    assert approx(sweat_heat_removal_J(1.0), LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G)


def test_sweat_heat_removal_scales_and_handles_zero():
    assert approx(sweat_heat_removal_J(2.5), 2.5 * 2426.0)
    assert approx(sweat_heat_removal_J(0.0), 0.0)


def test_sweat_heat_removal_rejects_negative_mass_or_bad_latent():
    assert raises(lambda: sweat_heat_removal_J(-1.0))
    assert raises(lambda: sweat_heat_removal_J(1.0, 0.0))
    assert raises(lambda: sweat_heat_removal_J(1.0, -2426.0))
    assert raises(lambda: sweat_heat_removal_J(float('nan')))


def test_latent_heat_source_names_the_audited_constant():
    assert isinstance(LATENT_HEAT_SOURCE, str) and LATENT_HEAT_SOURCE
    assert '2426' in LATENT_HEAT_SOURCE
    assert 'thermoregulation' in LATENT_HEAT_SOURCE


# ------------------------------------------------------------------------- ABSTAIN contract

def test_abstain_note_is_constant_and_documents_unknown_capacity():
    note = ABSTAIN_note()
    assert isinstance(note, str) and note
    assert note == ABSTAIN_note()
    assert 'ABSTAIN' in note
    assert 'unknown' in note.lower()
    assert 'no temperature' in note.lower()


# ------------------------------------------------------------------- reuse / no transport

def test_material_property_is_the_shared_region_mass_type():
    prop = _cp()
    assert isinstance(prop, region_mass_v1.MaterialProperty)
    assert prop.as_cp_J_per_kgK().units == 'J/(kg*K)'


def test_module_exposes_no_transport_helpers():
    # conduction/perfusion/spatial transport are out of scope; this module must not grow those symbols
    from bodytwin.geometry import thermal_capacity_v1 as mod
    forbidden = {'conduction', 'perfusion', 'blood_flow', 'conductance', 'diffusion',
                 'gradient', 'transport', 'spatial'}
    names = {n.lower() for n in dir(mod)}
    assert not (forbidden & names), sorted(forbidden & names)


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for fn in fns:
        fn()
        print('PASS', fn.__name__)
    print(f'{len(fns)} tests passed')
