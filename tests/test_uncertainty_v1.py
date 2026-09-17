#!/usr/bin/env python3
"""Acceptance tests for the uncertainty helpers (uncertainty_v1).

Run:
    OMP_NUM_THREADS=1 PYTHONPATH=src python tests/test_uncertainty_v1.py
    OMP_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q tests/test_uncertainty_v1.py

Tolerances were fixed BEFORE the run: exact algebraic identities (quadrature, relative
conversion, product/power propagation) are compared at RTOL = 1e-12 relative. All numbers
are small synthetic values; no physiological claims are made or tested.
"""
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

from bodytwin.geometry.uncertainty_v1 import (  # noqa: E402
    KINDS,
    PROVENANCE,
    Uncertainty,
    UncertaintyError,
    combine_standard,
    degrees_of_freedom_note,
    propagate_power,
    propagate_product,
    relative_of,
)
from bodytwin.geometry import region_mass_v1  # noqa: E402  (provenance vocabulary cross-check only)

RTOL = 1e-12  # fixed before the run for exact float64 / analytic formulas


def approx(a, b, rtol=RTOL):
    return math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=0.0)


def raises(fn):
    try:
        fn()
    except UncertaintyError:
        return True
    return False


def _abs(value, units='mm', provenance='measured', source='synthetic'):
    return Uncertainty(value, units, 'absolute', provenance, source)


def _rel(value, provenance='measured', source='synthetic'):
    return Uncertainty(value, '1', 'relative', provenance, source)


# --------------------------------------------------------------------------- quadrature

def test_quadrature_three_four_five():
    u = combine_standard(_abs(3.0, 'mm'), _abs(4.0, 'mm'))
    assert approx(u, 5.0), u
    assert type(u) is float


def test_quadrature_single_and_zero():
    assert approx(combine_standard(_abs(2.5, 'kg')), 2.5)
    assert approx(combine_standard(_abs(0.0, 'kg'), _abs(0.0, 'kg')), 0.0)


def test_combine_requires_operands():
    assert raises(lambda: combine_standard())
    assert raises(lambda: combine_standard(_abs(1.0), 1.0))


def test_mixing_kinds_rejected():
    assert raises(lambda: combine_standard(_abs(1.0, 'mm'), _rel(0.1)))


def test_mixing_absolute_units_rejected():
    assert raises(lambda: combine_standard(_abs(1.0, 'mm'), _abs(1.0, 'kg')))


# ----------------------------------------------------------------------- relative_of

def test_relative_of_round_trip():
    value = -40.0
    abs_u = _abs(2.0, 'mm', provenance='calibrated', source='cal-7')
    rel = relative_of(value, abs_u)
    assert rel.kind == 'relative'
    assert rel.units == '1'
    assert rel.provenance == 'calibrated'
    assert rel.source == 'cal-7'
    assert approx(rel.value, 2.0 / 40.0)
    # round-trip: rel * |value| recovers the absolute standard uncertainty
    assert approx(rel.value * abs(value), abs_u.value)


def test_relative_of_rejects_zero_and_nonfinite_value():
    for bad in (0.0, float('nan'), float('inf'), -float('inf')):
        assert raises(lambda b=bad: relative_of(b, _abs(1.0))), bad


def test_relative_of_rejects_relative_input():
    assert raises(lambda: relative_of(10.0, _rel(0.1)))
    assert raises(lambda: relative_of(10.0, 'not-an-uncertainty'))


# ------------------------------------------------------------------ propagate_product

def test_product_propagation_matches_hand_computation():
    # synthetic rho * A * L with independent relative standard uncertainties
    factors = [(1.05, 0.02), (12.0, 0.03), (200.0, 0.04)]
    product, u_rel = propagate_product(factors)
    assert approx(product, 1.05 * 12.0 * 200.0)
    assert approx(u_rel, math.sqrt(0.02 ** 2 + 0.03 ** 2 + 0.04 ** 2))
    assert approx(u_rel, math.sqrt(0.0029))


def test_product_propagation_zero_uncertainty():
    product, u_rel = propagate_product([(2.0, 0.0), (3.0, 0.0)])
    assert approx(product, 6.0)
    assert approx(u_rel, 0.0)


def test_product_propagation_rejects_bad_factors():
    assert raises(lambda: propagate_product([]))
    assert raises(lambda: propagate_product([(1.0,)]))
    assert raises(lambda: propagate_product([(1.0, float('nan'))]))
    assert raises(lambda: propagate_product([(float('inf'), 0.1)]))
    assert raises(lambda: propagate_product([(1.0, -0.1)]))
    assert raises(lambda: propagate_product([(1.0, 0.1), (float('inf'), 0.1)]))


# -------------------------------------------------------------------- propagate_power

def test_propagate_power_doubles_and_triples():
    u = 0.05
    assert approx(propagate_power(0.30, 2, u), 0.10)   # tau ~ L**2: doubled
    assert approx(propagate_power(0.30, 3, u), 0.15)   # tripled
    assert approx(propagate_power(0.30, -2, u), 0.10)  # magnitude of exponent
    assert approx(propagate_power(0.30, 1, u), u)
    assert approx(propagate_power(0.30, 0, u), 0.0)


def test_propagate_power_rejects_bad_input():
    assert raises(lambda: propagate_power(float('nan'), 2, 0.05))
    assert raises(lambda: propagate_power(0.3, float('inf'), 0.05))
    assert raises(lambda: propagate_power(0.3, 2, -0.05))
    assert raises(lambda: propagate_power(0.3, 2, float('nan')))


# --------------------------------------------------------------- construction guards

def test_negative_zero_nonfinite_rejected():
    assert raises(lambda: _abs(-1.0))
    assert raises(lambda: _rel(float('nan')))
    assert raises(lambda: _abs(float('inf')))
    assert raises(lambda: _abs(-float('inf')))
    assert raises(lambda: _abs(True))  # bool is not a legitimate uncertainty value


def test_zero_uncertainty_is_allowed():
    # >= 0 is required; exactly zero is a legitimate (degenerate) standard uncertainty
    assert _abs(0.0).value == 0.0
    assert _rel(0.0).value == 0.0


def test_missing_or_empty_provenance_rejected():
    for bad in ('', None, 'guessed', 'MEASURED'):
        assert raises(lambda b=bad: Uncertainty(1.0, 'mm', 'absolute', b)), bad
    assert raises(lambda: _abs(1.0, provenance='assumption', source=123))


def test_kind_validated_and_relative_units_dimensionless():
    assert raises(lambda: Uncertainty(1.0, 'mm', 'Absolute', 'measured'))
    assert raises(lambda: Uncertainty(1.0, 'mm', '1', 'measured'))
    assert raises(lambda: Uncertainty(1.0, 'mm', 'relative', 'measured'))  # wrong units
    assert raises(lambda: Uncertainty(1.0, '', 'absolute', 'measured'))    # empty units
    # the legal relative form
    assert _rel(0.1).units == '1'


# --------------------------------------------------------------------- provenance sync

def test_provenance_vocabulary_matches_region_mass_v1():
    assert tuple(PROVENANCE) == tuple(region_mass_v1.PROVENANCE)
    assert tuple(KINDS) == ('absolute', 'relative')
    # reuse an actual MaterialProperty to source the provenance vocabulary
    prop = region_mass_v1.MaterialProperty('density_synthetic', 1.05, 'kg/L', 'literal_cited')
    u = Uncertainty(0.05, prop.units, 'absolute', prop.provenance)
    assert u.provenance in PROVENANCE


def test_degrees_of_freedom_note_states_the_limits():
    note = degrees_of_freedom_note()
    assert isinstance(note, str) and note
    assert note == degrees_of_freedom_note()
    assert 'GUM' in note
    assert 'first-order' in note.lower()
    assert 'independent' in note.lower()


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for fn in fns:
        fn()
        print('PASS', fn.__name__)
    print(f'{len(fns)} tests passed')
