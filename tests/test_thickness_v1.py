#!/usr/bin/env python3
"""Acceptance tests for the thickness-distribution descriptor.

Run:
    OMP_NUM_THREADS=1 PYTHONPATH=src python tests/test_thickness_v1.py
    OMP_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q tests/test_thickness_v1.py

Tolerance is fixed BEFORE the run: exact algebra on a pure standard-library module is
compared at 1e-9 relative. All numbers here are synthetic; the only physiological content is
the cited poroelastic law tau = L^2/(H_A*k), which the test recomputes independently.
"""
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

from bodytwin.geometry import thickness_v1  # noqa: E402
from bodytwin.geometry.thickness_v1 import (  # noqa: E402
    PROVENANCE,
    ThicknessDistribution,
    ThicknessError,
    relaxation_time_relative_uncertainty,
    relaxation_time_s,
)

RTOL = 1e-9  # fixed before the run for exact algebra


def approx(actual, expected):
    return abs(actual - expected) <= RTOL * abs(expected)


def make_distribution(samples=(1.0, 2.0, 3.0, 4.0), **kwargs):
    params = {'method': 'synthetic thickness sweep', 'provenance': 'assumption'}
    params.update(kwargs)
    return ThicknessDistribution(samples, **params)


# --------------------------------------------------------------------------- distribution stats

def test_distribution_stats_hand_checkable():
    d = make_distribution([1, 2, 3, 4])
    assert d.n == 4
    assert d.minimum_mm == 1.0
    assert d.maximum_mm == 4.0
    assert approx(d.mean_mm, 2.5)
    assert approx(d.median_mm, 2.5)
    # population variance = 5/4 -> std = sqrt(1.25); CV = std / 2.5.
    expected_cv = math.sqrt(1.25) / 2.5
    assert approx(d.coefficient_of_variation, expected_cv)
    assert d.units == 'mm'
    assert tuple(d.samples_mm) == (1.0, 2.0, 3.0, 4.0)


def test_median_is_actual_middle_for_odd_n():
    d = make_distribution([5, 1, 3])
    assert approx(d.mean_mm, 3.0)
    assert approx(d.median_mm, 3.0)
    assert d.minimum_mm == 1.0 and d.maximum_mm == 5.0


def test_single_sample_is_a_distribution_with_zero_cv():
    d = make_distribution([2.0])
    assert d.n == 1
    assert approx(d.coefficient_of_variation, 0.0)
    assert d.summary()['n'] == 1


def test_summary_carries_units_and_enum_tags():
    d = make_distribution([1.0, 3.0], provenance='literal_cited', source='paper:table-2')
    s = d.summary()
    assert s['units'] == 'mm'
    assert s['kind'] == 'thickness_distribution'
    assert s['n'] == 2
    assert s['provenance'] in PROVENANCE
    assert s['source'] == 'paper:table-2'
    for key in ('minimum_mm', 'maximum_mm', 'mean_mm', 'median_mm',
                'coefficient_of_variation'):
        assert key in s


# --------------------------------------------------------------------------- validation

def test_units_m_rejected():
    for bad_units in ('m', 'cm', 'mm2', 'inch', '', 'MM'):
        try:
            make_distribution(units=bad_units)
            assert False, f'units={bad_units!r} must raise'
        except ThicknessError:
            pass


def test_empty_sample_rejected():
    for empty in ([], (), ()):
        try:
            make_distribution(empty)
            assert False, 'empty samples must raise'
        except ThicknessError:
            pass


def test_zero_and_negative_samples_rejected():
    for bad in ([0.0], [-1.0], [1.0, 0.0], [1.0, -0.5], [0]):
        try:
            make_distribution(bad)
            assert False, f'samples={bad!r} must raise'
        except ThicknessError:
            pass


def test_nonfinite_samples_rejected():
    for bad in ([float('nan')], [float('inf')], [1.0, float('-inf')]):
        try:
            make_distribution(bad)
            assert False, f'samples={bad!r} must raise'
        except ThicknessError:
            pass


def test_string_sample_rejected_not_coerced():
    try:
        make_distribution('1234')
        assert False, 'a string must not be accepted as a sample sequence'
    except ThicknessError:
        pass


def test_method_and_provenance_validated():
    try:
        make_distribution(method='')
        assert False, 'empty method must raise'
    except ThicknessError:
        pass
    try:
        make_distribution(provenance='guessed')
        assert False, 'unknown provenance must raise'
    except ThicknessError:
        pass
    assert PROVENANCE == ('measured', 'calibrated', 'literal_cited', 'assumption', 'UNKNOWN')


def test_thickness_error_is_value_error():
    assert issubclass(ThicknessError, ValueError)
    try:
        make_distribution([0.0])
    except ValueError:
        pass  # caught as a ValueError, as required by the API


# --------------------------------------------------------------------------- tau law

def test_tau_analytic_value_independent():
    L_mm = 1.0
    H_A = 0.7e6          # Pa
    k = 1e-15            # m^4/(N*s)
    # independent recomputation in the test, from the SI form tau = L_m^2 / (H_A * k)
    L_m = L_mm * 1.0e-3
    expected = (L_m ** 2) / (H_A * k)
    got = relaxation_time_s(L_mm, H_A, k)
    assert approx(got, expected)
    assert approx(got, 1.0e-6 / 0.7e-9)  # == 1428.571428... s


def test_tau_scaling_exponent_is_two():
    H_A = 0.7e6
    k = 1e-15
    tau_1 = relaxation_time_s(1.0, H_A, k)
    tau_2 = relaxation_time_s(2.0, H_A, k)   # L doubled
    tau_half = relaxation_time_s(0.5, H_A, k)
    assert approx(tau_2 / tau_1, 4.0)         # exponent 2
    assert approx(tau_1 / tau_half, 4.0)      # halving L quarters tau
    assert approx(tau_2, 4.0 * tau_1)


def test_tau_requires_positive_inputs():
    cases = [
        (0.0, 0.7e6, 1e-15),
        (-1.0, 0.7e6, 1e-15),
        (1.0, 0.0, 1e-15),
        (1.0, -0.7e6, 1e-15),
        (1.0, 0.7e6, 0.0),
        (1.0, 0.7e6, -1e-15),
        (float('nan'), 0.7e6, 1e-15),
        (1.0, float('inf'), 1e-15),
    ]
    for args in cases:
        try:
            relaxation_time_s(*args)
            assert False, f'{args!r} must raise'
        except ThicknessError:
            pass


def test_relative_uncertainty_doubles():
    assert approx(relaxation_time_relative_uncertainty(0.05), 0.10)
    assert approx(relaxation_time_relative_uncertainty(0.2), 0.4)
    assert relaxation_time_relative_uncertainty(0.0) == 0.0
    for bad in (-0.01, float('nan'), float('inf')):
        try:
            relaxation_time_relative_uncertainty(bad)
            assert False, f'{bad!r} must raise'
        except ThicknessError:
            pass


# --------------------------------------------------------------------------- scope guard

def test_no_physiological_claims_beyond_cited_law():
    doc = (thickness_v1.__doc__ or '').lower()
    assert 'no physiological claim' in doc
    assert 'l^2' in doc and 'h_a' in doc  # cites the law, does not extend it
    names = [n.upper() for n in dir(thickness_v1) if not n.startswith('__')]
    for banned in ('RANGE', 'MEASURED', 'PHYSIOLOGICAL', 'CLINICAL', 'GAIT', 'PATIENT'):
        assert not any(banned in n for n in names), f'unexpected empirical constant: {banned}'


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for fn in fns:
        fn()
        print('PASS', fn.__name__)
    print(f'{len(fns)} tests passed')
