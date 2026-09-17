"""Thickness distributions and the sensitivity of the cartilage relaxation law.

Adds a ``thickness_distribution`` descriptor and exposes the sensitivity of the cited
cartilage relaxation law.

Location: src/bodytwin/geometry/thickness_v1.py
It is additive, changes no existing symbol, and writes nothing.
Pure standard library; no I/O and no third-party imports.

What this provides
------------------
* ``ThicknessDistribution`` -- a ``thickness_distribution`` descriptor as a *distribution*
  (a sample of thicknesses), never a single scalar. Carries units, method, provenance and an
  optional source. Only ``mm`` is accepted; the caller converts first (raising beats guessing).
* ``relaxation_time_s`` -- the cited poroelastic law ``tau = L^2 / (H_A * k)`` in SI seconds,
  taking thickness in mm and converting to metres internally.
* ``relaxation_time_relative_uncertainty`` -- the exponent-2 sensitivity of that law: a
  relative thickness error ``dL/L`` propagates as ``dtau/tau = 2 * dL/L``.

No physiological claim is made here. The analytic helpers are arithmetic that reproduces the
law already computed in the read-only repo at
``src/bodytwin/cells/musculoskeletal/cartilage_poroelastic_relaxation.py``
(``tau(L_m, H_A, k) -> L_m**2 / (H_A * k)``); they do not assert that a particular tissue or
patient obeys the law, and no tissue constants are embedded.
"""
from __future__ import annotations

import math
import statistics
from dataclasses import dataclass
from typing import Sequence

__all__ = [
    'PROVENANCE',
    'MM_PER_M',
    'ThicknessError',
    'ThicknessDistribution',
    'relaxation_time_s',
    'relaxation_time_relative_uncertainty',
]

# Mirrors region_mass_v1.PROVENANCE (the five tags used by the I2 payload contract).
PROVENANCE = ('measured', 'calibrated', 'literal_cited', 'assumption', 'UNKNOWN')

# Fixed before use: 1 mm expressed in metres. Thickness is declared in mm; the law is SI.
MM_PER_M = 1.0e-3

# The only accepted thickness unit in the region payload (length convention: mm).
UNITS = 'mm'


class ThicknessError(ValueError):
    """A thickness distribution or relaxation-law input violates the I2 contract."""


def _positive_float(value, what):
    """Return ``float(value)`` if finite and strictly positive, else raise ThicknessError."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise ThicknessError(f'{what}: finite positive number required, got {value!r}')
    if not math.isfinite(v) or v <= 0.0:
        raise ThicknessError(f'{what}: finite positive number required, got {value!r}')
    return v


def _nonnegative_float(value, what):
    """Return ``float(value)`` if finite and non-negative, else raise ThicknessError."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise ThicknessError(f'{what}: finite non-negative number required, got {value!r}')
    if not math.isfinite(v) or v < 0.0:
        raise ThicknessError(f'{what}: finite non-negative number required, got {value!r}')
    return v


@dataclass(frozen=True)
class ThicknessDistribution:
    """A thickness sample (mm) for the ``thickness_distribution`` descriptor.

    ``samples_mm`` is an arbitrary-length sequence of finite, strictly positive millimetres.
    A single sample is allowed but the type is always a *distribution*: callers read the
    statistics instead of collapsing it to a scalar.
    """

    samples_mm: Sequence[float]
    method: str
    provenance: str
    units: str = UNITS
    source: str | None = None

    def __post_init__(self):
        if self.units != UNITS:
            raise ThicknessError(
                f'thickness units must be {UNITS!r}, got {self.units!r}; convert to mm before '
                'constructing the distribution')
        if type(self.method) is not str or not self.method:
            raise ThicknessError('method must be a nonempty string')
        if self.provenance not in PROVENANCE:
            raise ThicknessError(f'provenance {self.provenance!r} not in {PROVENANCE}')
        if self.source is not None and type(self.source) is not str:
            raise ThicknessError('source must be a string or None')
        if isinstance(self.samples_mm, (str, bytes, bytearray)):
            raise ThicknessError('samples_mm must be a sequence of numbers, not a string')
        try:
            samples = tuple(self.samples_mm)
        except TypeError:
            raise ThicknessError('samples_mm must be a sequence of finite positive numbers')
        if not samples:
            raise ThicknessError('samples_mm must contain at least one sample')
        checked = tuple(
            _positive_float(s, f'samples_mm[{i}]') for i, s in enumerate(samples))
        object.__setattr__(self, 'samples_mm', checked)

    @property
    def n(self):
        """Number of samples."""
        return len(self.samples_mm)

    @property
    def minimum_mm(self):
        return min(self.samples_mm)

    @property
    def maximum_mm(self):
        return max(self.samples_mm)

    @property
    def mean_mm(self):
        return math.fsum(self.samples_mm) / self.n

    @property
    def median_mm(self):
        return float(statistics.median(self.samples_mm))

    @property
    def std_mm(self):
        """Population standard deviation (mm); the denominator is ``n``, not ``n - 1``."""
        m = self.mean_mm
        return math.sqrt(math.fsum((x - m) ** 2 for x in self.samples_mm) / self.n)

    @property
    def coefficient_of_variation(self):
        """Population std / mean (dimensionless). Samples are positive so the mean is > 0."""
        return self.std_mm / self.mean_mm

    def summary(self):
        """Serialisable summary of the distribution, always carrying ``units``."""
        return {
            'kind': 'thickness_distribution',
            'units': self.units,
            'n': self.n,
            'minimum_mm': self.minimum_mm,
            'maximum_mm': self.maximum_mm,
            'mean_mm': self.mean_mm,
            'median_mm': self.median_mm,
            'coefficient_of_variation': self.coefficient_of_variation,
            'std_mm': self.std_mm,
            'method': self.method,
            'provenance': self.provenance,
            'source': self.source,
            'samples_mm': list(self.samples_mm),
        }


def relaxation_time_s(mean_thickness_mm, aggregate_modulus_pa, permeability_m4_per_Ns):
    """Poroelastic relaxation time ``tau = L^2 / (H_A * k)`` in SI seconds.

    Reproduces the repository law
    ``src/bodytwin/cells/musculoskeletal/cartilage_poroelastic_relaxation.py``::

        tau(L_m, H_A, k) = L_m ** 2 / (H_A * k)

    ``mean_thickness_mm`` is the characteristic length ``L`` in **mm** and is converted to
    metres (``* 1e-3``) before squaring; ``aggregate_modulus_pa`` is ``H_A`` in Pa and
    ``permeability_m4_per_Ns`` is ``k`` in m^4/(N*s). All three must be finite and positive.
    The result is in seconds because the inputs are SI once the mm->m conversion is applied.

    This is pure arithmetic; it makes no claim that any tissue follows the law.
    """
    length_mm = _positive_float(mean_thickness_mm, 'mean_thickness_mm')
    h_a = _positive_float(aggregate_modulus_pa, 'aggregate_modulus_pa')
    k = _positive_float(permeability_m4_per_Ns, 'permeability_m4_per_Ns')
    length_m = length_mm * MM_PER_M
    return (length_m ** 2) / (h_a * k)


def relaxation_time_relative_uncertainty(relative_thickness_uncertainty):
    """Propagate a relative thickness uncertainty to tau: ``dtau/tau = 2 * dL/L``.

    Because ``tau`` depends on ``L`` as ``L**2``, the relative error in ``tau`` is twice the
    relative error in ``L``. ``relative_thickness_uncertainty`` must be finite and
    non-negative (0 is valid and maps to 0).
    """
    rel = _nonnegative_float(relative_thickness_uncertainty,
                             'relative_thickness_uncertainty')
    return 2.0 * rel
