"""First-order uncertainty vocabulary and arithmetic for region geometry quantities.

Location: src/bodytwin/geometry/uncertainty_v1.py
It is additive. It changes no existing symbol, imports nothing at import time
beyond the standard library, performs no I/O and writes nothing.

Purpose
-------
The region payload has no per-field uncertainty band of its own. This module is a small, pure *uncertainty vocabulary and arithmetic layer* that a future
payload assembler can attach to the existing quantities (density, area, volume, mass, ...).

Scope and honesty
-----------------
These helpers implement only **first-order, independent-uncertainty** combination rules
(quadrature and the product/power propagation shortcuts). They are deliberately NOT a full
GUM (JCGM 100:2008) treatment: no coverage factors, no effective degrees of freedom, no
correlations, no non-linear terms. Provenance strings mirror ``region_mass_v1.PROVENANCE``
(``('measured', 'calibrated', 'literal_cited', 'assumption', 'UNKNOWN')``) so the vocabulary
stays consistent with the audited material registry; the test suite cross-checks that the
two tuples are identical.

Units convention
----------------
An ``absolute`` uncertainty carries the units of the quantity it qualifies (caller declared).
A ``relative`` uncertainty is dimensionless and MUST declare ``units == '1'``; it is the
standard uncertainty divided by the magnitude of the value it qualifies.

No physiological claims are made here; all arithmetic is numeric and synthetic.
"""
from __future__ import annotations

import math
import numbers
from dataclasses import dataclass

__all__ = [
    'UncertaintyError',
    'Uncertainty',
    'KINDS',
    'PROVENANCE',
    'combine_standard',
    'relative_of',
    'propagate_product',
    'propagate_power',
    'degrees_of_freedom_note',
]

KINDS = ('absolute', 'relative')
PROVENANCE = ('measured', 'calibrated', 'literal_cited', 'assumption', 'UNKNOWN')

_DEGREES_OF_FREEDOM_NOTE = (
    'First-order, independent-uncertainty rules only (quadrature / product / power '
    'shortcuts); NOT a full GUM treatment (no coverage factor, degrees of freedom, '
    'correlations or non-linear terms).'
)


class UncertaintyError(ValueError):
    """An uncertainty value, unit, kind, provenance or operand violates the contract."""


def _is_real(x):
    """True for a finite real number, explicitly excluding bool and non-finite floats."""
    return (isinstance(x, numbers.Real) and not isinstance(x, bool)
            and math.isfinite(float(x)))


def _require_real(x, what):
    if not _is_real(x):
        raise UncertaintyError(f'{what} must be a finite real number, got {x!r}')
    return float(x)


def _require_nonneg_real(x, what):
    v = _require_real(x, what)
    if v < 0.0:
        raise UncertaintyError(f'{what} must be finite and >= 0, got {x!r}')
    return v


@dataclass(frozen=True)
class Uncertainty:
    """A declared standard uncertainty for a quantity.

    ``value`` is a non-negative standard uncertainty. ``kind`` is ``'absolute'`` (``units``
    names the quantity's units) or ``'relative'`` (``units`` is ``'1'``). ``provenance`` is
    the vocabulary shared with ``region_mass_v1``. ``source`` optionally names where the
    value came from (a citation or a note); it is never required.
    """

    value: float
    units: str
    kind: str
    provenance: str
    source: str | None = None

    def __post_init__(self):
        _require_nonneg_real(self.value, 'Uncertainty.value')
        if type(self.units) is not str or not self.units:
            raise UncertaintyError('Uncertainty.units must be a nonempty string')
        if self.kind not in KINDS:
            raise UncertaintyError(
                f'Uncertainty.kind {self.kind!r} not in {KINDS}')
        if self.kind == 'relative' and self.units != '1':
            raise UncertaintyError(
                f"relative uncertainty is dimensionless and must declare units '1', "
                f'got {self.units!r}')
        if self.provenance not in PROVENANCE:
            raise UncertaintyError(
                f'Uncertainty.provenance {self.provenance!r} not in {PROVENANCE}')
        if self.source is not None and type(self.source) is not str:
            raise UncertaintyError('Uncertainty.source must be a string or None')


def combine_standard(*uncertainties):
    """Combine independent standard uncertainties in quadrature.

    All operands must share the same ``kind``; absolute operands must additionally share the
    same ``units`` (relative operands are all dimensionless). Returns ``sqrt(sum(u_i^2))``.
    """
    if not uncertainties:
        raise UncertaintyError('combine_standard needs at least one Uncertainty')
    for u in uncertainties:
        if not isinstance(u, Uncertainty):
            raise UncertaintyError(
                f'combine_standard operands must be Uncertainty, got {type(u).__name__}')
    kinds = {u.kind for u in uncertainties}
    if len(kinds) != 1:
        raise UncertaintyError(
            f'cannot combine mixed kinds {sorted(kinds)}; all must be absolute or all relative')
    kind = kinds.pop()
    if kind == 'absolute':
        units = {u.units for u in uncertainties}
        if len(units) != 1:
            raise UncertaintyError(
                f'cannot combine absolute uncertainties with mixed units {sorted(units)}')
    total = math.sqrt(math.fsum(float(u.value) ** 2 for u in uncertainties))
    return total


def relative_of(value, uncertainty):
    """Return ``uncertainty`` expressed as a relative standard uncertainty of ``value``.

    ``uncertainty`` must be an absolute ``Uncertainty``. ``value`` must be finite and
    non-zero (division by ``abs(value)``).
    """
    v = _require_real(value, 'value')
    if v == 0.0:
        raise UncertaintyError('relative_of rejects value == 0 (division by zero)')
    if not isinstance(uncertainty, Uncertainty):
        raise UncertaintyError('uncertainty must be an Uncertainty')
    if uncertainty.kind != 'absolute':
        raise UncertaintyError(
            'relative_of expects an absolute uncertainty; '
            f'got kind {uncertainty.kind!r}')
    return Uncertainty(
        value=abs(float(uncertainty.value) / abs(v)),
        units='1',
        kind='relative',
        provenance=uncertainty.provenance,
        source=uncertainty.source,
    )


def propagate_product(factors):
    """Propagate relative standard uncertainties through a product of factors.

    ``factors`` is a non-empty iterable of ``(value, relative_standard_uncertainty)`` pairs
    where each relative uncertainty is dimensionless (units ``'1'``). Returns
    ``(product, u_rel(product))`` with ``u_rel = sqrt(sum(u_i^2))`` for independent factors.
    """
    items = list(factors)
    if not items:
        raise UncertaintyError('propagate_product needs at least one (value, u_rel) factor')
    product = 1.0
    sum_sq = 0.0
    for i, pair in enumerate(items):
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            raise UncertaintyError(
                f'factor {i} must be a (value, relative_uncertainty) pair, got {pair!r}')
        value = _require_real(pair[0], f'factor {i} value')
        u_rel = _require_nonneg_real(pair[1], f'factor {i} relative uncertainty')
        product *= value
        sum_sq += u_rel * u_rel
    product = float(product)
    if not math.isfinite(product):
        raise UncertaintyError('product is not finite')
    return product, math.sqrt(sum_sq)


def propagate_power(value, exponent, relative_uncertainty):
    """Return the relative standard uncertainty of ``value ** exponent``.

    For a pure power the first-order relative uncertainty is ``|exponent| * u_rel``. This
    encodes the ``tau ~ L**2`` law explicitly: exponent 2 doubles the relative
    uncertainty, exponent 3 triples it.
    """
    _require_real(value, 'value')
    exponent = _require_real(exponent, 'exponent')
    u_rel = _require_nonneg_real(relative_uncertainty, 'relative_uncertainty')
    result = abs(exponent) * u_rel
    if not math.isfinite(result):
        raise UncertaintyError('propagated relative uncertainty is not finite')
    return float(result)


def degrees_of_freedom_note():
    """Short constant caveat about the limits of these helpers."""
    return _DEGREES_OF_FREEDOM_NOTE
