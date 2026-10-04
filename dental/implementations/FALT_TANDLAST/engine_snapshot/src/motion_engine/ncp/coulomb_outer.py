"""Global rational energy enclosure of every dissipative 3D contact solution.

Frozen sparse J, diagonal positive mass, hard gap box, normal compliance>=0.
Includes the exact circular Coulomb law. No existence/uniqueness or sharpness
claim; the returned status is distinct from a complete sharp mode union.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256
from math import isqrt
from operator import index
import json
from time import perf_counter
from .observable_cert import _q
from .gap_box import rational_json, dot

def sqrt_upper(x):
    """Dyadic upper bound, checked by squaring; no float in the proof."""
    x = Q(x)
    if x < 0:
        raise ValueError('negative square')
    den = 2 ** 40
    scaled = x.numerator * den * den
    floor = isqrt(scaled // x.denominator)
    if floor * floor * x.denominator < scaled:
        floor += 1
    result = Q(floor, den)
    if result * result < x:
        raise ArithmeticError('sqrt upper')
    return result

@dataclass(frozen=True)
class SparseCoulombBox:
    rows: tuple
    inverse_mass: tuple
    free_velocity: tuple
    initial_velocity: tuple
    gaps: tuple
    step: Q
    mu: tuple
    compliance: tuple

    @classmethod
    def make(cls, rows, inverse_mass, free_velocity, gaps, step, mu, *, compliance=None, initial_velocity=None):
        rows = tuple((tuple(((index(k), _q(v)) for (k, v) in row)) for row in rows))
        im = tuple((_q(v) for v in inverse_mass))
        vf = tuple((_q(v) for v in free_velocity))
        vi = tuple((_q(v) for v in (initial_velocity if initial_velocity is not None else [0] * len(im))))
        g = tuple((tuple((_q(v) for v in pair)) for pair in gaps))
        h = _q(step)
        mus = tuple((_q(v) for v in mu))
        d = tuple((_q(v) for v in (compliance if compliance is not None else [0] * len(g))))
        if not g or not im or len(rows) != 3 * len(g) or (len(vf) != len(im)) or (len(vi) != len(im)) or (len(mus) != len(g)) or (len(d) != len(g)) or (h <= 0) or any((v <= 0 for v in im)) or any((v < 0 for v in mus + d)) or any((len(pair) != 2 or pair[0] > pair[1] for pair in g)) or any((len(set((k for (k, v) in row))) != len(row) or any((k < 0 or k >= len(im) for (k, v) in row)) for row in rows)):
            raise ValueError('invalid sparse contact box')
        return cls(rows, im, vf, vi, g, h, mus, d)

    def digest(self):
        return sha256(json.dumps(rational_json(self.__dict__), sort_keys=True).encode()).hexdigest()

    def validate(self):
        """Reestablish theorem assumptions, including for a replaced dataclass."""
        n = len(self.gaps)
        nv = len(self.inverse_mass)
        if not n or not nv or len(self.rows) != 3 * n or (len(self.free_velocity) != nv) or (len(self.initial_velocity) != nv) or (len(self.mu) != n) or (len(self.compliance) != n) or (self.step <= 0) or any((v <= 0 for v in self.inverse_mass)) or any((v < 0 for v in self.mu + self.compliance)) or any((len(g) != 2 or g[0] > g[1] for g in self.gaps)) or any((len(set((k for (k, v) in row))) != len(row) or any((not isinstance(k, int) or k < 0 or k >= nv for (k, v) in row)) for row in self.rows)):
            raise ValueError('invalid theorem assumptions')
        scalars = list(self.inverse_mass + self.free_velocity + self.initial_velocity + self.mu + self.compliance) + list((v for pair in self.gaps for v in pair)) + [self.step] + [v for row in self.rows for (k, v) in row]
        if any((not isinstance(v, Q) for v in scalars)):
            raise TypeError('construct rational models with SparseCoulombBox.make')

def energy_enclosure(pb, escape):
    start = perf_counter()
    pb.validate()
    e = tuple((_q(v) for v in escape))
    if len(e) != len(pb.inverse_mass):
        raise ValueError('escape dimension')
    projections = [sum((v * e[k] for (k, v) in row), Q(0)) for row in pb.rows]
    margins = [projections[3 * i] - mu * (abs(projections[3 * i + 1]) + abs(projections[3 * i + 2])) for (i, mu) in enumerate(pb.mu)]
    alpha = min(margins)
    out = dict(status='UNCERTAIN', all_solutions_covered=False, existence='UNPROVED', bounds={}, input_sha256=pb.digest(), escape=e, margins=margins, alpha=alpha)
    if alpha <= 0:
        out.update(reason='NO_STRICT_DISSIPATIVE_SEPARATOR', seconds=perf_counter() - start)
        return out
    B = max([Q(0)] + [-lo / pb.step for (lo, hi) in pb.gaps])
    k = B / alpha
    e2 = sum((ei * ei / im for (ei, im) in zip(e, pb.inverse_mass)), Q(0))
    z = [v - k * ei for (v, ei) in zip(pb.free_velocity, e)]
    z2 = sum((zi * zi / im for (zi, im) in zip(z, pb.inverse_mass)), Q(0))
    en = sqrt_upper(e2)
    zn = sqrt_upper(z2)
    cap = en * zn / alpha
    center = [(v + k * ei) / 2 for (v, ei) in zip(pb.free_velocity, e)]
    radius = zn / 2
    bounds = {}
    for (j, (c, im)) in enumerate(zip(center, pb.inverse_mass)):
        delta = radius * sqrt_upper(im)
        lo = c - delta
        hi = c + delta
        bounds[f'velocity:{j}'] = dict(lo=lo, hi=hi)
        bounds[f'displacement:{j}'] = dict(lo=pb.step * lo, hi=pb.step * hi)
        bounds[f'contact_force:{j}'] = dict(lo=(lo - pb.free_velocity[j]) / im / pb.step, hi=(hi - pb.free_velocity[j]) / im / pb.step)
        bounds[f'net_force:{j}'] = dict(lo=(lo - pb.initial_velocity[j]) / im / pb.step, hi=(hi - pb.initial_velocity[j]) / im / pb.step)
    for i in range(len(pb.gaps)):
        bounds[f'normal_force:{i}'] = dict(lo=Q(0), hi=cap / pb.step)
    out.update(status='OUTER_ENCLOSURE', all_solutions_covered=True, reason='GLOBAL_DISSIPATIVE_ENERGY_THEOREM', bounds=bounds, B=B, k=k, escape_norm_squared=e2, shifted_velocity_norm_squared=z2, escape_norm_upper=en, shifted_velocity_norm_upper=zn, total_normal_impulse_upper=cap, velocity_center=center, velocity_radius_mass_norm_upper=radius, seconds=perf_counter() - start)
    return out

def verify_energy_enclosure(pb, answer):
    """Solver-free rational replay of separator, square bounds and all outputs."""
    try:
        pb.validate()
        if answer['input_sha256'] != pb.digest() or answer['status'] != 'OUTER_ENCLOSURE' or (not answer['all_solutions_covered']) or (answer['existence'] != 'UNPROVED'):
            return False
        e = [_q(v) for v in answer['escape']]
        if len(e) != len(pb.inverse_mass):
            return False
        proj = [sum((v * e[j] for (j, v) in row), Q(0)) for row in pb.rows]
        margins = [proj[3 * i] - mu * (abs(proj[3 * i + 1]) + abs(proj[3 * i + 2])) for (i, mu) in enumerate(pb.mu)]
        a = min(margins)
        if a <= 0 or a != _q(answer['alpha']) or margins != [_q(v) for v in answer['margins']]:
            return False
        B = max([Q(0)] + [-lo / pb.step for (lo, hi) in pb.gaps])
        k = B / a
        e2 = sum((v * v / im for (v, im) in zip(e, pb.inverse_mass)), Q(0))
        z = [v - k * ei for (v, ei) in zip(pb.free_velocity, e)]
        z2 = sum((v * v / im for (v, im) in zip(z, pb.inverse_mass)), Q(0))
        en = _q(answer['escape_norm_upper'])
        zn = _q(answer['shifted_velocity_norm_upper'])
        if en < 0 or zn < 0 or en * en < e2 or (zn * zn < z2):
            return False
        for (name, value) in [('B', B), ('k', k), ('escape_norm_squared', e2), ('shifted_velocity_norm_squared', z2), ('total_normal_impulse_upper', en * zn / a), ('velocity_radius_mass_norm_upper', zn / 2)]:
            if _q(answer[name]) != value:
                return False
        centers = [(v + k * ei) / 2 for (v, ei) in zip(pb.free_velocity, e)]
        if [_q(v) for v in answer['velocity_center']] != centers:
            return False
        expected_names = {f'{name}:{i}' for name in ['velocity', 'displacement', 'contact_force', 'net_force'] for i in range(len(centers))} | {f'normal_force:{i}' for i in range(len(pb.gaps))}
        if set(answer['bounds']) != expected_names:
            return False
        for (i, (c, im)) in enumerate(zip(centers, pb.inverse_mass)):
            vb = answer['bounds'][f'velocity:{i}']
            lo = _q(vb['lo'])
            hi = _q(vb['hi'])
            if not lo <= c <= hi or (c - lo) ** 2 < (zn / 2) ** 2 * im or (hi - c) ** 2 < (zn / 2) ** 2 * im:
                return False
            for (name, pair) in [('displacement', (pb.step * lo, pb.step * hi)), ('contact_force', ((lo - pb.free_velocity[i]) / im / pb.step, (hi - pb.free_velocity[i]) / im / pb.step)), ('net_force', ((lo - pb.initial_velocity[i]) / im / pb.step, (hi - pb.initial_velocity[i]) / im / pb.step))]:
                bound = answer['bounds'][f'{name}:{i}']
                if (_q(bound['lo']), _q(bound['hi'])) != pair:
                    return False
        return all(((_q(answer['bounds'][f'normal_force:{i}']['lo']), _q(answer['bounds'][f'normal_force:{i}']['hi'])) == (0, en * zn / a / pb.step) for i in range(len(pb.gaps))))
    except (KeyError, TypeError, ValueError, ArithmeticError):
        return False
