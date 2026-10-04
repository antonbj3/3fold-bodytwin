"""Rational virtual-work budgets with global circular-Coulomb coverage.

The caller supplies virtual fields, not presumed physical solutions. This
module never optimizes or samples slip directions. Every field is checked on
every contact, including contacts across the requested cluster boundary.
"""
from fractions import Fraction as Q
from time import perf_counter
from .observable_cert import _q
from .coulomb_outer import sqrt_upper, verify_energy_enclosure

def _field(pb, values):
    e = tuple((_q(v) for v in values))
    if len(e) != len(pb.inverse_mass):
        raise ValueError('virtual field dimension')
    proj = [sum((v * e[j] for (j, v) in row), Q(0)) for row in pb.rows]
    margins = tuple((proj[3 * i] - mu * (abs(proj[3 * i + 1]) + abs(proj[3 * i + 2])) for (i, mu) in enumerate(pb.mu)))
    return (e, margins)

def _mass_dot(pb, x, y):
    return sum((a * b / im for (a, b, im) in zip(x, y, pb.inverse_mass)), Q(0))

def _derive(pb, seed, energy_fields, budget_fields, clusters):
    pb.validate()
    if not verify_energy_enclosure(pb, seed):
        raise ValueError('invalid global reserve')
    B = tuple((max(Q(0), -lo / pb.step) for (lo, hi) in pb.gaps))
    balls = [(tuple((_q(x) for x in seed['velocity_center'])), _q(seed['velocity_radius_mass_norm_upper']))]
    energies = []
    for values in energy_fields:
        (e, a) = _field(pb, values)
        if any((ai < b for (ai, b) in zip(a, B))):
            raise ValueError('field misses contact gap budget')
        z = tuple((v - x for (v, x) in zip(pb.free_velocity, e)))
        c = tuple(((v + x) / 2 for (v, x) in zip(pb.free_velocity, e)))
        R = sqrt_upper(_mass_dot(pb, z, z)) / 2
        balls.append((c, R))
        energies.append(dict(field=e, margins=a, center=c, radius=R))
    fields = [seed['escape']] + [x['field'] for x in energies] + list(budget_fields)
    cap = [_q(seed['total_normal_impulse_upper'])] * len(pb.gaps)
    budgets = []
    for values in fields:
        (f, a) = _field(pb, values)
        if any((ai < 0 for ai in a)):
            raise ValueError('negative exterior contact margin')
        norm = sqrt_upper(_mass_dot(pb, f, f))
        supports = [_mass_dot(pb, f, tuple((ci - vf for (ci, vf) in zip(c, pb.free_velocity)))) + R * norm for (c, R) in balls]
        S = min(supports)
        if S < 0:
            raise ValueError('inconsistent nonnegative budget')
        budgets.append(dict(field=f, margins=a, norm_upper=norm, supports=supports, upper=S))
        cap = [min(old, S / ai) if ai > 0 else old for (old, ai) in zip(cap, a)]
    bounds = {}
    for (j, im) in enumerate(pb.inverse_mass):
        delta = sqrt_upper(im)
        lo = max((c[j] - R * delta for (c, R) in balls))
        hi = min((c[j] + R * delta for (c, R) in balls))
        if lo > hi:
            raise ValueError('empty ball intersection projection')
        bounds[f'velocity:{j}'] = dict(lo=lo, hi=hi)
        bounds[f'displacement:{j}'] = dict(lo=pb.step * lo, hi=pb.step * hi)
        bounds[f'contact_force:{j}'] = dict(lo=(lo - pb.free_velocity[j]) / im / pb.step, hi=(hi - pb.free_velocity[j]) / im / pb.step)
        bounds[f'net_force:{j}'] = dict(lo=(lo - pb.initial_velocity[j]) / im / pb.step, hi=(hi - pb.initial_velocity[j]) / im / pb.step)
    for (i, value) in enumerate(cap):
        bounds[f'normal_force:{i}'] = dict(lo=Q(0), hi=value / pb.step)
    cluster_bounds = {}
    normalized = []
    for (label, indices) in clusters:
        ids = tuple(indices)
        if not isinstance(label, str) or not label or label in cluster_bounds or (not ids) or (len(set(ids)) != len(ids)) or any((type(i) is not int or i < 0 or i >= len(cap) for i in ids)):
            raise ValueError('invalid cluster contact indices')
        upper = sum((cap[i] for i in ids), Q(0))
        for b in budgets:
            alpha = min((b['margins'][i] for i in ids))
            if alpha > 0:
                upper = min(upper, b['upper'] / alpha)
        cluster_bounds[label] = dict(lo=Q(0), hi=upper / pb.step, unit='N')
        normalized.append((label, ids))
    return dict(status='OUTER_ENCLOSURE', all_solutions_covered=True, existence='UNPROVED', sharp=False, input_sha256=pb.digest(), global_reserve=seed, energy_fields=[x['field'] for x in energies], budget_fields=[tuple((_q(v) for v in f)) for f in budget_fields], clusters=normalized, contact_gap_budgets=B, energy_proofs=energies, budget_proofs=budgets, bounds=bounds, cluster_force_bounds=cluster_bounds)

def virtual_budget_enclosure(pb, global_reserve, *, energy_fields=(), budget_fields=(), clusters=()):
    """Intersect sound global balls and globally valid local work budgets.

Clusters specify sets of contact indices; their force is the sum of normal
force magnitudes, not a signed resultant. Returned marginal projections do
not prove simultaneous attainment or existence for any/all gaps.
"""
    start = perf_counter()
    energy_fields = tuple((tuple(x) for x in energy_fields))
    budget_fields = tuple((tuple(x) for x in budget_fields))
    answer = _derive(pb, global_reserve, energy_fields, budget_fields, tuple(clusters))
    answer['seconds'] = perf_counter() - start
    return answer

def verify_virtual_budget(pb, answer):
    """Solver-free replay; all proof values, status and outputs are bound."""
    from .gap_box import rational_json
    try:
        expected = _derive(pb, answer['global_reserve'], answer['energy_fields'], answer['budget_fields'], answer['clusters'])
        actual = {k: v for (k, v) in answer.items() if k != 'seconds'}
        return rational_json(expected) == rational_json(actual)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False
