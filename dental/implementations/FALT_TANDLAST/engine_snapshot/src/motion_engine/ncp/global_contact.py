"""Exact sufficient global verdict for a freely rotating non-spinning rigid step.

Parallel vertical normals, contact points (x,y,-height), unit translational
mass, diagonal positive rotational mass, zero initial gaps, implicit Coulomb
friction. This is an optional exact-rational guard, not a numerical solver.
UNKNOWN is returned outside the proved regime. No finite compliance or
material selector is implicit in this API.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256
import json

@dataclass(frozen=True)
class GlobalContactVerdict:
    status: str
    reason: str
    velocity: tuple = ()
    impulse_witness: tuple = ()
    input_sha256: str = ''
    scope: str = 'rigid implicit Coulomb; global velocity, not individual impulses'

def _q(x):
    if isinstance(x, bool) or not isinstance(x, (int, str, Q)):
        raise ValueError('exact rational inputs required; floats are not certificates')
    return Q(x)

def global_rigid_translation_verdict(contact_xy, height, inverse_mass, free_velocity, mu):
    """Certify all contact modes, including opening, sticking and other roots.

    Input ordering is (vx,vy,vz,wx,wy,wz), and each impulse is (p,tx,ty).
    Units must be coherent; length/mass normalization is part of the model.
    Positive certificates carry an exact existence witness. The supported
    free state leads to pure +x translation, with rotations allowed.
    """
    try:
        xy = tuple((tuple((_q(t) for t in row)) for row in contact_xy))
        h = _q(height)
        mi = tuple((_q(t) for t in inverse_mass))
        vf = tuple((_q(t) for t in free_velocity))
        mu = _q(mu)
        if len(xy) < 4 or any((len(row) != 2 for row in xy)):
            raise ValueError('at least four planar contact points required')
        if len(mi) != 6 or len(vf) != 6 or any((t <= 0 for t in mi)):
            raise ValueError('six positive inverse masses and six velocities required')
        if mi[:3] != (1, 1, 1) or h < 0 or mu < 0:
            raise ValueError('unit translational mass, height>=0 and mu>=0 required')
        P = -vf[2]
        speed = vf[0] - mu * P
        if P <= 0 or speed <= 0:
            raise ValueError('positive normal load and reference slip required')
        if vf != (speed + mu * P, 0, -P, 0, -mi[4] * mu * h * P, 0):
            raise ValueError('free state outside non-spinning regime')
        X = max((abs(x) for (x, y) in xy))
        Y = max((abs(y) for (x, y) in xy))
        if X <= 0 or Y <= 0 or any(((x, y) not in xy for x in (-X, X) for y in (-Y, Y))):
            raise ValueError('four enclosing rectangle corners required')
        if sum((x for (x, y) in xy)) != 0 or sum((y for (x, y) in xy)) != 0:
            raise ValueError('uniform existence witness requires centered contacts')
        if min(X, Y) < 2 * mu * h:
            raise ValueError('global no-tilt inequality not certified')
        if 1 / mi[5] <= (mu * h) ** 2 / 4:
            raise ValueError('tangential work form not positive definite')
        n = len(xy)
        imp = tuple((t for _ in xy for t in (P / n, -mu * P / n, Q(0))))
        velocity = (speed, Q(0), Q(0), Q(0), Q(0), Q(0))
        J = []
        for (x, y) in xy:
            J.extend(((0, 0, 1, y, -x, 0), (1, 0, 0, 0, -h, -y), (0, 1, 0, h, 0, x)))
        computed = tuple((vf[k] + mi[k] * sum((row[k] * a for (row, a) in zip(J, imp))) for k in range(6)))
        if computed != velocity:
            raise ValueError('exact existence witness failed')
        encoded = json.dumps([[[str(t) for t in row] for row in xy], str(h), list(map(str, mi)), list(map(str, vf)), str(mu)], separators=(',', ':'))
        return GlobalContactVerdict('UNIQUE', 'global normal work and positive tangential work form', velocity, imp, sha256(encoded.encode()).hexdigest())
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as exc:
        return GlobalContactVerdict('UNKNOWN', str(exc))
