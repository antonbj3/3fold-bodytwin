"""Opt-in SI contact port. Requires field_engine.contact_port_v1 when imported.

Scene.b is velocity, not length. This adapter never treats b as a gap and never
claims that generic NCP output is bounded by its two gap-corner solutions.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as Q
import json
from field_engine.contact_port_v1 import Gap, Unit, Source, Assurance, ContactCandidate, Branch, BranchSet, Witness, Status, canonical, digest, exact as _port_exact

def exact(value):
    """Keep the field input policy while removing fixed-width Rational storage."""
    q = _port_exact(value)
    return Q(int(q.numerator), int(q.denominator))

@dataclass(frozen=True)
class SolverGapInput:
    gap: Gap
    dt_s: Q
    normal_offset_m_s: tuple[Q, Q]
    sigma_m_s: Q | None
    convention: str

    def __post_init__(self):
        if not isinstance(self.gap, Gap) or self.gap.unit != Unit.M:
            raise ValueError('solver requires typed SI gap')
        dt = exact(self.dt_s)
        if dt <= 0:
            raise ValueError('positive dt required')
        g = replace(self.gap, lower=exact(self.gap.lower), upper=exact(self.gap.upper), sigma=None if self.gap.sigma is None else exact(self.gap.sigma))
        object.__setattr__(self, 'gap', g)
        if self.convention == 'RIGID_SPECULATIVE':
            bounds = (max(g.lower, 0) / dt, max(g.upper, 0) / dt)
        elif self.convention == 'CLOTH_SIGNED':
            bounds = (g.lower / dt, g.upper / dt)
        else:
            raise ValueError('unsupported gap convention')
        expected = None if g.sigma is None else g.sigma / dt
        if self.sigma_m_s is not None:
            object.__setattr__(self, 'sigma_m_s', exact(self.sigma_m_s))
        if tuple(map(exact, self.normal_offset_m_s)) != bounds or self.sigma_m_s != expected:
            raise ValueError('gap, offset and sigma disagree')
        object.__setattr__(self, 'dt_s', dt)
        object.__setattr__(self, 'normal_offset_m_s', bounds)

def solver_gap(gap, *, dt_s, convention):
    """Known maps: rigid speculative max(g,0)/dt; cloth signed g/dt.

    Sigma/dt is the pre-clipping scale. For clipped rigid gaps it is not a
    Gaussian standard deviation; the source law and conditioning remain in gap.
    """
    dt = exact(dt_s)
    if dt <= 0:
        raise ValueError('positive dt required')
    g = replace(gap, lower=exact(gap.lower), upper=exact(gap.upper), sigma=None if gap.sigma is None else exact(gap.sigma)).to(Unit.M)
    if convention == 'RIGID_SPECULATIVE':
        (lo, hi) = (max(g.lower, 0), max(g.upper, 0))
    elif convention == 'CLOTH_SIGNED':
        (lo, hi) = (g.lower, g.upper)
    else:
        raise ValueError('explicit gap-to-velocity convention required')
    return SolverGapInput(g, dt, (lo / dt, hi / dt), None if g.sigma is None else g.sigma / dt, convention)

def from_motion_contact(contact, *, source, participants, frame):
    if 'gap' not in contact or 'b' in contact:
        raise ValueError('contact gap [m] required; b is velocity')
    g = exact(contact['gap'])
    return ContactCandidate(str(contact.get('id', 'motion')), participants, frame, Gap(g, g, Unit.M, source, assurance=Assurance.POINT))

def from_cloth_distance(distance_m, *, thickness_m, kind, source):
    d = exact(distance_m)
    h = exact(thickness_m)
    if h < 0:
        raise ValueError('nonnegative thickness required')
    if kind == 'SUPPORT':
        g = d - h / 2
    elif kind in ('VERTEX_TRIANGLE', 'EDGE_EDGE'):
        if d < 0:
            raise ValueError('nonnegative Euclidean self-contact distance required')
        g = d - h
    else:
        raise ValueError('explicit cloth contact kind required')
    return Gap(g, g, Unit.M, source, assurance=Assurance.POINT)

def from_mpm_cells(lower_cells, upper_cells, *, cell_size_m, source):
    dx = exact(cell_size_m)
    if dx <= 0:
        raise ValueError('positive physical cell size required')
    return Gap(exact(lower_cells) * dx, exact(upper_cells) * dx, Unit.M, source, assurance=Assurance.NUMERICAL, condition='represented cell band, not a continuum bound')

def solve_normal_gap(request, *, free_normal_m_s, inverse_mass_per_kg, mu=0):
    """Exact scalar frictionless NCP impulse enclosure (N s).

    License: lambda=max(0,-(v+offset(g))/G) is nonincreasing in g for G>0.
    The enclosure is conditional on the declared input gap event, independent
    of how its sigma was acquired. It is a model result, not a CCD guarantee.
    """
    if not isinstance(request, SolverGapInput) or request.gap.blocked:
        raise ValueError('resolved geometry completion required')
    if exact(mu) != 0:
        raise ValueError('coupled/frictional NCP needs its own coverage oracle')
    G = exact(inverse_mass_per_kg)
    v = exact(free_normal_m_s)
    if G <= 0:
        raise ValueError('positive normal inverse mass required')
    (lo, hi) = request.normal_offset_m_s
    return (max(Q(0), -(v + hi) / G), max(Q(0), -(v + lo) / G))

def point_scene(scene, request, *, contact_index, realized_gap_m):
    """Install an explicit admissible realization in unchanged ncp_ref.Scene.

    Existing normal offset must be zero: composing two bias conventions without
    a declared model is rejected. The interval request remains with the caller.
    """
    import numpy as np
    from dataclasses import replace
    g = exact(realized_gap_m)
    if request.gap.blocked:
        raise ValueError('unknown completion requires a geometry oracle')
    if not request.gap.lower <= g <= request.gap.upper:
        raise ValueError('realization outside gap band')
    if exact(float(scene.dt)) != request.dt_s:
        raise ValueError('Scene dt must match request exactly')
    if type(contact_index) is not int or not 0 <= contact_index < scene.n_c:
        raise ValueError('contact index')
    off = np.zeros(3 * scene.n_c) if scene.b_offset is None else np.asarray(scene.b_offset).copy()
    if off[3 * contact_index] != 0:
        raise ValueError('existing normal offset requires explicit composition')
    realized = max(g, 0) if request.convention == 'RIGID_SPECULATIVE' else g
    off[3 * contact_index] = float(realized / request.dt_s)
    return replace(scene, b_offset=off)

def from_global_verdict(verdict, *, source):
    data = {k: getattr(verdict, k) for k in ('status', 'reason', 'velocity', 'impulse_witness', 'input_sha256', 'scope')}
    w = Witness('global_contact_verdict', canonical(data), digest(data), source)
    bs = (Branch('rigid:0', canonical({'velocity': data['velocity']}), (w,)),) if data['velocity'] else ()
    complete = data['status'] in ('ENTYDIG', 'UNIQUE') and bool(bs)
    return BranchSet(bs, complete, data['scope'], w if complete else None, data['reason'])

def from_reachability(reach, *, source):
    """Keep full terminal identities, multiplicities, all edges and sticky UNKNOWN."""
    data = dict(status=reach.status, layers=[list(layer.items()) for layer in reach.layers], edges=reach.edges, reason=reach.reason)
    w = Witness('branch_set_reachability', canonical(data), digest(data), source)
    bs = tuple((Branch('terminal:' + str(i), canonical({'state': s, 'multiplicity': m}), (w,)) for (i, (s, m)) in enumerate(reach.layers[-1].items())))
    complete = reach.status == 'COMPLETE'
    return BranchSet(bs, complete, source.scope, w if complete else None, reach.reason)

def identical_normal_impulse_difference(port, first, second, *, free_normal_m_s, inverse_mass_per_kg, dt_s):
    """A small licensed relational query; generic affine gaps are not routed."""
    dt = exact(dt_s)
    G = exact(inverse_mass_per_kg)
    exact(free_normal_m_s)
    if dt <= 0 or G <= 0:
        raise ValueError('positive mass/time scale required')
    candidates = {c.id: c for c in port.contacts}
    aff = {a.candidate_id: a for a in port.affine_gaps}
    if first not in candidates or second not in candidates:
        raise ValueError('unbound contacts')
    (a, b) = (aff.get(first), aff.get(second))
    if a is None or b is None:
        raise ValueError('joint gap relation required')
    if a.constant_m != b.constant_m or dict(a.coefficients) != dict(b.coefficients):
        raise ValueError('no license for unequal affine gap laws')
    if candidates[first].gap.blocked or candidates[second].gap.blocked:
        raise ValueError('unknown geometry completion')
    return (Q(0), Q(0))
