"""Exact shared-velocity elimination for compliant parallel normal contacts.

Only J_i=1, one translational DOF and strictly positive explicit compliance.
The monotonicity license and unique scalar root cover the entire gap box.
Bounds are sharp marginal projections; joint states still share one velocity.
"""
from fractions import Fraction as Q
from time import perf_counter
from .gap_box import dot, rational_json

def _license(pb):
    if pb.law != 'normal' or len(pb.inverse_mass) != 1 or any((row != (Q(1),) for row in pb.J)) or (not pb.gaps) or (len(pb.J) != len(pb.gaps)) or (len(pb.compliance) != len(pb.gaps)) or (pb.step <= 0) or (pb.inverse_mass[0] <= 0) or (len(pb.free_velocity) != 1) or (len(pb.initial_velocity) != 1) or any((len(g) != 2 or g[0] > g[1] for g in pb.gaps)) or any((d <= 0 for d in pb.compliance)):
        raise ValueError('SHARED_VELOCITY_LICENSE: one DOF, J_i=1, d_i>0')

def scalar_point(pb, gaps):
    """Sort thresholds once; exact water filling, no numerical root tolerance."""
    _license(pb)
    if len(gaps) != len(pb.gaps):
        raise ValueError('gap dimension')
    thresholds = [-Q(g) / pb.step for g in gaps]
    order = sorted(range(len(gaps)), key=lambda i: thresholds[i], reverse=True)
    s = Q(0)
    t = Q(0)
    im = pb.inverse_mass[0]
    v = pb.free_velocity[0]
    for i in order:
        if v >= thresholds[i]:
            break
        s += 1 / pb.compliance[i]
        t += thresholds[i] / pb.compliance[i]
        v = (pb.free_velocity[0] + im * t) / (1 + im * s)
    p = [max(Q(0), (a - v) / d) for (a, d) in zip(thresholds, pb.compliance)]
    if v != pb.free_velocity[0] + im * sum(p):
        raise ArithmeticError('water-filling root')
    return p + list(gaps)

def _point_ok(pb, x):
    n = len(pb.gaps)
    if len(x) != 2 * n:
        return False
    (p, g) = (x[:n], x[n:])
    if any((not lo <= gi <= hi for (gi, (lo, hi)) in zip(g, pb.gaps))):
        return False
    v = pb.free_velocity[0] + pb.inverse_mass[0] * sum(p)
    w = [v + gi / pb.step + d * pi for (gi, d, pi) in zip(g, pb.compliance, p)]
    return all((pi >= 0 and wi >= 0 and (pi * wi == 0) for (pi, wi) in zip(p, w)))

def shared_gap_enclosure(pb):
    start = perf_counter()
    _license(pb)
    low = [a for (a, b) in pb.gaps]
    high = [b for (a, b) in pb.gaps]
    points = {}

    def save(g):
        key = tuple(g)
        if key not in points:
            points[key] = scalar_point(pb, g)
        return points[key]
    vmin = save(high)
    vmax = save(low)
    bounds = {}
    for (name, (c, off)) in pb.readouts().items():
        if name.startswith('normal_force:'):
            i = int(name.split(':')[1])
            gl = low.copy()
            gh = high.copy()
            gl[i] = high[i]
            gh[i] = low[i]
            xlo = save(gl)
            xhi = save(gh)
        else:
            xlo = vmin
            xhi = vmax
        bounds[name] = {side: off + dot(c, x) for (side, x) in [('lo', xlo), ('hi', xhi)]}
        bounds[name].update(lo_witness=xlo, hi_witness=xhi)
    cuts = sorted(set((-g / pb.step for pair in pb.gaps for g in pair)))
    out = dict(status='ENTYDIG' if all((b['lo'] == b['hi'] for b in bounds.values())) else "MULTIPLE", complete=True, input_sha256=pb.digest(), bounds=bounds, representation='shared_velocity_projection', velocity_cuts=cuts, projected_cells=len(cuts) + 1, point_queries=len(points))
    out['seconds'] = perf_counter() - start
    return out

def verify_shared_gap(pb, answer):
    """Replay license, prescribed extreme-gap witnesses and exact KKT only."""
    try:
        _license(pb)
        if not answer['complete'] or answer['input_sha256'] != pb.digest():
            return False
        cuts = sorted(set((-g / pb.step for pair in pb.gaps for g in pair)))
        if [Q(g) for g in answer['velocity_cuts']] != cuts or answer['projected_cells'] != len(cuts) + 1:
            return False
        if set(answer['bounds']) != set(pb.readouts()):
            return False
        low = [a for (a, b) in pb.gaps]
        high = [b for (a, b) in pb.gaps]
        for (name, (c, off)) in pb.readouts().items():
            for side in ['lo', 'hi']:
                g = (high if side == 'lo' else low).copy()
                if name.startswith('normal_force:'):
                    i = int(name.split(':')[1])
                    g = (low if side == 'lo' else high).copy()
                    g[i] = high[i] if side == 'lo' else low[i]
                b = answer['bounds'][name]
                x = [Q(t) for t in b[side + '_witness']]
                if x[len(pb.gaps):] != g or not _point_ok(pb, x) or off + dot(c, x) != Q(b[side]):
                    return False
        status = 'ENTYDIG' if all((Q(b['lo']) == Q(b['hi']) for b in answer['bounds'].values())) else "MULTIPLE"
        return answer['status'] == status
    except (ValueError, TypeError, KeyError, ArithmeticError):
        return False
