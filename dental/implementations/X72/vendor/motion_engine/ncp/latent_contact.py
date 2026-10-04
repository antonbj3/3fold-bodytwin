"""Exact one-latent-variable normal-LCP family, opt-in and bounded in size.

This is conventional parametric complementarity with a solver-free certificate.
SPD implies a unique impulse at each latent value, not a constant family output.
Physical input bounds, frozen contact frames and the constitutive law are caller
preconditions. This module neither derives scan bounds nor certifies Coulomb.
"""
from fractions import Fraction as Q
from itertools import product
from hashlib import sha256
import json, math

def exact(x):
    if isinstance(x, bool) or (isinstance(x, float) and (not math.isfinite(x))):
        raise ValueError('finite rational input required')
    return Q(x)

def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x

def digest(x):
    return sha256(json.dumps(encode(x), sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()

def _spd(G):
    n = len(G)
    if any((G[i][j] != G[j][i] for i in range(n) for j in range(n))):
        return False
    L = [[Q(0)] * n for _ in range(n)]
    d = []
    for i in range(n):
        L[i][i] = Q(1)
        for j in range(i):
            L[i][j] = (G[i][j] - sum((L[i][k] * d[k] * L[j][k] for k in range(j)))) / d[j]
        v = G[i][i] - sum((L[i][k] ** 2 * d[k] for k in range(i)))
        if v <= 0:
            return False
        d.append(v)
    return True

def make_family(G, q0, q1, box, *, step, source_id, source_sha256, unit_system='SI'):
    n = len(q0)
    if n < 1 or n > 8 or len(q1) != n or (len(G) != n) or any((len(r) != n for r in G)):
        raise ValueError('nonempty square family, dimension budget 8')
    G = tuple((tuple((exact(v) for v in r)) for r in G))
    q0 = tuple(map(exact, q0))
    q1 = tuple(map(exact, q1))
    if len(box) != 2:
        raise ValueError('latent box')
    box = tuple(map(exact, box))
    step = exact(step)
    if box[0] > box[1] or step <= 0:
        raise ValueError('ordered box, positive step')
    if not isinstance(source_id, str) or not source_id.strip() or (not isinstance(source_sha256, str)) or (len(source_sha256) != 64) or any((x not in '0123456789abcdef' for x in source_sha256)):
        raise ValueError('bound source identity required')
    if unit_system not in ('SI', 'EXTERNAL_UNSPECIFIED'):
        raise ValueError('explicit unit system required')
    if not _spd(G):
        raise ValueError('SPD_NORMAL_LCP_LICENSE')
    return dict(G=G, q0=q0, q1=q1, box=box, step=step, source_id=source_id, source_sha256=source_sha256, unit_system=unit_system, law='FRICTIONLESS_NORMAL', physical_status='UNKNOWN')

def _model(pb):
    if pb.get('law') != 'FRICTIONLESS_NORMAL' or pb.get('physical_status') != 'UNKNOWN':
        raise ValueError('model scope')
    return make_family(**{k: pb[k] for k in ('G', 'q0', 'q1', 'box', 'step', 'source_id', 'source_sha256', 'unit_system')})

def _solve(A, b):
    n = len(b)
    rows = [list(r) + [v] for (r, v) in zip(A, b)]
    for k in range(n):
        j = next((i for i in range(k, n) if rows[i][k]))
        (rows[k], rows[j]) = (rows[j], rows[k])
        f = rows[k][k]
        rows[k] = [v / f for v in rows[k]]
        for i in range(n):
            if i != k:
                f = rows[i][k]
                rows[i] = [x - f * y for (x, y) in zip(rows[i], rows[k])]
    return [r[-1] for r in rows]

def _affine_slack(pb, p0, p1):
    n = len(p0)
    G = pb['G']
    w0 = [pb['q0'][i] + sum((G[i][j] * p0[j] for j in range(n))) for i in range(n)]
    w1 = [pb['q1'][i] + sum((G[i][j] * p1[j] for j in range(n))) for i in range(n)]
    return (w0, w1)

def _readouts(pb, cells, readouts):
    n = len(pb['q0'])
    result = {}
    for (name, c) in readouts.items():
        c = tuple(map(exact, c))
        if not isinstance(name, str) or not name or len(c) != n:
            raise ValueError('dimensioned linear readout required')
        values = []
        for (k, cell) in enumerate(cells):
            for z in (cell['lo'], cell['hi']):
                p = [a + b * z for (a, b) in zip(cell['p0'], cell['p1'])]
                values.append((sum((x * y for (x, y) in zip(c, p))) / pb['step'], k, z))
        lo = min(values)
        hi = max(values)
        result[name] = dict(lo=lo[0], hi=hi[0], lo_witness=dict(cell=lo[1], z=lo[2]), hi_witness=dict(cell=hi[1], z=hi[2]))
    return result

def default_readouts(pb):
    n = len(pb['q0'])
    out = {f'normal_force:{i}': tuple((Q(int(j == i)) for j in range(n))) for i in range(n)}
    out['total_normal_force'] = tuple((Q(1) for _ in range(n)))
    if n >= 2:
        out['force_difference:0:1'] = tuple((Q(1 if j == 0 else -1 if j == 1 else 0) for j in range(n)))
    return out

def certify_family(problem, *, readouts=None):
    pb = _model(problem)
    n = len(pb['q0'])
    modes = []
    readouts = default_readouts(pb) if readouts is None else {k: tuple(map(exact, v)) for (k, v) in readouts.items()}
    if not readouts:
        raise ValueError('nonempty readout set')
    for bits in product((0, 1), repeat=n):
        ids = [i for i in range(n) if bits[i]]
        p0 = [Q(0)] * n
        p1 = p0.copy()
        if ids:
            A = [[pb['G'][i][j] for j in ids] for i in ids]
            a = _solve(A, [-pb['q0'][i] for i in ids])
            b = _solve(A, [-pb['q1'][i] for i in ids])
            for (i, x, y) in zip(ids, a, b):
                p0[i] = x
                p1[i] = y
        (w0, w1) = _affine_slack(pb, p0, p1)
        (lo, hi) = pb['box']
        for (a, b) in zip(p0 + w0, p1 + w1):
            if b > 0:
                lo = max(lo, -a / b)
            elif b < 0:
                hi = min(hi, -a / b)
            elif a < 0:
                (lo, hi) = (Q(1), Q(0))
                break
        if lo <= hi:
            modes.append(dict(lo=lo, hi=hi, p0=p0, p1=p1))
    cuts = sorted(set(pb['box']) | {m[k] for m in modes for k in ('lo', 'hi')})
    cells = []
    domains = list(zip(cuts, cuts[1:])) if len(cuts) > 1 else [(cuts[0], cuts[0])]
    for (lo, hi) in domains:
        mid = (lo + hi) / 2
        m = next((m for m in modes if m['lo'] <= mid <= m['hi']))
        cells.append(dict(lo=lo, hi=hi, p0=m['p0'], p1=m['p1']))
    bounds = _readouts(pb, cells, readouts)
    return dict(schema='threefold-normal-family/v1', status='FAMILY_ENCLOSURE', complete=True, uniqueness='UNIQUE_PER_LATENT_VALUE', physical_status='UNKNOWN', input_sha256=digest(pb), readouts=readouts, cells=cells, bounds=bounds, unit='N' if pb['unit_system'] == 'SI' else 'EXTERNAL_UNSPECIFIED', variation='CONSTANT' if all((b['lo'] == b['hi'] for b in bounds.values())) else 'VARIABLE')

def verify_family(problem, answer, *, readouts=None):
    """Replay cover and affine KKT; never calls a solver or certify_family."""
    try:
        pb = _model(problem)
        n = len(pb['q0'])
        requested = default_readouts(pb) if readouts is None else {k: tuple(map(exact, v)) for (k, v) in readouts.items()}
        if not requested:
            return False
        expected_unit = 'N' if pb['unit_system'] == 'SI' else 'EXTERNAL_UNSPECIFIED'
        if answer['schema'] != 'threefold-normal-family/v1' or answer['status'] != 'FAMILY_ENCLOSURE' or answer['complete'] is not True or (answer['uniqueness'] != 'UNIQUE_PER_LATENT_VALUE') or (answer['physical_status'] != 'UNKNOWN') or (answer['input_sha256'] != digest(pb)) or (answer['unit'] != expected_unit):
            return False
        if {k: tuple(map(exact, v)) for (k, v) in answer['readouts'].items()} != requested:
            return False
        cells = answer['cells']
        if not cells or len(cells) > 2 ** n + 1:
            return False
        rebuilt = []
        last = pb['box'][0]
        for c in cells:
            (lo, hi) = (exact(c['lo']), exact(c['hi']))
            p0 = list(map(exact, c['p0']))
            p1 = list(map(exact, c['p1']))
            if len(p0) != n or len(p1) != n or lo != last or (lo > hi) or (hi > pb['box'][1]):
                return False
            if lo == hi and pb['box'][0] != pb['box'][1]:
                return False
            (w0, w1) = _affine_slack(pb, p0, p1)
            if any((a + b * z < 0 for z in (lo, hi) for (a, b) in zip(p0 + w0, p1 + w1))):
                return False
            for (a, b, c0, d) in zip(p0, p1, w0, w1):
                if lo == hi:
                    if (a + b * lo) * (c0 + d * lo) != 0:
                        return False
                elif a * c0 != 0 or a * d + b * c0 != 0 or b * d != 0:
                    return False
            rebuilt.append(dict(lo=lo, hi=hi, p0=p0, p1=p1))
            last = hi
        if last != pb['box'][1]:
            return False
        bounds = _readouts(pb, rebuilt, requested)
        if encode(bounds) != encode(answer['bounds']):
            return False
        variation = 'CONSTANT' if all((b['lo'] == b['hi'] for b in bounds.values())) else 'VARIABLE'
        return answer['variation'] == variation
    except (ValueError, TypeError, KeyError, ArithmeticError, StopIteration, IndexError):
        return False

def port_family(port, G, free_velocity, *, step):
    """PORT v1 affine gap model, exactly one shared latent cause, SI transport.

    Gap bands must be hard MODEL_BOUND; statistical sigma alone is not enough.
    This creates a new normal-law query and does not promote upstream branches.
    """
    if port.get('schema') != 'threefold-contact-port/v1':
        raise ValueError('PORT v1 required')
    lat = port.get('latent_box', [])
    contacts = port.get('contacts', [])
    affine = port.get('affine_gaps', [])
    if len(lat) != 1 or not contacts or len(contacts) != len(affine) or (len(free_velocity) != len(contacts)):
        raise ValueError('one latent cause and complete affine contacts required')
    (latent, lo, hi) = lat[0]
    table = {a['candidate_id']: a for a in affine}
    if len(table) != len(affine) or len({c['id'] for c in contacts}) != len(contacts):
        raise ValueError('unique contact identities')
    h = exact(step)
    if h <= 0:
        raise ValueError('positive step')
    q0 = []
    q1 = []
    for (c, v) in zip(contacts, free_velocity):
        g = c['gap']
        a = table[c['id']]
        if g['unit'] != 'm' or g['assurance'] != 'MODEL_BOUND' or g.get('blocked', False):
            raise ValueError('hard SI gap bound required')
        coeff = a['coefficients']
        if len(coeff) != 1 or coeff[0][0] != latent:
            raise ValueError('shared latent identity required')
        g0 = exact(a['constant_m'])
        g1 = exact(coeff[0][1])
        ext = [g0 + g1 * exact(z) for z in (lo, hi)]
        if min(ext) < exact(g['lower']) or max(ext) > exact(g['upper']):
            raise ValueError('affine gap escapes declared band')
        q0.append(exact(v) + g0 / h)
        q1.append(g1 / h)
    return make_family(G, q0, q1, (lo, hi), step=h, source_id=port['producer'] + ':PORT-v1', source_sha256=digest(port))
