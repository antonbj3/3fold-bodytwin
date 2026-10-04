"""Complete univariate equilibrium sets for EXACT rational polynomials.

Global coverage: a Cauchy tail bound plus a replayable interval partition.
Local existence/uniqueness: exact rational interval Newton, strict inclusion.
No physical/continuum guarantee is inferred from the polynomial certificate.
"""
from fractions import Fraction as Q

def rational(v):
    if isinstance(v, bool) or not isinstance(v, (int, str, Q)):
        raise TypeError('exact int/string/Fraction required; floats are not certified inputs')
    return Q(v)

def add(a, b):
    return (a[0] + b[0], a[1] + b[1])

def neg(a):
    return (-a[1], -a[0])

def sub(a, b):
    return add(a, neg(b))

def mul(a, b):
    p = [x * y for x in a for y in b]
    return (min(p), max(p))

def div(a, b):
    if b[0] <= 0 <= b[1]:
        raise ZeroDivisionError('interval denominator crosses zero')
    return mul(a, (1 / b[1], 1 / b[0]))

def value(c, x):
    r = Q(0)
    for a in reversed(c):
        r = r * x + a
    return r

def enclosure(c, x):
    r = (Q(0), Q(0))
    for a in reversed(c):
        r = add(mul(r, x), (a, a))
    return r

def derivative(c):
    return tuple((i * c[i] for i in range(1, len(c)))) or (Q(0),)

def nonzero(i):
    return i[0] > 0 or i[1] < 0

def encoded(i):
    return [str(x) for x in i]

def decoded(i):
    return tuple((rational(x) for x in i))

def newton_box(c, bracket, tolerance):
    (lo, hi) = bracket
    d = derivative(c)
    while hi - lo > tolerance / 16:
        m = (lo + hi) / 2
        if value(c, m) == 0:
            lo = hi = m
            break
        if value(c, lo) * value(c, m) <= 0:
            hi = m
        else:
            lo = m
    if value(c, lo) == 0:
        hi = lo
    elif value(c, hi) == 0:
        lo = hi
    midpoint = (lo + hi) / 2
    radius = max(hi - lo, tolerance / 16)
    box = (midpoint - radius, midpoint + radius)
    slope = enclosure(d, box)
    if not nonzero(slope):
        return None
    n = sub((midpoint, midpoint), div((value(c, midpoint),) * 2, slope))
    if not box[0] < n[0] <= n[1] < box[1]:
        return None
    return {'box': encoded(box), 'bracket': encoded((lo, hi)), 'newton': encoded(n), 'derivative': encoded(slope)}

def certify(coefficients, *, tolerance='1/1000000000000', max_tiles=20000, max_depth=100):
    c = tuple((rational(x) for x in coefficients))
    while len(c) > 1 and c[-1] == 0:
        c = c[:-1]
    tol = rational(tolerance)
    if tol <= 0 or max_tiles < 1 or max_depth < 1:
        raise ValueError('invalid budget')
    out = {'schema': 1, 'scope': 'EXACT_RATIONAL_UNIVARIATE_POLYNOMIAL', 'coefficients': [str(x) for x in c], 'tolerance': str(tol), 'roots': [], 'tiles': [], 'status': 'OSAKER', 'complete': False}
    if c == (Q(0),):
        out['reason'] = 'CONTINUUM_OF_ROOTS'
        return out
    if len(c) == 1:
        out.update(status='GRENMANGD', complete=True, reason='NO_ROOTS', bound='0')
        return out
    bound = 1 + max((abs(x / c[-1]) for x in c[:-1]))
    out['bound'] = str(bound)
    stack = [(-bound, bound, 0)]
    d = derivative(c)
    processed = 0
    while stack:
        (lo, hi, depth) = stack.pop()
        processed += 1
        tile = {'interval': encoded((lo, hi))}
        fv = enclosure(c, (lo, hi))
        if nonzero(fv):
            tile['kind'] = 'EXCLUDED'
        elif nonzero(enclosure(d, (lo, hi))) and value(c, lo) * value(c, hi) <= 0:
            witness = newton_box(c, (lo, hi), tol)
            if witness is None:
                tile.update(kind='UNKNOWN', reason='NEWTON_NOT_STRICT')
            else:
                box = decoded(witness['box'])
                index = None
                for (j, other) in enumerate(out['roots']):
                    ob = decoded(other['box'])
                    union = (min(box[0], ob[0]), max(box[1], ob[1]))
                    if max(box[0], ob[0]) <= min(box[1], ob[1]) and nonzero(enclosure(d, union)):
                        index = j
                        break
                if index is None:
                    index = len(out['roots'])
                    out['roots'].append(witness)
                tile.update(kind='ROOT', root=index, bracket=witness['bracket'])
        elif depth >= max_depth or processed >= max_tiles:
            tile.update(kind='UNKNOWN', reason='PARTITION_BUDGET_OR_SINGULARITY')
        else:
            m = (lo + hi) / 2
            stack.extend([(m, hi, depth + 1), (lo, m, depth + 1)])
            continue
        out['tiles'].append(tile)
    out['processed_intervals'] = processed
    out['complete'] = all((t['kind'] != 'UNKNOWN' for t in out['tiles']))
    count = len(out['roots'])
    if out['complete']:
        out['status'] = 'ENTYDIG' if count == 1 else 'GRENMANGD'
        out['reason'] = 'COMPLETE_GLOBAL_COVER'
    else:
        out['reason'] = 'REMAINDER_UNRESOLVED'
    return out

def verify(cert, *, coefficients=None):
    """Replay every coverage tile and Newton inequality; never trust complete=True."""
    try:
        if cert.get('schema') != 1 or cert.get('scope') != 'EXACT_RATIONAL_UNIVARIATE_POLYNOMIAL':
            return False
        c = tuple((rational(x) for x in cert['coefficients']))
        if coefficients is not None:
            expected = tuple((rational(x) for x in coefficients))
            while len(expected) > 1 and expected[-1] == 0:
                expected = expected[:-1]
            if c != expected:
                return False
        if not c:
            return False
        if c == (Q(0),):
            return cert['status'] == 'OSAKER' and cert['complete'] is False and (cert['reason'] == 'CONTINUUM_OF_ROOTS') and (cert['roots'] == []) and (cert['tiles'] == [])
        if c[-1] == 0:
            return False
        if len(c) == 1:
            return cert['complete'] is True and cert['status'] == 'GRENMANGD' and (cert['roots'] == []) and (cert['tiles'] == []) and (cert['bound'] == '0')
        bound = 1 + max((abs(x / c[-1]) for x in c[:-1]))
        if rational(cert['bound']) != bound:
            return False
        tol = rational(cert['tolerance'])
        if tol <= 0:
            return False
        d = derivative(c)
        roots = cert['roots']
        for r in roots:
            (box, bracket) = (decoded(r['box']), decoded(r['bracket']))
            if not (box[0] < box[1] and box[1] - box[0] <= tol and (box[0] <= bracket[0] <= bracket[1] <= box[1])):
                return False
            if value(c, bracket[0]) * value(c, bracket[1]) > 0:
                return False
            m = sum(box) / 2
            slope = enclosure(d, box)
            if not nonzero(slope):
                return False
            n = sub((m, m), div((value(c, m),) * 2, slope))
            if decoded(r['newton']) != n or decoded(r['derivative']) != slope:
                return False
            if not box[0] < n[0] <= n[1] < box[1]:
                return False
        for i in range(len(roots)):
            for j in range(i):
                (a, b) = (decoded(roots[i]['box']), decoded(roots[j]['box']))
                if max(a[0], b[0]) <= min(a[1], b[1]):
                    return False
        previous = -bound
        referred = set()
        unknown = False
        for tile in sorted(cert['tiles'], key=lambda t: decoded(t['interval'])[0]):
            (lo, hi) = decoded(tile['interval'])
            if lo != previous or not lo < hi or hi > bound:
                return False
            previous = hi
            if tile['kind'] == 'EXCLUDED':
                if not nonzero(enclosure(c, (lo, hi))):
                    return False
            elif tile['kind'] == 'ROOT':
                index = tile['root']
                if isinstance(index, bool) or not isinstance(index, int) or (not 0 <= index < len(roots)):
                    return False
                if not nonzero(enclosure(d, (lo, hi))) or value(c, lo) * value(c, hi) > 0:
                    return False
                (a, b) = decoded(tile['bracket'])
                (rlo, rhi) = decoded(roots[index]['box'])
                if not lo <= a <= b <= hi or not rlo <= a <= b <= rhi:
                    return False
                if value(c, a) * value(c, b) > 0:
                    return False
                referred.add(index)
            elif tile['kind'] == 'UNKNOWN':
                unknown = True
            else:
                return False
        if previous != bound or referred != set(range(len(roots))):
            return False
        complete = not unknown
        status = 'OSAKER' if unknown else 'ENTYDIG' if len(roots) == 1 else 'GRENMANGD'
        return cert['complete'] is complete and cert['status'] == status
    except (KeyError, TypeError, ValueError, ZeroDivisionError, IndexError):
        return False
