"""Fixed-direction swept separation over a certified affine solution family.

For each latent cell, endpoint vertices are affine in the SAME latent variable.
The prescribed interpolated path is multiaffine in time, latent value and
primitive barycentric weights. A scalar separation test has a corner license;
Euclidean distance has no corner license here. UNKNOWN never accepts a step.
"""
from itertools import product
from fractions import Fraction as Q
from .latent_contact import exact, digest, verify_family

def _tensor(x):
    if len(x) != 2 or any((len(t) != 4 for t in x)) or any((len(p) != 3 for t in x for p in t)):
        raise ValueError('2x4x3 endpoint coefficients')
    return tuple((tuple((tuple(map(exact, p)) for p in t)) for t in x))

def _input(problem, family, geometry, complete):
    if complete is not True or not geometry or (not verify_family(problem, family)):
        raise ValueError('complete family and nonempty complete declared pair list required')
    pairs = []
    for pair in geometry:
        if pair['kind'] not in ('PT', 'EE') or len(pair['cells']) != len(family['cells']):
            raise ValueError('primitive/cell coverage')
        if not isinstance(pair['id'], str) or not pair['id']:
            raise ValueError('pair identity')
        skin = exact(pair.get('skin', 0))
        if skin < 0:
            raise ValueError('nonnegative combined skin')
        cells = []
        for g in pair['cells']:
            cells.append(dict(base=_tensor(g['base']), slope=_tensor(g['slope'])))
        pairs.append(dict(id=pair['id'], kind=pair['kind'], skin=skin, cells=cells))
    if len({p['id'] for p in pairs}) != len(pairs):
        raise ValueError('unique pair identity')
    return dict(model_sha256=digest(problem), family_sha256=digest(family), geometry=pairs, complete=True)

def _projected_values(g, domain, kind, n):
    groups = ((0,), (1, 2, 3)) if kind == 'PT' else ((0, 1), (2, 3))
    return [sum((n[k] * (g['base'][t][i][k] - g['base'][t][j][k] + z * (g['slope'][t][i][k] - g['slope'][t][j][k])) for k in range(3))) for (t, z, i, j) in product((0, 1), domain, *groups)]

def _separates(g, domain, kind, skin, n):
    if len(n) != 3 or not any(n):
        return None
    vals = _projected_values(g, domain, kind, n)
    margin = min(vals)
    norm2 = sum((x * x for x in n))
    return margin if margin > 0 and margin * margin > skin * skin * norm2 else None

def certify_sweep(problem, family, geometry, *, complete=False, directions=None):
    try:
        pb = _input(problem, family, geometry, complete)
    except (ValueError, TypeError, KeyError):
        return dict(status='UNKNOWN', complete=False, reason='UNLICENSED_OR_INCOMPLETE_INPUT')
    dirs = [tuple(map(exact, n)) for n in directions] if directions is not None else [tuple((Q(s if k == a else 0) for k in range(3))) for a in range(3) for s in (-1, 1)]
    certs = []
    for pair in pb['geometry']:
        for (i, (g, c)) in enumerate(zip(pair['cells'], family['cells'])):
            domain = tuple(map(exact, (c['lo'], c['hi'])))
            found = None
            for n in dirs:
                m = _separates(g, domain, pair['kind'], pair['skin'], n)
                if m is not None:
                    found = dict(pair_id=pair['id'], cell=i, direction=n, margin=m)
                    break
            if found is None:
                return dict(status='UNKNOWN', complete=False, reason='NO_FIXED_SEPARATING_DIRECTION')
            certs.append(found)
    return dict(schema='threefold-family-sweep/v1', status='SAFE_PRESCRIBED_LINEAR_PATH_FAMILY', complete=True, input_sha256=digest(pb), certificates=certs, physical_status='UNKNOWN', coverage='DECLARED_PAIR_LIST_ONLY')

def verify_sweep(problem, family, geometry, answer, *, complete=False):
    try:
        pb = _input(problem, family, geometry, complete)
        if answer['schema'] != 'threefold-family-sweep/v1' or answer['status'] != 'SAFE_PRESCRIBED_LINEAR_PATH_FAMILY' or answer['complete'] is not True or (answer['physical_status'] != 'UNKNOWN') or (answer['coverage'] != 'DECLARED_PAIR_LIST_ONLY') or (answer['input_sha256'] != digest(pb)):
            return False
        pairs = {p['id']: p for p in pb['geometry']}
        certs = answer['certificates']
        expected = {(p['id'], i) for p in pb['geometry'] for i in range(len(family['cells']))}
        if not expected or len(certs) != len(expected):
            return False
        seen = set()
        for c in certs:
            key = (c['pair_id'], c['cell'])
            if type(c['cell']) is not int or key not in expected or key in seen:
                return False
            p = pairs[key[0]]
            cell = family['cells'][key[1]]
            n = tuple(map(exact, c['direction']))
            m = _separates(p['cells'][key[1]], tuple(map(exact, (cell['lo'], cell['hi']))), p['kind'], p['skin'], n)
            if m is None or m != exact(c['margin']):
                return False
            seen.add(key)
        return seen == expected
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False
