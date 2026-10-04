"""Exact frozen-frame gap-box enclosure for normal and PLANAR Coulomb contact.

Finite mode closures cover all solutions of the declared single-step law.
HiGHS discovers LP bases; rational primal/dual and Farkas replay certifies them.
No sampling, caller completeness flag, or numerical rank certifies coverage.
3D Coulomb, changing frames, unbounded outputs or exhausted budgets -> UNCERTAIN.
Certificates apply to the rational model, not to scanner coverage or anatomy.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from time import perf_counter
from .observable_cert import _q, _rref, certify_observable, verify_observable

def dot(a, b):
    return sum((x * y for (x, y) in zip(a, b)), Q(0))

def solve_equalities(E, rhs, n):
    """Exact affine nullspace, using the existing observable certificate kernel."""
    if not E:
        return ([Q(0)] * n, [[Q(i == j) for j in range(n)] for i in range(n)])
    cert = certify_observable(E, rhs, [[Q(0)] * n])
    if not verify_observable(E, rhs, [[Q(0)] * n], cert):
        raise ArithmeticError('observable replay failed')
    if cert['status'] == 'INCOMPATIBLE':
        return (None, cert)
    (R, T, piv) = _rref(E)
    b = [dot(row, rhs) for row in T]
    x = [Q(0)] * n
    for (i, j) in enumerate(piv):
        x[j] = b[i]
    free = [j for j in range(n) if j not in piv]
    N = [[Q(i == j) for j in free] for i in range(n)]
    for (k, j) in enumerate(free):
        for (i, p) in enumerate(piv):
            N[p][k] = -R[i][j]
    return (x, N)

def unique_solution(E, b, n):
    (x, N) = solve_equalities(E, b, n)
    if x is None or (N and N[0]):
        raise ArithmeticError('basis not unique')
    return x

@dataclass(frozen=True)
class GapProblem:
    J: tuple
    inverse_mass: tuple
    free_velocity: tuple
    initial_velocity: tuple
    gaps: tuple
    step: Q
    mu: tuple
    compliance: tuple
    law: str

    @classmethod
    def make(cls, J, inverse_mass, free_velocity, gaps, step, *, law='normal', mu=None, compliance=None, initial_velocity=None):
        if law not in ('normal', 'planar_coulomb'):
            raise ValueError('unsupported law: 3D slip directions need a continuous certified oracle')
        J = tuple((tuple((_q(t) for t in row)) for row in J))
        im = tuple((_q(t) for t in inverse_mass))
        vf = tuple((_q(t) for t in free_velocity))
        vi = tuple((_q(t) for t in (initial_velocity if initial_velocity is not None else [0] * len(im))))
        gaps = tuple((tuple((_q(t) for t in row)) for row in gaps))
        h = _q(step)
        n = len(gaps)
        mus = tuple((_q(t) for t in (mu if mu is not None else [0] * n)))
        ds = tuple((_q(t) for t in (compliance if compliance is not None else [0] * n)))
        stride = 1 if law == 'normal' else 2
        if not n or not im or len(J) != stride * n or any((len(row) != len(im) for row in J)) or (len(vf) != len(im)) or (len(vi) != len(im)) or (len(mus) != n) or (len(ds) != n) or any((len(row) != 2 or row[0] > row[1] for row in gaps)) or (h <= 0) or any((t <= 0 for t in im)) or any((t < 0 for t in mus + ds)):
            raise ValueError('invalid dimensions, box, mass, friction, compliance or step')
        if law == 'normal' and any(mus):
            raise ValueError('normal-only law requires zero friction')
        return cls(J, im, vf, vi, gaps, h, mus, ds, law)

    @property
    def stride(self):
        return 1 if self.law == 'normal' else 2

    @property
    def modes(self):
        return ('open', 'closed') if self.law == 'normal' else ('open', 'stick', 'slide+', 'slide-')

    def digest(self):
        return sha256(json.dumps(self.as_dict(), sort_keys=True).encode()).hexdigest()

    def as_dict(self):
        return {k: rational_json(v) for (k, v) in self.__dict__.items()}

    def system(self, prefix):
        """H x<=r, E x=e for x=(all impulses, shared gap coordinates)."""
        (n, m) = (len(self.gaps), len(self.J))
        d = m + n
        if len(prefix) > n or any((mode not in self.modes for mode in prefix)):
            raise ValueError('invalid mode prefix')
        A = [[sum((self.J[i][k] * self.inverse_mass[k] * self.J[j][k] for k in range(len(self.inverse_mass)))) for j in range(m)] for i in range(m)]
        for (i, c) in enumerate(self.compliance):
            A[self.stride * i][self.stride * i] += c
        b = [dot(row, self.free_velocity) for row in self.J]
        (H, r, E, e) = ([], [], [], [])

        def ineq(row, rhs):
            H.append(row)
            r.append(rhs)

        def equality(row, rhs):
            E.append(row)
            e.append(rhs)

        def unit(k, value=Q(1)):
            row = [Q(0)] * d
            row[k] = value
            return row
        for (i, (lo, hi)) in enumerate(self.gaps):
            ni = self.stride * i
            ineq(unit(m + i), hi)
            ineq(unit(m + i, Q(-1)), -lo)
            ineq(unit(ni, Q(-1)), Q(0))
            wn = A[ni][:] + [Q(0)] * n
            wn[m + i] = 1 / self.step
            ineq([-t for t in wn], b[ni])
            if self.stride == 2:
                ti = ni + 1
                for sign in (-1, 1):
                    row = unit(ti, Q(sign))
                    row[ni] = -self.mu[i]
                    ineq(row, Q(0))
            if i >= len(prefix):
                continue
            mode = prefix[i]
            if mode == 'open':
                equality(unit(ni), Q(0))
                if self.stride == 2:
                    equality(unit(ti), Q(0))
            else:
                equality(wn, -b[ni])
                if self.stride == 2:
                    ut = A[ti][:] + [Q(0)] * n
                    if mode == 'stick':
                        equality(ut, -b[ti])
                    else:
                        s = Q(1 if mode == 'slide+' else -1)
                        row = unit(ti)
                        row[ni] = s * self.mu[i]
                        equality(row, Q(0))
                        ineq([-s * t for t in ut], s * b[ti])
        return (H, r, E, e)

    def readouts(self):
        """All generalized motion and body contact/net forces, plus each normal force.

        Translational DOFs use m/s,m,N; rotational DOFs use rad/s,rad,N m.
        Body grouping/orientation is the caller's explicit J contract.
        """
        (m, n) = (len(self.J), len(self.gaps))
        out = {}
        for (k, im) in enumerate(self.inverse_mass):
            c = [im * row[k] for row in self.J] + [Q(0)] * n
            out[f'velocity:{k}'] = (c, self.free_velocity[k])
            out[f'displacement:{k}'] = ([self.step * t for t in c], self.step * self.free_velocity[k])
            force = [row[k] / self.step for row in self.J] + [Q(0)] * n
            out[f'contact_force:{k}'] = (force, Q(0))
            out[f'net_force:{k}'] = (force, (self.free_velocity[k] - self.initial_velocity[k]) / im / self.step)
        for i in range(n):
            c = [Q(0)] * (m + n)
            c[self.stride * i] = 1 / self.step
            out[f'normal_force:{i}'] = (c, Q(0))
        return out

def rational_json(obj):
    if isinstance(obj, Q):
        return str(obj)
    if isinstance(obj, dict):
        return {k: rational_json(v) for (k, v) in obj.items()}
    if isinstance(obj, (tuple, list)):
        return [rational_json(v) for v in obj]
    return obj

def reduce_system(H, r, E, e):
    (x0, N) = solve_equalities(E, e, len(H[0]))
    if x0 is None:
        return {'kind': 'equality_infeasible', 'certificate': N}
    B = [[dot(row, col) for col in zip(*N)] for row in H]
    rhs = [ri - dot(row, x0) for (row, ri) in zip(H, r)]
    return {'x0': x0, 'N': N, 'B': B, 'r': rhs}

def replay_lp(B, r, c, cert):
    try:
        if cert['kind'] == 'farkas':
            y = [_q(t) for t in cert['y']]
            return len(y) == len(B) and all((t >= 0 for t in y)) and all((dot(col, y) == 0 for col in zip(*B))) and (dot(y, r) < 0)
        x = [_q(t) for t in cert['x']]
        if len(x) != len(c) or any((dot(row, x) > rhs for (row, rhs) in zip(B, r))):
            return False
        if cert['kind'] == 'feasible':
            return True
        y = [_q(t) for t in cert['y']]
        return cert['kind'] == 'optimum' and len(y) == len(B) and all((t >= 0 for t in y)) and ([dot(col, y) for col in zip(*B)] == [-t for t in c]) and (dot(c, x) == -dot(y, r))
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return False

def certified_lp(B, r, c, stats):
    """Discover a simplex basis; never certify float residuals or solver flags."""
    import numpy as np
    from scipy.optimize import linprog
    n = len(c)
    stats['lp_calls'] += 1
    if not n:
        bad = next((i for (i, rhs) in enumerate(r) if rhs < 0), None)
        if bad is not None:
            y = [Q(0)] * len(B)
            y[bad] = 1
            return {'kind': 'farkas', 'y': y}
        return {'kind': 'optimum', 'x': [], 'y': [Q(0)] * len(B)}
    t = perf_counter()
    H = np.array(B, dtype=float).reshape(len(B), n)
    rr = np.array(r, dtype=float)
    scale = np.maximum(np.max(np.abs(H), axis=1), np.abs(rr))
    scale = np.maximum(scale, 1e-12)
    result = linprog(np.array(c, dtype=float), A_ub=H / scale[:, None], b_ub=rr / scale, bounds=[(None, None)] * n, method='highs', options={'threads': 4})
    stats['discovery_seconds'] += perf_counter() - t
    if result.status == 2:
        P = [row + [-Q(1)] for row in B] + [[Q(0)] * n + [-Q(1)]]
        pr = r + [Q(0)]
        pc = [Q(0)] * n + [Q(1)]
        sub = certified_lp(P, pr, pc, stats)
        if sub['kind'] != 'optimum' or dot(pc, sub['x']) <= 0:
            raise ArithmeticError('infeasibility discovery not replayed')
        cert = {'kind': 'farkas', 'y': sub['y'][:-1]}
    elif result.status == 0:
        active = [i for (i, slack) in enumerate(result.ineqlin.residual) if abs(slack) <= 1e-07]
        rows = [B[i] for i in active]
        vals = [r[i] for i in active]
        (x0, N) = solve_equalities(rows, vals, n)
        if x0 is None:
            raise ArithmeticError('primal active equations inconsistent')
        piv = _rref(rows)[2] if rows else []
        free = [j for j in range(n) if j not in piv]
        x = None
        for denominator in (10 ** 6, 10 ** 9, 10 ** 12):
            z = [Q(float(result.x[j])).limit_denominator(denominator) for j in free]
            trial = [a + dot(row, z) for (a, row) in zip(x0, N)]
            if all((dot(row, trial) <= rhs for (row, rhs) in zip(B, r))):
                x = trial
                break
        if x is None:
            raise ArithmeticError('exact primal face witness rejected')
        support = [i for (i, y) in enumerate(result.ineqlin.marginals) if abs(y) > 1e-09]
        if support:
            D = [[B[i][j] for i in support] for j in range(n)]
            (vals, null) = solve_equalities(D, [-t for t in c], len(support))
            if vals is None:
                raise ArithmeticError('dual basis inconsistent')
            y = [Q(0)] * len(B)
            for (i, val) in zip(support, vals):
                y[i] = val
        else:
            y = [Q(0)] * len(B)
        cert = {'kind': 'optimum', 'x': x, 'y': y}
    else:
        raise ArithmeticError(f'LP discovery failed or unbounded: {result.status}')
    t = perf_counter()
    if not replay_lp(B, r, c, cert):
        raise ArithmeticError('exact LP witness rejected')
    stats['replay_seconds'] += perf_counter() - t
    return cert

def lift(reduced, z):
    return [a + dot(row, z) for (a, row) in zip(reduced['x0'], reduced['N'])]

def enclose_gap_box(problem, *, max_nodes=2048, exhaustive=False):
    """Complete mode union and sharp marginal readout intervals, or UNCERTAIN.

    `exhaustive=True` is the equally informed standard full-enumeration control.
    Every stored prefix/leaf is bound to the problem hash; replay rebuilds its
    equations rather than trusting claimed coefficients or a complete flag.
    """
    start = perf_counter()
    stats = dict(visited_nodes=0, lp_calls=0, endpoint_solutions=0, discovery_seconds=0.0, replay_seconds=0.0)
    out = dict(status='UNCERTAIN', reason='', input_sha256=problem.digest(), complete=False, records=[], bounds={}, stats=stats)
    if len(problem.J) + len(problem.gaps) > 128:
        out['reason'] = 'EXACT_ARITHMETIC_DIMENSION_BUDGET'
        stats['seconds'] = perf_counter() - start
        return out
    readouts = problem.readouts()
    stack = product(tuple(reversed(problem.modes)), repeat=len(problem.gaps)) if exhaustive else [()]
    cache = {}
    try:
        while True:
            if exhaustive:
                prefix = next(stack, None)
                if prefix is None:
                    break
            else:
                if not stack:
                    break
                prefix = stack.pop()
            if stats['visited_nodes'] >= max_nodes:
                raise ArithmeticError('PREFIX_BUDGET')
            stats['visited_nodes'] += 1
            (H, r, E, e) = problem.system(prefix)
            red = reduce_system(H, r, E, e)
            rec = dict(prefix=prefix)
            if red.get('kind') == 'equality_infeasible':
                rec.update(kind='equality_infeasible', certificate=red['certificate'])
                out['records'].append(rec)
                continue
            (B, rhs) = (red['B'], red['r'])
            nz = len(red['N'][0])
            feasible = certified_lp(B, rhs, [Q(0)] * nz, stats)
            if feasible['kind'] == 'farkas':
                rec.update(kind='pruned', certificate=feasible)
                out['records'].append(rec)
                continue
            if len(prefix) < len(problem.gaps):
                for mode in reversed(problem.modes):
                    stack.append(prefix + (mode,))
                continue
            rec.update(kind='leaf', feasible=feasible, endpoints={})
            for (name, (c, offset)) in readouts.items():
                cr = [dot(c, col) for col in zip(*red['N'])]
                base = offset + dot(c, red['x0'])
                for (side, sign) in [('lo', Q(1)), ('hi', Q(-1))]:
                    obj = [sign * t for t in cr]
                    key = tuple(obj)
                    if key not in cache:
                        cache[key] = None
                    cert = cache[key]
                    if cert is None:
                        cert = certified_lp(B, rhs, obj, stats)
                        if cert['kind'] != 'optimum':
                            raise ArithmeticError('feasible leaf lost')
                        stats['endpoint_solutions'] += 1
                        cache[key] = cert
                    value = base + dot(cr, cert['x'])
                    rec['endpoints'][name + ':' + side] = cert
                    bound = out['bounds'].setdefault(name, {})
                    if side not in bound or (value < bound[side] if side == 'lo' else value > bound[side]):
                        bound[side] = value
                        bound[side + '_witness'] = {'prefix': prefix, 'x': lift(red, cert['x'])}
            cache.clear()
            out['records'].append(rec)
        if not out['bounds']:
            raise ArithmeticError('EMPTY_MODEL: no existence witness')
        out['complete'] = True
        out['status'] = 'ENTYDIG' if all((b['lo'] == b['hi'] for b in out['bounds'].values())) else "MULTIPLE"
        out['reason'] = 'EXACT_FINITE_MODE_COVER'
        out['contact_candidates'] = contact_candidates(problem, out['bounds'])
    except (ArithmeticError, ValueError, OverflowError) as exc:
        out.update(status='UNCERTAIN', reason=str(exc), complete=False, bounds={})
    stats['seconds'] = perf_counter() - start
    return out

def contact_candidates(problem, bounds):
    """Sharp possible/necessary load, plus conservative measurement candidates.

    A varying load is a reason to investigate geometry, not proof that measuring
    that gap alone resolves force fibres or material uncertainty. Use paired
    single-gap witnesses for causal measurement claims (consumer PORT example).
    """
    return [dict(candidate=i, gap_interval=gap, normal_force_interval=(bounds[f'normal_force:{i}']['lo'], bounds[f'normal_force:{i}']['hi']), can_bear_load=bounds[f'normal_force:{i}']['hi'] > 0, must_bear_load=bounds[f'normal_force:{i}']['lo'] > 0, investigate_gap=gap[0] < gap[1] and bounds[f'normal_force:{i}']['lo'] < bounds[f'normal_force:{i}']['hi']) for (i, gap) in enumerate(problem.gaps)]

def verify_gap_enclosure(problem, answer):
    """Solver-free exact replay, including tree coverage and all claimed bounds."""
    try:
        if answer['input_sha256'] != problem.digest() or not answer['complete'] or answer['status'] not in ('ENTYDIG', "MULTIPLE"):
            return False
        terminal = {}
        bounds = {}
        for rec in answer['records']:
            prefix = tuple(rec['prefix'])
            if prefix in terminal:
                return False
            terminal[prefix] = rec['kind']
            (H, r, E, e) = problem.system(prefix)
            red = reduce_system(H, r, E, e)
            if rec['kind'] == 'equality_infeasible':
                if not verify_observable(E, e, [[Q(0)] * len(H[0])], rec['certificate']):
                    return False
                if rec['certificate']['status'] != 'INCOMPATIBLE':
                    return False
                continue
            (B, rhs) = (red['B'], red['r'])
            nz = len(red['N'][0])
            if rec['kind'] == 'pruned':
                if rec['certificate']['kind'] != 'farkas' or not replay_lp(B, rhs, [Q(0)] * nz, rec['certificate']):
                    return False
                continue
            if rec['kind'] != 'leaf' or len(prefix) != len(problem.gaps):
                return False
            if not replay_lp(B, rhs, [Q(0)] * nz, rec['feasible']):
                return False
            for (name, (c, off)) in problem.readouts().items():
                cr = [dot(c, col) for col in zip(*red['N'])]
                base = off + dot(c, red['x0'])
                for (side, sign) in [('lo', Q(1)), ('hi', Q(-1))]:
                    cert = rec['endpoints'][name + ':' + side]
                    if cert['kind'] != 'optimum' or not replay_lp(B, rhs, [sign * t for t in cr], cert):
                        return False
                    value = base + dot(cr, [_q(t) for t in cert['x']])
                    b = bounds.setdefault(name, {})
                    if side not in b or (value < b[side] if side == 'lo' else value > b[side]):
                        b[side] = value

        def covered(prefix):
            if prefix in terminal:
                return True
            if len(prefix) == len(problem.gaps):
                return False
            return all((covered(prefix + (mode,)) for mode in problem.modes))
        if not bounds or not covered(()) or set(bounds) != set(answer['bounds']):
            return False
        if any((p[:k] in terminal for p in terminal for k in range(len(p)))):
            return False
        for (name, b) in bounds.items():
            for side in ('lo', 'hi'):
                if _q(answer['bounds'][name][side]) != b[side]:
                    return False
                wit = answer['bounds'][name][side + '_witness']
                x = [_q(t) for t in wit['x']]
                (H, r, E, e) = problem.system(tuple(wit['prefix']))
                if len(wit['prefix']) != len(problem.gaps) or len(x) != len(H[0]):
                    return False
                if any((dot(row, x) > rhs for (row, rhs) in zip(H, r))) or any((dot(row, x) != rhs for (row, rhs) in zip(E, e))):
                    return False
                (c, off) = problem.readouts()[name]
                if off + dot(c, x) != b[side]:
                    return False
        status = 'ENTYDIG' if all((b['lo'] == b['hi'] for b in bounds.values())) else "MULTIPLE"
        if 'contact_candidates' in answer and rational_json(answer['contact_candidates']) != rational_json(contact_candidates(problem, bounds)):
            return False
        return status == answer['status']
    except (KeyError, TypeError, ValueError, ArithmeticError, ZeroDivisionError):
        return False

def solve_gap_box(J, inverse_mass, free_velocity, gaps, step, **kwargs):
    """Convenience fail-closed API. `max_nodes` controls work, never validity."""
    max_nodes = kwargs.pop('max_nodes', 2048)
    try:
        pb = GapProblem.make(J, inverse_mass, free_velocity, gaps, step, **kwargs)
        return enclose_gap_box(pb, max_nodes=max_nodes)
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as exc:
        return {'status': 'UNCERTAIN', 'reason': str(exc), 'complete': False, 'bounds': {}, 'records': []}
