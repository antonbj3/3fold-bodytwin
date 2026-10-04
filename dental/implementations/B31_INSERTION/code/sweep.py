"""Continuous affine triangle sweeps. Digital coordinates are exact IEEE dyadics.
Outward-rounded projection bounds certify separation. LP status never certifies.
"""
import itertools, time, warnings
from fractions import Fraction as Q
import numpy as np
from scipy.optimize import linprog
from rtree import index
DOWN = -np.inf
UP = np.inf

def rational(x):
    return Q(float(x))

def nextdown(x):
    return np.nextafter(x, DOWN)

def nextup(x):
    return np.nextafter(x, UP)

def dot_interval(v, n):
    p = v * n
    lo = nextdown(p)
    hi = nextup(p)
    return (nextdown(nextdown(lo[..., 0] + lo[..., 1]) + lo[..., 2]), nextup(nextup(hi[..., 0] + hi[..., 1]) + hi[..., 2]))

def separator_batch(A, B, D):
    """A and B: N x 3 x 3; any float axis is valid for projection proof."""
    n = len(A)
    if n == 0:
        return (np.empty(0, bool), np.empty(0))
    ea = np.roll(A, -1, axis=1) - A
    eb = np.roll(B, -1, axis=1) - B
    axes = [np.cross(ea[:, 0], ea[:, 1]), np.cross(eb[:, 0], eb[:, 1])]
    axes.extend((np.cross(ea[:, k], D) for k in range(3)))
    axes.extend((np.cross(ea[:, j], eb[:, k]) for j in range(3) for k in range(3)))
    axes.extend((np.cross(np.broadcast_to(D, (n, 3)), eb[:, k]) for k in range(3)))
    sep = np.zeros(n, bool)
    gap = np.zeros(n)
    for axis in axes:
        mag = np.abs(axis).max(1)
        axis = np.divide(axis, mag[:, None], out=np.zeros_like(axis), where=mag[:, None] > 0)
        (al, ah) = dot_interval(A, axis[:, None])
        (bl, bh) = dot_interval(B, axis[:, None])
        (dl, dh) = dot_interval(np.broadcast_to(D, (n, 3)), axis)
        alo = nextdown(al.min(1) + np.minimum(0, dl))
        ahi = nextup(ah.max(1) + np.maximum(0, dh))
        g = np.maximum(nextdown(bl.min(1) - ahi), nextdown(alo - bh.max(1)))
        good = (g > 0) & np.isfinite(g)
        sep |= good
        gap = np.maximum(gap, np.where(good, g, 0))
    return (sep, gap)

def rref_solve(M, b):
    a = [list(row) + [v] for (row, v) in zip(M, b)]
    m = len(a)
    n = len(M[0])
    piv = []
    r = 0
    for c in range(n):
        j = next((i for i in range(r, m) if a[i][c]), None)
        if j is None:
            continue
        (a[r], a[j]) = (a[j], a[r])
        v = a[r][c]
        a[r] = [x / v for x in a[r]]
        for i in range(m):
            if i != r and a[i][c]:
                f = a[i][c]
                a[i] = [x - f * y for (x, y) in zip(a[i], a[r])]
        piv.append(c)
        r += 1
        if r == m:
            break
    if any((all((x == 0 for x in row[:n])) and row[n] != 0 for row in a)):
        return None
    x = [Q(0)] * n
    for (i, c) in enumerate(piv):
        x[c] = a[i][-1]
    return x

def primal_matrix(A, B, D):
    M = [[rational(x) for x in list(A[:, j]) + list(-B[:, j]) + [D[j], 0.0]] for j in range(3)]
    M += [[Q(int(k < 3)) for k in range(8)], [Q(int(3 <= k < 6)) for k in range(8)], [Q(int(k >= 6)) for k in range(8)]]
    return (M, [Q(0)] * 3 + [Q(1)] * 3)

def validate_witness(A, B, D, w):
    try:
        x = list(map(Q, w['weights']))
        (M, b) = primal_matrix(A, B, D)
        valid = len(x) == 8 and min(x) >= 0 and all((sum((v * t for (v, t) in zip(row, x))) == rhs for (row, rhs) in zip(M, b)))
        if not valid:
            return False
        point = [float(sum((rational(A[k, j]) * x[k] for k in range(3))) + rational(D[j]) * x[6]) for j in range(3)]
        return Q(w['s_exact']) == x[6] and w['s_float'] == float(x[6]) and (w['point_mm'] == point)
    except (ValueError, KeyError, TypeError):
        return False

def lp_witness(A, B, D, maximize_s=False):
    (M, b) = primal_matrix(A, B, D)
    Mf = np.array(M, float)
    bf = np.array(b, float)
    c = np.zeros(8)
    c[6] = -1 if maximize_s else 0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sol = linprog(c, A_eq=Mf, b_eq=bf, bounds=(0, None), method='highs-ds', options={'threads': 4, 'primal_feasibility_tolerance': 1e-09, 'dual_feasibility_tolerance': 1e-09})
    if not sol.success:
        return (None, int(sol.status))
    support = np.flatnonzero(sol.x > 1e-10).tolist()
    candidates = [support, np.flatnonzero(sol.x > 0).tolist()]
    for size in range(max(len(support), 1), 7):
        candidates.extend((sorted(set(support) | set(q)) for q in itertools.combinations([i for i in range(8) if i not in support], size - len(support))))
    for supp in candidates:
        if not supp or len(supp) > 8:
            continue
        x = rref_solve([[row[k] for k in supp] for row in M], b)
        if x is None or min(x) < 0:
            continue
        full = [Q(0)] * 8
        for (k, v) in zip(supp, x):
            full[k] = v
        w = {'weights': [str(v) for v in full], 's_exact': str(full[6]), 's_float': float(full[6]), 'point_mm': [float(sum((rational(A[k, j]) * full[k] for k in range(3))) + rational(D[j]) * full[6]) for j in range(3)], 'arithmetic': 'exact rational on source float64 dyadics', 'LP_status': int(sol.status)}
        if validate_witness(A, B, D, w):
            return (w, 0)
    return (None, 0)

def independent_lp(A, B, D):
    """Five variables, implicit third barycentric weights. Independent constraint layout."""
    E = np.column_stack((A[0] - A[2], A[1] - A[2], B[2] - B[0], B[2] - B[1], D))
    rhs = B[2] - A[2]
    U = np.array([[1, 1, 0, 0, 0], [0, 0, 1, 1, 0]], float)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        s = linprog(np.zeros(5), A_eq=E, b_eq=rhs, A_ub=U, b_ub=np.ones(2), bounds=[(0, None)] * 4 + [(0, 1)], method='highs', options={'threads': 4, 'primal_feasibility_tolerance': 1e-09, 'dual_feasibility_tolerance': 1e-09})
    return dict(feasible=bool(s.success), status=int(s.status), max_residual=float(np.max(np.abs(E @ s.x - rhs))) if s.success else None, reference='SciPy HiGHS barycentric feasibility; float status paired with exact witness or interval separation')

class Obstacle:

    def __init__(self, tri):
        self.tri = np.asarray(tri, float)
        prop = index.Property()
        prop.dimension = 3
        self.lo = self.tri.min(1)
        self.hi = self.tri.max(1)
        self.tree = index.Index(((i, tuple(np.r_[a, b]), None) for (i, (a, b)) in enumerate(zip(self.lo, self.hi))), properties=prop) if len(tri) else None

def swept_bounds(A, D):
    lo = A.min(1)
    hi = A.max(1)
    return (nextdown(lo + np.minimum(0, D)), nextup(hi + np.maximum(0, D)))

def check_surface(tri, obstacle, D, face_ids=None, roles=None, stop_on_hit=True, budget_s=180):
    tick = time.perf_counter()
    tri = np.asarray(tri, float)
    D = np.asarray(D, float)
    face_ids = np.arange(len(tri)) if face_ids is None else np.asarray(face_ids)
    stats = dict(source_faces=len(tri), obstacle_faces=len(obstacle.tri), candidate_pairs=0, separated_pairs=0, LP_queries=0, unknown_pairs=0, processed_faces=0)
    if not len(obstacle.tri):
        return dict(status='UNKNOWN', reason='missing obstacle surface', coverage=stats, seconds=time.perf_counter() - tick)
    (lo, hi) = swept_bounds(tri, D)
    pending = []
    witnesses = []
    for (i, (a, b)) in enumerate(zip(lo, hi)):
        ids = list(obstacle.tree.intersection(tuple(np.r_[a, b])))
        stats['processed_faces'] += 1
        if ids:
            stats['candidate_pairs'] += len(ids)
            As = np.broadcast_to(tri[i], (len(ids), 3, 3))
            Bs = obstacle.tri[ids]
            (sep, g) = separator_batch(As, Bs, D)
            stats['separated_pairs'] += int(sep.sum())
            for j in np.flatnonzero(~sep):
                stats['LP_queries'] += 1
                (w, stat) = lp_witness(tri[i], Bs[j], D, maximize_s=True)
                if w:
                    w.update(source_face_id=int(face_ids[i]), obstacle_face_id=int(ids[j]), source_role=int(roles[i]) if roles is not None else None, source_triangle_mm=tri[i].tolist(), obstacle_triangle_mm=Bs[j].tolist(), translation_mm=D.tolist(), independent_control=independent_lp(tri[i], Bs[j], D))
                    witnesses.append(w)
                    if stop_on_hit:
                        return dict(status='COLLISION', witness=w, coverage=stats, seconds=time.perf_counter() - tick, full_enumeration=False)
                else:
                    stats['unknown_pairs'] += 1
                    pending.append([int(face_ids[i]), int(ids[j]), stat])
        if time.perf_counter() - tick > budget_s:
            return dict(status='COLLISION' if witnesses else 'UNKNOWN', reason='time budget, uncovered pairs abstain', witness=witnesses[0] if witnesses else None, coverage=stats, unresolved=pending[:20], seconds=time.perf_counter() - tick, full_enumeration=False)
    status = 'COLLISION' if witnesses else 'UNKNOWN' if pending else 'SURFACE_CLEAR_CERTIFIED'
    return dict(status=status, witness=witnesses[0] if witnesses else None, coverage=stats, seconds=time.perf_counter() - tick, full_enumeration=True, unresolved=pending[:20], certificate='All swept triangle pairs have positive outward-rounded separating projection or conservative AABB separation' if status == 'SURFACE_CLEAR_CERTIFIED' else None, solid_clearance='UNKNOWN; open surfaces/initial containment not established')
