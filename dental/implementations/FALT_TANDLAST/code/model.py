"""Dental consumer of archived rational kernels; no engine changes.

Rigid axial bodies with arbitrary positive finite tooth support stiffness,
normal-only contact, a force-controlled closure, declared contact incidence.
See QUANTITY_LICENSE.md for strict positivity and closed hull proofs.
"""
from dental_release.paths import expand as _release_expand
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import json, sys, time
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'engine_snapshot' / 'src'))
from motion_engine.ncp.gap_box import solve_equalities, certified_lp, lift, rational_json, dot
X21 = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/X21'))
X2 = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/X2'))
X54 = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/X54'))
SUPPORT = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace/notes/occlusion_support_intervals/SUPPORT_INTERVALS.json'))
FDI = tuple((q * 10 + t for q in range(1, 5) for t in range(1, 9)))
REGIONS = ('I_R', 'C_R', 'P_R', 'M_R', 'I_L', 'C_L', 'P_L', 'M_L')
TOTAL = (Q(205), Q(1963))

def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()

def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rational_json(value), indent=2, ensure_ascii=False, sort_keys=True, allow_nan=False) + '\n')

def stopped():
    if (ROOT / 'STOP').exists() or (ROOT.parents[1] / 'lanes/SOL_FALT_STOP').exists():
        raise RuntimeError('STOP_REQUESTED')

def region(t):
    side = 'R' if t // 10 in (1, 4) else 'L'
    n = t % 10
    return ('I' if n <= 2 else 'C' if n == 3 else 'P' if n <= 5 else 'M') + '_' + side

def edges_from_x21(case, refined=True):
    p = X21 / 'raw/refined_certificates' / (case + '.json')
    if not refined or not p.exists():
        p = X21 / 'raw/certificates' / (case + '.json')
    rows = json.loads(p.read_text())['pairs']
    edges = []
    for r in rows:
        (lo, hi) = r['minimum_gap_enclosure_mm']
        if lo is None or hi is None or lo > hi:
            raise ValueError('Invalid source gap enclosure')
        edges.append((int(r['upper_fdi']), int(r['lower_fdi']), Q(lo), Q(hi)))
    if not edges or len({e[:2] for e in edges}) != len(edges):
        raise ValueError('Nonempty unique tooth-pair graph required')
    return (edges, p)

def classify(edges, teeth):
    """Exact ALL-gap/ALL-positive-support quantifiers for one arch subset."""
    teeth = frozenset(teeth)
    incident = [i for (i, e) in enumerate(edges) if e[0] in teeth or e[1] in teeth]
    outside = [i for i in range(len(edges)) if i not in incident]
    if not incident:
        (cls, share, force, margin) = ('NEVER', (Q(0), Q(0)), (Q(0), Q(0)), None)
    elif not outside:
        (cls, share, force, margin) = ('MUST', (Q(100), Q(100)), TOTAL, None)
    else:
        a = min((edges[i][3] for i in incident))
        b = min((edges[i][2] for i in outside))
        margin = b - a
        cls = 'MUST' if margin >= 0 else 'CAN'
        (share, force) = ((Q(0), Q(100)), (Q(0), TOTAL[1]))
    return dict(classification=cls, can_bear=bool(incident), must_bear=cls == 'MUST', never_bears=cls == 'NEVER', force_hull_N=force, share_hull_pp=share, infimum_attained=cls != 'MUST' or not outside, necessity_margin_mm=margin, incident_count=len(incident), outside_count=len(outside), license='AXIAL_POSITIVE_SUPPORT_V1', physical_classification='UNCERTAIN', graph_coverage='DECLARED_PAIRS_ONLY', zero_force_physical_interpretation='UNKNOWN' if not incident else None)

def exact_point(edges, gaps, T, compliance):
    """Float QP discovery; archived exact equality/LP kernels certify KKT.

    The convex quadratic is .5 lambda^T B^T C B lambda + g^T lambda,
    lambda>=0, sum lambda=T. All inputs/replays are exact rationals.
    """
    import numpy as np
    from scipy.optimize import minimize
    start = time.perf_counter()
    n = len(edges)
    teeth = sorted({v for e in edges for v in e[:2]})
    idx = {t: i for (i, t) in enumerate(teeth)}
    cr = [Q(compliance[t]) for t in teeth]
    B = np.zeros((len(teeth), n))
    for (i, e) in enumerate(edges):
        B[idx[e[0]], i] = B[idx[e[1]], i] = 1
    A = B.T @ (np.array(cr, dtype=float)[:, None] * B)
    gr = list(map(Q, gaps))
    Tf = float(T)
    objective_scale = max(float(np.ptp(np.array(gr, dtype=float))), Tf * np.max(A), 1e-05)
    g = np.array(gr, dtype=float) / objective_scale
    AA = Tf * A / objective_scale
    res = minimize(lambda x: 0.5 * x @ AA @ x + g @ x, np.ones(n) / n, jac=lambda x: AA @ x + g, constraints={'type': 'eq', 'fun': lambda x: x.sum() - 1, 'jac': lambda x: np.ones(n)}, bounds=[(0, 1)] * n, method='SLSQP', options={'ftol': 1e-12, 'maxiter': 1500})
    discovery_s = time.perf_counter() - start
    if not res.success:
        raise ArithmeticError('QP discovery failed: ' + res.message)
    stats = dict(lp_calls=0, discovery_seconds=0.0, replay_seconds=0.0)
    Ar = [[sum((cr[k] for (k, t) in enumerate(teeth) if t in edges[i][:2] and t in edges[j][:2])) for j in range(n)] for i in range(n)]
    gradient = A @ (res.x * Tf) + np.array(gr, dtype=float)
    d0 = float(np.min(gradient))
    active = [i for i in range(n) if res.x[i] > 1e-08 or abs(gradient[i] - d0) < 1e-07]
    m = len(active)
    E = [[Ar[i][j] for j in active] + [Q(-1)] for i in active]
    E += [[Q(1)] * m + [Q(0)]]
    rhs = [-gr[i] for i in active] + [Q(T)]
    (x0, N) = solve_equalities(E, rhs, m + 1)
    if x0 is None:
        raise ArithmeticError('Exact QP equations incompatible')
    (H, h) = ([], [])
    for j in range(m):
        H.append([Q(-int(j == k)) for k in range(m)] + [Q(0)])
        h.append(Q(0))
    for i in range(n):
        H.append([-Ar[i][j] for j in active] + [Q(1)])
        h.append(gr[i])
    reduced_B = [[dot(row, col) for col in zip(*N)] for row in H]
    reduced_r = [b - dot(row, x0) for (row, b) in zip(H, h)]
    nz = len(N[0])
    cert = certified_lp(reduced_B, reduced_r, [Q(0)] * nz, stats)
    if cert['kind'] != 'optimum':
        raise ArithmeticError('QP active-set exact replay rejected')
    x = lift(dict(x0=x0, N=N), cert['x'])
    f = [Q(0)] * n
    for (i, val) in zip(active, x[:-1]):
        f[i] = val
    tooth_force = {t: sum((f[i] for (i, e) in enumerate(edges) if t in e[:2])) for t in teeth}
    closure = x[-1]
    w = [gr[i] - closure + sum((cr[k] * tooth_force[t] for (k, t) in enumerate(teeth) if t in edges[i][:2])) for i in range(n)]
    if sum(f) != T or any((fi < 0 or wi < 0 or fi * wi for (fi, wi) in zip(f, w))):
        raise ArithmeticError('Exact KKT rejected')
    return dict(edge_force_N=f, tooth_force_N=tooth_force, gap_mm=gr, closure_mm=closure, support_compliance_mm_per_N={t: cr[idx[t]] for t in teeth}, total_N=Q(T), exact_KKT=True, discovery_s=discovery_s, total_s=time.perf_counter() - start, lp_stats=stats, QP_objective='CONVEX_WITH_UNIQUE_TOOTH_FORCE; edge force may be nonunique')

def replay_point(edges, ans):
    try:
        f = list(map(Q, ans['edge_force_N']))
        g = list(map(Q, ans['gap_mm']))
        F = {int(k): Q(v) for (k, v) in ans['tooth_force_N'].items()}
        C = {int(k): Q(v) for (k, v) in ans['support_compliance_mm_per_N'].items()}
        (d, T) = (Q(ans['closure_mm']), Q(ans['total_N']))
        teeth = {t for e in edges for t in e[:2]}
        if len(f) != len(edges) or len(g) != len(edges) or set(F) != teeth or (set(C) != teeth):
            return False
        if any((c <= 0 for c in C.values())) or sum(f) != T or T <= 0:
            return False
        if any((F[t] != sum((f[i] for (i, e) in enumerate(edges) if t in e[:2])) for t in teeth)):
            return False
        w = [g[i] - d + C[u] * F[u] + C[v] * F[v] for (i, (u, v, _, _)) in enumerate(edges)]
        return all((edges[i][2] <= g[i] <= edges[i][3] and fi >= 0 and (wi >= 0) and (fi * wi == 0) for (i, (fi, wi)) in enumerate(zip(f, w))))
    except (KeyError, ValueError, TypeError, ZeroDivisionError):
        return False
