import numpy as np
import time
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack, vstack, diags
from .mesh import field_basis

def solve(t, r, regional=True, extra=None):
    start = time.perf_counter()
    (B, meta) = field_basis(t, regional)
    N = B.shape[1]
    ids = np.flatnonzero(t['query_good'])
    m = len(ids)
    qmap = np.full(len(t['xy']), -1, int)
    qmap[ids] = np.arange(m)
    hitkeep = t['query_good'][t['hit_queries']]
    q = t['hit_queries'][hitkeep]
    hmat = t['H'][hitkeep] @ B
    g = t['ceiling'][q] - t['hit_z'][hitkeep]
    ref = t['reference_gap'][ids]
    contact = (ref >= 0) & (ref <= 0.1)
    lower = np.where(contact, 0.0, 0.100001)
    L = csr_matrix((np.ones(len(q)), (np.arange(len(q)), qmap[q])), shape=(len(q), m))
    ownerH = t['H'][t['owner'][ids[contact]]] @ B
    upperN = int(contact.sum())
    upperId = np.flatnonzero(contact)
    U = csr_matrix((np.ones(upperN), (np.arange(upperN), upperId)), shape=(upperN, m))
    A = hstack([-hmat, -L, csr_matrix((len(q), upperN))], format='csr')
    rhs = g - lower[qmap[q]]
    Aupper = hstack([ownerH, csr_matrix((upperN, m)), -diags(np.ones(upperN))], format='csr')
    bupper = 0.1 - t['original_gap'][ids[contact]]
    A = vstack([A, Aupper], format='csr')
    rhs = np.r_[rhs, bupper]
    if extra is not None and extra[0].shape[0]:
        (X, xb) = extra
        A = vstack([A, hstack([X @ B, csr_matrix((len(xb), m + upperN))], format='csr')], format='csr')
        rhs = np.r_[rhs, xb]
    qw = t['query_weights'][ids]
    cfit = np.r_[np.zeros(N), qw, qw[contact]]
    cost = np.asarray(t['vertex_area_weights'] @ B).ravel()
    crem = np.r_[cost, np.zeros(m + upperN)]
    bounds = [(0.0, r['relief_cap_mm'])] * N + [(0.0, None)] * (m + upperN)
    options = dict(dual_feasibility_tolerance=1e-08, primal_feasibility_tolerance=1e-08)
    s1 = linprog(cfit, A_ub=A, b_ub=rhs, bounds=bounds, method='highs', options=options)
    if not s1.success:
        raise RuntimeError('Phase1 LP failed: ' + s1.message)
    s2 = linprog(crem, A_ub=vstack([A, csr_matrix(cfit[None, :])], format='csr'), b_ub=np.r_[rhs, s1.fun + 1e-08], bounds=bounds, method='highs', options=options)
    if not s2.success:
        raise RuntimeError('Phase2 LP failed: ' + s2.message)
    d = np.asarray(B @ s2.x[:N]).ravel()
    violation = float(np.max(A @ s2.x - rhs, initial=0.0))
    primal_dual_gap = float(abs(s2.fun - (s2.ineqlin.marginals @ np.r_[rhs, s1.fun + 1e-08] + s2.lower.marginals @ np.zeros(len(s2.x)) + s2.upper.marginals[:N] @ np.full(N, r['relief_cap_mm']))))
    return (d, dict(**meta, phase1_weighted_band_violation_mm3=s1.fun, phase2_area_integral_mm3=s2.fun, inequality_max_violation_mm=violation, dual_objective_gap_mm3=primal_dual_gap, seconds=time.perf_counter() - start, minimum_scope='Floating LP optimum of frozen nodal band surrogate then integral displacement; not a global binary-mask optimum', upper_contact_owner_constraint='Sufficient fixed original owner; actual upper envelope recomputed after solve', source_info='Same mask/grid/pose/cap as both numerical controls'))
