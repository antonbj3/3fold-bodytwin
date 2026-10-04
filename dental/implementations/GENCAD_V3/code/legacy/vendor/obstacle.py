"""Continuous projected antagonist obstacle on each roof triangle intersection."""
import numpy as np
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import splu
from scipy.optimize import minimize, linprog
from continuous import cross, inside, plane, broad_phase

def constraints(U, xy, faces, clearance):
    C = np.c_[xy, np.zeros(len(xy))][faces]
    pairs = broad_phase(U, C)
    up = plane(U)
    rr = []
    cc = []
    vv = []
    bb = []
    count = 0
    for first in range(0, len(pairs), 2048):
        ids = pairs[first:first + 2048]
        u = U[ids[:, 0], :, :2]
        l = C[ids[:, 1], :, :2]
        eu = np.roll(u, -1, axis=1) - u
        el = np.roll(l, -1, axis=1) - l
        delta = l[:, None, :, :] - u[:, :, None, :]
        den = cross(eu[:, :, None, :], el[:, None, :, :])
        ok = np.abs(den) > 1e-14
        t = np.divide(cross(delta, el[:, None, :, :]), den, out=np.zeros_like(den), where=ok)
        s = np.divide(cross(delta, eu[:, :, None, :]), den, out=np.zeros_like(den), where=ok)
        edgeok = ok & (t >= 0) & (t <= 1) & (s >= 0) & (s <= 1)
        inter = (u[:, :, None, :] + t[:, :, :, None] * eu[:, :, None, :]).reshape(-1, 9, 2)
        points = np.concatenate([u, l, inter], axis=1)
        valid = np.c_[inside(u, l), inside(l, u), edgeok.reshape(-1, 9)]
        (ii, jj) = np.where(valid)
        p = points[ii, jj]
        lt = l[ii]
        source = ids[ii, 0]
        candidate = ids[ii, 1]
        (a, b, c) = (lt[:, 0], lt[:, 1], lt[:, 2])
        det = (b[:, 1] - c[:, 1]) * (a[:, 0] - c[:, 0]) + (c[:, 0] - b[:, 0]) * (a[:, 1] - c[:, 1])
        w1 = ((b[:, 1] - c[:, 1]) * (p[:, 0] - c[:, 0]) + (c[:, 0] - b[:, 0]) * (p[:, 1] - c[:, 1])) / det
        w2 = ((c[:, 1] - a[:, 1]) * (p[:, 0] - c[:, 0]) + (a[:, 0] - c[:, 0]) * (p[:, 1] - c[:, 1])) / det
        weights = np.c_[w1, w2, 1 - w1 - w2]
        height = np.sum(np.c_[p, np.ones(len(p))] * up[source], axis=1) - clearance
        rr.extend(np.repeat(np.arange(count, count + len(p)), 3))
        cc.extend(faces[candidate].ravel())
        vv.extend(weights.ravel())
        bb.extend(height)
        count += len(p)
    A = coo_matrix((vv, (rr, cc)), shape=(count, len(xy))).tocsr()
    b = np.array(bb)
    return (A, b)

def hessian(faces, n, regularizer):
    edges = np.r_[faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    m = len(edges)
    L = coo_matrix((np.r_[np.ones(m), -np.ones(m)], (np.r_[np.arange(m), np.arange(m)], np.r_[edges[:, 0], edges[:, 1]])), shape=(m, n)).tocsr()
    return eye(n, format='csc') + regularizer * (L.T @ L).tocsc()

def qp(A, b, H, prior, maxiter=2000, gtol=1e-09):
    factor = splu(H)
    residual = A @ prior - b

    def fun(lam):
        v = A.T @ lam
        delta = factor.solve(v)
        return (0.5 * v @ delta - lam @ residual, A @ delta - residual)
    sol = minimize(fun, np.zeros(len(b)), jac=True, method='L-BFGS-B', bounds=[(0, None)] * len(b), options=dict(maxiter=maxiter, gtol=gtol, ftol=1e-12, maxcor=10))
    z = prior - factor.solve(A.T @ sol.x)
    return (z, dict(success=bool(sol.success), iterations=sol.nit, maximum_constraint_violation_mm=float(max(0, np.max(A @ z - b))), message=sol.message))

def lp_control(A, b, prior):
    from scipy.sparse import hstack, vstack
    n = len(prior)
    I = eye(n, format='csr')
    Z = coo_matrix((len(b), n))
    constraints = vstack([hstack([A, Z]), hstack([I, -I]), hstack([-I, -I])], format='csr')
    rhs = np.r_[b, prior, -prior]
    result = linprog(np.r_[np.zeros(n), np.ones(n)], A_ub=constraints, b_ub=rhs, bounds=[(None, None)] * n + [(0, None)] * n, method='highs', options={'threads': 1})
    return (result.x[:n] if result.success else None, dict(success=bool(result.success), message=result.message))
