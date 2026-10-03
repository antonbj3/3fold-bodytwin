"""Exact null-space formulation for L1 contact bands.

The contact observable has at most two free coordinates.  Box feasibility and
minimum effort still require the remaining null-space coordinates; they are
retained here instead of silently approximating the projected set.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.optimize import linprog
import clarabel


def _settings():
    s = clarabel.DefaultSettings()
    s.verbose = False
    s.max_iter = 120
    s.tol_gap_abs = 1e-8
    s.tol_feas = 1e-8
    return s


def contact_geometry(A, b, F0, cj, med, lat):
    """Return unconstrained min-effort solution, null basis, and contact rank.

    All coordinates are normalized activations. ``med`` and ``lat`` are
    ``(muscle coefficient, constant offset)`` pairs in newtons.
    """
    A = np.asarray(A, float); b = np.asarray(b, float)
    F0 = np.asarray(F0, float); cj = np.asarray(cj, float)
    if np.any(F0 <= 0):
        raise ValueError('F0 must be positive')
    B = A * F0[None, :]
    Q, R, piv = la.qr(B.T, mode='full', pivoting=True, check_finite=False)
    rank = np.linalg.matrix_rank(B, tol=1e-10*np.linalg.svd(B, compute_uv=False)[0])
    if rank != len(b):
        raise ValueError('A must have independent rows')
    particular = Q[:, :rank] @ la.solve_triangular(R[:rank, :rank].T, b[piv[:rank]], lower=True, check_finite=False)
    N = Q[:, rank:]
    if np.linalg.norm(B @ particular - b) > 1e-7 * max(1., np.linalg.norm(b)):
        raise ValueError('inconsistent equilibrium')
    C = np.stack([cj, np.asarray(med[0], float), np.asarray(lat[0], float)]) * F0[None, :]
    S = np.linalg.svd(C @ N, compute_uv=False)
    contact_rank = int(np.sum(S > (S[0] if len(S) else 0)*1e-10))
    return particular, N, contact_rank


def contact_band_fast(A, b, F0, cj, C0, med, lat, eps=None):
    """Return total-contact extrema, in N, with verified force witnesses.

    ``eps=None`` solves the LP. A nonnegative eps solves the E-band with
    E=sum((f/F0)**2) and the exact constrained E_min.  The function name is a
    compatibility entry point; speed is an empirical question.
    """
    p, N, rank = contact_geometry(A, b, F0, cj, med, lat)
    F0 = np.asarray(F0, float); cj = np.asarray(cj, float)
    k = N.shape[1]; bw_scale = max(1., abs(float(C0)), np.linalg.norm(np.asarray(b, float)))
    # The inequalities retain all nuisance coordinates.  This is the exact
    # extended formulation of the 2-D contact projection.
    rows = [-N, N]
    rhs = [p, 1-p]
    for v, offset in (med, lat):
        v = np.asarray(v, float)*F0/bw_scale
        rows.append(-(v @ N)[None, :])
        rhs.append(np.array([(float(offset)+v @ p*bw_scale)/bw_scale]))
    G = np.vstack(rows); h = np.concatenate(rhs)
    c = (cj*F0) @ N
    cp = float(C0 + (cj*F0) @ p)
    if eps is None:
        ends = []
        for sign in (1., -1.):
            sol = linprog(sign*c, A_ub=G, b_ub=h, bounds=[(None,None)]*k,
                          method='highs', options={'primal_feasibility_tolerance':1e-9,
                                                   'dual_feasibility_tolerance':1e-9})
            if not sol.success:
                raise ValueError(f'LP failed: {sol.message}')
            ends.append((cp+c@sol.x, (p+N@sol.x)*F0))
        emin = None
    else:
        if eps < 0: raise ValueError('eps must be nonnegative')
        Acone = sp.csc_matrix(G)
        cones = [clarabel.NonnegativeConeT(len(h))]
        base = clarabel.DefaultSolver(2*sp.eye(k,format='csc'), np.zeros(k), Acone, h,
                                     cones, _settings()).solve()
        if base.status != clarabel.SolverStatus.Solved:
            raise ValueError(f'minimum effort failed: {base.status}')
        emin = float(p@p + np.asarray(base.x)@np.asarray(base.x))
        radius2 = (1+eps)*emin - float(p@p)
        if radius2 < -1e-8: raise ValueError('negative energy radius')
        soc = sp.vstack([sp.csc_matrix((1,k)), -sp.eye(k,format='csc')],format='csc')
        Ac = sp.vstack([Acone,soc],format='csc')
        hc = np.r_[h,np.sqrt(max(0.,radius2)),np.zeros(k)]
        cc = [*cones,clarabel.SecondOrderConeT(k+1)]
        ends=[]
        for sign in (1.,-1.):
            sol=clarabel.DefaultSolver(sp.csc_matrix((k,k)),sign*c/bw_scale,Ac,hc,cc,_settings()).solve()
            if sol.status not in (clarabel.SolverStatus.Solved,clarabel.SolverStatus.AlmostSolved):
                raise ValueError(f'band failed: {sol.status}')
            z=np.asarray(sol.x)
            ends.append((cp+c@z,(p+N@z)*F0))
    for value, f in ends:
        x=f/F0
        if (np.linalg.norm(np.asarray(A)@f-np.asarray(b)) > 2e-5*max(1.,np.linalg.norm(b))
                or x.min() < -2e-5 or x.max() > 1+2e-5
                or any(offset+np.asarray(v)@f < -2e-4*bw_scale for v,offset in (med,lat))
                or (emin is not None and x@x > (1+eps)*emin+2e-5*max(1.,emin))):
            raise ValueError('solver returned an invalid force witness')
    return {'lower_N':float(ends[0][0]), 'upper_N':float(ends[1][0]),
            'lower_f':ends[0][1], 'upper_f':ends[1][1],
            'emin':emin, 'contact_rank':rank, 'nullity':k}
