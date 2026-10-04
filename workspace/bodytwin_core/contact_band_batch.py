"""Batched L1 contact bands with verified LP active-set reuse.

``trials`` is a sequence of L1 trial keys or loaded npz-like mappings.  The
returned arrays follow the in-range mask.  Bitwise-identical problems reuse endpoints. Adjacent LP frames may reuse an
active set only after primal-dual KKT verification; all others are solved.
This is an exact solver shortcut, not an unverified pose approximation.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def _problem(z, i):
    f0 = np.asarray(z['F0'], float)
    aa = np.asarray(z['A'][i], float) * f0[None, :]
    bb = np.asarray(z['b'][i], float)
    bw = float(z['bwN'])
    cj = np.asarray(z['cj'][i], float) * f0
    mk = np.asarray(z['Mkx'][i], float) * f0
    c0 = float(z['C0'][i]); mx = float(z['mx'][i]); d = float(z['d'])
    med = .5*cj + mk/d; lat = .5*cj - mk/d
    m0 = .5*c0 - mx/d; l0 = .5*c0 + mx/d
    scale = np.maximum(np.linalg.norm(aa, axis=1), 1.)
    rows = np.vstack((aa / scale[:, None], med / bw, lat / bw))
    lower = np.r_[bb/scale, -m0/bw, -l0/bw]
    upper = np.r_[bb/scale, np.inf, np.inf]
    return rows, lower, upper, cj/bw, c0, f0



def _batch_problems(z, frames):
    """Vectorize all frame geometry and contact-row preprocessing in NumPy."""
    f0 = np.asarray(z['F0'], float)
    aa = np.asarray(z['A'][frames], float) * f0[None, None, :]
    bb = np.asarray(z['b'][frames], float)
    bw = float(z['bwN']); d = float(z['d'])
    cj = np.asarray(z['cj'][frames], float) * f0[None, :]
    mk = np.asarray(z['Mkx'][frames], float) * f0[None, :]
    c0 = np.asarray(z['C0'][frames], float)
    mx = np.asarray(z['mx'][frames], float)
    scale = np.maximum(np.linalg.norm(aa, axis=2), 1.)
    eq = aa / scale[:, :, None]
    med = .5*cj + mk/d; lat = .5*cj - mk/d
    rows = np.concatenate((eq, med[:,None,:]/bw, lat[:,None,:]/bw), axis=1)
    rhs = bb/scale
    lower = np.column_stack((rhs, -(.5*c0-mx/d)/bw, -(.5*c0+mx/d)/bw))
    upper = lower.copy(); upper[:,6:] = np.inf
    return [(rows[j], lower[j], upper[j], cj[j]/bw, float(c0[j]), f0)
            for j in range(len(frames))]


def _fingerprint(problem, eps):
    h = hashlib.sha256()
    for arr in problem:
        a = np.asarray(arr, dtype='<f8')
        h.update(np.asarray(a.shape, dtype='<i8').tobytes())
        h.update(a.tobytes())
    h.update(np.asarray([np.nan if eps is None else eps], dtype='<f8').tobytes())
    return h.digest()


class _HighsPair:
    """Persistent primal simplex LP pair; changed frame coefficients retain bases."""
    def __init__(self, n):
        import highspy
        self.highspy = highspy
        self.solvers = []
        self.last_rows = None
        self.n = n
        for sign in (1., -1.):
            s = highspy.Highs()
            s.setOptionValue('output_flag', False)
            s.setOptionValue('solver', 'simplex')
            s.setOptionValue('threads', 1)
            s.setOptionValue('primal_feasibility_tolerance', 1e-9)
            s.setOptionValue('dual_feasibility_tolerance', 1e-9)
            s.addCols(n, np.zeros(n), np.zeros(n), np.ones(n), 0,
                      np.zeros(n+1, np.int32), np.array([], np.int32), np.array([], float))
            starts = np.arange(0, 8*n+1, n, dtype=np.int32)
            s.addRows(8, np.zeros(8), np.zeros(8), 8*n, starts,
                      np.tile(np.arange(n, dtype=np.int32), 8), np.zeros(8*n))
            self.solvers.append((sign, s))

    def solve(self, problem):
        rows, lower, upper, cost, c0, f0 = problem
        if rows.shape != (8, self.n):
            raise ValueError('expected six equilibrium rows and 166 muscles')
        indices = np.arange(self.n, dtype=np.int32)
        ends = []
        for sign, s in self.solvers:
            for r in range(8):
                s.changeRowBounds(r, float(lower[r]), float(upper[r]))
            # HiGHS has no bulk coefficient-update API; this cost is measured.
            for r in range(8):
                for j in range(self.n):
                    if self.last_rows is None or rows[r, j] != self.last_rows[r, j]:
                        s.changeCoeff(r, j, float(rows[r, j]))
            s.changeColsCost(self.n, indices, np.asarray(sign*cost, float))
            s.run()
            status = s.getModelStatus()
            if status != self.highspy.HighsModelStatus.kOptimal:
                raise ValueError(f'HiGHS failed: {status}')
            x = np.asarray(s.getSolution().col_value)
            if (np.max(np.abs(rows[:6]@x-lower[:6])) > 2e-5 or
                    np.min(rows[6:]@x-lower[6:]) < -2e-5 or
                    x.min() < -2e-5 or x.max() > 1+2e-5):
                raise ValueError('HiGHS returned invalid witness')
            ends.append((float(cost@x), x))
        self.last_rows = rows.copy()
        return ends




def _active_endpoint(problem, previous_x, sign):
    """Transfer an LP active set and certify the new KKT conditions."""
    rows, lower, upper, cost, c0, f0 = problem
    old = np.asarray(previous_x)
    at_lo = old <= 1e-8; at_hi = old >= 1-1e-8
    free = ~(at_lo|at_hi)
    active_contact = rows[6:]@old-lower[6:] <= 1e-8
    active_rows = np.r_[np.arange(6), 6+np.flatnonzero(active_contact)]
    if int(free.sum()) != len(active_rows):
        return None
    mat = rows[np.ix_(active_rows, np.flatnonzero(free))]
    rhs = lower[active_rows]-rows[np.ix_(active_rows,np.flatnonzero(at_hi))].sum(axis=1)
    try:
        xf = np.linalg.solve(mat,rhs)
        dual = np.linalg.solve(mat.T,sign*cost[free])
    except np.linalg.LinAlgError:
        return None
    x=np.zeros(len(cost));x[at_hi]=1.;x[free]=xf
    if (not np.isfinite(x).all() or x.min() < -1e-8 or x.max() > 1+1e-8 or
            np.max(abs(rows[:6]@x-lower[:6])) > 1e-8 or
            np.min(rows[6:]@x-lower[6:]) < -1e-8 or
            np.any(dual[6:] < -1e-8)):
        return None
    reduced=sign*cost-rows[active_rows].T@dual
    if (np.any(reduced[at_lo] < -1e-8) or
            np.any(reduced[at_hi] > 1e-8) or
            np.max(abs(reduced[free])) > 1e-8):
        return None
    # A feasible primal and feasible dual with this active set coincide.
    return float(cost@x),x


def _eps_direct(problem, eps, bw):
    """Direct activation-space Clarabel QP/SOCP, including constrained E_min."""
    import clarabel
    import scipy.sparse as sp
    rows, lower, upper, cost, c0, f0 = problem
    n = len(cost)
    eq = sp.csc_matrix(rows[:6])
    au = sp.vstack((sp.csc_matrix(-rows[6:]), sp.eye(n,format='csc'),
                    -sp.eye(n,format='csc')),format='csc')
    A = sp.vstack((eq,au),format='csc')
    b = np.r_[lower[:6],-lower[6:],np.ones(n),np.zeros(n)]
    cones = [clarabel.ZeroConeT(6),clarabel.NonnegativeConeT(2+2*n)]
    settings = clarabel.DefaultSettings()
    settings.verbose=False;settings.max_iter=120
    settings.tol_gap_abs=1e-8;settings.tol_feas=1e-8
    base=clarabel.DefaultSolver(2*sp.eye(n,format='csc'),np.zeros(n),A,b,cones,settings).solve()
    if base.status != clarabel.SolverStatus.Solved:
        raise ValueError(f'minimum effort failed: {base.status}')
    xb=np.asarray(base.x);emin=float(xb@xb)
    soc=sp.vstack((sp.csc_matrix((1,n)),-sp.eye(n,format='csc')),format='csc')
    Ac=sp.vstack((A,soc),format='csc')
    bc=np.r_[b,np.sqrt((1+eps)*emin),np.zeros(n)]
    cc=[*cones,clarabel.SecondOrderConeT(n+1)]
    ends=[]
    for sign in (1.,-1.):
        sol=clarabel.DefaultSolver(sp.csc_matrix((n,n)),sign*cost,Ac,bc,cc,settings).solve()
        if sol.status not in (clarabel.SolverStatus.Solved,clarabel.SolverStatus.AlmostSolved):
            raise ValueError(f'epsilon endpoint failed: {sol.status}')
        x=np.asarray(sol.x)
        if (np.max(abs(rows[:6]@x-lower[:6]))>2e-5 or
                np.min(rows[6:]@x-lower[6:]) < -2e-4 or
                x.min() < -2e-5 or x.max() > 1+2e-5 or
                x@x > (1+eps)*emin+2e-5*max(1.,emin)):
            raise ValueError('Clarabel returned invalid witness')
        ends.append(float(c0/bw+cost@x))
    return tuple(ends)


def contact_band_batch(trials, eps=None, tol=.01):
    """Compute L1 contact bands for each trial's ``inr`` frames.

    ``eps=None`` selects LP. ``tol`` is in body weights; only exact identity
    or a KKT-certified active set is accepted, hence any nonnegative tolerance
    yields the same endpoints. Epsilon bands use direct Clarabel QP/SOCP.
    """
    if tol < 0 or (eps is not None and eps < 0):
        raise ValueError('tol and eps must be nonnegative')
    if isinstance(trials, (str, Path)) or hasattr(trials, 'files'):
        trials = [trials]
    locations_path = ROOT/'results/CX-INVERSEOC/locations.json'
    locations = json.loads(locations_path.read_text())['trials'] if locations_path.exists() else {}
    out = {}; cache = {}; counts = {'frames':0, 'identity_hits':0, 'exact_solves':0,
                                    'highs_basis_reuses':0, 'highs_failures':0, 'active_set_hits':0}
    start = time.perf_counter()
    for item in trials:
        key = item if isinstance(item, str) else (item.stem if isinstance(item, Path) else f'array_{len(out)}')
        z = (np.load(ROOT/'results/L1/prep'/f'{key}.npz') if isinstance(item, str)
             else np.load(item) if isinstance(item, Path) else item)
        frames = (np.array([r['frame'] for r in locations[key]['rows']], dtype=int)
                  if key in locations else np.flatnonzero(z['inr']))
        lo = np.full(len(frames), np.nan); hi = lo.copy()
        bw = float(z['bwN']); persistent = None; previous_ends = None
        problems = _batch_problems(z, frames)
        for j, i in enumerate(frames):
            p = problems[j]
            digest = _fingerprint(p, eps)
            counts['frames'] += 1
            if digest in cache:
                lo[j], hi[j] = cache[digest]
                counts['identity_hits'] += 1
                continue
            if eps is None:
                ends = None
                if previous_ends is not None:
                    candidate = [_active_endpoint(p, previous_ends[k][1], sign)
                                 for k, sign in enumerate((1., -1.))]
                    if all(v is not None for v in candidate):
                        ends = candidate
                        counts['active_set_hits'] += 1
                if ends is None:
                    try:
                        if persistent is None:
                            persistent = _HighsPair(len(z['F0']))
                        elif j > 0:
                            counts['highs_basis_reuses'] += 1
                        ends = persistent.solve(p)
                    except (ImportError, ValueError):
                        counts['highs_failures'] += 1
                        persistent = None
                        ends = _scipy_lp_with_witness(p)
                lo[j], hi[j] = [p[4]/bw + v[0] for v in ends]
                previous_ends = ends
            else:
                lo[j], hi[j] = _eps_direct(p, eps, bw)
            cache[digest] = (float(lo[j]), float(hi[j]))
            counts['exact_solves'] += 1
        out[key] = {'frame':frames, 'lower_BW':lo, 'upper_BW':hi, 'bw_N':bw}
    out['_meta'] = {**counts, 'elapsed_s':time.perf_counter()-start,
                    'certificate':'bitwise identity or primal-dual KKT for LP; exact Clarabel fallback otherwise'}
    return out



def _scipy_lp_with_witness(p):
    from scipy.optimize import linprog
    rows, lower, upper, cost, c0, f0 = p
    ans = []
    for sign in (1., -1.):
        sol = linprog(sign*cost, A_eq=rows[:6], b_eq=lower[:6],
                      A_ub=-rows[6:], b_ub=-lower[6:],
                      bounds=[(0., 1.)]*len(cost), method='highs')
        if not sol.success:
            raise ValueError(sol.message)
        ans.append((float(cost@sol.x), np.asarray(sol.x)))
    return ans


def _scipy_lp(p, bw):
    from scipy.optimize import linprog
    rows, lower, upper, cost, c0, f0 = p
    ans = []
    for sign in (1., -1.):
        sol = linprog(sign*cost, A_eq=rows[:6], b_eq=lower[:6],
                      A_ub=-rows[6:], b_ub=-lower[6:],
                      bounds=[(0., 1.)]*len(cost), method='highs')
        if not sol.success:
            raise ValueError(sol.message)
        ans.append(c0/bw + cost@sol.x)
    return tuple(map(float, ans))
