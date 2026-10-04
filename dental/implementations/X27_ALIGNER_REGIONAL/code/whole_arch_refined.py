"""R8 same floating QP: whitened primal and active KKT dual refinement."""
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize, LinearConstraint
from whole_arch import build

def solve_floating(s, active, method='dual'):
    e = np.array(s['arch']['teeth'][str(active)]['buccal_unit'])
    d = -0.2 * e
    b = np.array([p['normal'] @ d * (p['fdi'] == active) - p['gap'] for p in s['patch']])
    K = s['K']
    A = s['A']
    W = s['W']
    C = s['C']
    if method == 'dual':
        res = minimize(lambda l: 0.5 * l @ W @ l - b @ l, np.zeros(len(b)), jac=lambda l: W @ l - b, method='SLSQP', bounds=[(0, None)] * len(b), constraints=[LinearConstraint(C, np.zeros(len(C)), np.zeros(len(C)))], options=dict(ftol=1e-13, maxiter=1500))
        if not res.success:
            raise RuntimeError('R8 initial dual: ' + res.message)
        mask = res.x > 1e-07
        seen = set()
        for iteration in range(100):
            key = tuple(np.where(mask)[0])
            if key in seen:
                raise RuntimeError('R8 active set cycle')
            seen.add(key)
            j = np.where(mask)[0]
            cj = C[:, j]
            mat = np.block([[W[np.ix_(j, j)], cj.T], [cj, np.zeros((len(C), len(C)))]])
            solution = np.linalg.lstsq(mat, np.r_[b[j], np.zeros(len(C))], rcond=1e-12)[0]
            lam = np.zeros(len(b))
            lam[j] = solution[:len(j)]
            rigid = solution[len(j):]
            if lam.min() < -1e-10:
                mask[np.argmin(lam)] = False
                continue
            gap = W @ lam + C.T @ rigid - b
            if gap.min() < -1e-09:
                mask[np.argmin(gap)] = True
                continue
            if np.max(np.abs(C @ lam)) > 1e-09:
                raise RuntimeError('R8 dual equilibrium residual')
            break
        else:
            raise RuntimeError('R8 active set maximum iterations')
        x = s['P'] @ A.T @ lam + s['N'] @ rigid
        iterations = int(res.nit) + iteration + 1
    else:
        (vals, V) = eigh(K)
        positive = vals > vals[-1] * 1e-10
        np_ = int(positive.sum())
        T = np.column_stack((V[:, positive] / np.sqrt(vals[positive]), s['N']))
        B = A @ T
        start = np.zeros(len(K))
        for k in range(5):
            start[6 + k::11] = np.maximum(b[k::5], 0)
        v0 = np.r_[np.sqrt(vals[positive]) * (V[:, positive].T @ start), s['N'].T @ start]
        res = minimize(lambda v: 0.5 * v[:np_] @ v[:np_], v0, jac=lambda v: np.r_[v[:np_], np.zeros(s['null_modes'])], method='SLSQP', constraints=[LinearConstraint(B, b, np.full(len(b), np.inf))], options=dict(ftol=1e-13, maxiter=1500))
        if not res.success:
            raise RuntimeError('R8 whitened primal: ' + res.message)
        x = T @ res.x
        mask = np.abs(A @ x - b) < 1e-07
        lam = np.zeros(len(b))
        lam[mask] = np.linalg.lstsq(A[mask].T, K @ x, rcond=None)[0]
        iterations = int(res.nit)
    gap = A @ x - b
    wrenches = {str(i): np.zeros(6) for i in s['ids']}
    regional = []
    Ftotal = np.zeros(3)
    Mtotal = np.zeros(3)
    for (p, l, g) in zip(s['patch'], lam, gap):
        F = -l * p['normal']
        M = np.cross(p['lever'], F)
        wrenches[str(p['fdi'])] += np.r_[F, M]
        Ftotal += F
        Mtotal += np.cross(p['point'], F)
        regional.append(dict(fdi=p['fdi'], reaction_N=float(l), gap_mm=float(g), force_N=F.tolist()))
    return dict(active_fdi=active, wrenches={k: v.tolist() for (k, v) in wrenches.items()}, regional=regional, energy_Nmm=float(0.5 * x @ K @ x), null_modes=s['null_modes'], solver_iterations=iterations, minimum_gap_mm=float(gap.min()), minimum_multiplier_N=float(lam.min()), force_balance_N=float(np.linalg.norm(Ftotal)), torque_balance_Nmm=float(np.linalg.norm(Mtotal)), physical_status='UNKNOWN_UNMATCHED_GAP_GEOMETRY_CONSTITUTIVE', external_housing_ground_supports=0, moment_origin='Undeformed labelled crown centroid', resolution_level='PER_TOOTH', time_scale='SIMULTANEOUS')
