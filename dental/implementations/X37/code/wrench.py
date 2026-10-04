"""Generalized multi-tooth load, reusing X27 R8 energy and its two solvers.
No edits to predecessor. Region/membrane closures retain predecessor limitations.
"""
from dental_release.paths import expand as _release_expand
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize, LinearConstraint
PARENT = Path(_release_expand('@DENTAL_INPUT_ROOT@/artifacts/LANE_X27_ALIGNER_REGIONAL/code'))
sys.path.insert(0, str(PARENT))
from whole_arch import build

def solve_load(s, b, method):
    (K, A, W, C) = (s['K'], s['A'], s['W'], s['C'])
    if method == 'dual':
        opt = minimize(lambda l: 0.5 * l @ W @ l - b @ l, np.zeros(len(b)), jac=lambda l: W @ l - b, method='SLSQP', bounds=[(0, None)] * len(b), constraints=[LinearConstraint(C, np.zeros(len(C)), np.zeros(len(C)))], options={'ftol': 1e-13, 'maxiter': 1500})
        if not opt.success:
            raise RuntimeError('dual: ' + opt.message)
        mask = opt.x > 1e-07
        seen = set()
        for it in range(100):
            key = tuple(np.where(mask)[0])
            if key in seen:
                raise RuntimeError('dual active set cycle')
            seen.add(key)
            ids = np.where(mask)[0]
            cj = C[:, ids]
            kk = np.block([[W[np.ix_(ids, ids)], cj.T], [cj, np.zeros((len(C), len(C)))]])
            sol = np.linalg.lstsq(kk, np.r_[b[ids], np.zeros(len(C))], rcond=1e-12)[0]
            lam = np.zeros(len(b))
            lam[ids] = sol[:len(ids)]
            rigid = sol[len(ids):]
            if lam.min() < -1e-10:
                mask[np.argmin(lam)] = False
                continue
            gap = W @ lam + C.T @ rigid - b
            if gap.min() < -1e-09:
                mask[np.argmin(gap)] = True
                continue
            if np.max(np.abs(C @ lam)) > 1e-09:
                raise RuntimeError('dual equilibrium')
            break
        else:
            raise RuntimeError('dual maximum iteration')
        x = s['P'] @ A.T @ lam + s['N'] @ rigid
        iterations = int(opt.nit) + it + 1
    else:
        (vals, V) = eigh(K)
        pos = vals > vals[-1] * 1e-10
        npos = int(pos.sum())
        T = np.column_stack([V[:, pos] / np.sqrt(vals[pos]), s['N']])
        B = A @ T
        start = np.zeros(len(K))
        for k in range(5):
            start[6 + k::11] = np.maximum(b[k::5], 0)
        v0 = np.r_[np.sqrt(vals[pos]) * (V[:, pos].T @ start), s['N'].T @ start]
        opt = minimize(lambda v: 0.5 * v[:npos] @ v[:npos], v0, jac=lambda v: np.r_[v[:npos], np.zeros(s['null_modes'])], method='SLSQP', constraints=[LinearConstraint(B, b, np.full(len(b), np.inf))], options={'ftol': 1e-13, 'maxiter': 1500})
        if not opt.success:
            raise RuntimeError('primal: ' + opt.message)
        x = T @ opt.x
        active = np.abs(A @ x - b) < 1e-07
        lam = np.zeros(len(b))
        lam[active] = np.linalg.lstsq(A[active].T, K @ x, rcond=None)[0]
        iterations = int(opt.nit)
    gap = A @ x - b
    forces = {str(i): np.zeros(6) for i in s['ids']}
    F = np.zeros(3)
    M = np.zeros(3)
    for (patch, l) in zip(s['patch'], lam):
        f = -l * patch['normal']
        m = np.cross(patch['lever'], f)
        forces[str(patch['fdi'])] += np.r_[f, m]
        F += f
        M += np.cross(patch['point'], f)
    return {'wrenches': {k: v.tolist() for (k, v) in forces.items()}, 'iterations': iterations, 'minimum_gap_mm': float(gap.min()), 'minimum_multiplier_N': float(lam.min()), 'force_balance_N': float(np.linalg.norm(F)), 'moment_balance_Nmm': float(np.linalg.norm(M)), 'complementarity_Nmm': float(np.max(np.abs(lam * gap))), 'null_modes': s['null_modes']}

def gate(a, b):
    df = max((np.max(np.abs(np.array(a['wrenches'][k])[:3] - np.array(b['wrenches'][k])[:3])) for k in a['wrenches']))
    dm = max((np.max(np.abs(np.array(a['wrenches'][k])[3:] - np.array(b['wrenches'][k])[3:])) for k in a['wrenches']))
    numerical = all((r['minimum_gap_mm'] >= -1e-07 and r['minimum_multiplier_N'] >= -1e-08 and (r['force_balance_N'] <= 1e-06) and (r['moment_balance_Nmm'] <= 1e-05) for r in [a, b])) and df <= 1e-05 and (dm <= 1e-05)
    return {'pass': bool(numerical), 'max_force_component_difference_N': float(df), 'max_moment_component_difference_Nmm': float(dm)}
