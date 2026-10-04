"""Whole-arch floating shell reduction. Internal energy only, six rigid modes."""
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize, LinearConstraint
from contact_model import REGIONS, thickness, skew

def build(arch, records):
    teeth = arch['teeth']
    order = [i for i in list(range(28, 20, -1)) + list(range(11, 19)) if str(i) in teeth]
    n = len(order)
    nd = 11 * n
    K = np.zeros((nd, nd))
    A = np.zeros((5 * n, nd))
    patch = []
    z = np.array(arch['frame']['occlusal_unit'])
    z /= np.linalg.norm(z)
    centers = [np.array(teeth[str(i)]['centroid_mm']) for i in order]
    heights = [teeth[str(i)]['crown_height_mm'] for i in order]
    mean_h = np.mean([thickness(r, records, 'regional', 1) for r in REGIONS])
    for (j, i) in enumerate(order):
        H = heights[j]
        c = centers[j]
        e = np.array(teeth[str(i)]['buccal_unit'])
        e -= z * (e @ z)
        e /= np.linalg.norm(e)
        if H <= 1:
            raise ValueError('Invalid source crown height')
        spans = [np.linalg.norm(c - centers[k]) for k in [j - 1, j + 1] if 0 <= k < n]
        width = min(spans) / 2
        if width < 1:
            raise ValueError('Submillimetre source neighbour span')
        a = 0.35 * H
        points = [a * e - 0.3 * H * z, a * e + 0.2 * H * z, 0.5 * H * z, -a * e + 0.2 * H * z, -a * e - 0.3 * H * z]
        normals = [e - 0.2 * z, e, z, -e, -e - 0.2 * z]
        for (k, (region, r, normal)) in enumerate(zip(REGIONS, points, normals)):
            normal = normal / np.linalg.norm(normal)
            signature = np.r_[normal, np.cross(r, normal)]
            h = thickness(region, records, 'regional', i % 10)
            kp = 2746 * (width / 5) * h ** 3 / (4 * (H / 2) ** 3)
            v = np.zeros(nd)
            v[11 * j:11 * j + 6] = -signature
            v[11 * j + 6 + k] = 1
            K += kp * np.outer(v, v)
            A[5 * j + k, 11 * j + 6 + k] = 1
            gap = next((r['median_mm'] for r in records if r['arm'] == 'TS' and r['quantity'] == 'passive_gap' and (r['region'] == region)))
            patch.append(dict(fdi=i, normal=normal, point=c + r, lever=r, gap=gap, signature=signature))
    for j in range(n - 1):
        mid = (centers[j] + centers[j + 1]) / 2
        span = np.linalg.norm(centers[j + 1] - centers[j])
        width = (heights[j] + heights[j + 1]) / 2
        kt = 2746 * width * mean_h ** 3 / (4 * span ** 3)
        kr = 2746 * width * mean_h ** 3 / (12 * span)
        D = np.zeros((3, nd))
        D[:, 11 * j:11 * j + 3] = np.eye(3)
        D[:, 11 * j + 3:11 * j + 6] = -skew(mid - centers[j])
        D[:, 11 * (j + 1):11 * (j + 1) + 3] = -np.eye(3)
        D[:, 11 * (j + 1) + 3:11 * (j + 1) + 6] = skew(mid - centers[j + 1])
        K += kt * D.T @ D
        T = np.zeros((3, nd))
        T[:, 11 * j + 3:11 * j + 6] = np.eye(3)
        T[:, 11 * (j + 1) + 3:11 * (j + 1) + 6] = -np.eye(3)
        K += kr * T.T @ T
    (vals, V) = eigh(K)
    tol = vals[-1] * 1e-10
    positive = vals > tol
    N = V[:, ~positive]
    P = V[:, positive] / vals[positive] @ V[:, positive].T
    W = A @ P @ A.T
    C = N.T @ A.T
    return dict(arch=arch, ids=order, K=K, A=A, W=W, C=C, N=N, P=P, patch=patch, z=z, null_modes=N.shape[1], minimum_eigenvalue=float(vals[0]))

def solve_floating(s, active, method='dual'):
    e = np.array(s['arch']['teeth'][str(active)]['buccal_unit'])
    d = -0.2 * e
    b = np.array([np.dot(p['normal'], d) * (p['fdi'] == active) - p['gap'] for p in s['patch']])
    K = s['K']
    A = s['A']
    W = s['W']
    C = s['C']
    if method == 'dual':
        res = minimize(lambda l: 0.5 * l @ W @ l - b @ l, np.zeros(len(b)), jac=lambda l: W @ l - b, method='SLSQP', bounds=[(0, None)] * len(b), constraints=[LinearConstraint(C, np.zeros(len(C)), np.zeros(len(C)))], options=dict(ftol=1e-13, maxiter=1500))
        if not res.success:
            raise RuntimeError('Floating dual: ' + res.message)
        lam = res.x
        positive = lam > 1e-07
        rigid = np.linalg.lstsq(C[:, positive].T, (b - W @ lam)[positive], rcond=None)[0]
        x = s['P'] @ A.T @ lam + s['N'] @ rigid
    else:
        start = np.zeros(len(K))
        start[6::11] = np.maximum(b[::5], 0)
        start[7::11] = np.maximum(b[1::5], 0)
        start[8::11] = np.maximum(b[2::5], 0)
        start[9::11] = np.maximum(b[3::5], 0)
        start[10::11] = np.maximum(b[4::5], 0)
        res = minimize(lambda x: 0.5 * x @ K @ x, start, jac=lambda x: K @ x, method='SLSQP', constraints=[LinearConstraint(A, b, np.full(len(b), np.inf))], options=dict(ftol=1e-13, maxiter=1500))
        if not res.success:
            raise RuntimeError('Floating primal: ' + res.message)
        x = res.x
        mask = np.abs(A @ x - b) < 1e-06
        lam = np.zeros(len(b))
        lam[mask] = np.linalg.lstsq(A[mask].T, K @ x, rcond=None)[0]
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
    return dict(active_fdi=active, wrenches={k: v.tolist() for (k, v) in wrenches.items()}, regional=regional, energy_Nmm=float(0.5 * x @ K @ x), null_modes=s['null_modes'], solver_iterations=int(res.nit), minimum_gap_mm=float(gap.min()), force_balance_N=float(np.linalg.norm(Ftotal)), torque_balance_Nmm=float(np.linalg.norm(Mtotal)), physical_status='UNKNOWN_UNMATCHED_GAP_GEOMETRY_CONSTITUTIVE', external_housing_ground_supports=0, moment_origin='Undeformed labelled crown centroid', resolution_level='PER_TOOTH', time_scale='SIMULTANEOUS')
