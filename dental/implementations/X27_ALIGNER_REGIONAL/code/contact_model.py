"""Regional compliant-patch model. A research reduction, not continuum shell FE.

Hard contact at fixed normals; poses and membrane/anchor compliance are explicit
unmeasured closures. All wrenches are in the source arch frame, N and N mm.
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import nnls, minimize, LinearConstraint
REGIONS = ['Buccogingival', 'Buccal', 'Incisal/occlusal', 'Palatal', 'Palatogingival']

def skew(v):
    (a, b, c) = v
    return np.array([[0, -c, b], [c, 0, -a], [-b, a, 0.0]])

def neighbours(fdi, teeth):
    (q, p) = divmod(fdi, 10)
    prev = fdi - 1 if p > 1 else (q + 1 if q % 2 else q - 1) * 10 + 1
    return [k for k in [prev, fdi + 1] if str(k) in teeth and k % 10 <= 8]

def thickness(region, records, kind, position, quantile='median_mm'):
    if kind == 'nominal_contact':
        return 0.75
    key = ('Anterior teeth' if position <= 3 else 'Posterior teeth') if kind == 'tooth_type_only' else region
    return next((r[quantile] for r in records if r['arm'] == 'TS' and r['quantity'] == 'finished_thickness' and (r['region'] == key)))

def make_system(arch, active, records, kind='regional', quantile='median_mm', E=2746.0):
    teeth = arch['teeth']
    ids = [active] + neighbours(active, teeth)
    if len(ids) < 3:
        raise ValueError('Two immediate source-labelled neighbours required')
    if any((teeth[str(i)]['crown_height_mm'] <= 1 for i in ids)):
        raise ValueError('Crown height <=1 mm; geometry reduction refused')
    n = len(ids)
    nd = 11 * n
    K = np.zeros((nd, nd))
    ground = np.zeros_like(K)
    A = np.zeros((5 * n, nd))
    patch = []
    z = np.array(arch['frame']['occlusal_unit'])
    z /= np.linalg.norm(z)
    centers = [np.array(teeth[str(i)]['centroid_mm']) for i in ids]
    heights = [teeth[str(i)]['crown_height_mm'] for i in ids]
    mean_h = np.mean([thickness(r, records, kind, active % 10, quantile) for r in REGIONS])
    for (j, i) in enumerate(ids):
        base = 11 * j
        H = heights[j]
        c = centers[j]
        e = np.array(teeth[str(i)]['buccal_unit'])
        e -= z * (e @ z)
        e /= np.linalg.norm(e)
        spans = [np.linalg.norm(c - centers[k]) for k in range(n) if k != j]
        width = float(min(spans)) / 2
        if width < 1:
            raise ValueError('Submillimetre span: not acceptable tooth geometry')
        k_anchor = E * H * mean_h ** 3 / (4 * (2 * width) ** 3)
        g = np.diag([k_anchor] * 3 + [k_anchor * H ** 2 / 12] * 3)
        ground[base:base + 6, base:base + 6] += g
        K[base:base + 6, base:base + 6] += g
        a = 0.35 * H
        points = [a * e - 0.3 * H * z, a * e + 0.2 * H * z, 0.5 * H * z, -a * e + 0.2 * H * z, -a * e - 0.3 * H * z]
        normals = [e - 0.2 * z, e, z, -e, -e - 0.2 * z]
        for (k, (region, r, normal)) in enumerate(zip(REGIONS, points, normals)):
            normal = normal / np.linalg.norm(normal)
            signature = np.r_[normal, np.cross(r, normal)]
            h = thickness(region, records, kind, i % 10, quantile)
            kp = E * (width / 5) * h ** 3 / (4 * (H / 2) ** 3)
            v = np.zeros(nd)
            v[base:base + 6] = -signature
            v[base + 6 + k] = 1
            K += kp * np.outer(v, v)
            A[5 * j + k, base + 6 + k] = 1
            gap = 0.0 if kind == 'nominal_contact' else next((t[quantile] for t in records if t['arm'] == 'TS' and t['quantity'] == 'passive_gap' and (t['region'] == region)))
            patch.append(dict(fdi=i, housing_index=j, region=region, normal=normal.tolist(), point_mm=(c + r).tolist(), lever_mm=r.tolist(), signature=signature.tolist(), thickness_mm=h, gap_mm=gap, k_patch_N_per_mm=kp, resolution_level='PER_SURFACE_REGION', state_status='CONDITIONAL_REGION_SCENARIO_NOT_A_MEASURED_TOOTH_MAP'))
    for j in range(1, n):
        mid = (centers[0] + centers[j]) / 2
        span = np.linalg.norm(centers[j] - centers[0])
        width = (heights[0] + heights[j]) / 2
        kt = E * width * mean_h ** 3 / (4 * span ** 3)
        kr = E * width * mean_h ** 3 / (12 * span)
        D = np.zeros((3, nd))
        D[:, :3] = np.eye(3)
        D[:, 3:6] = -skew(mid - centers[0])
        D[:, 11 * j:11 * j + 3] = -np.eye(3)
        D[:, 11 * j + 3:11 * j + 6] = skew(mid - centers[j])
        K += kt * D.T @ D
        T = np.zeros((3, nd))
        T[:, 3:6] = np.eye(3)
        T[:, 11 * j + 3:11 * j + 6] = -np.eye(3)
        K += kr * T.T @ T
    factor = cho_factor(K, lower=True)
    W = A @ cho_solve(factor, A.T)
    return dict(K=K, ground=ground, A=A, factor=factor, W=W, patch=patch, ids=ids, centers=centers, z=z, kind=kind, E_MPa=E, active=active, arch=arch, quantile=quantile)

def solve(system, activation=0.2, offsets=None, seating_N=0.0, method='dual', relaxation_factor=1.0, tooth_pose=None):
    n = len(system['ids'])
    offsets = np.zeros(n) if offsets is None else np.asarray(offsets)
    K = system['K'] * relaxation_factor
    A = system['A']
    factor = cho_factor(K, lower=True)
    f = np.zeros(len(K))
    for j in range(n):
        f[11 * j:11 * j + 3] = -system['z'] * seating_N / n
    tooth = system['arch']['teeth'][str(system['active'])]
    activation_vec = -activation * np.array(tooth['buccal_unit'])
    if tooth_pose is None:
        b = np.array([np.dot(p['normal'], activation_vec) * (p['fdi'] == system['active']) - p['gap_mm'] - np.dot(p['normal'], system['z']) * offsets[p['housing_index']] for p in system['patch']])
    else:
        b = np.array([np.dot(p['signature'], tooth_pose) * (p['fdi'] == system['active']) - p['gap_mm'] - np.dot(p['normal'], system['z']) * offsets[p['housing_index']] for p in system['patch']])
    x0 = cho_solve(factor, f)
    W = A @ cho_solve(factor, A.T)
    be = b - A @ x0
    if method == 'dual':
        L = np.linalg.cholesky(W)
        lam = nnls(L.T, solve_triangular(L, be, lower=True), maxiter=500)[0]
        x = cho_solve(factor, f + A.T @ lam)
    else:
        L = np.linalg.cholesky(K)
        C = solve_triangular(L, A.T, lower=True).T
        fy = solve_triangular(L, f, lower=True)
        res = minimize(lambda y: 0.5 * y @ y - fy @ y, fy, jac=lambda y: y - fy, method='SLSQP', constraints=[LinearConstraint(C, b, np.full(len(b), np.inf))], options=dict(ftol=1e-12, maxiter=500))
        if not res.success:
            raise RuntimeError('Primal control: ' + res.message)
        x = solve_triangular(L.T, res.x, lower=False)
        active = np.abs(A @ x - b) < 1e-07
        lam = np.zeros(len(b))
        lam[active] = np.linalg.lstsq(A[active].T, K @ x - f, rcond=None)[0]
    gap = A @ x - b
    wrenches = {str(i): np.zeros(6) for i in system['ids']}
    shell_force = np.zeros(3)
    shell_torque = np.zeros(3)
    regional = []
    for (k, (p, l)) in enumerate(zip(system['patch'], lam)):
        F = -l * np.array(p['normal'])
        M = np.cross(p['lever_mm'], F)
        wrenches[str(p['fdi'])] += np.r_[F, M]
        shell_force -= F
        shell_torque -= np.cross(p['point_mm'], F)
        regional.append(dict(fdi=p['fdi'], region=p['region'], reaction_N=float(l), gap_mm=float(gap[k]), thickness_mm=p['thickness_mm'], force_N=F.tolist(), moment_Nmm=M.tolist(), contact=bool(l > 1e-08)))
    supports = -(system['ground'] * relaxation_factor) @ x
    for (j, c) in enumerate(system['centers']):
        force = supports[11 * j:11 * j + 3] + f[11 * j:11 * j + 3]
        shell_force += force
        shell_torque += np.cross(c, force) + supports[11 * j + 3:11 * j + 6]
    Fa = wrenches[str(system['active'])][:3]
    e = np.array(tooth['buccal_unit'])
    return dict(kind=system['kind'], activation_mm=activation, seating_N=seating_N, offsets_mm=offsets.tolist(), active_fdi=system['active'], wrenches={k: v.tolist() for (k, v) in wrenches.items()}, regional=regional, active_buccolingual_N=float(abs(Fa @ e)), active_force_norm_N=float(np.linalg.norm(Fa)), active_axial_N=float(Fa @ system['z']), moment_origin='Undeformed source mesh crown centroid; linearized virtual work', residuals=dict(penetration_mm=max(0.0, -float(gap.min())), complementarity_Nmm=float(np.max(np.abs(lam * gap))), stationarity_N=float(np.linalg.norm(K @ x - f - A.T @ lam, np.inf)), reaction_negativity_N=max(0.0, -float(lam.min())), global_force_balance_N=float(np.linalg.norm(shell_force)), global_torque_balance_Nmm=float(np.linalg.norm(shell_torque))), resolution_level='PER_TOOTH', contact_resolution_level='PER_SURFACE_REGION', time_scale='SIMULTANEOUS', physical_status='UNKNOWN_UNMATCHED_SPECIMEN_AND_CLOSURES', housing_state=x.tolist())

def insertion(system, path='simultaneous', seating_N=0):
    n = len(system['ids'])
    trace = []
    if path == 'simultaneous':
        schedules = [np.full(n, s) for s in [0.6, 0.4, 0.2, 0.1, 0]]
    else:
        state = np.full(n, 0.6)
        schedules = []
        order = list(range(n)) if path == 'front_first' else list(range(n - 1, -1, -1))
        for j in order:
            for s in [0.4, 0.2, 0.1, 0]:
                state = state.copy()
                state[j] = s
                schedules.append(state)
    for (stage, s) in enumerate(schedules):
        r = solve(system, offsets=s, seating_N=seating_N)
        r['stage'] = stage
        r['path'] = path
        trace.append(r)
    return trace
