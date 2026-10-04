import numpy as np
import time, resource
from pathlib import Path
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import nnls, root
from . import loader, ncp

def labels_adapter(lam, G, b, mu, tol_lam=1e-10, tol_rel=1e-06):
    l = np.asarray(lam).reshape(-1, 3)
    scale = max(np.max(np.abs(l)), 1e-30)
    ln = l[:, 0]
    lt = np.linalg.norm(l[:, 1:], axis=1)
    return np.where(ln <= tol_lam + tol_rel * scale, 'open', np.where(lt >= mu * ln * (1 - tol_rel) - tol_lam, 'slip', 'stick'))

def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, Path):
        return str(x)
    return x

def operator(geom, support_scale=None, crown=True, rotation=True, mu=0.2):
    """Use P1 rebuild for incidence and for W. Minv here means K^-1, not mass."""
    patches = geom['patches']
    teeth = [(jaw, int(t)) for jaw in ('upper', 'lower') for t in sorted(geom['centers'][jaw], key=int)]
    lookup = {k: i for (i, k) in enumerate(teeth)}
    bodies = []
    for (i, (jaw, t)) in enumerate(teeth):
        scale = 1 if support_scale is None else support_scale.get((jaw, t), 1.0)
        c = np.array([1 / 75.0, 1 / 75.0, 1 / 750.0, 1 / 7500.0, 1 / 7500.0, 1 / 7500.0]) / scale
        if not rotation:
            c[3:] = 0
        bodies.append(dict(id=i, dof=6, minv_block=np.diag(c).tolist(), block_dim=6))
    contacts = []
    for (i, p) in enumerate(patches):
        a = lookup['upper', p['upper_fdi']]
        b = lookup['lower', p['lower_fdi']]
        contacts.append(dict(id=i, a=a, b=b, normal=p['basis'][0], t1=p['basis'][1], t2=p['basis'][2], lever_a=(np.asarray(p['upper_point_mm']) - geom['centers']['upper'][str(p['upper_fdi'])]).tolist(), lever_b=(np.asarray(p['lower_point_mm']) - geom['centers']['lower'][str(p['lower_fdi'])]).tolist()))
    scene = dict(model='body_contacts', n_c=len(patches), bodies=bodies, contacts=contacts, eta=[0.0] * len(patches), b=[0.0] * (3 * len(patches)), mu=[mu] * len(patches))
    (_, _, _, J, C) = loader.rebuild(scene)
    radius = np.sqrt(np.array([p['area_mm2'] for p in patches]) / np.pi)
    E = 97000.0
    nu = 0.3
    shear = E / (2 * (1 + nu))
    cn = 2 * (1 - nu ** 2) / (2 * radius * E)
    ct = 2 * (2 - nu) / (8 * radius * shear)
    ce = np.column_stack([cn, ct, ct]).ravel() if crown else np.zeros(3 * len(patches))
    JJ = np.column_stack([J, np.eye(3 * len(patches))])
    CC = np.zeros((JJ.shape[1], JJ.shape[1]))
    CC[:len(C), :len(C)] = C
    CC[len(C):, len(C):] = np.diag(ce)
    s = dict(model='operator', n_c=len(patches), operator=dict(J=JJ.tolist(), Minv=CC.tolist()), eta=[0.0] * len(patches), b=[0.0] * (3 * len(patches)), mu=[mu] * len(patches))
    (W, _, mus, _, _) = loader.rebuild(s)
    return (W, J, C, ce, mus, teeth)

def load_vector(geom, closure=0.05, lateral=0.0, offset=0.0):
    b = []
    for p in geom['patches']:
        Dd = np.asarray(p['basis'])
        disp = Dd @ np.array([-lateral, 0.0, -closure])
        disp[0] += p['gap_mm'] + offset * Dd[0, 2]
        b.extend(disp)
    return np.array(b)

def normal_nnls(W, b):
    """Independent frictionless QP: L L.T=A, min||L.T*x+L^-1*b||²."""
    inds = np.arange(0, len(b), 3)
    A = W[np.ix_(inds, inds)]
    bb = b[inds]
    L = cholesky(A, lower=True)
    rhs = solve_triangular(L, -bb, lower=True)
    (x, _) = nnls(L.T, rhs, maxiter=10000, atol=1e-12)
    lam = np.zeros(len(b))
    lam[inds] = x
    return lam

def solve(geom, closure=0.05, lateral=0.0, offset=0.0, mu=0.2, support_scale=None, crown=True, rotation=True, controls=False, iters=2000):
    st = time.perf_counter()
    cpu = time.process_time()
    (W, J, C, ce, mus, teeth) = operator(geom, support_scale, crown, rotation, mu)
    b = load_vector(geom, closure, lateral, offset)
    prep = time.perf_counter() - st
    tick = time.perf_counter()
    warm = normal_nnls(W, b)
    (lam, hist, labels) = ncp.solve_ncp_pgs(W, b, mus, iters=iters, tol=1e-07, warm=warm)
    elapsed = time.perf_counter() - tick
    L = lam.reshape(-1, 3)
    gap = (W @ lam + b).reshape(-1, 3)
    active = L[:, 0] > 1e-07
    rho = ncp._rho_of(W, len(mus))
    res = ncp.natural_residual(lam, W, b, mus, rho)
    nc = float(max(0, -gap[:, 0].min()))
    comp = float(np.max(np.abs(L[:, 0] * gap[:, 0])))
    cone = float(max(0, np.max(np.linalg.norm(L[:, 1:], axis=1) - mu * L[:, 0])))
    forces = np.array([np.asarray(p['basis']).T @ f for (p, f) in zip(geom['patches'], L)])
    normalforce = {jaw: {} for jaw in ('upper', 'lower')}
    verticalforce = {jaw: {} for jaw in ('upper', 'lower')}
    vectors = {jaw: {} for jaw in ('upper', 'lower')}
    for (jaw, t) in teeth:
        ids = np.array([p[jaw + '_fdi'] == t for p in geom['patches']])
        normalforce[jaw][str(t)] = float(L[ids, 0].sum())
        verticalforce[jaw][str(t)] = float(forces[ids, 2].sum())
        vectors[jaw][str(t)] = ((1 if jaw == 'upper' else -1) * forces[ids].sum(0)).tolist()
    shares = {jaw: {t: 100 * f / max(sum(normalforce[jaw].values()), 1e-30) for (t, f) in normalforce[jaw].items()} for jaw in ('upper', 'lower')}
    wrench = J.T @ lam
    q = C @ wrench
    rigid_q = q.reshape(-1, 6)
    rows = [dict(**p, lambda_N=f, vector_on_upper_N=force, post_gap_mm=d, state=str(s), mean_normal_traction_MPa=float(f[0] / p['area_mm2']), traction_scope='regional equivalent-area mean, not peak pressure') for (p, f, force, d, s) in zip(geom['patches'], L, forces, gap, labels)]
    equal = float(np.linalg.norm(sum((np.asarray(v) for jaw in vectors.values() for v in jaw.values()), np.zeros(3))))
    gates = dict(natural_residual=res <= 1e-06, nonpenetration=nc <= 1e-06, complementarity=comp <= 1e-05, coulomb_cone=cone <= 1e-07, equal_opposite=equal <= 1e-08)
    ctl = {}
    if controls:
        t = time.perf_counter()
        zero = np.zeros_like(mus)
        (a, h, s) = ncp.solve_ncp_pgs(W, b, zero, iters=2000, tol=1e-07, warm=warm)
        ctl['frictionless_NNLS'] = dict(max_difference_N=float(np.max(np.abs(a - warm))), natural_residual_N=ncp.natural_residual(warm, W, b, zero), wall_s=time.perf_counter() - t, pass_gate=bool(np.max(np.abs(a - warm)) <= 0.0001))

        def mapping(x):
            u = (W @ x + b).reshape(-1, 3)
            xx = x.reshape(-1, 3)
            return (xx - ncp.proj_cone(xx - rho[:, None] * (u + ncp.desaxce(u, mus)), mus)).ravel()
        t = time.perf_counter()
        guess = lam + np.random.default_rng(54).normal(0, 0.01, len(lam)) * max(1, np.max(np.abs(lam)))
        rr = root(mapping, guess, method='krylov', options=dict(fatol=1e-07, maxiter=100))
        diff = float(np.max(np.abs(rr.x - lam)))
        rrres = float(np.max(np.abs(mapping(rr.x))))
        ctl['independent_scipy_root'] = dict(max_difference_N=diff, natural_residual_N=rrres, success=bool(rr.success), wall_s=time.perf_counter() - t, pass_gate=bool(diff <= 0.0001 and rrres <= 1e-06))
        bad = lam.copy()
        bad[int(np.argmax(np.abs(bad)))] += max(1, 0.2 * np.max(np.abs(lam)))
        badres = float(np.max(np.abs(mapping(bad))))
        ctl['injected_force_error'] = dict(natural_residual_N=badres, rejected=bool(badres > 1e-06))
    out = dict(case=geom['case'], closure_mm=closure, lateral_mm=lateral, gap_offset_scenario_mm=offset, mu=mu, normal_force_N=normalforce, vertical_force_N=verticalforce, force_shares_pp=shares, net_force_on_teeth_N=vectors, total_normal_N=float(L[:, 0].sum()), total_vertical_N=float(forces[:, 2].sum()), argmax={jaw: max(shares[jaw], key=shares[jaw].get) for jaw in shares}, rows=rows, support_motion=[dict(jaw=jaw, fdi=t, translation_mm=rigid_q[i, :3], rotation_rad=rigid_q[i, 3:], wrench_N_Nmm=wrench.reshape(-1, 6)[i]) for (i, (jaw, t)) in enumerate(teeth)], numerics=dict(natural_residual_N=res, minimum_normal_gap_mm=float(gap[:, 0].min()), normal_nonpenetration_violation_mm=nc, normal_complementarity_Nmm=comp, cone_violation_N=cone, equal_opposite_N=equal, iterations=len(hist), gates=gates, all_pass=all(gates.values()), W_min_eigenvalue_mm_per_N=float(np.linalg.eigvalsh(W)[0])), controls=ctl, cost=dict(operator_s=prep, solve_s=elapsed, total_s=time.perf_counter() - st, cpu_s=time.process_time() - cpu, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024), resolution='PER_SURFACE_REGION', claim_type=['information_link', 'capability'], physical_status='CONDITIONAL_SIMULATION; support/load/preload/FDI/sensor UNKNOWN', branches={str(k): int(np.sum(labels == k)) for k in np.unique(labels)}, rigorous_error_enclosure='MISSING: sampled geometry, reduced elasticity and floating NCP residual are not a continuum or parameter-box enclosure')
    return (clean(out), dict(W=W, b=b, lambda_N=lam, J=J, C=C, crown_C=ce, mu=mus, residual_history=hist))
ncp.classify = labels_adapter
