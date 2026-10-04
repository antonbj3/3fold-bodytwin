from common import *
from mechanics import *
from coupled_mpc import SurfaceTieElastic
from whole import three_mf
from scipy.sparse import coo_matrix
from scipy.optimize import minimize, LinearConstraint
import time, trimesh

def run():
    freeze(ROOT / 'PREREG_R5.json', dict(round='R5', claim_type='capability', capability='Re-optimize whole crown with a conforming virtual-work surface tie rather than nearest die-node coupling', obstacle='R4 nearest-node tie p95 displacement 0.299mm is larger than nominal cement film and can alter local stiffness', changed_operation='Project slave crown interface nodes to closest die triangles; barycentric MPC u_slave=sum w_i u_die_i; exact reduced equilibrium Q^T K Q and unchanged local shape adjoint', consumer='K39 whole-crown shape optimizer with physically declared compliant support', external_referent=dict(kind='published_dataset', locator=['doi:10.5281/zenodo.10597292', str(X1 / 'inputs/geometry/D1_model.npz')], compared_quantity='Independent STS tooth-derived crown/die surfaces; no measured optimized fracture endpoint', refutes_us=True), metrics=dict(interface_projection_p95_mm_max=0.06, partition_of_unity_error_max=1e-12, adjoint_relative_error_max=0.0001, equilibrium_relative_residual_max=1e-08, intaglio_identity_error=0.0, volume_relative_increase_max=1e-08, min_tetra_volume_ratio=0.5, both_load_force_ratio_min=1.02), criterion_scope='0.06mm is half the source 0.12mm grid step, a numerical interface gate, not cement measurement', optimizer_maxiter=20, strongest_equally_informed_control='Frozen R4 design under new identical barycentric tie; R3 rigid design likewise retained', falsifiers=['interface projection p95 >.06 fails geometric port', 'partition sum or virtual work mismatch', 'any load gain <=1.02 fails useful shape gain', 'absolute scale/milling/contact remain UNKNOWN without lab'], full_cost=dict(preparation='Triangle projection and sparse MPC, no remeshing', fit='Two sparse primal/adjoint solves per iteration', validation='FD + geometry + all prior faults retained', queries='Counted', fallback='Numerically closed prototype only; no physical calibration')))
    start = time.perf_counter()
    a = np.load(X1 / 'inputs/geometry/D1_model.npz')
    raw = np.load(DATA / 'WHOLE_RAW.npz')
    cV = a['cV']
    cT = a['cT']
    dV = a['dV']
    dT = a['dT']
    nc = len(cV)
    ne = len(cT)
    V = np.vstack([cV, dV])
    T = np.vstack([cT, dT + nc])
    modes = np.concatenate([raw['modes'], np.zeros((len(dV), 3, 4))])
    slaves = raw['fixed_nodes']
    grid = read(X1 / 'inputs/geometry/D1_grid.json')
    df = boundary(dT)
    dm = trimesh.Trimesh(dV, df, process=False)
    closest = []
    dist = []
    fi = []
    for ids in np.array_split(np.arange(len(slaves)), max(1, int(np.ceil(len(slaves) / 128)))):
        (pp, dd, ff) = trimesh.proximity.closest_point(dm, cV[slaves[ids]])
        closest.append(pp)
        dist.append(dd)
        fi.append(ff)
    closest = np.concatenate(closest)
    dist = np.concatenate(dist)
    fi = np.concatenate(fi)
    tri = df[fi]
    bary = trimesh.triangles.points_to_barycentric(dV[tri], closest)
    active = np.setdiff1d(np.arange(len(V)), slaves)
    mapping = np.full(len(V), -1, int)
    mapping[active] = np.arange(len(active))
    qr = list(active)
    qc = list(mapping[active])
    qv = list(np.ones(len(active)))
    qr.extend(np.repeat(slaves, 3))
    qc.extend(mapping[nc + tri].ravel())
    qv.extend(bary.ravel())
    Q = coo_matrix((qv, (qr, qc)), shape=(len(V), len(active))).tocsr()
    rowerr = float(np.max(abs(np.asarray(Q.sum(1)).ravel() - 1)))
    fixed = mapping[nc + np.flatnonzero(dV[:, 2] < grid['z_m'] - 3 + 0.03)]
    E = np.r_[np.full(ne, 201400.0), np.full(len(dT), 18000.0)]
    saved = np.load(X1 / 'inputs/fe/L005_lo_M1_k6_thin_3Y.npz')
    bn = np.unique(raw['faces'])
    top = bn[(cV[bn, 2] > grid['z_m'] + 1.0) & ~np.isin(bn, slaves)]
    models = {}
    baseJ = {}
    prior = {}
    r3 = np.asarray(read(ROOT / 'raw/R3.json')['coefficients_mm'])
    r4 = np.asarray(read(ROOT / 'raw/R4.json')['coefficients_mm'])
    for angle in [0, 30]:
        center = saved['lc__axial' if angle == 0 else 'lc__offaxis30'][0, :3]
        w = np.exp(-np.sum((cV[top, :2] - center[:2]) ** 2, axis=1) / (2 * 0.6 ** 2))
        w /= w.sum()
        load = np.zeros_like(V)
        load[top] = w[:, None] * np.array([np.sin(np.radians(angle)), 0, -np.cos(np.radians(angle))])
        m = SurfaceTieElastic(V, T, fixed, load, modes, ne, Q, E)
        k = str(angle)
        models[k] = m
        baseJ[k] = m.evaluate(np.zeros(4))['J']
        prior[k] = dict(R3_ratio=float(np.exp(baseJ[k] - m.evaluate(r3)['J'])), R4_ratio=float(np.exp(baseJ[k] - m.evaluate(r4)['J'])))
    (inv, B, v0, sign) = geometry(cV, cT)
    dv = []
    for j in range(4):
        dp = np.concatenate([np.zeros((ne, 4, 1)), raw['modes'][cT, :, j]], axis=2)
        dv.append(v0 * np.einsum('eij,eji->e', inv, dp))
    dv = np.array(dv).T
    M = raw['modes'][:, 2, :]
    lin = np.vstack([np.c_[M, np.zeros(len(M))], np.c_[dv, np.zeros(ne)], np.r_[dv.sum(0), 0][None, :]])

    def risk(y):
        return np.array([y[4] - (m.evaluate(y[:4])['J'] - baseJ[k]) for (k, m) in models.items()])

    def jac(y):
        return np.array([np.r_[-m.evaluate(y[:4])['gradient'], 1] for m in models.values()])
    sol = minimize(lambda y: (y[4], np.r_[np.zeros(4), 1]), np.zeros(5), jac=True, method='SLSQP', constraints=[LinearConstraint(lin, np.r_[np.full(len(M), -0.2), -0.5 * v0, -np.inf], np.r_[np.full(len(M), 0.2), np.full(ne, np.inf), 0]), dict(type='ineq', fun=risk, jac=jac)], bounds=[(-0.2, 0.2)] * 4 + [(None, None)], options=dict(maxiter=20, ftol=1e-08))
    contrasts = {}
    gc = []
    for (k, m) in models.items():
        rr = m.evaluate(sol.x[:4])
        contrasts[k] = dict(**prior[k], optimized_ratio=float(np.exp(baseJ[k] - rr['J'])), equilibrium_relative_residual=rr['equilibrium_relative_residual'])
        gc.append(dict(angle_deg=int(k), **gradient_check(m, sol.x[:4])))
    final = models['0'].evaluate(sol.x[:4])
    new = final['vertices'][:nc]
    vol = geometry(new, cT)[2]
    dest = ROOT / 'exports/whole_STS_surface_tie'
    dest.mkdir(parents=True, exist_ok=True)
    mesh = trimesh.Trimesh(new, raw['faces'], process=True)
    if mesh.volume < 0:
        mesh.invert()
    mesh.export(dest / 'optimized.stl')
    three_mf(mesh, dest / 'optimized.3mf')
    identity = float(np.max(abs(new[slaves] - cV[slaves])))
    out = dict(round='R5', claim_type='capability', optimizer_success=bool(sol.success), iterations=int(sol.nit), coefficients_mm=sol.x[:4], FEM_queries={k: m.calls for (k, m) in models.items()}, contrasts=contrasts, gradient_checks=gc, interface_projection_mm=dict(median=float(np.median(dist)), p95=float(np.quantile(dist, 0.95)), max=float(dist.max()), resolution='PER_POINT'), partition_of_unity_error=rowerr, intaglio_identity_error_mm=identity, volume_relative_change=float(vol.sum() / v0.sum() - 1), min_tetra_volume_ratio=float(np.min(vol / v0)), calibrated_force05_N=None, physical_Pareto='UNKNOWN', gates=dict(interface_p95=float(np.quantile(dist, 0.95)) <= 0.06, partition=rowerr <= 1e-12, adjoint=all((g['pass_gate'] for g in gc)), equilibrium=all((v['equilibrium_relative_residual'] <= 1e-08 for v in contrasts.values())), volume=vol.sum() <= v0.sum() * (1 + 1e-08), intaglio=identity == 0, min_tetra_volume=float(np.min(vol / v0)) >= 0.5 - 1e-08, both_load_gain=all((v['optimized_ratio'] > 1.02 for v in contrasts.values())), export=bool(mesh.is_watertight and mesh.volume > 0), physical_manufacture=False), seconds=time.perf_counter() - start, external_referent=read(ROOT / 'PREREG_R5.json')['external_referent'], rigorous_design_box_enclosure='MISSING')
    dump(ROOT / 'raw/R5.json', out)
    np.savez_compressed(DATA / 'SURFACE_TIE_RAW.npz', vertices=new, faces=raw['faces'], interface_slaves=slaves, interface_triangles=tri, barycentric_weights=bary, projection_points=closest, projection_distances=dist, coefficients_mm=sol.x[:4])
    freeze(ROOT / 'FROZEN_SURFACE_TIE_PREDICTIONS.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R5.json'), contrasts=contrasts, ratio_operational_windows={k: [v['optimized_ratio'] * np.exp(-0.3), v['optimized_ratio'] * np.exp(0.3)] for (k, v) in contrasts.items()}, exports=inventory(dest), raw_sha256=sha(DATA / 'SURFACE_TIE_RAW.npz'), measurement='NOT_RUN', calibrated_force05_N=None, time_scale='HANDOVER'))
    state('R5_SURFACE_TIE_FROZEN', out['gates'], 'Freeze laboratory comparison and verify reproducible package; physical force/CAM/antagonist measurements remain missing')
    print(clean(dict(contrasts=contrasts, interface=out['interface_projection_mm'], gates=out['gates'], seconds=out['seconds'])))
if __name__ == '__main__':
    run()
