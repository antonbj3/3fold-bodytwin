from common import *
from mechanics import *
from coupled import CoupledElastic
from whole import three_mf
from scipy.spatial import cKDTree
from scipy.optimize import minimize, LinearConstraint
import sys, time, trimesh

def prereg():
    freeze(ROOT / 'PREREG_R4.json', dict(round='R4', claim_type='capability', capability='Find a whole-crown shape improving the worst of axial and 30deg offaxis conditional fracture response on a deformable 18GPa die', obstacle='R3 rigid support and one fixed load may manufacture an apparent shape gain', changed_operation='Two-domain elastic stiffness, crown-to-die nearest-node kinematic tie, two loads and minimax epigraph shape optimization at fixed crown volume', consumer='K39 prospective matched support/load contrast', external_referent=dict(kind='published_code', locator=[str(X1 / 'vendor/crown_design_fe.py'), str(X1 / 'inputs/geometry/D1_model.npz')], compared_quantity='Source crown/die geometry and declared support/material/load scenarios, not a physical optimized-crown force measurement', refutes_us=True), metrics=dict(adjoint_relative_error_max=0.0001, equilibrium_relative_residual_max=1e-08, volume_relative_increase_max=1e-08, intaglio_identity_error=0.0, min_tetra_volume_ratio=0.5, conditional_gain_both_loads_min=1.02, prospective_ratio_log_window=0.3, valid_fault_rejection=True), mechanical_closures=['Nearest-node tie distances retained, not native surface tie', 'C3D4 instead of C3D10; discretization error not enclosed', '18GPa nu=.30 die from X1; not a new measurement', 'Geometry-derived fixed force patches, load does not update with contact', 'm=2.89 borrowed from G0.8 force law; m2/m4 sensitivity is conditional'], optimizer_maxiter=20, full_cost=dict(preparation='Existing whole mesh and raw fixed-node/mode arrays', fit='Two sparse equilibrium/adjoint pairs per iteration', validation='FD, both load contrasts, original controls with valid-coordinate injection', questions='FEM queries recorded', fallback='Any failed gate retained; absolute quantile stays UNKNOWN'), strongest_equally_informed_control='R3 frozen design evaluated under both new support/load scenarios, original shape at exactly same volume', falsifiers=['either load ratio <=1.02 fails useful gain', 'gradient mismatch >1e-4', 'changed intaglio or volume excess', 'valid antagonist/connector/film/wall fault survives unchanged predicate']))

def run():
    prereg()
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
    fixedc = raw['fixed_nodes']
    grid = read(X1 / 'inputs/geometry/D1_grid.json')
    (dist, near) = cKDTree(dV).query(cV[fixedc])
    node_map = np.arange(len(V))
    node_map[fixedc] = nc + near
    fixed = nc + np.flatnonzero(dV[:, 2] < grid['z_m'] - 3 + 0.03)
    E = np.r_[np.full(ne, 201400.0), np.full(len(dT), 18000.0)]
    saved = np.load(X1 / 'inputs/fe/L005_lo_M1_k6_thin_3Y.npz')
    bn = np.unique(raw['faces'])
    top = bn[(cV[bn, 2] > grid['z_m'] + 1.0) & ~np.isin(bn, fixedc)]
    models = {}
    baselines = {}
    frozen = {}
    r3coef = np.asarray(read(ROOT / 'raw/R3.json')['coefficients_mm'])
    for angle in [0, 30]:
        center = saved['lc__axial' if angle == 0 else 'lc__offaxis30'][0, :3]
        w = np.exp(-np.sum((cV[top, :2] - center[:2]) ** 2, axis=1) / (2 * 0.6 ** 2))
        w /= w.sum()
        load = np.zeros_like(V)
        direction = np.array([np.sin(np.radians(angle)), 0, -np.cos(np.radians(angle))])
        load[top] = w[:, None] * direction
        model = CoupledElastic(V, T, fixed, load, modes, ne, node_map, E)
        models[str(angle)] = model
        baselines[str(angle)] = model.evaluate(np.zeros(4))['J']
        frozen[str(angle)] = np.exp(baselines[str(angle)] - model.evaluate(r3coef)['J'])
    (inv, B, v0, sign) = geometry(cV, cT)
    dv = []
    for j in range(4):
        dp = np.concatenate([np.zeros((ne, 4, 1)), raw['modes'][cT, :, j]], axis=2)
        dv.append(v0 * np.einsum('eij,eji->e', inv, dp))
    dv = np.array(dv).T
    M = raw['modes'][:, 2, :]
    volgrad = dv.sum(0)
    lin = np.vstack([np.c_[M, np.zeros(len(M))], np.c_[dv, np.zeros(ne)], np.r_[volgrad, 0][None, :]])
    lo = np.r_[np.full(len(M), -0.2), -0.5 * v0, -np.inf]
    hi = np.r_[np.full(len(M), 0.2), np.full(ne, np.inf), 0]

    def risk(y):
        return np.array([y[4] - (m.evaluate(y[:4])['J'] - baselines[k]) for (k, m) in models.items()])

    def riskjac(y):
        return np.array([np.r_[-m.evaluate(y[:4])['gradient'], 1] for m in models.values()])
    sol = minimize(lambda y: (y[4], np.r_[np.zeros(4), 1.0]), np.zeros(5), jac=True, method='SLSQP', constraints=[LinearConstraint(lin, lo, hi), dict(type='ineq', fun=risk, jac=riskjac)], bounds=[(-0.2, 0.2)] * 4 + [(None, None)], options=dict(maxiter=20, ftol=1e-08))
    contrasts = {}
    gc = []
    for (k, m) in models.items():
        r = m.evaluate(sol.x[:4])
        contrasts[k] = dict(R3_frozen_ratio_on_compliant_support=float(frozen[k]), robust_optimized_ratio=float(np.exp(baselines[k] - r['J'])), equilibrium_relative_residual=r['equilibrium_relative_residual'])
        gc.append(dict(angle_deg=int(k), **gradient_check(m, sol.x[:4])))
    final = models['0'].evaluate(sol.x[:4])
    new = final['vertices'][:nc]
    dest = ROOT / 'exports/whole_STS_robust'
    dest.mkdir(parents=True, exist_ok=True)
    mesh = trimesh.Trimesh(new, raw['faces'], process=True)
    if mesh.volume < 0:
        mesh.invert()
    mesh.export(dest / 'optimized.stl')
    three_mf(mesh, dest / 'optimized.3mf')
    vol = geometry(new, cT)[2]
    identity = float(np.max(abs(new[fixedc] - cV[fixedc])))
    r1 = read(ROOT / 'raw/R1.json')
    out = dict(round='R4', claim_type='capability', coefficients_mm=sol.x[:4], optimizer_success=bool(sol.success), optimizer_message=str(sol.message), iterations=int(sol.nit), FEM_queries={k: m.calls for (k, m) in models.items()}, contrasts=contrasts, gradient_checks=gc, interface_tie_distance_mm=dict(median=float(np.median(dist)), p95=float(np.quantile(dist, 0.95)), max=float(dist.max()), resolution='PER_POINT', validity='Numerical mapping mismatch, not measured cement thickness'), intaglio_identity_error_mm=identity, volume_relative_change=float(vol.sum() / v0.sum() - 1), min_tetra_volume_ratio=float(np.min(vol / v0)), calibrated_force05_N=None, physical_Pareto='UNKNOWN', gates=dict(adjoint=all((g['pass_gate'] for g in gc)), equilibrium=all((v['equilibrium_relative_residual'] <= 1e-08 for v in contrasts.values())), volume=vol.sum() <= v0.sum() * (1 + 1e-08), intaglio=identity == 0, min_tetra_volume=float(np.min(vol / v0)) >= 0.5 - 1e-08, both_load_gain=all((v['robust_optimized_ratio'] > 1.02 for v in contrasts.values())), export=bool(mesh.is_watertight and mesh.volume > 0), physical_manufacture=False), external_referent=read(ROOT / 'PREREG_R4.json')['external_referent'], seconds=time.perf_counter() - start, rigorous_design_box_enclosure='MISSING')
    dump(ROOT / 'raw/R4.json', out)
    np.savez_compressed(DATA / 'ROBUST_RAW.npz', vertices=new, faces=raw['faces'], coefficients_mm=sol.x[:4], load0=models['0'].load, load30=models['30'].load, node_map=node_map, tie_distances=dist)
    freeze(ROOT / 'FROZEN_ROBUST_PREDICTIONS.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R4.json'), contrasts=contrasts, ratio_operational_windows={k: [v['robust_optimized_ratio'] * np.exp(-0.3), v['robust_optimized_ratio'] * np.exp(0.3)] for (k, v) in contrasts.items()}, exports=inventory(dest), raw_sha256=sha(DATA / 'ROBUST_RAW.npz'), measurement='NOT_RUN', physical_quantile_N=None, time_scale='HANDOVER', resolution='PER_TOOTH'))
    state('R4_ROBUST_WHOLE_FROZEN', out['gates'], 'Valid-coordinate rejecting faults, lab comparison, reproducibility and handoff')
    print(clean(dict(contrasts=contrasts, gates=out['gates'], seconds=out['seconds'])))
    return out
if __name__ == '__main__':
    run()
