from common import *
from native_port import model, deck_contract, q_matrix
import sys, trimesh
from scipy.sparse import coo_matrix, kron, eye, bmat, csc_matrix
from scipy.sparse.linalg import splu
sys.path.insert(0, str(X60 / 'code'))
from mechanics import geometry, material, stress_tensor

def interpolate10(V4, V10, T4, T10):
    n = len(V4)
    rr = list(range(n))
    cc = list(range(n))
    vv = list(np.ones(n))
    ends = {}
    pairs = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
    for (k, (i, j)) in enumerate(pairs):
        for (a, b, mid) in zip(T4[:, i], T4[:, j], T10[:, k + 4]):
            ends[int(mid)] = (int(a), int(b))
    for (mid, (a, b)) in sorted(ends.items()):
        rr.extend([mid, mid])
        cc.extend([a, b])
        vv.extend([0.5, 0.5])
    L = coo_matrix((vv, (rr, cc)), shape=(len(V10), n)).tocsr()
    assert np.max(abs(L @ V4 - V10)) < 1e-06
    return L

def assemble(V, T, E, nu):
    (inv, B, vol, sign) = geometry(V, T)
    C = np.array([material(e, v) for (e, v) in zip(E, nu)])
    CB = np.einsum('eij,ejk->eik', C, B)
    ke = vol[:, None, None] * np.einsum('eji,ejk->eik', B, CB)
    ed = (3 * T[:, :, None] + np.arange(3)).reshape(-1, 12)
    K = coo_matrix((ke.ravel(), (np.repeat(ed, 12, axis=1).ravel(), np.tile(ed, (1, 12)).ravel())), shape=(3 * len(V),) * 2).tocsc()
    return (K, B, C, vol, ed)

def reduce_solver(K, Q, fix):
    q = kron(Q, eye(3), format='csc')
    qr = (q.T @ K @ q).tocsc()
    active = np.arange(qr.shape[0])
    ix = np.setdiff1d(active, (3 * fix[:, None] + np.arange(3)).ravel())
    fac = splu(qr[ix][:, ix])

    def solve(f):
        fr = np.asarray(q.T @ f.ravel()).ravel()
        ur = np.zeros(qr.shape[0])
        ur[ix] = fac.solve(fr[ix])
        u = np.asarray(q @ ur).ravel()
        res = float(np.linalg.norm((qr @ ur - fr)[ix]) / np.linalg.norm(fr[ix]))
        return (u, res)
    return (solve, q)

def kkt_solver(K, slaves, masters, bary, fixed_original):
    n = K.shape[0] // 3
    C = coo_matrix((np.r_[np.ones(len(slaves)), -bary.ravel()], (np.r_[np.arange(len(slaves)), np.repeat(np.arange(len(slaves)), 3)], np.r_[slaves, masters.ravel()])), shape=(len(slaves), n)).tocsr()
    C = kron(C, eye(3), format='csc')
    free = np.setdiff1d(np.arange(K.shape[0]), (3 * fixed_original[:, None] + np.arange(3)).ravel())
    Cf = C[:, free]
    A = bmat([[K[free][:, free], Cf.T], [Cf, csc_matrix((C.shape[0],) * 2)]], format='csc')
    fac = splu(A)

    def solve(f):
        z = fac.solve(np.r_[f.ravel()[free], np.zeros(C.shape[0])])
        u = np.zeros(K.shape[0])
        u[free] = z[:len(free)]
        return (u, float(np.linalg.norm(A @ z - np.r_[f.ravel()[free], np.zeros(C.shape[0])]) / np.linalg.norm(f.ravel()[free])))
    return solve

def run():
    verify_inputs()
    t0 = time.perf_counter()
    (a, g, m) = model()
    deck = deck_contract()
    port = np.load(DATA / 'NATIVE_PORT.npz')
    slaves = port['slaves']
    allfaces = __import__('crown_design_fe').boundary_faces10(m.Td)
    allmesh = trimesh.Trimesh(m.Vd, allfaces[:, 2:5], process=False)
    P = []
    di = []
    fi = []
    for s in range(0, len(slaves), 128):
        (p, d, i) = trimesh.proximity.closest_point(allmesh, m.Vc[slaves[s:s + 128]])
        P.append(p)
        di.append(d)
        fi.append(i)
    P = np.concatenate(P)
    di = np.concatenate(di)
    fi = np.concatenate(fi)
    tri = allfaces[fi, 2:5]
    bary = trimesh.triangles.points_to_barycentric(m.Vd[tri], P)
    sdf = m.grid.sample(a['die'], m.Vc[slaves], cval=10.0)
    good = di <= 0.06
    corners = slaves < len(a['cV'])
    goodcorner = corners & good
    cv = a['cV']
    dv = a['dV']
    nc = len(cv)
    V = np.vstack([cv, dv])
    T = np.vstack([a['cT'], a['dT'] + nc])
    ne = len(a['cT'])
    Lc = interpolate10(cv, m.Vc, a['cT'], m.Tc)
    Ld = interpolate10(dv, m.Vd, a['dT'], m.Td)
    loads = [np.vstack([Lc.T @ f[:m.nc], Ld.T @ f[m.nc:]]) for f in deck['loads']]
    E = np.r_[np.full(ne, 201400.0), np.full(len(a['dT']), 18000.0)]
    nu = np.r_[np.full(ne, 0.3), np.full(len(a['dT']), 0.31)]
    (K, B, C, vol, ed) = assemble(V, T, E, nu)
    fixed = nc + m.fix_nodes[m.fix_nodes < len(dv)]
    sets = {'native_all': (slaves[corners], port['master_nodes'][corners, :3] + nc, port['barycentric'][corners]), 'all_master_all': (slaves[corners], tri[corners] + nc, bary[corners]), 'coverage_qualified': (slaves[goodcorner], tri[goodcorner] + nc, bary[goodcorner])}
    contrasts = {}
    solutions = {}
    solves = 0
    qq = {}
    for (name, (ss, mt, w)) in sets.items():
        (Q, map_) = q_matrix(len(V), ss, mt, w)
        (fun, q) = reduce_solver(K, Q, map_[fixed])
        qq[name] = (Q, map_)
        rows = []
        us = []
        for (angle, f) in zip([0, 30], loads):
            (u, res) = fun(f)
            solves += 1
            energy = float(0.5 * np.dot(f.ravel(), u))
            us.append(u)
            strain = np.einsum('eij,ej->ei', B, u[ed])
            stress = np.einsum('eij,ej->ei', C, strain)
            sp = np.maximum(np.linalg.eigvalsh(stress_tensor(stress))[:, -1], 0)
            haz = float(vol[:ne] @ sp[:ne] ** 2.89)
            rows.append({'angle_deg': angle, 'strain_energy_N_mm_per_unit_load': energy, 'residual': res, 'native_total_force_N': f.sum(0), 'conditional_volume_hazard': haz, 'hazard_scope': 'PHENOMENOLOGICAL independent flaws m2.89 from Prott force distribution; no strength inference', 'n_slave_corner_nodes': len(ss), 'unity_error': float(np.max(abs(np.asarray(Q.sum(1)).ravel() - 1)))})
        contrasts[name] = rows
        solutions[name] = us
    (ss, mt, w) = sets['coverage_qualified']
    full = kkt_solver(K, ss, mt, w, fixed)
    check = []
    for (i, f) in enumerate(loads):
        (uf, res) = full(f)
        solves += 1
        ur = solutions['coverage_qualified'][i]
        err = float(np.linalg.norm(uf - ur) / np.linalg.norm(ur))
        check.append({'angle_deg': [0, 30][i], 'displacement_relative_L2_error': err, 'KKT_residual': res, 'energy_difference_N_mm': float(0.5 * np.dot(f.ravel(), uf - ur)), 'pass': bool(err <= 1e-07 and res <= 1e-08)})
    faxis = loads[0]
    loaded = np.flatnonzero(np.linalg.norm(faxis[:nc], axis=1) > 0)
    xy = cv[loaded, :2]
    pair = np.unravel_index(np.argmax(np.sum((xy[:, None, :] - xy[None, :, :]) ** 2, axis=2)), (len(loaded),) * 2)
    nodes = loaded[list(pair)]
    fw = []
    energy = []
    (Q, map_) = qq['native_all']
    (fun, q) = reduce_solver(K, Q, map_[fixed])
    for node in nodes:
        ff = np.zeros_like(V)
        ff[node, 2] = -1.0
        (uw, _) = fun(ff)
        solves += 1
        fw.append(ff)
        energy.append(float(0.5 * np.dot(ff.ravel(), uw)))
    area = read(R / 'raw/R1_NATIVE_PORT.json')['native_support_area_mm2']
    summary = [np.r_[ff.sum(0), area] for ff in fw]
    ident = float(np.max(abs(summary[0] - summary[1])))
    diff = abs(energy[1] - energy[0]) / max(energy)
    suff = {'summary_a': summary[0], 'summary_b': summary[1], 'identity_error': ident, 'support_area_mm2': area, 'node_ids': nodes, 'coordinates_mm': cv[nodes], 'strain_energy_N_mm': energy, 'relative_downstream_difference': diff, 'summary_sufficient': False if ident == 0 and diff > 0.01 else 'UNREFUTED', 'minimum_extension': 'Keep local load vector and node/facet-owner transfer Q; total load plus area loses the necessary association. No proof of globally minimal dimension.', 'origin': 'our_own_fixture on fixed anatomical FE model; not physical validation'}
    out = {'round': 'R2_SUPPORT_COVERAGE', 'claim_type': 'capability', 'unrestricted_projection_mm': {'max': float(di.max()), 'p95': float(np.quantile(di, 0.95)), 'median': float(np.median(di)), 'resolution': 'PER_POINT'}, 'die_SDF_samples_mm': {'max_abs': float(np.max(abs(sdf))), 'p95_abs': float(np.quantile(abs(sdf), 0.95)), 'worst_native_node_signed_value': float(sdf[np.argmax(port['distances_mm'])]), 'scope': 'Trilinear SDF values, not exact Euclidean distance; no physical metrology'}, 'native_coverage': {'required_nodes': len(slaves), 'qualified_nodes': int(good.sum()), 'missing_nodes': int((~good).sum()), 'missing_fraction': float((~good).mean()), 'required_corner_nodes': int(corners.sum()), 'qualified_corner_nodes': int(goodcorner.sum())}, 'original_model_preserved': {'geometry': True, 'material': True, 'load_port': True, 'support_sets_for_alternative': False}, 'contrasts_P1_conditional': contrasts, 'same_Q_full_FE_control': check, 'sufficiency': suff, 'FE_solve_calls': solves, 'gate': {'unrestricted_all_projection': bool(di.max() <= 0.06), 'qualified_projection': bool(di[good].max() <= 0.06), 'full_native_coverage': bool(good.all()), 'B04_native_set_parity': False, 'same_Q_FE_control': all((x['pass'] for x in check)), 'summary_counterexample': bool(ident == 0 and diff > 0.01), 'optimization_allowed': False}, 'rigorous_enclosure': 'MISSING; no linear sensitivity or physical force claim', 'cost': cost(t0)}
    np.savez_compressed(DATA / 'SUPPORT_COVERAGE.npz', unrestricted_projection_points=P, unrestricted_distances_mm=di, qualified_mask=good, sdf_die_at_native_slaves=sdf, master_face=fi, master_nodes=tri, barycentric=bary, P1_loads=np.array(loads), P1_vertices=V, P1_tetra=T, P1_native_displacements=np.array(solutions['native_all']), P1_qualified_displacements=np.array(solutions['coverage_qualified']), sufficiency_loads=np.array(fw), sufficiency_node_ids=nodes)
    dump(R / 'raw/R2_SUPPORT_COVERAGE.json', out)
    state('R2_SUPPORT_COVERAGE_DECIDED', out['gate'], 'Native solver parity, external source contrast and changed interface construction (not shape optimization)')
    print(json.dumps(clean(out), indent=2))
if __name__ == '__main__':
    run()
