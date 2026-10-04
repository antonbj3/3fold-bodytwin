"""P1 tetrahedral elastic roof; contact-map influence is a conditional simulation."""
import time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from geometry import *

def constitutive(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[np.arange(3), np.arange(3)] += 2 * mu
    D[np.arange(3, 6), np.arange(3, 6)] = mu
    return D

def mesh(xy, z, faces, t=1.2, layers=3):
    n = len(xy)
    xyz = np.concatenate([np.c_[xy, z - t + k * t / layers] for k in range(layers + 1)])
    tet = []
    for k in range(layers):
        for (a, b, c) in np.sort(faces, axis=1):
            (a, b, c) = (a + k * n, b + k * n, c + k * n)
            tet.extend([[a, b, c, c + n], [a, b, b + n, c + n], [a, a + n, b + n, c + n]])
    return (xyz, np.asarray(tet, int))

def operators(xyz, tet, E, nu):
    X = np.concatenate([np.ones((len(tet), 4, 1)), xyz[tet]], axis=2)
    inv = np.linalg.inv(X)
    grad = inv[:, 1:, :].transpose(0, 2, 1)
    vol = np.abs(np.linalg.det(X)) / 6
    B = np.zeros((len(tet), 6, 12))
    for i in range(4):
        (gx, gy, gz) = grad[:, i, :].T
        B[:, 0, 3 * i] = gx
        B[:, 1, 3 * i + 1] = gy
        B[:, 2, 3 * i + 2] = gz
        B[:, 3, 3 * i] = gy
        B[:, 3, 3 * i + 1] = gx
        B[:, 4, 3 * i + 1] = gz
        B[:, 4, 3 * i + 2] = gy
        B[:, 5, 3 * i] = gz
        B[:, 5, 3 * i + 2] = gx
    Dm = constitutive(E, nu)
    ke = np.einsum('eai,ab,ebj,e->eij', B, Dm, B, vol, optimize=True)
    dofs = (3 * tet[:, :, None] + np.arange(3)).reshape(-1, 12)
    K = coo_matrix((ke.ravel(), (np.repeat(dofs, 12, axis=1).ravel(), np.tile(dofs, (1, 12)).ravel())), shape=(len(xyz) * 3,) * 2).tocsr()
    return (K, B, Dm, vol, dofs)

def solve_cases(xy, z, faces, uv, loads, contract):
    t = time.perf_counter()
    (xyz, tet) = mesh(xy, z, faces, contract['thickness_mm'])
    (K, B, Dm, vol, dofs) = operators(xyz, tet, contract['E_MPa'], contract['nu'])
    n = len(xy)
    support = np.flatnonzero(np.sum(uv ** 2, axis=1) >= 0.7 ** 2)
    fixed = (support[:, None] * 3 + np.arange(3)).ravel()
    free = np.setdiff1d(np.arange(len(xyz) * 3), fixed)
    factor = splu(K[free][:, free].tocsc())
    out = {}
    for (name, load) in loads.items():
        if load is None:
            out[name] = dict(status='UNKNOWN_NO_CONTACT')
            continue
        f = np.zeros(len(xyz) * 3)
        f[3 * (len(xyz) - n + np.arange(n)) + 2] = -load
        u = np.zeros_like(f)
        u[free] = factor.solve(f[free])
        reaction = K @ u - f
        strain = np.einsum('eai,ei->ea', B, u[dofs])
        sig = strain @ Dm.T
        S = np.zeros((len(tet), 3, 3))
        S[:, 0, 0] = sig[:, 0]
        S[:, 1, 1] = sig[:, 1]
        S[:, 2, 2] = sig[:, 2]
        S[:, 0, 1] = S[:, 1, 0] = sig[:, 3]
        S[:, 1, 2] = S[:, 2, 1] = sig[:, 4]
        S[:, 0, 2] = S[:, 2, 0] = sig[:, 5]
        principal = np.linalg.eigvalsh(S)
        tensile = np.maximum(principal[:, -1], 0)
        von = np.sqrt(0.5 * ((sig[:, 0] - sig[:, 1]) ** 2 + (sig[:, 1] - sig[:, 2]) ** 2 + (sig[:, 2] - sig[:, 0]) ** 2) + 3 * np.sum(sig[:, 3:] ** 2, axis=1))
        balance = np.sum(reaction[fixed].reshape(-1, 3), axis=0) + np.sum(f.reshape(-1, 3), axis=0)
        out[name] = dict(status='SIMULATED', maximum_tensile_MPa=float(tensile.max()), p95_tensile_MPa=float(np.quantile(tensile, 0.95)), maximum_von_Mises_MPa=float(von.max()), elastic_energy_Nmm=float(0.5 * u @ (K @ u)), equilibrium_residual_N=float(np.linalg.norm(reaction[free])), force_balance_N=float(np.linalg.norm(balance)), total_load_N=float(load.sum()), displacement_max_mm=float(np.linalg.norm(u.reshape(-1, 3), axis=1).max()))
    return (out, dict(wall_seconds=time.perf_counter() - t, nodes=len(xyz), tetrahedra=len(tet), dofs=len(xyz) * 3))

def normalize(mask, weights, F):
    w = mask * weights
    return w * F / w.sum() if w.sum() > 0 else None

def run():
    pr = json.loads((H / 'PREREG_R1.json').read_text())
    preds = json.loads((H / 'raw/PREDICTIONS_R1.json').read_text())
    rows = json.loads((H / 'raw/RESULTS_R1_ROWS.json').read_text())
    selected = []
    for fdi in pr['tooth_fdi']:
        rr = next((r for r in rows if r.get('status') == 'SCORED' and r['fdi'] == fdi), None)
        if rr:
            selected.append(rr)
    results = []
    for rr in selected:
        (case, fdi) = (rr['case'], rr['fdi'])
        state('FE_CROSSLOAD', 'R1 contact gates decided; stress is conditional', 'Contact-map roof FE ' + str(case) + '/' + str(fdi))
        pp = next((r for r in preds if r['case'] == case and r['fdi'] == fdi))
        z = np.load(pp['file'])
        ref = np.load(D / (str(case) + '_' + str(fdi) + '_reference.npz'))
        zz = ref['original_z'].copy()
        good = np.isfinite(zz)
        if not good.all():
            zz[~good] = cKDTree(z['xy'][good]).query(z['xy'][~good])[1]
            zz[~good] = ref['original_z'][good][zz[~good].astype(int)]
        loads = {arm: normalize(ref['0.1_' + arm + '_contact'], z['weights'], 100.0) for arm in ['original', 'practice', 'informed']}
        rad = np.linalg.norm(z['xy'] - z['center'], axis=1)
        bottom = np.full(len(rad), np.inf)
        inside = rad <= 3
        bottom[inside] = 3 - np.sqrt(9 - rad[inside] ** 2)
        offset = float(np.max(zz[inside] - bottom[inside]))
        sphere_gap = offset + bottom - zz
        sphere = np.isfinite(sphere_gap) & (sphere_gap <= 0.1)
        loads['sphere'] = normalize(sphere, z['weights'], 100.0)
        (ans, cost) = solve_cases(z['xy'], zz, z['faces'], z['uv'], loads, pr['fe_contract'])
        original = ans['original']
        for (name, a) in ans.items():
            if a['status'] == 'SIMULATED' and original['status'] == 'SIMULATED':
                a['relative_peak_tensile_error'] = abs(a['maximum_tensile_MPa'] / max(original['maximum_tensile_MPa'], 1e-12) - 1)
        np.savez_compressed(D / (str(case) + '_' + str(fdi) + '_fe_loads.npz'), **{k: v for (k, v) in loads.items() if v is not None}, sphere_mask=sphere)
        results.append(dict(case=case, fdi=fdi, type=rr['type'], crossloaded_same_original_roof=True, answers=ans, cost=cost))
        print('FE', case, fdi, cost['wall_seconds'], flush=True)
    dump(H / 'raw/FE_RESULTS.json', dict(claim_type='capability', rows=results, force_status='UNKNOWN physical pressure; normalized 100N scenario', external_referent=dict(kind='our_own_fixture', locator='code/fe.py same-original-roof cross-load model; NOT external stress reference', compared_quantity='Conditional P1 roof tensile stress from different geometric contact masks', refutes_us=False), physical_stress_accuracy='UNKNOWN', constitutive_contract=pr['fe_contract']))
if __name__ == '__main__':
    run()
