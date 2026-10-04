from dental_release.paths import expand as _release_expand
from common import *
from elastic_primitives import geometry, strain_matrix, stress_tensor, material, boundary
from scipy.sparse import coo_matrix, eye, kron
from scipy.sparse.linalg import splu
from scipy.optimize import minimize, LinearConstraint
import time

class MultiLoad:

    def __init__(self, V, T, modes, nc, E, Q, fixed, loads, name):
        self.V = V
        self.T = T
        self.modes = modes
        self.nc = nc
        self.name = name
        self.C = np.array([material(e) for e in E])
        self.Q = kron(Q, eye(3), format='csc')
        self.edof = (3 * T[:, :, None] + np.arange(3)).reshape(-1, 12)
        self.rows = np.repeat(self.edof, 12, axis=1).ravel()
        self.cols = np.tile(self.edof, (1, 12)).ravel()
        self.free = np.setdiff1d(np.arange(self.Q.shape[1]), (3 * np.array(fixed)[:, None] + np.arange(3)).ravel())
        self.F = np.stack(loads, axis=-1).reshape(-1, 2)
        self.fr = np.asarray(self.Q.T @ self.F)
        self.sign = geometry(V, T)[3]
        self.calls = 0
        self.last = None

    def solve(self, x):
        if self.last is not None and np.array_equal(x, self.last[0]):
            return self.last[1]
        self.calls += 1
        V = self.V + np.einsum('nip,p->ni', self.modes, x)
        (inv, B, vol, sgn) = geometry(V, self.T)
        if np.any(sgn != self.sign):
            raise ValueError('inverted tetra')
        CB = np.einsum('eij,ejk->eik', self.C, B)
        ke = vol[:, None, None] * np.einsum('eji,ejk->eik', B, CB)
        K = coo_matrix((ke.ravel(), (self.rows, self.cols)), shape=(3 * len(V),) * 2).tocsc()
        Kr = (self.Q.T @ K @ self.Q).tocsc()
        ix = self.free
        fac = splu(Kr[ix][:, ix])
        ur = np.zeros((Kr.shape[0], 2))
        ur[ix] = fac.solve(self.fr[ix])
        u = np.asarray(self.Q @ ur)
        ue = u[self.edof]
        strain = np.einsum('eij,ejk->eik', B, ue)
        stress = np.einsum('eij,ejk->eik', self.C, strain)
        residual = float(np.linalg.norm((Kr @ ur - self.fr)[ix]) / np.linalg.norm(self.fr[ix]))
        out = dict(V=V, inv=inv, B=B, CB=CB, vol=vol, fac=fac, u=u, ue=ue, stress=stress, Kr=Kr, residual=residual)
        self.last = (np.array(x).copy(), out)
        return out

    def values_grad(self, x, angles, ms):
        a = self.solve(x)
        n = len(self.T)
        pairs = [(angle, m) for angle in angles for m in ms]
        weights = a['vol'].copy()
        weights[self.nc:] = 0
        J = []
        rhs = []
        local = []
        for (angle, m) in pairs:
            t = np.deg2rad(angle)
            c = np.array([np.cos(t), np.sin(t)])
            s = a['stress'] @ c
            ue = a['ue'] @ c
            (vals, vectors) = np.linalg.eigh(stress_tensor(s))
            sp = np.maximum(vals[:, -1], 0)
            v = vectors[:, :, -1]
            H = float(weights @ sp ** m)
            J.append(np.log(H) / m)
            spec = np.c_[v[:, 0] ** 2, v[:, 1] ** 2, v[:, 2] ** 2, 2 * v[:, 0] * v[:, 1], 2 * v[:, 1] * v[:, 2], 2 * v[:, 0] * v[:, 2]]
            d = weights[:, None] * (sp ** (m - 1))[:, None] * spec / H
            r = np.zeros(3 * len(self.V))
            np.add.at(r, self.edof.ravel(), np.einsum('ei,eij->ej', d, a['CB']).ravel())
            rhs.append(np.asarray(self.Q.T @ r).ravel())
            local.append((s, ue, sp, m, H, d))
        rr = np.stack(rhs, axis=-1)
        lr = np.zeros_like(rr)
        lr[self.free] = a['fac'].solve(rr[self.free])
        lam = np.asarray(self.Q @ lr)
        grad = np.empty((len(pairs), self.modes.shape[2]))
        for j in range(self.modes.shape[2]):
            dp = np.concatenate([np.zeros((n, 4, 1)), self.modes[self.T, :, j]], axis=2)
            di = -np.einsum('eij,ejk,ekl->eil', a['inv'], dp, a['inv'])
            db = strain_matrix(di[:, 1:, :].transpose(0, 2, 1))
            dv = a['vol'] * np.einsum('eij,eji->e', a['inv'], dp)
            dw = dv.copy()
            dw[self.nc:] = 0
            for (k, (s, ue, sp, m, H, d)) in enumerate(local):
                le = lam[self.edof, k]
                el = np.einsum('eij,ej->ei', a['B'], le)
                delu = np.einsum('eij,ej->ei', db, ue)
                dell = np.einsum('eij,ej->ei', db, le)
                ds = np.einsum('eij,ej->ei', self.C, delu)
                direct = dw @ sp ** m / (m * H) + np.sum(d * ds)
                dku = dv @ np.einsum('ei,ei->e', el, s) + np.sum(a['vol'][:, None] * (dell * s + el * ds))
                grad[k, j] = direct - dku
        return (np.array(J), grad)

    def fields(self, x):
        a = self.solve(x)
        return dict(stress_basis=a['stress'][:self.nc].copy(), vol=a['vol'][:self.nc].copy(), residual=a['residual'])

def build():
    z = np.load(X1 / 'inputs/geometry/D1_model.npz')
    raw = np.load(_release_expand('@DENTAL_WORK_ROOT@/X60/WHOLE_RAW.npz'))
    tie = np.load(_release_expand('@DENTAL_WORK_ROOT@/X60/SURFACE_TIE_RAW.npz'))
    cV = z['cV']
    cT = z['cT']
    dV = z['dV']
    dT = z['dT']
    nc = len(cV)
    ne = len(cT)
    modes0 = raw['modes']
    slaves = tie['interface_slaves']
    tri = tie['interface_triangles']
    bary = tie['barycentric_weights']
    V = np.vstack([cV, dV])
    T = np.vstack([cT, dT + nc])
    modes = np.concatenate([modes0, np.zeros((len(dV), 3, 4))])
    active = np.setdiff1d(np.arange(len(V)), slaves)
    mapping = np.full(len(V), -1, int)
    mapping[active] = np.arange(len(active))
    qr = np.r_[active, np.repeat(slaves, 3)]
    qc = np.r_[mapping[active], mapping[nc + tri].ravel()]
    qv = np.r_[np.ones(len(active)), bary.ravel()]
    Q = coo_matrix((qv, (qr, qc)), shape=(len(V), len(active))).tocsr()
    zm = read(X1 / 'inputs/geometry/D1_grid.json')['z_m']
    fixed = mapping[nc + np.flatnonzero(dV[:, 2] < zm - 3 + 0.03)]
    fz = raw['load']
    fx = np.zeros_like(fz)
    fx[:, 0] = -fz[:, 2]
    loads = [np.vstack([f, np.zeros_like(dV)]) for f in [fz, fx]]
    models = [MultiLoad(V, T, modes, ne, np.r_[np.full(ne, 201400.0), np.full(len(dT), E)], Q, fixed, loads, str(int(E))) for E in [8600.0, 18000.0]]
    models.append(MultiLoad(cV, cT, modes0, ne, np.full(ne, 201400.0), eye(nc, format='csr'), slaves, [fz, fx], 'rigid'))
    return (models, cV, cT, modes0, slaves, raw['faces'])

def run():
    start = time.perf_counter()
    p = read(ROOT / 'PREREG_B.json')
    (models, V, T, modes, slaves, faces) = build()
    angles = p['robustness']['angle_fit_grid']
    ms = p['robustness']['m_fit_grid']
    zero = np.zeros(4)
    base = {}
    fields = {}
    baseline = {}
    for model in models:
        print('base', model.name, flush=True)
        base[model.name] = model.values_grad(zero, angles, ms)[0]
        fields[model.name] = model.fields(zero)
    prior = {}
    for name in ['R3', 'R4', 'R5']:
        x = np.array(read(X60 / f'raw/{name}.json')['coefficients_mm'])
        prior[name] = {}
        for model in models:
            f = model.fields(x)
            prior[name][model.name] = f
    (inv, _, v0, _) = geometry(V, T)
    dv = []
    for j in range(4):
        dp = np.concatenate([np.zeros((len(T), 4, 1)), modes[T, :, j]], axis=2)
        dv.append(v0 * np.einsum('eij,eji->e', inv, dp))
    dv = np.array(dv).T
    M = modes[:, 2, :]
    lin = np.vstack([np.c_[M, np.zeros(len(M))], np.c_[dv, np.zeros(len(T))], np.r_[dv.sum(0), 0][None, :]])
    constraints = [LinearConstraint(lin, np.r_[np.full(len(M), -0.2), -0.5 * v0, -np.inf], np.r_[np.full(len(M), 0.2), np.full(len(T), np.inf), 0])]
    cache = {}
    history = []

    def evaluate(y):
        key = tuple(y[:4])
        if key not in cache:
            values = []
            grads = []
            for model in models:
                (v, g) = model.values_grad(y[:4], angles, ms)
                values.extend(v - base[model.name])
                grads.extend(g)
            cache.clear()
            cache[key] = (np.array(values), np.array(grads))
        return cache[key]

    def cfun(y):
        return y[4] - evaluate(y)[0]

    def cjac(y):
        return np.c_[-evaluate(y)[1], np.ones(len(models) * len(angles) * len(ms))]

    def cb(y):
        h = dict(iteration=len(history) + 1, coefficients_mm=y[:4].tolist(), worst_force_ratio=float(np.exp(-np.max(evaluate(y)[0]))))
        history.append(h)
        save(ROOT / 'raw/B_OPTIMIZER_CHECKPOINT.json', history)
        print(h, flush=True)
    sol = minimize(lambda y: (y[4], np.r_[np.zeros(4), 1.0]), np.zeros(5), jac=True, method='SLSQP', bounds=[(-0.2, 0.2)] * 4 + [(None, None)], constraints=constraints + [dict(type='ineq', fun=cfun, jac=cjac)], callback=cb, options=dict(maxiter=p['metrics']['optimizer_maxiter'], ftol=1e-08))
    final = {model.name: model.fields(sol.x[:4]) for model in models}
    model = models[0]
    (val, grad) = model.values_grad(sol.x[:4], [15], [6.0])
    fd = []
    h = 1e-05
    for j in range(4):
        step = np.zeros(4)
        step[j] = h
        vv = []
        for xx in [sol.x[:4] + step, sol.x[:4] - step]:
            f = model.fields(xx)
            s = f['stress_basis'] @ np.array([np.cos(np.pi / 12), np.sin(np.pi / 12)])
            sp = np.maximum(np.linalg.eigvalsh(stress_tensor(s))[:, -1], 0)
            vv.append(np.log(f['vol'] @ sp ** 6) / 6)
        fd.append((vv[0] - vv[1]) / (2 * h))
    err = float(np.linalg.norm(grad[0] - fd) / max(np.linalg.norm(fd), 1e-12))
    new = V + np.einsum('nip,p->ni', modes, sol.x[:4])
    vol = geometry(new, T)[2]
    arrays = dict(vertices_reference=V, vertices_optimized=new, tetra=T, faces=faces, modes=modes, slaves=slaves)
    for (label, ff) in [('reference', fields), ('optimized', final)] + list(prior.items()):
        for (support, f) in ff.items():
            for k in ['stress_basis', 'vol']:
                arrays[label + '__' + support + '__' + k] = f[k]
    np.savez_compressed(DATA / 'B_FIELDS.npz', **arrays)
    out = dict(claim_type='capability', round='B', coefficients_mm=sol.x[:4], optimizer_success=bool(sol.success), optimizer_message=sol.message, history=history, iterations=int(sol.nit), factorizations={m.name: m.calls for m in models}, gradient_check=dict(adjoint=grad[0], finite_difference=fd, relative_error=err, pass_gate=err <= 0.0001, rigorous_design_box_enclosure='MISSING'), residuals={k: v['residual'] for (k, v) in final.items()}, intaglio_coordinate_identity_error_mm=float(np.max(abs(new[slaves] - V[slaves]))), volume_relative_change=float(vol.sum() / v0.sum() - 1), min_tetra_volume_ratio=float(np.min(vol / v0)), max_vertex_shift_mm=float(np.max(np.linalg.norm(new - V, axis=1))), seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, field_path=DATA / 'B_FIELDS.npz', field_sha256=sha(DATA / 'B_FIELDS.npz'), physical_certification='UNKNOWN', external_referent=p['external_referent'])
    save(ROOT / 'raw/B_FIT.json', out)
    freeze(ROOT / 'FROZEN_PREDICTIONS_B.json', dict(claim_type='capability', coefficients_mm=sol.x[:4], prereg_sha256=sha(ROOT / 'PREREG_B.json'), field_sha256=out['field_sha256'], lab_measurement='NOT_RUN', prediction='Relative quantile ordering only under common volume flaw law and declared FE support; evaluator will enclose angle/m before any lab observation', absolute_force_quantile_N=None))
    state('B_FIELDS_FROZEN', dict(optimizer_success=bool(sol.success), adjoint_error=err), 'Certify continuous parameter box and evaluate inherited forms under same scenarios')
    print(json.dumps(clean(out)), flush=True)
if __name__ == '__main__':
    run()
