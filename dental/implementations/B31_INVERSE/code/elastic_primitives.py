"""Linear tetra elasticity, spectral Weibull functional, exact fixed-topology shape adjoint.

No physical crown calibration is implied. Fixed load and zero Dirichlet support
are explicit closures. Derivatives are local; a rigorous design-box enclosure is absent.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

def material(E=201400.0, nu=0.3):
    C = np.zeros((6, 6))
    C[:3, :3] = nu
    np.fill_diagonal(C[:3, :3], 1 - nu)
    C[3:, 3:] = np.eye(3) * (1 - 2 * nu) / 2
    return E / ((1 + nu) * (1 - 2 * nu)) * C

def strain_matrix(g):
    B = np.zeros((len(g), 6, 12))
    for i in range(4):
        (x, y, z) = g[:, i, :].T
        j = 3 * i
        B[:, 0, j] = x
        B[:, 1, j + 1] = y
        B[:, 2, j + 2] = z
        B[:, 3, j] = y
        B[:, 3, j + 1] = x
        B[:, 4, j + 1] = z
        B[:, 4, j + 2] = y
        B[:, 5, j] = z
        B[:, 5, j + 2] = x
    return B

def geometry(V, T):
    P = np.concatenate([np.ones((len(T), 4, 1)), V[T]], axis=2)
    det = np.linalg.det(P)
    vol = abs(det) / 6
    if np.min(vol) < 1e-12:
        raise ValueError('degenerate tetra')
    inv = np.linalg.inv(P)
    B = strain_matrix(inv[:, 1:, :].transpose(0, 2, 1))
    return (inv, B, vol, np.sign(det))

def stress_tensor(s):
    a = np.zeros((len(s), 3, 3))
    a[:, 0, 0] = s[:, 0]
    a[:, 1, 1] = s[:, 1]
    a[:, 2, 2] = s[:, 2]
    a[:, 0, 1] = a[:, 1, 0] = s[:, 3]
    a[:, 1, 2] = a[:, 2, 1] = s[:, 4]
    a[:, 0, 2] = a[:, 2, 0] = s[:, 5]
    return a

class Elastic:

    def __init__(self, V, T, fixed, load, m=2.89, modes=None):
        self.V = np.asarray(V, float)
        self.T = np.asarray(T, int)
        self.load = np.asarray(load, float).ravel()
        self.m = m
        self.C = material()
        self.calls = 0
        self.edof = (3 * self.T[:, :, None] + np.arange(3)).reshape(-1, 12)
        self.free = np.setdiff1d(np.arange(3 * len(V)), (3 * np.asarray(fixed)[:, None] + np.arange(3)).ravel())
        self.rows = np.repeat(self.edof, 12, axis=1).ravel()
        self.cols = np.tile(self.edof, (1, 12)).ravel()
        self.modes = np.zeros((len(V), 3, 0)) if modes is None else np.asarray(modes, float)
        self.sign = geometry(self.V, self.T)[3]
        self.last = None

    def evaluate(self, x, gradient=True):
        x = np.asarray(x, float)
        if self.last is not None and np.array_equal(x, self.last[0]) and (not gradient or self.last[1].get('gradient') is not None):
            return self.last[1]
        self.calls += 1
        V = self.V + np.einsum('nip,p->ni', self.modes, x)
        (inv, B, vol, sign) = geometry(V, self.T)
        if np.any(sign != self.sign):
            raise ValueError('tetra inverted')
        CB = np.einsum('ij,ejk->eik', self.C, B)
        ke = vol[:, None, None] * np.einsum('eji,ejk->eik', B, CB)
        K = coo_matrix((ke.ravel(), (self.rows, self.cols)), shape=(3 * len(V),) * 2).tocsc()
        factor = splu(K[self.free][:, self.free])
        u = np.zeros(3 * len(V))
        u[self.free] = factor.solve(self.load[self.free])
        ue = u[self.edof]
        strain = np.einsum('eij,ej->ei', B, ue)
        s = strain @ self.C
        (ev, v) = np.linalg.eigh(stress_tensor(s))
        sp = np.maximum(ev[:, -1], 0)
        vv = v[:, :, -1]
        H = float(vol @ sp ** self.m)
        if H <= 0:
            raise ValueError('zero tensile hazard')
        J = np.log(H) / self.m
        out = dict(J=float(J), hazard=H, volume_mm3=float(vol.sum()), max_tensile_MPa_per_N=float(sp.max()), equilibrium_relative_residual=float(np.linalg.norm((K @ u - self.load)[self.free]) / max(np.linalg.norm(self.load[self.free]), 1e-30)), gradient=None, vertices=V, stress_MPa_per_N=s, positive_principal=sp, element_volumes=vol, displacement_mm_per_N=u.reshape(-1, 3))
        if gradient:
            vs = np.c_[vv[:, 0] ** 2, vv[:, 1] ** 2, vv[:, 2] ** 2, 2 * vv[:, 0] * vv[:, 1], 2 * vv[:, 1] * vv[:, 2], 2 * vv[:, 0] * vv[:, 2]]
            dHs = vol[:, None] * self.m * np.power(sp, self.m - 1)[:, None] * vs
            rhs = np.zeros(len(u))
            np.add.at(rhs, self.edof.ravel(), np.einsum('ei,eij->ej', dHs, CB).ravel() / (self.m * H))
            lam = np.zeros(len(u))
            lam[self.free] = factor.solve(rhs[self.free])
            le = lam[self.edof]
            el = np.einsum('eij,ej->ei', B, le)
            grad = []
            for j in range(self.modes.shape[2]):
                dP = np.concatenate([np.zeros((len(self.T), 4, 1)), self.modes[self.T, :, j]], axis=2)
                dInv = -np.einsum('eij,ejk,ekl->eil', inv, dP, inv)
                dB = strain_matrix(dInv[:, 1:, :].transpose(0, 2, 1))
                dv = vol * np.einsum('eij,eji->e', inv, dP)
                deu = np.einsum('eij,ej->ei', dB, ue)
                delam = np.einsum('eij,ej->ei', dB, le)
                direct = (dv @ sp ** self.m + np.sum(dHs * (deu @ self.C))) / (self.m * H)
                dKU = dv @ np.einsum('ei,ei->e', el, s) + np.sum(vol[:, None] * (delam * s + el * (deu @ self.C)))
                grad.append(float(direct - dKU))
            out['gradient'] = np.array(grad)
        self.last = (x.copy(), out)
        return out

def roof_mesh(xy, faces, z, inner):
    n = len(xy)
    V = np.vstack([np.c_[xy, inner], np.c_[xy, z]])
    (a, b, c) = np.sort(np.asarray(faces), axis=1).T
    T = np.concatenate([np.c_[a, b, c, c + n], np.c_[a, b, b + n, c + n], np.c_[a, a + n, b + n, c + n]])
    return (V, T)

def roof_modes(xy):
    q = xy - xy.mean(0)
    q = q / np.maximum(np.max(abs(q), axis=0), 1e-08)
    (x, y) = q.T
    M = np.c_[np.ones(len(x)), x, y, x * y, x * x - y * y, np.exp(-4 * ((x - 0.4) ** 2 + (y - 0.2) ** 2)), np.exp(-4 * ((x + 0.4) ** 2 + (y + 0.2) ** 2))]
    return M

def roof_model(t, z, inner):
    xy = np.asarray(t['xy'])
    n = len(xy)
    (V, T) = roof_mesh(xy, t['faces'], z, inner)
    M = roof_modes(xy)
    modes = np.zeros((2 * n, 3, M.shape[1]))
    modes[n:, 2, :] = M
    g = np.asarray(t['ceiling']) - np.asarray(t['prior'])
    ok = np.isfinite(g)
    if ok.any():
        a = np.exp(-np.clip(g - np.nanmin(g), 0, 4) / 0.15)
        a[~ok] = 0
    else:
        a = np.exp(-4 * np.sum((xy / np.maximum(abs(xy).max(0), 1e-05)) ** 2, axis=1))
    w = np.asarray(t['weights']) * a
    w /= w.sum()
    load = np.zeros_like(V)
    load[n:, 2] = -w
    return (Elastic(V, T, np.arange(n), load, modes=modes), M)

def boundary(T):
    fs = np.concatenate([T[:, [0, 2, 1]], T[:, [0, 1, 3]], T[:, [0, 3, 2]], T[:, [1, 2, 3]]])
    key = np.sort(fs, axis=1)
    (_, i, c) = np.unique(key, axis=0, return_index=True, return_counts=True)
    return fs[i[c == 1]]

def gradient_check(model, x, step=1e-05):
    r = model.evaluate(x)
    ad = r['gradient']
    fd = []
    for i in range(len(x)):
        dx = np.zeros(len(x))
        dx[i] = step
        fd.append((model.evaluate(x + dx, False)['J'] - model.evaluate(x - dx, False)['J']) / (2 * step))
    fd = np.array(fd)
    err = float(np.linalg.norm(ad - fd) / max(np.linalg.norm(fd), 1e-12))
    return dict(adjoint=ad.tolist(), central_difference=fd.tolist(), relative_error=err, step_mm=step, pass_gate=err <= 0.0001, rigorous_design_box_enclosure='MISSING')
