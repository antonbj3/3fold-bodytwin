"""Two material elastic solids with explicit nearest-node tied interface.

The interface mapping is a numerical closure, with distances retained. It is
not the native C3D10 surface interpolation tie and does not certify cement.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from mechanics import geometry, strain_matrix, stress_tensor, material

class CoupledElastic:

    def __init__(self, V, T, fixed, load, modes, nc_elements, node_map, E, m=2.89):
        self.V = V
        self.T = T
        self.modes = modes
        self.m = m
        self.nc_elements = nc_elements
        self.calls = 0
        self.last = None
        self.C = np.array([material(e) for e in E])
        self.edof = (3 * node_map[T[:, :, None]] + np.arange(3)).reshape(-1, 12)
        n = len(V) * 3
        active = np.unique(self.edof)
        fix = (3 * np.asarray(fixed)[:, None] + np.arange(3)).ravel()
        self.free = np.setdiff1d(active, fix)
        self.load = np.zeros(n)
        np.add.at(self.load, (3 * node_map[:, None] + np.arange(3)).ravel(), load.ravel())
        self.rows = np.repeat(self.edof, 12, axis=1).ravel()
        self.cols = np.tile(self.edof, (1, 12)).ravel()
        self.sign = geometry(V, T)[3]

    def evaluate(self, x, gradient=True):
        x = np.asarray(x, float)
        if self.last is not None and np.array_equal(x, self.last[0]) and (not gradient or self.last[1]['gradient'] is not None):
            return self.last[1]
        self.calls += 1
        V = self.V + np.einsum('nip,p->ni', self.modes, x)
        (inv, B, vol, sign) = geometry(V, self.T)
        if np.any(sign != self.sign):
            raise ValueError('inverted coupled tetra')
        CB = np.einsum('eij,ejk->eik', self.C, B)
        ke = vol[:, None, None] * np.einsum('eji,ejk->eik', B, CB)
        K = coo_matrix((ke.ravel(), (self.rows, self.cols)), shape=(3 * len(V),) * 2).tocsc()
        fac = splu(K[self.free][:, self.free])
        u = np.zeros(3 * len(V))
        u[self.free] = fac.solve(self.load[self.free])
        ue = u[self.edof]
        eu = np.einsum('eij,ej->ei', B, ue)
        s = np.einsum('eij,ej->ei', self.C, eu)
        (vals, vec) = np.linalg.eigh(stress_tensor(s))
        sp = np.maximum(vals[:, -1], 0)
        v = vec[:, :, -1]
        weight = vol.copy()
        weight[self.nc_elements:] = 0
        H = float(weight @ sp ** self.m)
        J = float(np.log(H) / self.m)
        out = dict(J=J, hazard=H, volume_mm3=float(vol[:self.nc_elements].sum()), gradient=None, vertices=V, positive_principal=sp[:self.nc_elements], stress_MPa_per_N=s[:self.nc_elements], equilibrium_relative_residual=float(np.linalg.norm((K @ u - self.load)[self.free]) / np.linalg.norm(self.load[self.free])))
        if gradient:
            spec = np.c_[v[:, 0] ** 2, v[:, 1] ** 2, v[:, 2] ** 2, 2 * v[:, 0] * v[:, 1], 2 * v[:, 1] * v[:, 2], 2 * v[:, 0] * v[:, 2]]
            dHs = weight[:, None] * self.m * (sp ** (self.m - 1))[:, None] * spec
            rhs = np.zeros(len(u))
            np.add.at(rhs, self.edof.ravel(), (np.einsum('ei,eij->ej', dHs, CB) / (self.m * H)).ravel())
            lam = np.zeros(len(u))
            lam[self.free] = fac.solve(rhs[self.free])
            el = np.einsum('eij,ej->ei', B, lam[self.edof])
            gr = []
            for j in range(self.modes.shape[2]):
                dp = np.concatenate([np.zeros((len(self.T), 4, 1)), self.modes[self.T, :, j]], axis=2)
                di = -np.einsum('eij,ejk,ekl->eil', inv, dp, inv)
                db = strain_matrix(di[:, 1:, :].transpose(0, 2, 1))
                dv = vol * np.einsum('eij,eji->e', inv, dp)
                dw = dv.copy()
                dw[self.nc_elements:] = 0
                deu = np.einsum('eij,ej->ei', db, ue)
                delam = np.einsum('eij,ej->ei', db, lam[self.edof])
                ds = np.einsum('eij,ej->ei', self.C, deu)
                direct = (dw @ sp ** self.m + np.sum(dHs * ds)) / (self.m * H)
                dKU = dv @ np.einsum('ei,ei->e', el, s) + np.sum(vol[:, None] * (delam * s + el * ds))
                gr.append(float(direct - dKU))
            out['gradient'] = np.asarray(gr)
        self.last = (x.copy(), out)
        return out
