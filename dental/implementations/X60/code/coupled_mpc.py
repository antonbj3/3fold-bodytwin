"""Surface barycentric multipoint tie, preserving exact virtual work through Q^T K Q."""
from coupled import *
from scipy.sparse import kron, eye

class SurfaceTieElastic(CoupledElastic):

    def __init__(self, V, T, fixed, load, modes, nc_elements, Q, E, m=2.89):
        super().__init__(V, T, [], load, modes, nc_elements, np.arange(len(V)), E, m)
        self.Q = kron(Q, eye(3), format='csc')
        self.reduced_free = np.setdiff1d(np.arange(self.Q.shape[1]), (3 * np.asarray(fixed)[:, None] + np.arange(3)).ravel())

    def evaluate(self, x, gradient=True):
        x = np.asarray(x, float)
        if self.last is not None and np.array_equal(x, self.last[0]) and (not gradient or self.last[1]['gradient'] is not None):
            return self.last[1]
        self.calls += 1
        V = self.V + np.einsum('nip,p->ni', self.modes, x)
        (inv, B, vol, sign) = geometry(V, self.T)
        if np.any(sign != self.sign):
            raise ValueError('inverted surface-tied tetra')
        CB = np.einsum('eij,ejk->eik', self.C, B)
        ke = vol[:, None, None] * np.einsum('eji,ejk->eik', B, CB)
        K = coo_matrix((ke.ravel(), (self.rows, self.cols)), shape=(3 * len(V),) * 2).tocsc()
        Kr = (self.Q.T @ K @ self.Q).tocsc()
        fr = np.asarray(self.Q.T @ self.load).ravel()
        ix = self.reduced_free
        fac = splu(Kr[ix][:, ix])
        ur = np.zeros(Kr.shape[0])
        ur[ix] = fac.solve(fr[ix])
        u = np.asarray(self.Q @ ur).ravel()
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
        out = dict(J=J, hazard=H, volume_mm3=float(vol[:self.nc_elements].sum()), gradient=None, vertices=V, positive_principal=sp[:self.nc_elements], stress_MPa_per_N=s[:self.nc_elements], equilibrium_relative_residual=float(np.linalg.norm((Kr @ ur - fr)[ix]) / np.linalg.norm(fr[ix])))
        if gradient:
            spec = np.c_[v[:, 0] ** 2, v[:, 1] ** 2, v[:, 2] ** 2, 2 * v[:, 0] * v[:, 1], 2 * v[:, 1] * v[:, 2], 2 * v[:, 0] * v[:, 2]]
            dHs = weight[:, None] * self.m * (sp ** (self.m - 1))[:, None] * spec
            rhs = np.zeros(len(u))
            np.add.at(rhs, self.edof.ravel(), (np.einsum('ei,eij->ej', dHs, CB) / (self.m * H)).ravel())
            rr = np.asarray(self.Q.T @ rhs).ravel()
            lr = np.zeros(Kr.shape[0])
            lr[ix] = fac.solve(rr[ix])
            lam = np.asarray(self.Q @ lr).ravel()
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
