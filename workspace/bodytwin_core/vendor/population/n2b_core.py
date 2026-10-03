"""N2b core: N14's full-body recruitment system for many individuals at once (numpy or torch backend).

Individual = N14 base system with muscle end nodes moved (attachment perturbation, straight end segments only) and/or
strength factors (strength law).  Per step t and individual:  A(delta) = A_t + P_t dW(delta),  b = b_t (rigid-body demand
is geometry-independent), Ap = A(delta) * N_eff.  dW: end node e moved by R_se delta (segment frame), unit force at e is
sign*u (u = unit vector end -> neighbour), neighbour gets the opposite change; moments about segment origins.
Joint forces (force on reaction segment a):  F_j = U_j - G_j (f - f_J) - Rp_j dW f   (r = r_J + R^+(W f_J - W' f)).
Shoulder = GH contact force = sum of GH-edge actuator force rows on the humerus.
Min/max: batched primal-dual interior point (Mehrotra) for U* = min U s.t. Aa = b, 0 <= a <= U, then the K1/N2 tie-break
box QP (p = 2, U*(1+1e-6)) with the copied N2 solver."""
import numpy as np

try:
    import torch
except Exception:   # cloud: no torch
    torch = None


class BE:
    """minimal backend shim: numpy ('np') or torch ('cpu'/'cuda')"""
    def __init__(self, kind='np'):
        self.kind = kind
        self.np = kind == 'np'
        if not self.np:
            self.dev = kind
            import n2b_solver_torch as S
        else:
            import n2b_solver_np as S
        self.S = S

    def asarray(self, x):
        x = np.asarray(x)
        if self.np:
            return x
        if x.dtype.kind in 'iub':
            return torch.as_tensor(x, device=self.dev)
        return torch.as_tensor(np.ascontiguousarray(x), dtype=torch.float64, device=self.dev)

    def tonp(self, x):
        return x if self.np else x.detach().cpu().numpy()

    def zeros(self, shape):
        return np.zeros(shape) if self.np else torch.zeros(shape, dtype=torch.float64, device=self.dev)

    def cross(self, a, b):
        return np.cross(a, b) if self.np else torch.linalg.cross(a, b, dim=-1)

    def norm(self, x, axis=-1, keepdims=False):
        return np.linalg.norm(x, axis=axis, keepdims=keepdims) if self.np else torch.linalg.norm(x, dim=axis, keepdim=keepdims)

    def einsum(self, s, *a):
        return np.einsum(s, *a, optimize=True) if self.np else torch.einsum(s, *a)

    def solve(self, K, r):
        return np.linalg.solve(K, r[..., None])[..., 0] if self.np else torch.linalg.solve(K, r)

    def amax(self, x, axis=-1):
        return x.max(axis) if self.np else x.amax(axis)

    def amin(self, x, axis=-1):
        return x.min(axis) if self.np else x.amin(axis)

    def where(self, c, a, b):
        return np.where(c, a, b) if self.np else torch.where(c, a, b)

    def maximum(self, a, b):
        return np.maximum(a, b) if self.np else torch.maximum(a, b if torch.is_tensor(b) else torch.as_tensor(b, dtype=a.dtype, device=a.device))

    def minimum(self, a, b):
        return np.minimum(a, b) if self.np else torch.minimum(a, b if torch.is_tensor(b) else torch.as_tensor(b, dtype=a.dtype, device=a.device))

    def sync(self):
        if not self.np and self.dev == 'cuda':
            torch.cuda.synchronize()


class Problem:
    def __init__(self, sys_path, aux_path, be):
        self.be = be
        z = np.load(sys_path, allow_pickle=True); x = np.load(aux_path, allow_pickle=True)
        self.names = [str(s) for s in z['names']]; self.region = [str(s) for s in z['region']]
        self.T, self.m, self.n = z['A'].shape
        self.mm = z['m']; assert (self.mm == self.m).all()
        self.Nnp = z['N'].T.copy()            # (T, n)
        self.FJnp = z['FJ'].T.copy()          # (T, n)
        self.A = be.asarray(z['A']); self.b = be.asarray(z['b'])
        self.N = be.asarray(self.Nnp); self.FJ = be.asarray(self.FJnp)
        R6 = x['P'].shape[2]; self.R6 = R6
        Pz = np.concatenate([x['P'], np.zeros(x['P'].shape[:2] + (1,))], 2)       # extra zero column = "no row"
        self.P = be.asarray(Pz)
        self.QTnp = x['QT']; self.mQ = x['mQ']; self.Pnp = x['P']
        self.joints = [str(s) for s in x['joints']]
        Rpz = np.concatenate([x['Rp'], np.zeros(x['Rp'].shape[:3] + (1,))], 3)
        self.Rp = be.asarray(Rpz); self.G = be.asarray(x['G']); self.UJ = be.asarray(x['UJ'])
        self.GH = {s: (be.asarray(x['GH_' + s]), be.asarray(x['GHcols_' + s].astype(np.int64))) for s in ('R', 'L')}
        self.ident = x['ident']
        # ends
        self.E = len(x['e_col'])
        self.e_col = x['e_col'].astype(np.int64); self.e_grp = x['e_grp'].astype(np.int64); self.groups = [str(s) for s in x['groups']]
        self.e_which = np.array([str(s) for s in x['e_which']])
        re_, rn_ = x['e_row_e'].astype(np.int64), x['e_row_n'].astype(np.int64)
        idx = np.zeros((self.E, 12), np.int64)
        for k in range(6):
            idx[:, k] = np.where(re_ >= 0, re_ + k, R6); idx[:, 6 + k] = np.where(rn_ >= 0, rn_ + k, R6)
        self.e_idx = be.asarray(idx)
        self.e_sign = be.asarray(x['e_sign'])
        self.PE = be.asarray(x['PE']); self.PN = be.asarray(x['PN']); self.RE = be.asarray(x['RE'])
        self.SRE = be.asarray(x['SRE']); self.SRN = be.asarray(x['SRN'])
        o = np.where(self.e_which == 'origin')[0]; i = np.where(self.e_which == 'insertion')[0]
        assert len(set(self.e_col[o])) == len(o) and len(set(self.e_col[i])) == len(i)
        self.sets = [(be.asarray(o), be.asarray(self.e_col[o])), (be.asarray(i), be.asarray(self.e_col[i]))]
        self.e_col_t = be.asarray(self.e_col); self.e_grp_t = be.asarray(self.e_grp)
        self.G_n = len(self.groups)

    # ---------------- geometry perturbation ----------------
    def dW(self, t, dl):
        """dl (M, G, 3) local end offsets (m) -> dw (M, E, 12) = [dF_e, dM_e, dF_n, dM_n] per unit actuator force"""
        be = self.be
        d = dl[:, self.e_grp_t, :]                                     # (M,E,3)
        dg = be.einsum('eij,mej->mei', self.RE[t], d)
        pe, pn = self.PE[t], self.PN[t]
        u = pn - pe; u = u / be.norm(u, keepdims=True)
        pe2 = pe[None] + dg; u2 = pn[None] - pe2; u2 = u2 / be.norm(u2, keepdims=True)
        sg = self.e_sign[None, :, None]
        dF = sg * (u2 - u[None])
        dMe = be.cross(pe2 - self.SRE[t][None], sg * u2) - be.cross((pe - self.SRE[t])[None], sg * u[None])
        dFn = -dF
        dMn = be.cross((pn - self.SRN[t])[None].expand_as(dFn) if not be.np else np.broadcast_to((pn - self.SRN[t])[None], dFn.shape), dFn)
        if be.np:
            return np.concatenate([dF, dMe, dFn, dMn], -1)
        return torch.cat([dF, dMe, dFn, dMn], -1)

    def system(self, t, dl=None, sfac=None):
        """returns Ap (M,m,n), b (M,m), Neff (M,n), dw (M,E,12) or None"""
        be = self.be
        M = 1 if dl is None else dl.shape[0]
        A = self.A[t][None].repeat(M, 1, 1) if not be.np else np.repeat(self.A[t][None], M, 0)
        dw = None
        if dl is not None:
            dw = self.dW(t, dl)
            Pg = self.P[t][:, self.e_idx]                              # (m, E, 12)
            dA = be.einsum('mek,Mek->Mme', Pg, dw)                     # (M, m, E)
            for ends, cols in self.sets:
                A[:, :, cols] += dA[:, :, ends]
        Neff = self.N[t][None] * (1.0 if sfac is None else sfac)
        if be.np:
            Neff = np.broadcast_to(Neff, (M, self.n)).copy()
            b = np.repeat(self.b[t][None], M, 0)
        else:
            Neff = Neff.expand(M, self.n).clone()
            b = self.b[t][None].repeat(M, 1)
        return A * Neff[:, None, :], b, Neff, dw

    def joint_forces(self, t, f, dw):
        """f (M,n) actuator forces -> dict name -> (M,3) force on segment a; plus GH_R, GH_L"""
        be = self.be
        out = {}
        df = f - self.FJ[t][None]
        for j, nm in enumerate(self.joints):
            F = self.UJ[j, t][None] - be.einsum('kn,Mn->Mk', self.G[j, t], df)
            if dw is not None:
                Rg = self.Rp[j, t][:, self.e_idx]                      # (3,E,12)
                F = F - be.einsum('keq,Meq,Me->Mk', Rg, dw, f[:, self.e_col_t])
            out[nm] = F
        for s in ('R', 'L'):
            Wg, cols = self.GH[s]
            out['gh_' + s.lower()] = be.einsum('kc,Mc->Mk', Wg[t], f[:, cols])
        return out

    def out_of_range(self, t, f, dw):
        """|(I - U_k U_k^T) Q^T dW f| / |b| for a few individuals (numpy)"""
        f = self.be.tonp(f); dw = self.be.tonp(dw); idx = self.be.tonp(self.e_idx)
        M = f.shape[0]; x = np.zeros((M, self.R6 + 1))
        for q in range(12):
            np.add.at(x.T, idx[:, q], (dw[:, :, q] * f[:, self.e_col]).T)
        x = x[:, :self.R6]
        qx = x @ self.QTnp[t, :self.mQ[t]].T; px = x @ self.Pnp[t].T
        bn = np.linalg.norm(self.be.tonp(self.b[t]))
        return np.sqrt(np.maximum((qx ** 2).sum(1) - (px ** 2).sum(1), 0)) / bn


# ---------------- batched min/max LP (primal-dual interior point, Mehrotra) ----------------
def lp_minmax(be, Ap, b, tol=1e-10, maxit=150):
    """U* = min U s.t. Ap a = b, 0 <= a <= U, batched. Returns dict U, a, lam (normalized), iters, conv, lb (certified
    dual bound lam.b / sum (Ap^T lam)^+ on the normalized problem)."""
    B, m, n = Ap.shape
    sc = be.norm(b); sc = be.maximum(sc, 1e-300)
    Ah = Ap / sc[:, None, None]; bh = b / sc[:, None]
    At = Ah.transpose(0, 2, 1) if be.np else Ah.transpose(1, 2)
    # start: least-norm solution shifted positive
    G = Ah @ At
    trG = (G.diagonal(0, 1, 2) if be.np else torch.diagonal(G, dim1=1, dim2=2)).sum(-1)
    eye = np.eye(m) if be.np else torch.eye(m, dtype=Ap.dtype, device=Ap.device)
    y = be.solve(G + (1e-12 * trG / m)[:, None, None] * eye, bh)
    a = (At @ y[..., None])[..., 0]
    abs_a = a if be.np else a
    sh = be.amax(-a) ; sh = be.maximum(sh, 0.0)
    scale = (abs(a) if be.np else a.abs()).mean(-1)
    a = a + (sh + 0.5 * scale)[:, None]
    U = be.amax(a) * 1.5
    w = U[:, None] - a
    z1 = (np.ones_like(a) if be.np else torch.ones_like(a)) / n
    z2 = z1 * 1.0
    lam = y * 0.0
    ones = z1 * 0 + 1
    act = (np.ones(B, bool) if be.np else torch.ones(B, dtype=torch.bool, device=Ap.device))
    iters = (np.zeros(B, np.int64) if be.np else torch.zeros(B, dtype=torch.long, device=Ap.device))

    def mv(x): return (Ah @ x[..., None])[..., 0]
    def mtv(x): return (At @ x[..., None])[..., 0]

    def step_len(x, dx):
        r = be.where(dx < 0, -x / be.where(dx < 0, dx, -ones), ones * 1e300)
        return be.minimum(be.amin(r), 1.0)

    for it in range(maxit):
        rp = bh - mv(a)                                  # A da = rp
        r1 = mtv(lam) + z1 - z2                          # -A^T dl - dz1 + dz2 = -( -A^T l - z1 + z2 ) = r1
        r2 = z2.sum(-1) - 1.0                            # -sum dz2 = -(1 - sum z2)
        mu = ((a * z1).sum(-1) + (w * z2).sum(-1)) / (2 * n)
        gap = 2 * n * mu / be.maximum(U, 1e-300)
        conv = (be.norm(rp) <= tol) & (be.norm(r1) <= tol * (1 + be.norm(z1))) & (abs(r2) <= tol if be.np else r2.abs() <= tol) & (gap <= tol)
        act = act & ~conv
        if not bool(act.any()):
            break
        D = z1 / a + z2 / w; e = z2 / w; E = e.sum(-1); Di = 1.0 / D
        ADe = mv(Di * e)
        K = (Ah * Di[:, None, :]) @ At
        Kf = be.zeros((B, m + 1, m + 1))
        Kf[:, :m, :m] = K; Kf[:, :m, m] = ADe; Kf[:, m, :m] = ADe; Kf[:, m, m] = (e * Di * e).sum(-1) - E

        def direction(r4, r5):
            g = r1 + r4 / a - r5 / w
            rhs = be.zeros((B, m + 1))
            rhs[:, :m] = rp - mv(Di * g)
            rhs[:, m] = -r2 - (r5 / w).sum(-1) - (e * Di * g).sum(-1)   # sum dz2 = 1 - sum z2
            sol = be.solve(Kf, rhs)
            dl, dU = sol[:, :m], sol[:, m]
            da = Di * (g + mtv(dl) + e * dU[:, None])
            dz1 = (r4 - z1 * da) / a
            dw = dU[:, None] - da
            dz2 = (r5 - z2 * dw) / w
            return da, dU, dw, dl, dz1, dz2

        da, dU, dw_, dl, dz1, dz2 = direction(-a * z1, -w * z2)
        ap = be.minimum(step_len(a, da), step_len(w, dw_)); ad = be.minimum(step_len(z1, dz1), step_len(z2, dz2))
        mu_aff = (((a + ap[:, None] * da) * (z1 + ad[:, None] * dz1)).sum(-1) + ((w + ap[:, None] * dw_) * (z2 + ad[:, None] * dz2)).sum(-1)) / (2 * n)
        sig = (mu_aff / be.maximum(mu, 1e-300)) ** 3
        smu = (sig * mu)[:, None]
        da, dU, dw_, dl, dz1, dz2 = direction(smu - a * z1 - da * dz1, smu - w * z2 - dw_ * dz2)
        ap = 0.995 * be.minimum(step_len(a, da), step_len(w, dw_)); ad = 0.995 * be.minimum(step_len(z1, dz1), step_len(z2, dz2))
        ap = be.where(act, ap, ap * 0); ad = be.where(act, ad, ad * 0)
        a = a + ap[:, None] * da; U = U + ap * dU; w = w + ap[:, None] * dw_   # w as own variable (no cancellation)
        lam = lam + ad[:, None] * dl; z1 = z1 + ad[:, None] * dz1; z2 = z2 + ad[:, None] * dz2
        iters = iters + act
    s = mtv(lam); den = be.maximum(s, 0.0).sum(-1); lb = (lam * bh).sum(-1) / be.maximum(den, 1e-300)
    return dict(U=U, a=a, lam=lam, lb=lb, iters=iters, conv=~act, resid=be.norm(bh - mv(a)))


def qp_box_ipm(be, Ap, b, U, tol=1e-10, maxit=300):
    """tie-break box QP  min 1/2 sum a^2  s.t. Ap a = b, 0 <= a <= U (per problem U), batched primal-dual IPM
    (Mehrotra), Schur complement m x m.  N2b: replaces the semismooth-Newton box QP, which stalls at U ~ U*
    because fewer than m muscles stay interior (degenerate Hessian)."""
    B, m, n = Ap.shape
    sc = be.maximum(be.norm(b), 1e-300)
    Ah = Ap / sc[:, None, None]; bh = b / sc[:, None]
    At = Ah.transpose(0, 2, 1) if be.np else Ah.transpose(1, 2)
    def mv(x): return (Ah @ x[..., None])[..., 0]
    def mtv(x): return (At @ x[..., None])[..., 0]
    Uc = U[:, None]
    a = 0.5 * Uc + 0 * Ah[:, 0, :]
    w = Uc - a
    z1 = 0 * a + 1.0; z2 = 0 * a + 1.0
    lam = 0 * bh
    ones = 0 * a + 1.0
    act = (np.ones(B, bool) if be.np else torch.ones(B, dtype=torch.bool, device=Ap.device))
    iters = (np.zeros(B, np.int64) if be.np else torch.zeros(B, dtype=torch.long, device=Ap.device))
    def step_len(x, dx):
        r = be.where(dx < 0, -x / be.where(dx < 0, dx, -ones), ones * 1e300)
        return be.minimum(be.amin(r), 1.0)
    for it in range(maxit):
        rp = bh - mv(a)
        rd = a - mtv(lam) - z1 + z2
        mu = ((a * z1).sum(-1) + (w * z2).sum(-1)) / (2 * n)
        conv = (be.norm(rp) <= tol) & (be.norm(rd) <= tol * (1 + be.norm(a))) & (2 * n * mu <= tol * (1 + (a * a).sum(-1)))
        act = act & ~conv
        if not bool(act.any()):
            break
        D = 1.0 + z1 / a + z2 / w; Di = 1.0 / D
        K = (Ah * Di[:, None, :]) @ At
        def direction(r4, r5):
            g = -rd + r4 / a - r5 / w
            dl = be.solve(K, rp - mv(Di * g))
            da = Di * (g + mtv(dl))
            return da, -da, dl, (r4 - z1 * da) / a, (r5 + z2 * da) / w
        da, dw_, dl, dz1, dz2 = direction(-a * z1, -w * z2)
        ap = be.minimum(step_len(a, da), step_len(w, dw_)); ad = be.minimum(step_len(z1, dz1), step_len(z2, dz2))
        al = be.minimum(ap, ad)
        mu_aff = (((a + al[:, None] * da) * (z1 + al[:, None] * dz1)).sum(-1) + ((w + al[:, None] * dw_) * (z2 + al[:, None] * dz2)).sum(-1)) / (2 * n)
        smu = ((mu_aff / be.maximum(mu, 1e-300)) ** 3 * mu)[:, None]
        da, dw_, dl, dz1, dz2 = direction(smu - a * z1 - da * dz1, smu - w * z2 - dw_ * dz2)
        al = 0.995 * be.minimum(be.minimum(step_len(a, da), step_len(w, dw_)), be.minimum(step_len(z1, dz1), step_len(z2, dz2)))
        al = be.where(act, al, al * 0)
        a = a + al[:, None] * da; w = w + al[:, None] * dw_
        lam = lam + al[:, None] * dl; z1 = z1 + al[:, None] * dz1; z2 = z2 + al[:, None] * dz2
        iters = iters + act
    resid = be.norm(bh - mv(a))
    return dict(a=a, lam=lam, s=mtv(lam), resid=resid, iters=iters, conv=(~act) & (resid <= 1e-10))


def solve_minmax_lp(be, Ap, b):
    lp = lp_minmax(be, Ap, b)
    Ut = lp['U'] * (1 + 1e-6)
    r = qp_box_ipm(be, Ap, b, Ut)
    r['U_ipm'] = lp['U']; r['U_lb'] = lp['lb']; r['ipm_iters'] = lp['iters']; r['ipm_conv'] = lp['conv']
    return r


def solve(be, Ap, b, crit):
    if crit == 'p2':      # tol 1e-12 on |Ap a - b|/|b| (conv criterion 1e-10); 1e-13 stalls single GPU items at maxit
        return be.S.solve_poly(Ap, b, 2, tol=1e-12)
    if crit == 'p3':
        return be.S.solve_poly(Ap, b, 3, tol=1e-12)
    if crit == 'minmax':
        return solve_minmax_lp(be, Ap, b)
    if crit == 'minmax_bisect':
        return be.S.solve_minmax(Ap, b)
    raise ValueError(crit)
