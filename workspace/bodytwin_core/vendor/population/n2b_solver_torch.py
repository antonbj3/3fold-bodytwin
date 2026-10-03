"""N2b copy of results/N2/n2_solver.py (unchanged). N2 batched muscle recruitment (polynomial and min-max criteria) for many problems at once. PyTorch float64, CPU or CUDA.

Problem per batch item:  min sum_j a_j^p   s.t.  Ap a = b (m rows),  0 <= a <= U      (a = f / F0)
Dual (m-dim, concave): a(lam) = clip(((Ap^T lam)^+ / p)^(1/(p-1)), 0, U).
Solved with semismooth Newton on J(lam) = -lam.b - sum_j phi(s_j),  phi(s) = a^p - s a,  grad J = Ap a - b,
Hess J = Ap diag(da/ds) Ap^T, Armijo backtracking. Every problem is iterated independently (frozen once converged),
so the arithmetic per problem does not depend on the other problems in the batch except through kernel choice.

min/max: U* by bisection on the box QP (p=2) feasibility; each QP also yields the certified dual lower bound
lam.b / sum (Ap^T lam)^+ <= U*, and every converged QP a certified feasible upper bound. Then the K1 tie-break
min sum a^2 s.t. 0 <= a <= U*(1+1e-6).
"""
import torch


def _amap(s, p, U):
    sp = torch.clamp(s, min=0.0)
    if p == 2:
        a0 = sp / 2.0
    else:
        a0 = (sp / p) ** (1.0 / (p - 1.0))
    a = torch.minimum(a0, U[:, None])
    inner = (a0 > 0) & (a0 < U[:, None])
    if p == 2:
        w = torch.where(inner, torch.full_like(s, 0.5), torch.zeros_like(s))
    else:
        # curvature da/ds = a/((p-1)s) is unbounded as s -> 0+; evaluate it at s >= 1e-12 max(s) (see n2_np)
        sfl = torch.maximum(s, 1e-12 * sp.amax(-1, keepdim=True))
        w = torch.where(inner, (torch.minimum(sfl, p * U[:, None] ** (p - 1)) / p) ** (1.0 / (p - 1.0)) / ((p - 1.0) * sfl), torch.zeros_like(s))
    return a, w


def _J(lam, s, a, b, p):
    return -(lam * b).sum(-1) - (a ** p - s * a).sum(-1)


def solve_poly(Ap, b, p, U=None, tol=1e-14, maxit=200, lam0=None, stop_infeasible=False):
    """Ap (B,m,n), b (B,m). Returns dict a (B,n), lam (B,m), resid (B,) = |Ap a - b|/|b|, iters (B,), conv (B,)."""
    B, m, n = Ap.shape
    dt, dev = Ap.dtype, Ap.device
    U = torch.full((B,), float('inf'), dtype=dt, device=dev) if U is None else U
    sc = torch.linalg.norm(b, dim=-1).clamp_min(1e-300)
    bh = b / sc[:, None]; Ah = Ap / sc[:, None, None]          # same a, normalized rows
    G = Ah @ Ah.transpose(1, 2)
    reg0 = 1e-12 * torch.diagonal(G, dim1=1, dim2=2).sum(-1) / m
    eye = torch.eye(m, dtype=dt, device=dev)
    if lam0 is None:
        lam = 2.0 * torch.linalg.solve(G + reg0[:, None, None] * eye, bh)   # unconstrained p=2 optimum
        if p != 2:   # rescale so the start has the activity scale of the p=2 optimum: s = p a^(p-1)
            s2 = (Ah.transpose(1, 2) @ lam[..., None])[..., 0]
            atyp = torch.clamp(s2, min=0).sum(-1) / 2.0 / (s2 > 0).sum(-1).clamp_min(1)
            atyp = atyp.clamp_min(1e-12)
            lam = lam * (p * atyp ** (p - 1.0) / (2.0 * atyp))[:, None]
    else:
        lam = lam0.clone()
    mu = torch.zeros(B, dtype=dt, device=dev)   # Levenberg damping, raised when the line search fails
    s = (Ah.transpose(1, 2) @ lam[..., None])[..., 0]
    a, w = _amap(s, p, U)
    grad = (Ah @ a[..., None])[..., 0] - bh
    J = _J(lam, s, a, bh, p)
    active = torch.linalg.norm(grad, dim=-1) > tol
    iters = torch.zeros(B, dtype=torch.long, device=dev)
    for it in range(maxit):
        if not bool(active.any()):
            break
        idx = active.nonzero()[:, 0]
        Ai, bi, li, gi, wi, Ji, Ui = Ah[idx], bh[idx], lam[idx], grad[idx], w[idx], J[idx], U[idx]
        H = (Ai * wi[:, None, :]) @ Ai.transpose(1, 2)
        trH = torch.diagonal(H, dim1=1, dim2=2).sum(-1)
        reg = reg0[idx] + 1e-12 * trH + mu[idx] * (trH / m + reg0[idx])
        d = -torch.linalg.solve(H + reg[:, None, None] * eye, gi)
        slope = (gi * d).sum(-1)
        t = torch.ones(len(idx), dtype=dt, device=dev)
        done = torch.zeros(len(idx), dtype=torch.bool, device=dev)
        ln, sn, an, wn, Jn = li.clone(), None, None, None, Ji.clone()
        sn = torch.empty_like(s[idx]); an = torch.empty_like(sn); wn = torch.empty_like(sn)
        gnorm = torch.linalg.norm(gi, dim=-1)
        for _ in range(40):
            lt = li + t[:, None] * d
            st = (Ai.transpose(1, 2) @ lt[..., None])[..., 0]
            at, wt = _amap(st, p, Ui)
            Jt = _J(lt, st, at, bi, p)
            # Armijo on J; near the optimum J stops resolving (decrease ~ |grad|^2 < eps|J|), so a step that
            # strictly reduces |grad| is also accepted there.
            gt = torch.linalg.norm((Ai @ at[..., None])[..., 0] - bi, dim=-1)
            ok = (~done) & ((Jt <= Ji + 1e-4 * t * slope) | ((gnorm < 1e-6) & (gt < (1 - 1e-4 * t) * gnorm)))
            ln[ok], sn[ok], an[ok], wn[ok], Jn[ok] = lt[ok], st[ok], at[ok], wt[ok], Jt[ok]
            done |= ok
            if bool(done.all()):
                break
            t = torch.where(done, t, t * 0.5)
        # items whose line search failed keep lam (stalled); recompute their state
        if not bool(done.all()):
            nd = ~done
            sn[nd] = (Ai[nd].transpose(1, 2) @ li[nd][..., None])[..., 0]
            an[nd], wn[nd] = _amap(sn[nd], p, Ui[nd])
        gn = (Ai @ an[..., None])[..., 0] - bi
        lam[idx], s[idx], a[idx], w[idx], J[idx], grad[idx] = ln, sn, an, wn, Jn, gn
        iters[idx] += 1
        mu[idx] = torch.where(done, mu[idx] * 0.1, torch.clamp(mu[idx] * 10.0, min=1e-6))
        still = (torch.linalg.norm(gn, dim=-1) > tol) & (done | (mu[idx] < 1e8))
        if stop_infeasible:   # certificate lam.b > U sum(s+) proves b outside Ap[0,U]^n: stop iterating
            cert = (ln * bi).sum(-1) > Ui * torch.clamp(sn, min=0).sum(-1) * (1 + 1e-12)
            still = still & ~cert
        active[idx] = still
    resid = torch.linalg.norm(grad, dim=-1)
    return dict(a=a, lam=lam, s=s, resid=resid, iters=iters, conv=resid <= 1e-10, sc=sc)


def dual_lower_bound(Ap, b, lam_n):
    """Certified lower bound on U* from any dual vector (normalized problem): lam.b / sum (Ap^T lam)^+."""
    sc = torch.linalg.norm(b, dim=-1).clamp_min(1e-300)
    s = ((Ap / sc[:, None, None]).transpose(1, 2) @ lam_n[..., None])[..., 0]
    den = torch.clamp(s, min=0).sum(-1)
    num = (lam_n * b / sc[:, None]).sum(-1)
    return torch.where(den > 0, num / den, torch.where(num > 0, torch.full_like(num, float('inf')), torch.zeros_like(num)))


def solve_minmax(Ap, b, rel_gap=1e-9, max_bisect=60, qp_maxit=80):
    """Returns dict a (tie-break solution), Ustar_ub, Ustar_lb, feasible (cone), resid."""
    B = Ap.shape[0]
    base = solve_poly(Ap, b, 2)                      # no upper bound: feasible iff b in cone
    cone_ok = base['conv']
    hi = torch.where(cone_ok, base['a'].max(-1).values, torch.full_like(base['resid'], float('inf')))
    lo = torch.clamp(dual_lower_bound(Ap, b, base['lam']), min=0.0)
    lo = torch.where(cone_ok, torch.minimum(lo, hi), lo)
    lam_ws = base['lam'].clone()
    for k in range(max_bisect):
        act = cone_ok & ((hi - lo) > rel_gap * hi)
        if not bool(act.any()):
            break
        idx = act.nonzero()[:, 0]
        mid = 0.5 * (lo[idx] + hi[idx])
        r = solve_poly(Ap[idx], b[idx], 2, U=mid, maxit=qp_maxit, lam0=lam_ws[idx], stop_infeasible=True)
        lb = dual_lower_bound(Ap[idx], b[idx], r['lam'])
        feas = r['conv']
        hi[idx] = torch.where(feas, torch.minimum(hi[idx], r['a'].max(-1).values.clamp_max(mid)), hi[idx])
        # infeasible: dual certificate lb > mid proves it; otherwise non-convergence is taken as infeasible (counted)
        lo[idx] = torch.where(feas, torch.maximum(lo[idx], torch.minimum(lb, mid)), torch.maximum(lo[idx], mid))
        lam_ws[idx] = torch.where(feas[:, None], r['lam'], lam_ws[idx])
    Ut = torch.where(cone_ok, hi * (1 + 1e-6), torch.ones_like(hi))
    r = solve_poly(Ap, b, 2, U=Ut, lam0=lam_ws)
    return dict(a=r['a'], lam=r['lam'], s=r['s'], resid=r['resid'], conv=r['conv'] & cone_ok, U_ub=hi, U_lb=lo,
                cone_ok=cone_ok, U=Ut)


def kkt_report(Ap, b, out, p, U):
    """Explicit KKT check on the normalized problem. Returns per-problem primal residual, bound violation,
    complementarity violation and interior stationarity residual (relative)."""
    a, s = out['a'], out['s']
    U = U[:, None] if U.dim() == 1 else U
    viol = torch.clamp(-a, min=0).amax(-1) + torch.clamp(a - U, min=0).amax(-1)
    tolz = 1e-12
    at0 = a <= tolz; atU = a >= U * (1 - 1e-12); inner = ~(at0 | atU)
    pa = p * a ** (p - 1)
    comp0 = torch.where(at0, torch.clamp(s, min=0), torch.zeros_like(s)).amax(-1)       # a=0 needs s<=0
    compU = torch.where(atU & torch.isfinite(U), torch.clamp(p * U ** (p - 1) - s, min=0), torch.zeros_like(s)).amax(-1)
    stat = torch.where(inner, (pa - s).abs() / s.abs().clamp_min(1e-300), torch.zeros_like(s)).amax(-1)
    return dict(primal=out['resid'], bound_viol=viol, comp_lo=comp0, comp_hi=compU, stationarity=stat)
