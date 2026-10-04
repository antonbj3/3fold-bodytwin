"""Emulator error bounds at active-set switches (CLOUD-W4-A303-BOUND).

Original synthetic parametric convex QPs only. Run:  python3 emu_bound.py
Deps: numpy, scipy, quadprog.  Writes results.json.
"""
import json, math, sys, time
import numpy as np
import quadprog
from scipy.optimize import linprog, minimize
from scipy.interpolate import RBFInterpolator

TOL = 1e-9
STRATA = [(0.0, 0.005), (0.005, 0.02), (0.02, 0.05), (0.05, math.inf)]
N_PER_STRATUM = 100      # per instance
N_CAL = 1000             # per instance
N_UNIF = 500             # per instance (secondary)
NDIR, TSTEP, TMAX = 32, 0.001, 0.06
FD = 1e-6


# ---------------------------------------------------------------- families
class Muscle:
    """Planar arm, 8 muscles, 2 torque equalities, 0<=f<=N(q)."""
    name, n = "muscle", 8

    def __init__(self, seed):
        r = np.random.default_rng(seed)
        self.N0 = r.uniform(300, 1200, 8)
        self.phi = r.uniform(0, 2 * np.pi, 8)
        s = np.array([1, 1, -1, -1, 1, -1, 1, -1.0])            # flex/ext
        self.a_sh = s * r.uniform(0.01, 0.04, 8) * np.array([1, 1, 1, 1, 1, 1, 0, 0])
        self.a_el = np.r_[s[:6] * r.uniform(0.005, 0.03, 6), s[6:] * r.uniform(0.01, 0.04, 2)]
        self.b_el = r.uniform(0.005, 0.015, 8)
        self.c_el = r.uniform(-1, 1, 8)
        self.a_el[:2] = s[:2] * r.uniform(0.01, 0.03, 2)      # biarticular pair
        self.mmax = 1.0
        self.w = np.ones(8)
        qs = np.linspace(0, 1, 41)
        self.mmax = 0.9 * min(self._maxload(t) for t in qs)

    def _geom(self, th):
        q = 0.3 + 1.8 * th[:, 0]
        el = self.a_el + np.sign(self.a_el) * self.b_el * np.sin(q[:, None] + self.c_el)
        sh = self.a_sh * (1 + 0.2 * np.cos(q[:, None] + self.phi))
        R = np.stack([sh, el], 1)                               # (N,2,8)
        N = self.N0 * (1 + 0.3 * np.cos(q[:, None] - self.phi))
        l1, l2, q1 = 0.30, 0.35, 0.5
        J2 = l2 * np.cos(q1 + q)                                 # x-lever for vertical load
        J1 = l1 * np.cos(q1) + J2
        return q, R, N, J1, J2

    def _tau(self, th, J1, J2):
        m = 1.0 + (self.mmax - 1.0) * th[:, 1]
        g = 9.81
        return np.stack([(m + 1.0) * g * J1, (m + 0.5) * g * J2], 1)

    def _maxload(self, t1):
        th = np.array([[t1, 0.0]])
        q, R, N, J1, J2 = self._geom(th)
        # maximise m s.t. R f = [(m+1) g J1, (m+.5) g J2], 0<=f<=N
        g = 9.81
        A = np.c_[R[0], -g * np.array([J1[0], J2[0]])]
        b = g * np.array([J1[0], 0.5 * J2[0]])
        res = linprog(np.r_[np.zeros(8), -1], A_eq=A, b_eq=b,
                      bounds=[(0, N[0, i]) for i in range(8)] + [(0, None)])
        assert res.status == 0
        return res.x[-1]

    def data(self, th):
        th = np.atleast_2d(th)
        q, R, N, J1, J2 = self._geom(th)
        k = len(th)
        Q = np.zeros((k, 8, 8)); idx = np.arange(8)
        Q[:, idx, idx] = 2.0 / N ** 2
        G = np.broadcast_to(np.r_[-np.eye(8), np.eye(8)], (k, 16, 8)).copy()
        h = np.c_[np.zeros((k, 8)), N]
        return dict(Q=Q, c=np.zeros((k, 8)), A=R, b=self._tau(th, J1, J2), G=G, h=h)


class Generic:
    """Random strictly convex pQP, n=6, 11 inequalities, nonlinear data."""
    name, n = "generic", 6

    def __init__(self, seed):
        r = np.random.default_rng(seed)
        n = 6
        M = r.normal(size=(n, n)); self.Q0 = M.T @ M / n + 0.5 * np.eye(n)
        P1 = r.normal(size=(n, n)); self.P1 = P1.T @ P1 / n
        P2 = r.normal(size=(n, n)); self.P2 = P2.T @ P2 / n
        self.c0 = r.normal(size=n)
        self.C = 2.0 * r.normal(size=(n, 6))
        self.Gg = r.normal(size=(4, n))
        self.p0 = r.uniform(0.2, 0.6, n)
        self.hs = r.uniform(0.05, 0.4, 4)
        self.hf = r.normal(size=(4, 2))
        self.w = r.uniform(0.5, 1.5, n)

    def data(self, th):
        th = np.atleast_2d(th); k = len(th); n = 6
        t1, t2 = th[:, 0], th[:, 1]
        Q = (self.Q0 + (t1 ** 2)[:, None, None] * self.P1
             + np.sin(np.pi * t2)[:, None, None] * self.P2)
        phi = np.stack([t1, t2, t1 * t2, np.cos(2 * t1), np.sin(2 * t2), t1 ** 2], 1)
        c = self.c0 + phi @ self.C.T
        s = 1 + 0.5 * t1
        p = self.p0 * (s / self.p0.sum())[:, None] * 1.3          # strictly feasible point
        hg = p @ self.Gg.T + self.hs + 0.3 * np.sin(th @ self.hf.T * 3) ** 2
        G = np.broadcast_to(np.r_[-np.eye(n), -np.ones((1, n)), self.Gg], (k, 11, n)).copy()
        h = np.c_[np.zeros((k, n)), -s, hg]
        return dict(Q=Q, c=c, A=np.zeros((k, 0, n)), b=np.zeros((k, 0)), G=G, h=h)


def one(d, i):
    return {k: v[i] for k, v in d.items()}


# ---------------------------------------------------------------- exact solve
def kkt(Q, c, A, b, G, h, S):
    """Equality-constrained solve with inequalities S active. Returns x, nu_eq, lam_S."""
    AS = np.r_[A, G[S]]; bS = np.r_[b, h[S]]
    n, m = Q.shape[0], AS.shape[0]
    K = np.zeros((n + m, n + m)); K[:n, :n] = Q; K[:n, n:] = AS.T; K[n:, :n] = AS
    z = np.linalg.solve(K, np.r_[-c, bS])
    return z[:n], z[n:n + len(b)], z[n + len(b):]


def solve_exact(d):
    Q, c, A, b, G, h = d["Q"], d["c"], d["A"], d["b"], d["G"], d["h"]
    Cm = np.r_[A, -G].T
    bm = np.r_[b, -h]
    x0 = quadprog.solve_qp(Q, -c, Cm, bm, len(b))[0]
    slack = h - G @ x0
    S = np.where(slack < 1e-8 * (1 + np.abs(h)))[0]
    # polish: drop active constraints with negative multipliers (degenerate ties)
    for _ in range(len(S) + 1):
        x, nu, lam = kkt(Q, c, A, b, G, h, S)
        if len(lam) and lam.min() < -TOL:
            S = np.delete(S, np.argmin(lam)); continue
        break
    lam_full = np.zeros(len(h)); lam_full[S] = lam
    slack = h - G @ x
    stat = Q @ x + c + A.T @ nu + G.T @ lam_full
    sc = 1 + np.abs(c).max() + np.abs(Q).max()
    ok = (slack.min() > -1e-7 * (1 + np.abs(h).max()) and lam_full.min() > -1e-7 * sc
          and np.abs(stat).max() < 1e-8 * sc and np.abs(np.r_[0.0, A @ x - b]).max() < 1e-8 * (1 + np.abs(np.r_[0.0, b]).max())
          and np.abs(x - x0).max() < 1e-5 * (1 + np.abs(x0).max()))
    return dict(x=x, nu=nu, lam=lam_full, S=S, slack=slack, ok=bool(ok))


# ---------------------------------------------------------------- one-sided derivative
def data_dir(fam, th0, u):
    dp = fam.data(np.array([th0 + FD * u, th0 - FD * u]))
    return {k: (v[0] - v[1]) / (2 * FD) for k, v in dp.items()}


def one_sided(fam, th0, sol, u):
    """Directional derivative of x*, multipliers, slacks at th0 along unit u."""
    d = one(fam.data(th0[None]), 0); dd = data_dir(fam, th0, u)
    x, nu, lam, S = sol["x"], sol["nu"], sol["lam"], sol["S"]
    sc = 1 + np.abs(np.r_[0.0, lam]).max()
    Sp = np.array([i for i in S if lam[i] > 1e-9 * sc], int)     # strongly active
    W = np.array([i for i in range(len(d["h"])) if i not in Sp and sol["slack"][i] < 1e-9], int)
    Q = d["Q"]
    rhs_lin = dd["Q"] @ x + dd["c"] + dd["A"].T @ nu + dd["G"][Sp].T @ lam[Sp]
    Ae = np.r_[d["A"], d["G"][Sp]]
    be = np.r_[dd["b"] - dd["A"] @ x, dd["h"][Sp] - dd["G"][Sp] @ x]
    if len(W) == 0:
        n, m = Q.shape[0], len(be)
        K = np.zeros((n + m, n + m)); K[:n, :n] = Q; K[:n, n:] = Ae.T; K[n:, :n] = Ae
        z = np.linalg.solve(K, np.r_[-rhs_lin, be]); dx, dmu, dw = z[:n], z[n:], np.zeros(0)
    else:
        Gi = d["G"][W]; hi = dd["h"][W] - dd["G"][W] @ x
        Cm = np.r_[Ae, -Gi].T
        dx, _, _, _, mult, _ = quadprog.solve_qp(Q, -rhs_lin, Cm, np.r_[be, -hi], len(be))
        dmu = mult[:len(be)]; dw = mult[len(be):]
    dlam = np.zeros(len(d["h"]))
    dlam[Sp] = dmu[d["A"].shape[0]:]
    if len(W): dlam[W] = dw
    dslack = dd["h"] - dd["G"] @ x - d["G"] @ dx
    return dx, dlam, dslack


# ---------------------------------------------------------------- distance to switch
def in_region(fam, ths, S, nb):
    """Boolean per theta: frozen active set S still KKT-optimal."""
    d = fam.data(ths); k = len(ths)
    AS = np.concatenate([d["A"], d["G"][:, S]], 1); bS = np.c_[d["b"], d["h"][:, S]]
    n, m = d["Q"].shape[1], AS.shape[1]
    K = np.zeros((k, n + m, n + m)); K[:, :n, :n] = d["Q"]
    K[:, :n, n:] = AS.transpose(0, 2, 1); K[:, n:, :n] = AS
    z = np.linalg.solve(K, np.c_[-d["c"], bS][..., None])[..., 0]
    x, lam = z[:, :n], z[:, n + nb:]
    slack = d["h"] - np.einsum("kij,kj->ki", d["G"], x)
    mask = np.ones(slack.shape[1], bool); mask[S] = False
    okS = (lam > -1e-9 * (1 + np.abs(lam).max(1, keepdims=True))).all(1) if len(S) else np.ones(k, bool)
    okI = (slack[:, mask] > -1e-9 * (1 + np.abs(d["h"]).max(1, keepdims=True))).all(1)
    return okS & okI


ANG = np.linspace(0, 2 * np.pi, NDIR, endpoint=False)
DIRS = np.c_[np.cos(ANG), np.sin(ANG)]
TS = np.arange(1, int(round(TMAX / TSTEP)) + 1) * TSTEP


def dist_to_switch(fam, th, sol):
    nb = fam.data(th[None])["A"].shape[1]
    S = sol["S"]
    pts = th + (DIRS[:, None, :] * TS[None, :, None]).reshape(-1, 2)
    ok = in_region(fam, pts, S, nb).reshape(NDIR, len(TS))
    best = math.inf
    for j in range(NDIR):
        bad = np.where(~ok[j])[0]
        if not len(bad): continue
        k = bad[0]
        lo, hi = (TS[k - 1] if k else 0.0), TS[k]
        if lo >= best: continue
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if in_region(fam, (th + mid * DIRS[j])[None], S, nb)[0]: lo = mid
            else: hi = mid
        best = min(best, 0.5 * (lo + hi))
    return best


# ---------------------------------------------------------------- emulators
class Emulators:
    def __init__(self, fam):
        self.fam = fam
        g = np.linspace(0, 1, 11)
        self.anch = np.array([[a, b] for a in g for b in g])
        self.sols = [solve_exact(one(fam.data(t[None]), 0)) for t in self.anch]
        assert all(s["ok"] for s in self.sols)
        self.y = np.array([fam.w @ s["x"] for s in self.sols])
        self.rbf = RBFInterpolator(self.anch, self.y, kernel="thin_plate_spline", degree=1)

    def predict(self, th):
        i = int(np.argmin(((self.anch - th) ** 2).sum(1)))
        th0, s0 = self.anch[i], self.sols[i]
        v = th - th0; h = float(np.linalg.norm(v)); u = v / h
        dx, dlam, dslack = one_sided(self.fam, th0, s0, u)
        e1 = self.y[i] + h * (self.fam.w @ dx)
        d = one(self.fam.data(th[None]), 0)
        x2, _, _ = kkt(d["Q"], d["c"], d["A"], d["b"], d["G"], d["h"], s0["S"])
        b2 = self.fam.w @ x2
        b1 = float(self.rbf(th[None])[0])
        lin_sl = s0["slack"] + h * dslack; lin_lam = s0["lam"] + h * dlam
        inact = np.ones(len(lin_sl), bool); inact[s0["S"]] = False
        flag = bool((lin_sl[inact] < 0).any() or (lin_lam[s0["S"]] < 0).any())
        return dict(E1=e1, B1=b1, B2=b2, h=h, flag=flag, anchor=i)


# ---------------------------------------------------------------- pipeline
def make_family(kind, k):
    """Deterministic seed with rejection rule from ANALYSIS_PLAN.md."""
    seed = 1000 * (1 if kind == "muscle" else 2) + 10 * k
    while True:
        fam = Muscle(seed) if kind == "muscle" else Generic(seed)
        g = np.linspace(0, 1, 21)
        try:
            sets = set()
            for a in g:
                for b in g:
                    s = solve_exact(one(fam.data(np.array([[a, b]])), 0))
                    if not s["ok"]: raise ValueError("kkt")
                    sets.add(tuple(s["S"]))
        except ValueError:
            seed += 1; continue
        if len(sets) >= 4:
            fam.seed, fam.nsets = seed, len(sets)
            return fam
        seed += 1


def evaluate_point(fam, em, th):
    d = one(fam.data(th[None]), 0)
    sol = solve_exact(d)
    y = fam.w @ sol["x"]
    p = em.predict(th)
    p.update(y=float(y), ok=sol["ok"], theta=th.tolist())
    return p, sol


def stratum_of(dist):
    for j, (a, b) in enumerate(STRATA):
        if a <= dist < b: return j
    return len(STRATA) - 1


def main():
    t0 = time.time()
    out = dict(instances=[], cal=[], test=[], unif=[], indep=[], fdcheck=[])
    for kind in ["muscle", "generic"]:
        for k in range(3):
            fam = make_family(kind, k)
            em = Emulators(fam)
            inst = f"{kind}{k}"
            out["instances"].append(dict(id=inst, family=kind, seed=fam.seed, n_active_sets_21x21=fam.nsets,
                                         y_range=[float(em.y.min()), float(em.y.max())]))
            r = np.random.default_rng(fam.seed + 7)
            # calibration (uniform)
            for th in r.uniform(0, 1, (N_CAL, 2)):
                p, _ = evaluate_point(fam, em, th); p.update(inst=inst, family=kind); out["cal"].append(p)
            # uniform test (secondary)
            for th in r.uniform(0, 1, (N_UNIF, 2)):
                p, sol = evaluate_point(fam, em, th)
                p.update(inst=inst, family=kind, dist=float(dist_to_switch(fam, th, sol)))
                out["unif"].append(p)
            # stratified test
            cnt = [0] * len(STRATA); ncand = 0
            while min(cnt) < N_PER_STRATUM:
                th = r.uniform(0, 1, 2); ncand += 1
                d = one(fam.data(th[None]), 0); sol = solve_exact(d)
                dist = dist_to_switch(fam, th, sol); j = stratum_of(dist)
                if cnt[j] >= N_PER_STRATUM: continue
                cnt[j] += 1
                p, _ = evaluate_point(fam, em, th)
                p.update(inst=inst, family=kind, dist=float(dist), stratum=j)
                out["test"].append(p)
            out["instances"][-1]["n_candidates_stratified"] = ncand
            # independent check: SLSQP on 50 stratified test points
            tp = [p for p in out["test"] if p["inst"] == inst]
            for p in [tp[i] for i in r.choice(len(tp), 50, replace=False)]:
                d = one(fam.data(np.array([p["theta"]])), 0)
                cons = [dict(type="ineq", fun=lambda x, d=d: d["h"] - d["G"] @ x, jac=lambda x, d=d: -d["G"])]
                if d["A"].shape[0]:
                    cons.append(dict(type="eq", fun=lambda x, d=d: d["A"] @ x - d["b"], jac=lambda x, d=d: d["A"]))
                x0 = np.zeros(fam.n) + (d["h"][-1] if kind == "muscle" else 0.3)
                res = minimize(lambda x: 0.5 * x @ d["Q"] @ x + d["c"] @ x, x0, jac=lambda x: d["Q"] @ x + d["c"],
                               constraints=cons, method="SLSQP", options=dict(ftol=1e-14, maxiter=500))
                ys = float(fam.w @ res.x)
                out["indep"].append(dict(inst=inst, y_exact=p["y"], y_slsqp=ys, success=bool(res.success),
                                         rel_diff=abs(ys - p["y"]) / abs(p["y"])))
            # finite-difference check of the one-sided derivative at 30 anchors
            for i in r.choice(len(em.anch), 30, replace=False):
                u = DIRS[r.integers(NDIR)]; th0 = em.anch[i]
                dx, _, _ = one_sided(fam, th0, em.sols[i], u)
                eps = 1e-6
                ye = fam.w @ solve_exact(one(fam.data((th0 + eps * u)[None]), 0))["x"]
                fd = (ye - em.y[i]) / eps
                out["fdcheck"].append(dict(inst=inst, analytic=float(fam.w @ dx), fd=float(fd)))
            print(inst, "seed", fam.seed, "sets", fam.nsets, "cand", ncand, f"{time.time()-t0:.0f}s", flush=True)
    json.dump(out, open("raw_points.json", "w"))
    summarize(out)


# ---------------------------------------------------------------- stats
def wilson(k, n, z=1.959964):
    if n == 0: return [None, None]
    p = k / n; den = 1 + z * z / n
    ctr = (p + z * z / (2 * n)) / den; hw = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [ctr - hw, ctr + hw]


def conformal(cal, key, alpha=0.1):
    q = {}
    for fam in ["muscle", "generic"]:
        for fl in [False, True]:
            sc = np.array([abs(p[key] - p["y"]) / abs(p[key]) / (p["h"] if fl else p["h"] ** 2)
                           for p in cal if p["family"] == fam and p["flag"] == fl])
            n = len(sc)
            kq = math.ceil((n + 1) * (1 - alpha))
            q[(fam, fl)] = float(np.sort(sc)[kq - 1]) if 0 < kq <= n else math.inf
            q[(fam, fl, "n")] = n
    return q


def covered(p, key, q):
    s = p["h"] if p["flag"] else p["h"] ** 2
    return abs(p[key] - p["y"]) <= q[(p["family"], p["flag"])] * s * abs(p[key])


def summarize(out):
    res = dict(task="CLOUD-W4-A303-BOUND", units="relative error |y_hat-y|/|y| (dimensionless; y in N for muscle, "
               "arbitrary units for generic); distances in normalised parameter units",
               instances=out["instances"])
    assert all(p["ok"] for p in out["cal"] + out["test"] + out["unif"]), "KKT verification failed"
    res["kkt_verified_points"] = len(out["cal"]) + len(out["test"]) + len(out["unif"])
    em = {}
    for key in ["E1", "B1", "B2"]:
        q = conformal(out["cal"], key)
        e = {}
        for setname in ["test", "unif"]:
            pts = out[setname]
            rel = np.array([abs(p[key] - p["y"]) / abs(p["y"]) for p in pts])
            cov = np.array([covered(p, key, q) for p in pts])
            dist = np.array([p["dist"] for p in pts]); strat = np.array([stratum_of(x) for x in dist])
            fams = np.array([p["family"] for p in pts]); insts = np.array([p["inst"] for p in pts])
            blk = dict(n=len(pts), max_rel_err=float(rel.max()), p99_rel_err=float(np.quantile(rel, .99)),
                       median_rel_err=float(np.median(rel)), frac_rel_err_gt_3pct=float((rel > 0.03).mean()),
                       coverage=float(cov.mean()), strata=[])
            for j, (a, b) in enumerate(STRATA):
                m = strat == j
                sb = dict(stratum=f"S{j+1} [{a},{b})", n=int(m.sum()))
                if m.sum():
                    sb.update(coverage=float(cov[m].mean()), wilson95=wilson(int(cov[m].sum()), int(m.sum())),
                              max_rel_err=float(rel[m].max()), median_rel_err=float(np.median(rel[m])),
                              per_family={f: float(cov[m & (fams == f)].mean()) for f in ["muscle", "generic"] if (m & (fams == f)).any()},
                              per_family_max_rel_err={f: float(rel[m & (fams == f)].max()) for f in ["muscle", "generic"] if (m & (fams == f)).any()},
                              per_instance_coverage_range=[float(min(cov[m & (insts == i)].mean() for i in set(insts[m]))),
                                                           float(max(cov[m & (insts == i)].mean() for i in set(insts[m])))])
                blk["strata"].append(sb)
            if setname == "test":
                worst = np.argsort(-rel)[:5]
                blk["worst_cases"] = [dict(inst=pts[i]["inst"], theta=pts[i]["theta"], dist=pts[i]["dist"], h=pts[i]["h"],
                                           flag=pts[i]["flag"], y=pts[i]["y"], y_hat=pts[i][key], rel_err=float(rel[i]))
                                      for i in worst]
            e[setname] = blk
        e["conformal_q"] = {f"{a}|flag={b}": v for (a, b, *r), v in q.items() if not r}
        e["conformal_ncal"] = {f"{a}|flag={b}": v for (a, b, *r), v in q.items() if r}
        em[key] = e
    res["emulators"] = em
    ind = out["indep"]
    res["independent_check_slsqp"] = dict(n=len(ind), n_success=sum(p["success"] for p in ind),
                                          max_rel_diff=max(p["rel_diff"] for p in ind),
                                          median_rel_diff=float(np.median([p["rel_diff"] for p in ind])))
    fd = out["fdcheck"]
    res["fd_check_one_sided_derivative"] = dict(n=len(fd), max_abs_diff_rel=max(abs(p["analytic"] - p["fd"]) / (1 + abs(p["fd"])) for p in fd))
    # exploratory (not preregistered): split by whether anchor->query crosses a switch,
    # proxied by the frozen-set re-solve B2 being exact (rel err < 1e-9)
    xp = {}
    for setname in ["test", "unif"]:
        P = out[setname]
        e1 = np.array([abs(p["E1"] - p["y"]) / abs(p["y"]) for p in P])
        same = np.array([abs(p["B2"] - p["y"]) / abs(p["y"]) < 1e-9 for p in P])
        fl = np.array([p["flag"] for p in P])
        xp[setname] = dict(n_no_crossing=int(same.sum()), E1_max_no_crossing=float(e1[same].max()),
                           n_crossing=int((~same).sum()), E1_max_crossing=float(e1[~same].max()),
                           flag_recall_on_crossings=float(fl[~same].mean()), flag_rate_no_crossing=float(fl[same].mean()),
                           n_E1_gt_3pct=int((e1 > 0.03).sum()), n_E1_gt_3pct_crossing=int(((e1 > 0.03) & ~same).sum()))
    res["exploratory_crossing_split"] = xp
    t = em["E1"]["test"]
    ca = t["max_rel_err"] <= 0.03 and t["n"] >= 1000
    cb = all(0.85 <= s["coverage"] <= 0.95 for s in t["strata"])
    cb_strict = all(0.85 <= v <= 0.95 for s in t["strata"] for v in s["per_family"].values())
    res["criterion"] = dict(a_max_rel_err_le_3pct_on_ge1000=bool(ca), b_coverage_85_95_every_stratum=bool(cb),
                            b_strict_per_family=bool(cb_strict), PASS=bool(ca and cb))
    json.dump(res, open("results.json", "w"), indent=1)
    print(json.dumps(res["criterion"]))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--summarize":
        summarize(json.load(open("raw_points.json")))
    else:
        main()
