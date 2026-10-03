#!/usr/bin/env python3
"""Implicit KKT derivative parity on a seeded synthetic planar 2-body model.

Force solver (strongly convex, equality constrained QP):
    f*(p) = argmin_f 1/2 f^T W f   s.t.   A(p) f = b(p)
KKT:  [W  A^T][f]   [0]
      [A   0 ][l] = [b]
Implicit derivative w.r.t. a geometry parameter p_j (differentiate both rows):
    K [df; dl] = [ -dA^T l ;  db - dA f ]
dA/dp, db/dp are obtained with the complex-step method on the (analytic)
geometry, to machine precision; the solve itself is never differentiated.

Checks: (1) central finite differences at 5 step sizes of the full solver,
(2) complex-step through the whole solver (independent of the KKT formula).
Optional nonnegative-muscle mode ("nn"): exact active set by enumeration;
on a fixed active set the problem is again an equality-constrained QP.

Units: forces N, moments N*m, angles rad, attachment coordinate m.
Derivatives: N/rad (q1, q2) and N/m (attachment offset d).
Run:  python3 kkt_parity.py   (writes results.json)
"""
import itertools
import json
import sys

import numpy as np

G = 9.81
STEPS = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]    # multiplied by per-parameter scale
PSCALE = np.array([1.0, 1.0, 0.1])       # rad, rad, m
PNAMES = ["q1_rad", "q2_rad", "d_m"]
TOL_REL = 1e-4                           # preregistered threshold
COND_FLAG = 1e8                          # flag: cond(A) above this = (near) singular
SVAL_FLAG = 1e-10                        # flag: sigma_min(A)/sigma_max(A) below this
MARGIN_FLAG = 1e-9                       # nn mode: complementarity margin (N)
NMUS = 6


def rot(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s], [s, c]])


def cross2(r, f):
    return r[0] * f[1] - r[1] * f[0]


def make_case(seed, kind="regular"):
    rng = np.random.default_rng(seed)
    L1, L2 = rng.uniform(0.28, 0.36), rng.uniform(0.24, 0.30)
    c = dict(
        seed=seed, kind=kind, L1=L1, L2=L2,
        m1=rng.uniform(1.5, 3.0), m2=rng.uniform(1.0, 2.0),          # kg (varied mass)
        Fext=rng.uniform(-30, 30, 2), Mext=rng.uniform(-3, 3),        # N, N*m (varied wrench)
        s=rng.uniform(300, 1500, NMUS),                               # muscle strength N
        wR=1e-2,                                                      # reaction weight
        p0=np.array([rng.uniform(0.3, 1.2), rng.uniform(0.3, 1.4), 0.0]),
    )
    # muscle attachments: (origin body, local origin, insertion body, local insertion)
    # body 0 = ground (global coords), body 1/2 local frame along segment from proximal joint
    def off():
        return rng.uniform(0.02, 0.05) * rng.choice([-1, 1])
    o = []
    # mono-articular joint 0 (ground -> body1), agonist/antagonist
    o.append((0, np.array([rng.uniform(-0.12, -0.06), 0.04]), 1, np.array([rng.uniform(0.08, 0.15), 0.03])))
    o.append((0, np.array([rng.uniform(0.06, 0.12), -0.04]), 1, np.array([rng.uniform(0.08, 0.15), -0.03])))
    # mono-articular joint 1 (body1 -> body2)
    o.append((1, np.array([L1 * rng.uniform(0.4, 0.7), 0.035]), 2, np.array([rng.uniform(0.04, 0.08), 0.02])))
    o.append((1, np.array([L1 * rng.uniform(0.4, 0.7), -0.035]), 2, np.array([rng.uniform(0.04, 0.08), -0.02])))
    # bi-articular (ground -> body2)
    o.append((0, np.array([-0.05, 0.05]), 2, np.array([rng.uniform(0.05, 0.1), 0.03])))
    o.append((0, np.array([0.05, -0.05]), 2, np.array([rng.uniform(0.05, 0.1), -0.03])))
    c["mus"] = o
    if kind == "singular":
        # every joint-1 muscle and bi-articular muscle line runs through joint-1 centre:
        # body-2 moment balance about the joint cannot be generated -> rank(A) = 5
        c["mus"] = [o[0], o[1],
                    (1, np.array([0.5 * L1, 0.0]), 2, np.array([0.06, 0.0])),
                    (1, np.array([0.6 * L1, 0.0]), 2, np.array([0.07, 0.0])),
                    (1, np.array([0.3 * L1, 0.0]), 2, np.array([0.08, 0.0])),
                    (1, np.array([0.4 * L1, 0.0]), 2, np.array([0.09, 0.0]))]
        c["p0"][1] = 0.0          # straight chain: body1 and body2 collinear
    if kind == "near_singular":
        c["mus"] = [o[0], o[1],
                    (1, np.array([0.5 * L1, 1e-9]), 2, np.array([0.06, 0.0])),
                    (1, np.array([0.6 * L1, -1e-9]), 2, np.array([0.07, 0.0])),
                    (1, np.array([0.3 * L1, 0.0]), 2, np.array([0.08, 1e-9])),
                    (1, np.array([0.4 * L1, 0.0]), 2, np.array([0.09, -1e-9]))]
        c["p0"][1] = 0.0
    return c


def assemble(c, p):
    """A (6 x 10), b (6); works for complex p (complex-step)."""
    q1, q2, d = p
    R1, R2 = rot(q1), rot(q1 + q2)
    J0 = np.zeros(2)
    J1 = R1 @ np.array([c["L1"], 0.0])
    frames = {0: (np.zeros(2), np.eye(2)), 1: (J0, R1), 2: (J1, R2)}
    com = {1: R1 @ np.array([c["L1"] / 2, 0.0]), 2: J1 + R2 @ np.array([c["L2"] / 2, 0.0])}
    dt = np.result_type(p, float)
    A = np.zeros((6, NMUS + 4), dtype=dt)
    b = np.zeros(6, dtype=dt)

    def add(body, col, F, at, sign=1.0):
        if body == 0:
            return
        k = 3 * (body - 1)
        A[k:k + 2, col] += sign * F
        A[k + 2, col] += sign * cross2(at - com[body], F)

    for i, (bo, lo, bi, li) in enumerate(c["mus"]):
        lo = lo.copy().astype(dt)
        li = li.copy().astype(dt)
        if i == 0:
            li = li + np.array([d, 0.0])          # geometry parameter d: insertion offset
        po = frames[bo][0] + frames[bo][1] @ lo
        pi = frames[bi][0] + frames[bi][1] @ li
        u = (po - pi) / np.sqrt(np.sum((po - pi) ** 2))   # analytic norm (complex-safe)
        add(bi, i, u, pi, 1.0)
        add(bo, i, -u, po, 1.0)
    # reactions: R0 (ground on body1) columns 6,7; R1 (body1 on body2) columns 8,9
    for k in range(2):
        e = np.zeros(2); e[k] = 1.0
        add(1, NMUS + k, e, J0, 1.0)
        add(2, NMUS + 2 + k, e, J1, 1.0)
        add(1, NMUS + 2 + k, e, J1, -1.0)
    # loads: gravity on both bodies, external wrench at distal tip of body 2
    g = np.array([0.0, -G])
    tip = J1 + R2 @ np.array([c["L2"], 0.0])
    b[0:2] = -c["m1"] * g
    b[3:5] = -c["m2"] * g - c["Fext"]
    b[2] = 0.0
    b[5] = -(cross2(tip - com[2], c["Fext"]) + c["Mext"])
    return A, b


def weights(c):
    return np.concatenate([1.0 / c["s"] ** 2, np.full(4, c["wR"])])


def kkt_solve(W, A, b):
    n, m = A.shape[1], A.shape[0]
    K = np.zeros((n + m, n + m), dtype=A.dtype)
    K[:n, :n] = np.diag(W)
    K[:n, n:] = A.T
    K[n:, :n] = A
    rhs = np.concatenate([np.zeros(n, dtype=A.dtype), b])
    z = np.linalg.solve(K, rhs)
    return z[:n], z[n:], K


def solve(c, p, mode="eq"):
    """Return f (10,), lambda (6,), active set (tuple of muscles fixed at 0)."""
    A, b = assemble(c, p)
    W = weights(c)
    if mode == "eq":
        f, l, _ = kkt_solve(W, A, b)
        return f, l, ()
    best = None
    for r in range(NMUS + 1):
        for S in itertools.combinations(range(NMUS), r):
            keep = [j for j in range(NMUS + 4) if j not in S]
            As = A[:, keep]
            if np.linalg.matrix_rank(np.real(As)) < 6:
                continue
            fs, l, _ = kkt_solve(W[keep], As, b)
            f = np.zeros(NMUS + 4, dtype=fs.dtype); f[keep] = fs
            mu = np.real(A[:, list(S)].T @ l) if S else np.zeros(0)   # bound multipliers
            if np.all(np.real(f[:NMUS]) >= -1e-12) and np.all(mu >= -1e-12):
                if best is None:
                    best = (f, l, S, np.min(np.concatenate([np.real(f[[j for j in range(NMUS) if j not in S]]), mu])))
    if best is None:
        raise RuntimeError("no feasible active set")
    return best[0], best[1], best[2]


def solve_on_set(c, p, S):
    """Equality-constrained QP on a fixed active set (complex-safe)."""
    A, b = assemble(c, p)
    W = weights(c)
    keep = [j for j in range(NMUS + 4) if j not in S]
    fs, l, K = kkt_solve(W[keep], A[:, keep], b)
    f = np.zeros(NMUS + 4, dtype=fs.dtype); f[keep] = fs
    return f, l, K, keep


def implicit_jac(c, p, S=()):
    """df/dp via the KKT implicit function theorem. dA, db by complex step (h=1e-30)."""
    f, l, K, keep = solve_on_set(c, p, S)
    n = len(keep)
    J = np.zeros((NMUS + 4, 3))
    for j in range(3):
        pc = p.astype(complex); pc[j] += 1e-30j
        Ac, bc = assemble(c, pc)
        dA = np.imag(Ac[:, keep]) / 1e-30
        db = np.imag(bc) / 1e-30
        rhs = np.concatenate([-dA.T @ l, db - dA @ f[keep]])
        dz = np.linalg.solve(K, rhs)
        J[keep, j] = dz[:n]
    return J, K


def complex_step_jac(c, p, S=()):
    J = np.zeros((NMUS + 4, 3))
    for j in range(3):
        pc = p.astype(complex); pc[j] += 1e-30j
        f, _, _, _ = solve_on_set(c, pc, S)
        J[:, j] = np.imag(f) / 1e-30
    return J


def fd_jacs(c, p, mode):
    """Central FD of the full solver at each step; flags active-set change in stencil."""
    out = []
    _, _, S0 = solve(c, p, mode)
    for h in STEPS:
        J = np.zeros((NMUS + 4, 3)); changed = False; signchg = False
        fmax = 0.0
        for j in range(3):
            e = np.zeros(3); e[j] = h * PSCALE[j]
            fp, _, Sp = solve(c, p + e, mode)
            fm, _, Sm = solve(c, p - e, mode)
            changed |= (Sp != S0) or (Sm != S0)
            signchg |= bool(np.any(np.sign(fp[:NMUS]) != np.sign(fm[:NMUS])))
            fmax = max(fmax, np.max(np.abs(fp)), np.max(np.abs(fm)))
            J[:, j] = (fp - fm) / (2 * e[j])
        out.append(dict(h=h, J=J, active_set_change=bool(changed), muscle_sign_change=signchg, fmax=fmax))
    return out, S0


def relerr(J, Jref):
    """Norm-wise relative error, columns scaled by PSCALE (N per unit-scale step)."""
    D = (J - Jref) * PSCALE
    return float(np.linalg.norm(D) / np.linalg.norm(Jref * PSCALE))


def analyse(c, mode="eq"):
    p = c["p0"].copy()
    A, b = assemble(c, p)
    sv = np.linalg.svd(A, compute_uv=False)
    condA = float(sv[0] / sv[-1]) if sv[-1] > 0 else float("inf")
    rec = dict(seed=c["seed"], kind=c["kind"], mode=mode, m1_kg=c["m1"], m2_kg=c["m2"],
               Fext_N=c["Fext"].tolist(), Mext_Nm=c["Mext"], p0=p.tolist(),
               sigma_min_over_max_A=float(sv[-1] / sv[0]), cond_A=condA)
    flags = []
    if sv[-1] / sv[0] < SVAL_FLAG:
        flags.append("SINGULAR_A")
    elif condA > COND_FLAG:
        flags.append("ILL_CONDITIONED_A")
    if "SINGULAR_A" in flags:
        rec.update(flags=flags, status="FLAGGED_NOT_DIFFERENTIABLE")
        return rec
    try:
        fd, S0 = fd_jacs(c, p, mode)
    except (np.linalg.LinAlgError, RuntimeError) as ex:
        rec.update(flags=flags + ["SOLVE_FAILED:" + str(ex)], status="FLAGGED")
        return rec
    f0, l0, _ = solve(c, p, mode)
    Jimp, K = implicit_jac(c, p, S0)
    Jcs = complex_step_jac(c, p, S0)
    rec["active_set"] = list(S0)
    rec["f0_N"] = f0.tolist()
    rec["cond_K"] = float(np.linalg.cond(K))
    if mode == "nn":
        inact = [j for j in range(NMUS) if j not in S0]
        mu = A[:, list(S0)].T @ l0 if S0 else np.zeros(0)
        margin = float(np.min(np.concatenate([f0[inact], mu])))
        rec["complementarity_margin_N"] = margin
        if margin < MARGIN_FLAG:
            flags.append("WEAK_COMPLEMENTARITY")
    rows = []
    for k, r in enumerate(fd):
        row = dict(h_scaled=r["h"], relerr_vs_implicit=relerr(r["J"], Jimp),
                   relerr_vs_complexstep=relerr(r["J"], Jcs),
                   active_set_change=r["active_set_change"], muscle_sign_change=r["muscle_sign_change"])
        # step-to-step self consistency (does NOT use implicit result)
        if k + 1 < len(fd):
            row["selfdiff_next"] = relerr(r["J"], fd[k + 1]["J"])
        # budget: truncation ~ |D(h)-D(h/10)|; rounding ~ eps*fmax/h relative to |J|
        row["rounding_est"] = float(np.finfo(float).eps * r["fmax"] / r["h"]
                                    / (np.linalg.norm(Jimp * PSCALE) / np.sqrt(Jimp.size)))
        rows.append(row)
        if r["active_set_change"]:
            flags.append(f"ACTIVE_SET_CHANGE_h={r['h']:g}")
    for k, row in enumerate(rows):
        row["trunc_est"] = rows[k].get("selfdiff_next", float("nan"))
    # identified stable step: minimise FD self-difference among steps whose own and
    # next stencil have no active-set change; the chosen step is the finer of the pair
    cand = [k for k in range(len(rows) - 1)
            if not rows[k]["active_set_change"] and not rows[k + 1]["active_set_change"]]
    rec["fd"] = rows
    rec["relerr_implicit_vs_complexstep"] = relerr(Jimp, Jcs)
    rec["jac_implicit"] = Jimp.tolist()
    if not cand:
        flags.append("NO_STABLE_STEP")
        rec.update(flags=flags, status="FLAGGED")
        return rec
    k = min(cand, key=lambda k: rows[k]["selfdiff_next"])
    ks = k + 1
    rec["stable_step"] = rows[ks]["h_scaled"]
    rec["relerr_at_stable_step"] = rows[ks]["relerr_vs_implicit"]
    rec["pass"] = bool(rec["relerr_at_stable_step"] < TOL_REL)
    rec["flags"] = flags
    rec["status"] = "FLAGGED" if flags else "REGULAR"
    return rec


def tune_kink_case(seed):
    """nn mode: move q1 (by bisection) onto a point where the active set switches,
    so the solution map has a kink at p0 (not differentiable; must be flagged)."""
    c = make_case(seed)
    p = c["p0"].copy()
    def aset(q1):
        pp = p.copy(); pp[0] = q1
        try:
            return solve(c, pp, "nn")[2]
        except RuntimeError:
            return None
    S0 = aset(p[0])
    grid = p[0] + np.linspace(-0.5, 0.5, 201)
    sets = [aset(q) for q in grid]
    i0 = int(np.argmin(np.abs(grid - p[0])))
    idx = [i for i in range(len(grid) - 1)
           if sets[i] is not None and sets[i + 1] is not None and sets[i] != sets[i + 1]]
    if S0 is None or not idx:
        raise RuntimeError("no feasible active-set switch within q1 +/- 0.5 rad")
    i = min(idx, key=lambda i: abs(i - i0))
    lo, hi, Slo = grid[i], grid[i + 1], sets[i]
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if aset(mid) == Slo:
            lo = mid
        else:
            hi = mid
    c["p0"][0] = 0.5 * (lo + hi)
    c["kind"] = "active_set_kink"
    return c


def _r(x, n=3):
    return None if x is None else float("%.*e" % (n, x))


def compact(out):
    """Small summary written to results.json (full record via --raw)."""
    def case(r):
        d = dict(seed=r["seed"], kind=r["kind"], mode=r["mode"], status=r["status"],
                 flags=r.get("flags", []), m1_kg=round(r["m1_kg"], 4), m2_kg=round(r["m2_kg"], 4),
                 Fext_N=[round(v, 3) for v in r["Fext_N"]], Mext_Nm=round(r["Mext_Nm"], 4),
                 p0=[round(v, 6) for v in r["p0"]], cond_A=_r(r["cond_A"]))
        if "cond_K" in r:
            d["cond_K"] = _r(r["cond_K"])
        if "complementarity_margin_N" in r:
            d["complementarity_margin_N"] = _r(r["complementarity_margin_N"])
        if "fd" in r:
            d["relerr_fd_vs_implicit_by_step"] = [_r(x["relerr_vs_implicit"]) for x in r["fd"]]
            d["relerr_implicit_vs_complexstep"] = _r(r["relerr_implicit_vs_complexstep"])
        if "stable_step" in r:
            d["stable_step"] = r["stable_step"]
            d["relerr_at_stable_step"] = _r(r["relerr_at_stable_step"])
            d["pass"] = r["pass"]
        return d

    def budget(rs):
        rows = []
        for i, h in enumerate(STEPS):
            e = [r["fd"][i]["relerr_vs_implicit"] for r in rs]
            t = [r["fd"][i]["trunc_est"] for r in rs if i + 1 < len(STEPS)]
            ro = [r["fd"][i]["rounding_est"] for r in rs]
            rows.append(dict(h_scaled=h, relerr_median=_r(float(np.median(e))), relerr_max=_r(max(e)),
                             truncation_est_max=_r(max(t)) if t else None,
                             rounding_est_max=_r(max(ro)),
                             cases_with_muscle_sign_change_in_stencil=sum(r["fd"][i]["muscle_sign_change"] for r in rs)))
        return rows

    reg = out["regular_eq"]
    return dict(
        task="CLOUD-W4-A907-PARITY", data="synthetic only (seeded), no empirical data",
        units=dict(force="N", q1_q2="rad", d="m", jacobian="N/rad, N/m"),
        relerr_definition="||(J_fd-J_imp)*diag(PSCALE)||_F / ||J_imp*diag(PSCALE)||_F, PSCALE=[1 rad,1 rad,0.1 m]",
        steps_scaled=STEPS, tolerance=TOL_REL, criterion=out["criterion"],
        regular_eq_summary=dict(
            relerr_at_stable_step_min=_r(min(r["relerr_at_stable_step"] for r in reg)),
            relerr_at_stable_step_median=_r(float(np.median([r["relerr_at_stable_step"] for r in reg]))),
            relerr_at_stable_step_max=_r(max(r["relerr_at_stable_step"] for r in reg)),
            implicit_vs_complexstep_max=_r(max(r["relerr_implicit_vs_complexstep"] for r in reg)),
            cond_K_range=[_r(min(r["cond_K"] for r in reg)), _r(max(r["cond_K"] for r in reg))],
            cond_A_range=[_r(min(r["cond_A"] for r in reg)), _r(max(r["cond_A"] for r in reg))],
            stable_steps=sorted(set(r["stable_step"] for r in reg))),
        error_budget_regular_eq=budget(reg),
        regular_eq=[case(r) for r in reg],
        regular_nn_exploratory=[dict(seed=r["seed"], status=r["status"], flags=r.get("flags", []),
                                     relerr_at_stable_step=_r(r.get("relerr_at_stable_step")))
                                for r in out["regular_nn"]],
        failure_cases=[case(r) for r in out["failure_cases"]],
        private_A907_reported_value="2.6e-8 (context only, not reproduced)",
        private_A907_parity="UNKNOWN")


def main():
    out = dict(regular_eq=[], regular_nn=[], failure_cases=[])
    for s in range(20):
        c = make_case(s)
        out["regular_eq"].append(analyse(c, "eq"))
        out["regular_nn"].append(analyse(c, "nn"))
    out["failure_cases"].append(analyse(make_case(100, "singular"), "eq"))
    out["failure_cases"].append(analyse(make_case(101, "near_singular"), "eq"))
    # post-hoc harness correction (see RESULTS.md): kink located by bisection in q1
    # on nn-feasible base seeds; first run perturbed Mext on seeds 102/103
    for s in (0, 2, 6, 8):
        try:
            out["failure_cases"].append(analyse(tune_kink_case(s), "nn"))
        except RuntimeError as ex:
            out["failure_cases"].append(dict(seed=s, kind="active_set_kink", status="TUNE_FAILED:" + str(ex)))
    reg = out["regular_eq"]
    fails = out["failure_cases"]
    crit_a = sum(bool(r.get("pass")) for r in reg)
    crit_b = all(r.get("flags") for r in fails)
    out["criterion"] = dict(
        regular_pass_count=crit_a, regular_n=len(reg),
        all_singular_or_active_set_cases_flagged=bool(crit_b),
        nn_mode_pass_count=sum(bool(r.get("pass")) for r in out["regular_nn"]),
        PASS=bool(crit_a == 20 and crit_b),
        private_A907_parity="UNKNOWN")
    if "--raw" in sys.argv:
        json.dump(out, open("raw_results.json", "w"), indent=1)
    json.dump(compact(out), open("results.json", "w"), indent=1)
    return out


if __name__ == "__main__":
    o = main()  # usage: python3 kkt_parity.py [--raw]
    print(json.dumps(o["criterion"], indent=1))
    for grp in ("regular_eq", "regular_nn", "failure_cases"):
        for r in o[grp]:
            print(grp, r["seed"], r.get("kind"), r.get("status"), r.get("flags"),
                  "h*=%s" % r.get("stable_step"), "err=%.2e" % r.get("relerr_at_stable_step", np.nan),
                  "cs=%.1e" % r.get("relerr_implicit_vs_complexstep", np.nan))
