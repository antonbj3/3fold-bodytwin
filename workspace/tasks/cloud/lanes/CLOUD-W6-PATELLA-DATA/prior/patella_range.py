#!/usr/bin/env python3
"""Original synthetic planar knee-extensor countermodel (stdlib only).

Sagittal plane, femur-fixed frame, units mm and degrees; x anterior, y superior,
origin at the centre of a circular posterior femoral condyle (radius R).

Tibia: flexion theta rotates the tibia by -theta about the femoral frame; the
femoral condyle centre, seen from the tibia, moves posteriorly by s*R*theta
(s = rollback fraction; s=1 pure rolling, s=0 fixed hinge).
Tibial tuberosity T: fixed point in the tibia frame.
Patella apex P: point on a trochlear circle (radius Rp, centre (cx, cy)) whose
position is fixed by an inextensible patellar tendon, |P - T| = L0 (set at 0 deg).
Quadriceps: straight line from fixed proximal femoral point Q to P.

Quantities (moment arms positive = extension moment from tendon tension):
  r_geo : signed perpendicular distance from the tibia's instant centre of
          rotation (IC, computed from two poses) to the patellar-tendon line.
  r_vw  : virtual-work moment arm, +dL_PT/dtheta (per rad) with P held fixed (central FD);
          tension resists the tendon lengthening that flexion causes -> extension moment.
  r_eff_vw : effective quadriceps moment arm  +dL_Q/dtheta   (virtual work, whole chain).
  r_eff_eq : r_geo * F_PT/F_Q from tangential static equilibrium of a
             frictionless point patella on the trochlea (independent check).

This is a method check on synthetic geometry, NOT empirical validation.
"""
import json, math, sys

BASE = dict(R=22.0, s=0.5, Tx=32.0, Ty=-40.0, Rp=30.0, cx=8.0, cy=12.0, gamma0=35.0,
            Qx=10.0, Qy=250.0)
TH = [i * 0.5 for i in range(0, 241)]  # 0..120 deg, 0.5 deg grid
H = 1e-4  # deg, finite-difference step


def rot(v, a):
    c, s_ = math.cos(a), math.sin(a)
    return (c * v[0] - s_ * v[1], s_ * v[0] + c * v[1])


def tibia_to_femur(p, th, g):
    """Map tibia-frame point p to femur frame at flexion th (deg)."""
    t = math.radians(th)
    ct = (-g["s"] * g["R"] * t, g["R"])  # femoral centre in tibia frame
    d = (p[0] - ct[0], p[1] - ct[1])
    return rot(d, -t)


def T_at(th, g):
    # tuberosity given relative to the plateau point under the condyle centre at 0 deg
    return tibia_to_femur((g["Tx"], g["Ty"]), th, g)


def trochlea(beta, g):
    b = math.radians(beta)
    return (g["cx"] + g["Rp"] * math.cos(b), g["cy"] + g["Rp"] * math.sin(b))


def dist(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


def L0(g):
    return dist(trochlea(g["gamma0"], g), T_at(0.0, g))


def P_at(th, g, _cache={}):
    """Patella apex: trochlear angle beta with |P-T| = L0, beta <= gamma0 branch.
    Bisection on beta in [gamma0-200, gamma0+5]; the root nearest gamma0 moving distally."""
    T = T_at(th, g)
    L = L0(g)
    f = lambda b: dist(trochlea(b, g), T) - L
    # march distally from gamma0+5 until sign change
    b_hi = g["gamma0"] + 5.0
    fh = f(b_hi)
    b = b_hi
    step = 0.5
    while b > g["gamma0"] - 270:
        b2 = b - step
        f2 = f(b2)
        if (fh > 0) != (f2 > 0):
            lo, hi, flo = b2, b, f2
            for _ in range(80):
                m = 0.5 * (lo + hi)
                fm = f(m)
                if (fm > 0) == (flo > 0):
                    lo, flo = m, fm
                else:
                    hi = m
            return trochlea(0.5 * (lo + hi), g), 0.5 * (lo + hi)
        b, fh = b2, f2
    raise RuntimeError("patella constraint infeasible at theta=%g" % th)


def instant_centre(th, g):
    """IC of tibia w.r.t. femur from velocities of two tibial points (FD)."""
    pts = [(0.0, 0.0), (50.0, -80.0)]
    vs, xs = [], []
    for p in pts:
        a = tibia_to_femur(p, th + H, g)
        b = tibia_to_femur(p, th - H, g)
        vs.append(((a[0] - b[0]) / (2 * H), (a[1] - b[1]) / (2 * H)))
        xs.append(tibia_to_femur(p, th, g))
    w = -math.pi / 180.0  # angular velocity (rad/deg) of tibia in femur frame
    # v = w k x (x - C)  ->  C = x - (k x v)/w ... use v = w*(-(y-Cy), x-Cx)
    x, v = xs[0], vs[0]
    C = (x[0] - v[1] / w, x[1] + v[0] / w)
    # consistency with second point
    x2, v2 = xs[1], vs[1]
    C2 = (x2[0] - v2[1] / w, x2[1] + v2[0] / w)
    return C, dist(C, C2)


def r_geo(th, g):
    C, _ = instant_centre(th, g)
    P, _ = P_at(th, g)
    T = T_at(th, g)
    L = dist(P, T)
    u = ((P[0] - T[0]) / L, (P[1] - T[1]) / L)  # unit, tibia -> patella
    # moment of unit tension on tibia about C (z-component), extension = -z (tibia rotates -theta)
    mz = (T[0] - C[0]) * u[1] - (T[1] - C[1]) * u[0]
    return mz


def r_vw(th, g):
    P, _ = P_at(th, g)  # held fixed
    lp = dist(P, T_at(th + H, g))
    lm = dist(P, T_at(th - H, g))
    return (lp - lm) / (2 * H) * 180.0 / math.pi  # mm per rad


def LQ(th, g):
    return dist((g["Qx"], g["Qy"]), P_at(th, g)[0])


def r_eff_vw(th, g):
    return (LQ(th + H, g) - LQ(th - H, g)) / (2 * H) * 180.0 / math.pi


def r_eff_eq(th, g):
    P, beta = P_at(th, g)
    T = T_at(th, g)
    Q = (g["Qx"], g["Qy"])
    b = math.radians(beta)
    t = (-math.sin(b), math.cos(b))  # trochlear tangent
    uQ = ((Q[0] - P[0]) / dist(Q, P), (Q[1] - P[1]) / dist(Q, P))
    uT = ((T[0] - P[0]) / dist(T, P), (T[1] - P[1]) / dist(T, P))
    ratio = -(uQ[0] * t[0] + uQ[1] * t[1]) / (uT[0] * t[0] + uT[1] * t[1])  # F_PT/F_Q
    return r_geo(th, g) * ratio


def peak(xs, ys):
    i = max(range(len(ys)), key=lambda k: ys[k])
    x = xs[i]
    if 0 < i < len(ys) - 1:  # parabolic refinement
        y0, y1, y2 = ys[i - 1], ys[i], ys[i + 1]
        den = y0 - 2 * y1 + y2
        if den != 0:
            dx = 0.5 * (y0 - y2) / den
            x = xs[i] + dx * (xs[1] - xs[0])
            return x, y1 - 0.25 * (y0 - y2) * dx, "interior"
    return x, ys[i], "endpoint_0" if i == 0 else "endpoint_120"


def curves(g, grid=TH):
    return {k: [f(t, g) for t in grid] for k, f in
            (("r_geo", r_geo), ("r_vw", r_vw), ("r_eff_eq", r_eff_eq), ("r_eff_vw", r_eff_vw))}


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-12)


def compare(g):
    c = curves(g)
    out = {}
    for a, b in (("r_geo", "r_vw"), ("r_eff_eq", "r_eff_vw")):
        pa, pb = peak(TH, c[a]), peak(TH, c[b])
        pw = max(rel(x, y) for x, y in zip(c[a], c[b]) if abs(y) > 1.0)
        out[a + "_vs_" + b] = dict(
            peak_A=dict(theta_deg=round(pa[0], 3), value_mm=round(pa[1], 4), kind=pa[2]),
            peak_B=dict(theta_deg=round(pb[0], 3), value_mm=round(pb[1], 4), kind=pb[2]),
            peak_value_rel_diff=rel(pa[1], pb[1]),
            peak_angle_abs_diff_deg=abs(pa[0] - pb[0]),
            peak_angle_rel_diff_of_range=abs(pa[0] - pb[0]) / 120.0,
            max_pointwise_rel_diff=pw)
    ic_err = max(instant_centre(t, g)[1] for t in TH)
    return c, out, ic_err


def sweep():
    """Sensitivity: rollback s x tuberosity anterior offset Tx; classify r_geo peak."""
    grid = [i * 2.0 for i in range(61)]
    res = []
    for s in (0.0, 0.25, 0.5, 0.75, 1.0):
        for Tx in (15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0):
            g = dict(BASE, s=s, Tx=Tx)
            try:
                ys = [r_geo(t, g) for t in grid]
                pk = peak(grid, ys)
                vw = [r_vw(t, g) for t in grid]
                pv = peak(grid, vw)
                res.append(dict(s=s, Tx_mm=Tx, peak_theta_deg=round(pk[0], 2),
                                peak_mm=round(pk[1], 2), kind=pk[2],
                                margin_mm=round(pk[1] - max(ys[0], ys[-1]), 3),
                                vw_peak_theta_deg=round(pv[0], 2), vw_kind=pv[2],
                                max_rel_diff_geo_vw=max(rel(a, b) for a, b in zip(ys, vw)),
                                r0_mm=round(ys[0], 2), r120_mm=round(ys[-1], 2)))
            except RuntimeError as e:
                res.append(dict(s=s, Tx_mm=Tx, kind="infeasible", error=str(e)))
    return res


def slope0(g, th=0.0, d=0.5):
    return (r_geo(th + d, g) - r_geo(th, g)) / d


def boundary_s():
    """For each Tx, bisect rollback s at which dr/dtheta|0 = 0 (endpoint-0 -> interior)."""
    out = []
    for Tx in (15.0, 25.0, 35.0, 45.0):
        f = lambda s: slope0(dict(BASE, s=s, Tx=Tx))
        lo, hi = 0.0, 1.0
        if (f(lo) > 0) == (f(hi) > 0):
            out.append(dict(Tx_mm=Tx, s_star=None)); continue
        for _ in range(40):
            m = 0.5 * (lo + hi)
            if (f(m) > 0) == (f(lo) > 0):
                lo = m
            else:
                hi = m
        out.append(dict(Tx_mm=Tx, s_star=round(0.5 * (lo + hi), 4)))
    return out


def sweep_gamma():
    """Secondary sweep: patella start angle gamma0 (patella alta/baja proxy)."""
    grid = [i * 2.0 for i in range(61)]
    res = []
    for gm in (10.0, 20.0, 35.0, 50.0, 65.0):
        g = dict(BASE, gamma0=gm)
        try:
            ys = [r_geo(t, g) for t in grid]
            pk = peak(grid, ys)
            res.append(dict(gamma0_deg=gm, peak_theta_deg=round(pk[0], 2),
                            peak_mm=round(pk[1], 2), kind=pk[2]))
        except RuntimeError as e:
            res.append(dict(gamma0_deg=gm, kind="infeasible", error=str(e)))
    return res


def main():
    c, cmp_, ic_err = compare(BASE)
    sw = sweep()
    sg = sweep_gamma()
    bs = boundary_s()
    kinds = sorted({r["kind"] for r in sw + sg})
    feas = [r for r in sw if r["kind"] != "infeasible"]
    sweep_worst = max(r["max_rel_diff_geo_vw"] for r in feas)
    sweep_kind_agree = all(r["kind"] == r["vw_kind"] and abs(r["peak_theta_deg"] - r["vw_peak_theta_deg"]) < 1.2 for r in feas)
    math_ok = all(v["peak_value_rel_diff"] < 0.01 and v["peak_angle_rel_diff_of_range"] < 0.01
                  and v["max_pointwise_rel_diff"] < 0.01 for v in cmp_.values())
    math_ok = math_ok and sweep_worst < 0.01 and sweep_kind_agree
    robust_interior = [r for r in feas if r["kind"] == "interior" and r["margin_mm"] > 0.5]
    sep_ok = bool(robust_interior) and any(r["kind"].startswith("endpoint") for r in feas) \
        and any(b["s_star"] is not None for b in bs)
    table = [dict(theta_deg=t, **{k: round(c[k][i], 3) for k in c})
             for i, t in enumerate(TH) if t % 15 == 0]
    res = dict(
        task="CLOUD-W4-PATELLA-RANGE", units=dict(angle="deg knee flexion", moment_arm="mm"),
        baseline_params=BASE, finite_diff_step_deg=H, grid="0..120 deg step 0.5",
        instant_centre_two_point_max_discrepancy_mm=ic_err,
        baseline_table_every_15deg=table, checks=cmp_, sensitivity_s_Tx=sw,
        sensitivity_gamma0=sg, boundary_s_star_slope0_zero=bs,
        sweep_worst_rel_diff_geo_vw=sweep_worst, sweep_peak_kind_agree_geo_vw=sweep_kind_agree,
        n_robust_interior_margin_gt_0p5mm=len(robust_interior), peak_kinds_observed=kinds,
        criterion_math_extrema_vw_within_1pct=math_ok,
        criterion_sensitivity_separates_endpoint_interior=sep_ok,
        criterion_empirical_peak_reproduction="UNKNOWN",
        empirical_reference=dict(
            citation="Krevolin JL, Pandy MG, Pearce JC. Moment arm of the patellar tendon in the "
                     "human knee. J Biomech 2004;37(5):785-788",
            doi="10.1016/j.jbiomech.2003.09.010",
            direct_access="BLOCKED: doi.org, api.crossref.org, pubmed.ncbi.nlm.nih.gov, europepmc, "
                          "openalex, findanexpert.unimelb.edu.au all refused by egress proxy (HTTP 403)",
            claim_seen="peak patellar-tendon moment arm 4-6 cm across 6 cadaver knees, maximal near "
                       "45 deg flexion; moment arm taken about the finite screw axis of tibia rel. femur",
            claim_provenance="web-search result summary only (search listing incl. "
                             "https://pubmed.ncbi.nlm.nih.gov/15047009/); abstract text NOT fetched verbatim",
            digitized_curves="none available; abstract does not supply curves"))
    out = sys.argv[1] if len(sys.argv) > 1 else "results.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
        f.write("\n")
    print("math PASS" if math_ok else "math FAIL", "| separation PASS" if sep_ok else "| separation FAIL",
          "| empirical UNKNOWN ->", out)


if __name__ == "__main__":
    main()
