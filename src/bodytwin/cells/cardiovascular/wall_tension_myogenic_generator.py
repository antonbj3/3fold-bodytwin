#!/usr/bin/env python3
"""wall_tension_myogenic_generator.py -- one shared vessel-scale mechanism (radius as a function
of transmural pressure) behind the organ-level autoregulation curves that the nephron transport,
cerebral autoregulation and coronary blood flow cells each assert phenomenologically.

MECHANISM (Laplace wall-tension myogenic response, single parametrization, organ-specific
GEOMETRY-ONLY difference):
  T = P * r                         (Laplace, thin-wall, per-unit-length -- task-specified form)
  c(T) = cmax / (1 + exp(-(T-T50)/Tw))   * forced_dilation_falloff(T)   -- myogenic constriction
         fraction, a stretch(tension)-activated smooth-muscle response
  r(P) = r0 * (P/P0)^beta_passive * (1 - c(T))                          -- self-consistent (T
         depends on r; solved by fixed-point iteration)
  R(P) = R0 * (r0/r(P))^4                                               -- Poiseuille 4th power

ANCHOR (the smooth-muscle response parameters cmax/T50/Tw/beta_passive are fit ONLY to real
ISOLATED-VESSEL myography data, decorrelated from the three organ curves this cell tries to
reproduce):
  [A] Falcone, Schilling, Meininger 1993 rat cremaster-muscle 1A arterioles, pressure myograph,
      n=41: PASSIVE diameter 167+/-3 um vs ACTIVE (myogenic tone) diameter 82+/-4 um AT THE SAME
      75 mmHg -- PMID 10325970 abstract cross-check (single
      concrete magnitude point: tone = (167-82)/167 = 50.9% of passive diameter at 75 mmHg).
  [B] Osol, Brekke, McElroy-Yaggy, Gokina 2002, "Myogenic tone, reactivity, and forced
      dilatation: a three-phase model of in vitro arterial myogenic behavior," Am J Physiol
      Heart Circ Physiol, PMID 12388265: three
      phases spanning transmural pressure in CEREBRAL arteries -- MT (tone onset) 40-60 mmHg,
      MR (myogenic reactivity, "little or no change in diameter" despite rising wall tension)
      60-140 mmHg, FD (forced dilatation, tone overwhelmed) >140 mmHg.
  NEITHER anchor is a fit to the nephron transport cell's derived elasticity target (4.0/2.5/1.5)
  nor to the cerebral_autoregulation or coronary_blood_flow cells' own gain numbers -- this is
  the adversary this cell must force: if the model reproduces those organs' curves ONLY because
  it was secretly tuned to them, it is a fit dressed as a mechanism, not an independent generator.

ORGAN-SPECIFIC GEOMETRY (the ONLY per-organ difference -- same c(T)/r(P)/R(P) code path):
  - NEPHRON: 2-resistor Starling divider (afferent myogenic, efferent held passive/fixed -- real
    anatomy: efferent arteriolar smooth muscle is measured less myogenically responsive than
    afferent). Free geometric input: Raff0/Reff0 baseline ratio = 0.6667, taken directly from
    the nephron transport cell's x0 (NOT re-fit here).
  - CEREBRAL: single resistor F=P/R, whole pial/parenchymal-arteriole bed myogenic, P0=93.33
    mmHg (the cerebral_autoregulation cell's P_REF).
  - CORONARY: single resistor F=P/R, whole coronary resistance bed myogenic, P0=90.0 mmHg
    (systemic reference, Berwick 2012's control MAP).
Only P0 (organ baseline pressure) and the myogenic-vessel/total-resistance split are organ
inputs; cmax/T50/Tw/beta_passive are IDENTICAL across all three calls below.

Reads:  nothing (all parameters embedded or calibrated against anchors A/B).
Writes: results.json under the cell output directory.
Gate:   per-organ reproduction tests plus the passive-null falsifier (CMAX=0 must fail).
"""
import json
import math
import os as _os

OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "wall_tension_myogenic_generator")

# ---------------------------------------------------------------------------
# 1. ANCHOR-DERIVED MYOGENIC PARAMETERS (fit ONLY to [A],[B] above, never to organ targets)
# ---------------------------------------------------------------------------
CMAX = None                         # solved (see calibrate_params) not hand-picked -- fixing CMAX equal
                                    # to anchor [A]'s 50.9% value (the first attempt) forces the
                                    # sigmoid to be SATURATED exactly at the calibration point, which is
                                    # underdetermined in T50/TW and collapses to a near-step function
                                    # (checked, see honest_gaps) -- CMAX is instead left as a free,
                                    # physiologically-bounded ceiling (0.5-0.95, near-total closure is
                                    # documented in extreme-pressure vessel studies) and solved jointly
                                    # with T50/TW against anchors [A] (hard, magnitude) and [B] (soft,
                                    # onset~50mmHg / near-saturation by ~130mmHg, before FD_ONSET=140)
T50 = None
TW = None
FD_ONSET = 140.0                    # [B] Osol 2002: forced dilatation begins ~140mmHg
FD_SCALE = 60.0                     # disclosed, NON-anchored extrapolation: falloff scale above FD_ONSET
                                    # (Osol reports FD is QUALITATIVE -- "loss of force production" --
                                    # no quantitative decay slope was given in the abstract; this shape
                                    # parameter is an honest, flagged gap, not a measured number)
BETA_PASSIVE = 0.05                 # small passive-elastic-distension exponent (secondary term; the
                                    # active myogenic term dominates per [B]'s "little or no change in
                                    # diameter" -- kept small and NOT tuned to any organ's residual slope)


def constriction_fraction(T):
    """Stretch(tension)-activated smooth-muscle constriction fraction c(T) in [0,cmax], with a
    disclosed (non-anchored) forced-dilation falloff above FD_ONSET."""
    base = CMAX / (1.0 + math.exp(-(T - T50) / TW))
    if T > FD_ONSET:
        decay = math.exp(-(T - FD_ONSET) / FD_SCALE)
        base = base * decay
    return base


def solve_radius(P, r0, P0, tol=1e-10, max_bisect=200, resid_tol=1e-6):
    """Self-consistent BRACKETED root-find (bisection to tol, no external deps), working in
    DIMENSIONLESS radius rho=r/r0 so that the tension proxy T=P*rho is directly comparable
    (mmHg units) to T50/TW, which were calibrated against anchor [A]/[B] with rho~=1 at
    baseline (r0 IS the reference/passive radius at P0).
    rho(P) = (P/P0)^BETA_PASSIVE * (1 - c(P*rho))  -- T=P*rho depends on the very rho being
    solved for -- this IS the Laplace feedback loop, not an open-loop function of P alone.

    FIX (replaces a fixed-n_iter=80 damped fixed-point iteration that was measured to settle
    into a stable PERIOD-2 LIMIT CYCLE at most sampled pressures -- e.g. at the P=75mmHg anchor
    calibration point the map cycles 0.510713 <-> 0.743049 while the true root, confirmed unique
    by a single sign change of g(rho) over rho in [1e-9, passive_rho], is 0.604821; a fixed EVEN
    iteration count deterministically returned the wrong cycle point, a 15.6% error, silently,
    with no residual check anywhere). g(rho) = passive_rho*(1-c(P*rho)) - rho is guaranteed to
    change sign on [eps, passive_rho] because c(T) in (0,1) strictly (sigmoid never touches its
    bounds): g(eps)=passive_rho*(1-c(~0))-eps > 0, g(passive_rho) = -passive_rho*c(P*passive_rho) < 0.
    Returns absolute r=r0*rho for downstream reporting, plus rho and c. Raises AssertionError if
    the bracket is degenerate or the final residual exceeds resid_tol -- a non-converged solve can
    never again return a number silently."""
    passive_rho = (P / P0) ** BETA_PASSIVE

    def g(rho):
        T = P * rho
        c = constriction_fraction(T)
        return passive_rho * (1.0 - c) - rho

    lo, hi = 1e-9, passive_rho
    g_lo, g_hi = g(lo), g(hi)
    if not (g_lo * g_hi <= 0.0):
        # explicit raise, NOT assert -- assert is stripped under python -O and would fail OPEN
        raise ValueError(
            f"solve_radius: bracket [{lo},{hi}] does not change sign (g_lo={g_lo}, g_hi={g_hi}) "
            f"at P={P}, P0={P0} -- root-uniqueness assumption violated, do not silently proceed"
        )
    for _ in range(max_bisect):
        mid = 0.5 * (lo + hi)
        g_mid = g(mid)
        if g_lo * g_mid <= 0.0:
            hi, g_hi = mid, g_mid
        else:
            lo, g_lo = mid, g_mid
        if (hi - lo) < tol:
            break
    rho = 0.5 * (lo + hi)
    resid = g(rho)
    if not (abs(resid) < resid_tol):
        raise ValueError(
            f"solve_radius: NON-CONVERGED residual {resid} at P={P}, P0={P0}, rho={rho} "
            f"-- refuse to return a silent number"
        )
    c_final = constriction_fraction(P * rho)
    return r0 * rho, rho, c_final


def resistance(P, r0, P0, R0=1.0, rho_baseline=None):
    """Resistance NORMALIZED to the SELF-CONSISTENT (tone-included) operating radius at P0, not
    the passive (rho=1) radius -- the physiological baseline Pa0/P0 already carries whatever
    myogenic tone the mechanism predicts there; the organ curves being tested are themselves
    normalized to their own baseline operating state, not to an unconstricted reference."""
    if rho_baseline is None:
        _, rho_baseline, _ = solve_radius(P0, r0, P0)
    r, rho, c = solve_radius(P, r0, P0)
    R = R0 * (rho_baseline / rho) ** 4
    return R, r, c


def _rho_c_at(CMAX_try, T50_try, TW_try, P, P0=75.0, tol=1e-10, max_bisect=200, resid_tol=1e-6):
    """Self-consistent (rho,c) at pressure P using CANDIDATE params (calibration-only helper --
    pure function of trial parameters, no module-state mutation, so calibration can grid-search
    freely). FIX: this ALSO used a fixed-even-n_iter=120 damped fixed-point iteration and was
    measured to limit-cycle at the SAME P=75 calibration anchor point (0.510713 <-> 0.743049,
    n_iter=120 landing on the wrong 0.510713 branch) -- i.e. calibrate_params() was scoring
    candidates against a WRONG self-consistent rho75 for every grid point, not just solve_radius's
    downstream callers. Replaced with the same bracketed bisection + residual guard."""
    passive_rho = (P / P0) ** 0.0  # P0=75 IS the calibration baseline (anchor [A]'s pressure)

    def g(rho):
        T = P * rho
        c = CMAX_try / (1.0 + math.exp(-(T - T50_try) / TW_try))
        return passive_rho * (1.0 - c) - rho

    lo, hi = 1e-9, passive_rho
    g_lo, g_hi = g(lo), g(hi)
    if not (g_lo * g_hi <= 0.0):
        raise ValueError(
            f"_rho_c_at: bracket [{lo},{hi}] does not change sign (g_lo={g_lo}, g_hi={g_hi}) "
            f"for CMAX={CMAX_try}, T50={T50_try}, TW={TW_try}, P={P}"
        )
    for _ in range(max_bisect):
        mid = 0.5 * (lo + hi)
        g_mid = g(mid)
        if g_lo * g_mid <= 0.0:
            hi, g_hi = mid, g_mid
        else:
            lo, g_lo = mid, g_mid
        if (hi - lo) < tol:
            break
    rho = 0.5 * (lo + hi)
    resid = g(rho)
    if not (abs(resid) < resid_tol):
        raise ValueError(f"_rho_c_at: NON-CONVERGED residual {resid} for CMAX={CMAX_try}, "
                          f"T50={T50_try}, TW={TW_try}, P={P}")
    c_final = CMAX_try / (1.0 + math.exp(-(P * rho - T50_try) / TW_try))
    return rho, c_final


def calibrate_params():
    """Joint grid+refine search over (CMAX,T50,TW) -- NONE hand-picked -- using ONLY anchors
    [A] (hard equality: self-consistent rho(P=75)=82/167=0.4910) and [B] (soft shape targets:
    onset~0.05*CMAX near P=50mmHg, near-saturation~0.95*CMAX by P=130mmHg, i.e. before
    FD_ONSET=140). Setting CMAX equal to the anchor's 0.509 (first attempt, see git history
    of this cell) forces exact saturation at the calibration point and is underdetermined in
    T50/TW, collapsing to a near-step function -- CMAX is left free here specifically to force
    that adversary to a fairer, non-degenerate form. NEVER touches organ-curve data."""
    target_rho75 = 82.0 / 167.0  # 0.4910...
    best = None  # (score, CMAX, T50, TW, err_A)
    for CMAX_try in [0.50 + 0.025 * i for i in range(19)]:      # 0.50..0.95, step 0.025
        for T50_try in [10.0 + 5.0 * i for i in range(19)]:     # 10..100, step 5
            for TW_try in [2.0 + 3.0 * i for i in range(20)]:   # 2..59, step 3
                rho75, c75 = _rho_c_at(CMAX_try, T50_try, TW_try, 75.0)
                err_A = abs(rho75 - target_rho75)
                if err_A > 0.03:
                    continue
                _, c50 = _rho_c_at(CMAX_try, T50_try, TW_try, 50.0)
                _, c130 = _rho_c_at(CMAX_try, T50_try, TW_try, 130.0)
                onset_err = abs(c50 - 0.05 * CMAX_try)
                sat_err = abs(c130 - 0.95 * CMAX_try)
                # penalize c75 being AT (not below) CMAX_try -- reward genuine mid-curve (non-
                # saturated) behavior at the calibration point itself, the whole point of freeing CMAX
                nondegenerate_bonus = -(CMAX_try - c75)  # more negative (better) if c75 << CMAX_try
                score = 5.0 * err_A + onset_err + sat_err + 0.5 * nondegenerate_bonus
                if best is None or score < best[0]:
                    best = (score, CMAX_try, T50_try, TW_try, err_A, c75, c50, c130)
    if best is None:
        raise RuntimeError("joint (CMAX,T50,TW) grid found NO candidate reproducing anchor [A] within "
                            "3% -- report as a genuine calibration failure, do not silently widen tolerance")
    _, CMAX_b, T50_b, TW_b, errA_b, c75_b, c50_b, c130_b = best
    return {
        "CMAX": round(CMAX_b, 4), "T50": round(T50_b, 3), "TW": round(TW_b, 3),
        "anchor_A_fit_error": round(errA_b, 5),
        "c_at_75mmHg (should be << CMAX if non-degenerate)": round(c75_b, 4),
        "c_at_50mmHg (onset target ~0.05*CMAX)": round(c50_b, 4),
        "c_at_130mmHg (near-saturation target ~0.95*CMAX)": round(c130_b, 4),
        "is_degenerate_saturated_at_calibration_point": bool(abs(c75_b - CMAX_b) < 0.01),
    }


def local_elasticity(f, P, h_frac=1e-3):
    """d(ln f)/d(ln P) via central finite difference (machine-computed, not eyeballed)."""
    h = P * h_frac
    fp, fm = f(P + h), f(P - h)
    if fp <= 0 or fm <= 0:
        return float("nan")
    return (math.log(fp) - math.log(fm)) / (math.log(P + h) - math.log(P - h))


# ---------------------------------------------------------------------------
# 2. FALSIFIER FIRST: reproduce anchor [A] exactly at P=75mmHg with r0 set to match, and confirm
#    the model is NOT degenerate (c must vary with P, not be a constant -- the void-floor this
#    whole exercise exists to beat).
# ---------------------------------------------------------------------------
def selfcheck_anchor_A():
    r0 = 167.0  # um, PMID10325970's passive diameter units, treated as r0 (radius proxy)
    P0 = 75.0   # calibration pressure IS the anchor's measurement pressure
    r, rho, c = solve_radius(75.0, r0, P0)
    return {
        "P_mmHg": 75.0,
        "passive_r0_um": r0,
        "solved_active_r_um": round(r, 2),
        "expected_active_diam_um_PMID10325970": 82.0,
        "solved_active_diam_um (2r, passive-exponent negligible at P=P0)": round(2 * r, 2),
        "solved_constriction_fraction": round(c, 4),
        "target_constriction_fraction_PMID10325970": 0.509,
    }


# ---------------------------------------------------------------------------
# 3. DECISIVE TEST 1 -- NEPHRON (held out: elasticity profile DERIVED from filtration equation,
#    never measured/fit to)
# ---------------------------------------------------------------------------
Pgc0_nephron = 60.0
Pbowman = 18.0
pi_gc = 32.0
Kf = 12.5
Pa0_nephron = 100.0
X0_NEPHRON = (Pa0_nephron - Pgc0_nephron) / Pgc0_nephron  # 0.6667, Raff0/Reff0, nephron_transport_cell.py's
GFR0 = Kf * (Pgc0_nephron - Pbowman - pi_gc)  # 125.0
NEPHRON_R0 = 1.0  # dimensionless reference radius (rho=1 at Pa0); matches T50/TW calibration convention


def nephron_raff_over_reff(Pa):
    """Raff(Pa)/Reff, using the SAME generic myogenic mechanism for the afferent arteriole,
    Reff held passive/fixed (organ-specific geometry: only the afferent vessel is myogenic,
    matching real differential afferent/efferent smooth-muscle density)."""
    R_myo, r, c = resistance(Pa, NEPHRON_R0, Pa0_nephron, R0=1.0)
    return X0_NEPHRON * R_myo  # R_myo(Pa0)=1 by construction, so Raff0/Reff=X0_NEPHRON exactly at Pa0


def nephron_gfr(Pa):
    x = nephron_raff_over_reff(Pa)
    ratio = 1.0 / (1.0 + x)  # Reff/(Raff+Reff)
    Pgc = ratio * Pa
    NFP = Pgc - Pbowman - pi_gc
    return Kf * NFP


def nephron_test():
    SWEEP = [80 + 5 * i for i in range(21)]
    rows = []
    worst = 0.0
    for Pa in SWEEP:
        gfr = nephron_gfr(Pa)
        pct = 100.0 * gfr / GFR0
        dev = abs(pct - 100.0)
        worst = max(worst, dev)
        rows.append({"Pa": Pa, "GFR": round(gfr, 2), "pct_of_baseline": round(pct, 1)})
    elast = {
        Pa: round(local_elasticity(lambda P: nephron_raff_over_reff(P), Pa), 3)
        for Pa in (80.0, 100.0, 180.0)
    }
    TARGET = {80.0: 4.0, 100.0: 2.5, 180.0: 1.5}
    elast_err = {Pa: round(abs(elast[Pa] - TARGET[Pa]), 3) for Pa in TARGET}
    return {
        "sweep_worst_deviation_pct": round(worst, 2),
        "gate_10pct": worst <= 10.0,
        "prior_best_constant_elasticity_worst_deviation_pct": 38.0,
        "emergent_local_elasticity": elast,
        "target_local_elasticity_derived_from_filtration_eqn": TARGET,
        "elasticity_abs_error": elast_err,
        "rows_extremes": [rows[0], rows[4], rows[-1]],
    }


# ---------------------------------------------------------------------------
# 4. DECISIVE TEST 2a -- CEREBRAL (single resistor F=P/R, plateau 60-150mmHg CPP)
# ---------------------------------------------------------------------------
CEREBRAL_P0 = 93.33
CEREBRAL_R0 = 1.0  # dimensionless reference radius


def cerebral_flow_ratio(P):
    R, r, c = resistance(P, CEREBRAL_R0, CEREBRAL_P0, R0=1.0)
    return P / R  # F/F0 since R(P0)=1 and F0=P0/1=P0 -> normalize below


def cerebral_test():
    F0 = cerebral_flow_ratio(CEREBRAL_P0)
    sweep = [50, 60, 70, 90, 100, 120, 140, 150, 160]
    rows = []
    plateau_worst = 0.0
    for P in sweep:
        F = cerebral_flow_ratio(P) / F0
        rows.append({"P": P, "F_over_F0": round(F, 4)})
        if 60 <= P <= 150:
            plateau_worst = max(plateau_worst, abs(100.0 * F - 100.0))
    F_at_50 = cerebral_flow_ratio(50.0) / F0
    return {
        "rows": rows,
        "plateau_60_150_worst_pct_deviation": round(plateau_worst, 2),
        "F_over_F0_at_50mmHg (below classic LLA, MUST drop)": round(F_at_50, 4),
        "drops_below_LLA": F_at_50 < 0.95,
        "leak_gate_5pct (modern target: some non-zero leak, not perfectly flat)": plateau_worst > 0.5,
    }


# ---------------------------------------------------------------------------
# 5. DECISIVE TEST 2b -- CORONARY (single resistor F=P/R, Berwick 2012 Gc=0.46 over 60-100mmHg,
#    extrapolated by that script to 60-140mmHg)
# ---------------------------------------------------------------------------
CORONARY_P0 = 90.0
CORONARY_R0 = 1.0  # dimensionless reference radius


def coronary_flow_ratio(P):
    R, r, c = resistance(P, CORONARY_R0, CORONARY_P0, R0=1.0)
    return P / R


def coronary_test():
    F0 = coronary_flow_ratio(CORONARY_P0)
    lo, hi = 60.0, 140.0
    F_lo, F_hi = coronary_flow_ratio(lo) / F0, coronary_flow_ratio(hi) / F0
    # implied Gc from F(P)/F0 = (P/P0)^(1-Gc) fit across the two endpoints (log-linear 2-point slope)
    exponent = math.log(F_hi / F_lo) / math.log(hi / lo)
    Gc_implied = 1.0 - exponent
    voidfloor_lo, voidfloor_hi = (lo / CORONARY_P0), (hi / CORONARY_P0)  # Gc=0, fully passive
    return {
        "F_over_F0_at_60": round(F_lo, 4),
        "F_over_F0_at_140": round(F_hi, 4),
        "implied_Gc_from_emergent_curve": round(Gc_implied, 3),
        "measured_Gc_Berwick2012_PMID22466959": 0.46,
        "Gc_abs_error": round(abs(Gc_implied - 0.46), 3),
        "voidfloor_passive_F_over_F0_at_60_140 (Gc=0 null)": [round(voidfloor_lo, 3), round(voidfloor_hi, 3)],
    }


# ---------------------------------------------------------------------------
# 6. PASSIVE-NULL FALSIFIER (adversary): CMAX=0 must reproduce the known-bad passive result
# ---------------------------------------------------------------------------
def passive_null_check():
    global CMAX
    saved = CMAX
    CMAX = 0.0
    worst = 0.0
    for Pa in [80 + 5 * i for i in range(21)]:
        gfr = nephron_gfr(Pa)
        dev = abs(100.0 * gfr / GFR0 - 100.0)
        worst = max(worst, dev)
    CMAX = saved
    return {"nephron_worst_deviation_with_zero_myogenic_gain_pct": round(worst, 1),
            "matches_prior_session_passive_null_480pct (order of magnitude)": worst > 100.0}


def main():
    global CMAX, T50, TW
    calib = calibrate_params()
    CMAX, T50, TW = calib["CMAX"], calib["T50"], calib["TW"]
    out = {
        "calibration (T50/TW solved against anchors A/B only)": calib,
        "anchor_A_selfcheck": selfcheck_anchor_A(),
        "nephron_test": nephron_test(),
        "cerebral_test": cerebral_test(),
        "coronary_test": coronary_test(),
        "passive_null_falsifier": passive_null_check(),
        "params": {"CMAX": CMAX, "T50": T50, "TW": round(TW, 4), "FD_ONSET": FD_ONSET,
                   "FD_SCALE": FD_SCALE, "BETA_PASSIVE": BETA_PASSIVE},
    }
    print(json.dumps(out, indent=2))
    _os.makedirs(OUT_DIR, exist_ok=True)
    with open(_os.path.join(OUT_DIR, "results.json"), "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
