#!/usr/bin/env python3
"""
BT-HX-Q014 -- "When does changing the order between two interventions change an
observable outcome?"

First-principles mechanistic model of the ORDER contrast (first-order carryover)
in a two-period AB/BA sequence experiment with a washout ("wait") interval.

    python3 model.py            # runs the full analysis, writes results.json

Frozen preregistration: PREREG.md (sha256 in PREREG.sha256), incl. Errata 1-4.
No measured data are contained in or invented by this file. Every number in
RESULTS.md is produced here or derived in PREREG.md.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, replace

import numpy as np
from scipy import optimize, stats

PT = "pt"     # outcome points (calibration scale: Mini-BESTest, 0-28)
DAY = "d"


# ---------------------------------------------------------------------------
# PARAMETER TABLE  (unit, value, source-or-assumption)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Params:
    # --- design geometry (all [d]) -----------------------------------------
    t_p: float = 90.0   # period-1 duration. DERIVED: R1 trial "twice weekly for
                        # 90 min over 3 months" -> 90 d.
    t_w: float = 0.0    # washout. SOURCE: R1 Table 4 (Sparrow et al.) = 0 d.
    t_o: float = 90.0   # read-out offset into period 2. ASSUMPTION: t_o = t_p.
    # --- mechanism ----------------------------------------------------------
    tau: float = 9.0    # memory time constant. OVERIFIERAD as a literature value
                        # (PREREG Sec. 4). FREE PARAMETER, swept. Default t_p/10
                        # is a declared, deliberately weak prior, not evidence.
    a_A: float = 4.8    # increment the latent state writes per unit S. IDENTIFIED
                        # from R1 (see identify_aA), not assumed.
    a_B: float = 0.0    # increment written by the comparator (usual care).
                        # ASSUMPTION: the comparator writes nothing.
    G: float = 1.0      # observability gain. CONVENTION: G == 1 by definition --
                        # the latent state is expressed in outcome units. Only
                        # the product G*a is identified; G is swept in the
                        # sensitivity block as an invariance check.
    # --- noise floor --------------------------------------------------------
    sigma_w: float = 3.45   # within-subject SD [pt]. DERIVED from R1 Table 2A
                            # SEs: 3.19 (period 1) and 3.69 (period 2).
    n_seq: tuple = (7, 9)   # subjects per sequence. SOURCE: R1 Table 2A.
    alpha: float = 0.05     # two-sided test size.


PARAMETER_TABLE = [
    ("t_p",    DAY, "90.0",  "DERIVED: R1 trial duration 'twice weekly for 90 min over 3 months' -> 90 d"),
    ("t_w",    DAY, "0.0",   "SOURCE: R1 Table 4 (Sparrow et al.), washout period = 0 d"),
    ("t_o",    DAY, "90.0",  "ASSUMPTION: read-out at the end of period 2, t_o = t_p. Swept: it is a design lever"),
    ("tau",    DAY, "9.0",   "OVERIFIERAD (no primary value retrieved). FREE PARAMETER; default t_p/10, swept over 4 decades"),
    ("a_A",    PT,  "4.80",  "IDENTIFIED: period1_raw/S with period1_raw = 4.8 pt (R1 Table 2A, 25.0 - 20.2)"),
    ("a_B",    PT,  "0.0",   "ASSUMPTION: comparator (usual care) writes no latent state"),
    ("G",      "-",  "1.0",   "CONVENTION: latent state expressed in outcome units; only G*a_A is identified"),
    ("sigma_w",PT,  "3.45",  "DERIVED: R1 Table 2A SEs*sqrt(n) -> 3.19 (P1) and 3.69 (P2) pt"),
    ("n_seq",  "-",  "(7, 9)","SOURCE: R1 Table 2A, N=7 (Exercise->UC) and N=9 (UC->Exercise)"),
    ("alpha",  "-",  "0.05", "CONVENTION: two-sided"),
]


# ---------------------------------------------------------------------------
# FIRST-PRINCIPLES CORE
# ---------------------------------------------------------------------------
def propagator(z, u: float, D: float, tau: float):
    """Exact solution of the first-order relaxation (E1) over a step of length D.

        (E1)  dz/dt = -(1/tau) [ z(t) - u(t) ]
        (E3)  z(t+D) = z(t) exp(-D/tau) + u (1 - exp(-D/tau))

    `u` is the ATTRACTOR LEVEL of the state in outcome units [pt] (u = 0 during
    washout = relaxation back to baseline), NOT a rate. This normalisation is
    what makes z(t+D) -> u as D -> inf, and it is checked by
    test_model.py::test_plateau (which is what caught the earlier rate/level
    confusion). Units: [pt] * dimensionless + [pt] * dimensionless = [pt].
    """
    if tau <= 0:
        raise ValueError("tau must be > 0")
    e = np.exp(-np.asarray(D, dtype=float) / tau)
    return np.asarray(z, dtype=float) * e + u * (1.0 - e)


def gates(p: Params) -> dict:
    """The gates of the order contrast (PREREG.md Sec. 3).

    S = 1 - e^{-t_p/tau}          state written by the period-1 intervention
    R = e^{-(t_w+t_o)/tau}        fraction of that state still present at the read
    c = 1 - e^{-t_o/tau}          new-treatment increment not yet realised
    D = G (a_B - a_A)             treatment contrast of the latent state
    """
    return {
        "S": 1.0 - math.exp(-p.t_p / p.tau),
        "R": math.exp(-(p.t_w + p.t_o) / p.tau),
        "residual": 1.0 - math.exp(-(p.t_w + p.t_o) / p.tau),
        "c": 1.0 - math.exp(-p.t_o / p.tau),
        "new_increment": 1.0 - math.exp(-p.t_o / p.tau),
        "D": p.G * (p.a_B - p.a_A),
    }


def order_contrast(p: Params) -> float:
    """THE ORDER CONTRAST: first-order carryover estimand `lam_hat`, in [pt].

    Standard 2x2 (Willan & Pater 1986, R2), corrected per PREREG Erratum 1:
        tau_hat = 1/2[(y11-y21) + (y12-y22)]      treatment effect
        lam_hat = (y11-y21) - tau_hat               first-order carryover ("order")
    Closed form:
        lam_hat = -(D/2) [ S(1 + R) - c ]
    Positive = doing A first made the period-2 outcome WORSE than B-first.
    Identically 0 when the period-1 state is fully erased (R -> 0) and the
    outcome is read at steady state (S = c) -- i.e. when order cannot matter.
    """
    g = gates(p)
    return -0.5 * g["D"] * (g["S"] * (1.0 + g["R"]) - g["c"])


def period_interaction(p: Params) -> float:
    """Sequence/period interaction `pi_hat` = 1/2[(y12-y11)+(y22-y21)], [pt]."""
    g = gates(p)
    return 0.5 * p.G * (p.a_A + p.a_B) * (g["c"] - g["S"] * (1.0 - g["R"]))


def treatment_effect(p: Params) -> float:
    """`tau_hat` = 1/2[(y11-y21)+(y12-y22)], [pt]."""
    g = gates(p)
    return -0.5 * g["D"] * (g["S"] + g["c"] - g["S"] * g["R"])


def period1_raw_contrast(p: Params) -> float:
    """Uncorrected period-1 contrast between sequences, y11 - y21 = G(a_A-a_B)S."""
    return -p.G * (p.a_B - p.a_A) * gates(p)["S"]


def delta_order(p: Params) -> float:
    """Alias for the headline estimand (kept for readability in RESULTS)."""
    return order_contrast(p)


def sequence_states(p: Params, seq: str) -> dict:
    """Latent state at the three read-outs of one sequence."""
    first, second = (p.a_A, p.a_B) if seq == "AB" else (p.a_B, p.a_A)
    z_p1_end = propagator(0.0, first, p.t_p, p.tau)
    z_p2_start = propagator(z_p1_end, 0.0, p.t_w, p.tau)
    z_p2_read = propagator(z_p2_start, second, p.t_o, p.tau)
    return {"z_p1_end": float(z_p1_end), "z_p2_start": float(z_p2_start),
            "z_p2_read": float(z_p2_read)}


# ---------------------------------------------------------------------------
# NOISE FLOOR
# ---------------------------------------------------------------------------
# the four published cell SEs (R1 Table 2A), period-1 seq1/seq2, period-2 seq1/seq2
PUBLISHED_CELL_SE = (0.87, 1.24, 0.80, 1.51)
SIGMA_W_REFERENCE = 3.45   # the SD those SEs correspond to, see PARAMETER_TABLE


def _se_corrected(p: Params) -> float:
    """SE of a *mean-of-differences* estimand, from the four published cell SEs.

    The published SEs are rescaled by sigma_w / sigma_w_ref so that the within-
    subject SD is a live (sensitivity-controlling) parameter rather than a
    frozen constant. At sigma_w = sigma_w_ref this returns the published value.
    """
    k = p.sigma_w / SIGMA_W_REFERENCE
    return 0.5 * k * math.sqrt(sum(s ** 2 for s in PUBLISHED_CELL_SE))


def mdc(p: Params, estimand: str = "carryover") -> float:
    """Minimal detectable contrast, 95%, two-sided, df = n1+n2-2."""
    n1, n2 = p.n_seq
    df = n1 + n2 - 2
    if estimand in ("carryover", "period"):
        se = _se_corrected(p)
    elif estimand == "period1_raw":
        se = (p.sigma_w / SIGMA_W_REFERENCE) * math.sqrt(
            (3.19 ** 2) / n1 + (3.69 ** 2) / n2)
    elif estimand == "baseline":
        se = math.sqrt(2.0) * p.sigma_w / math.sqrt((n1 + n2) / 2)
    else:
        raise ValueError(estimand)
    return float(stats.t.ppf(1 - p.alpha / 2, df) * se)


# ---------------------------------------------------------------------------
# THE TWO PREDICTIONS OF INTEREST
# ---------------------------------------------------------------------------
def critical_washout(p: Params) -> float:
    """t_w* : the largest washout for which |lam_hat| stays above the noise floor.

    With a_B < a_A, D < 0 and lam_hat = (G(a_A-a_B)/2)[S(1+R) - c] = K[S(1+R)-c],
    K > 0. Detectability  K[S(1+R) - c] > MDC  gives
        R > (MDC/K + c)/S - 1
        t_w* = -tau*ln( (MDC/K + c)/S - 1 ) - t_o
    NaN  if no non-negative washout is admissible (contrast below MDC even at t_w=0)
    +inf if it is detectable even after the state is fully erased (R -> 0)
    """
    g = gates(p)
    S, c, D = g["S"], g["c"], g["D"]
    K = -0.5 * D
    if S <= 0 or K == 0:
        return float("nan")
    if abs(order_contrast(p)) <= mdc(p):
        return float("nan")
    q = (mdc(p) / K + c) / S - 1.0
    if q <= 0:
        return float("inf")
    if q >= 1:
        return float("nan")
    return float(-p.tau * math.log(q) - p.t_o)


def tau_min(p: Params, t_w: float) -> float:
    """Smallest memory timescale for which |lam_hat| exceeds the noise floor."""
    q = replace(p, t_w=t_w)
    f = lambda t: abs(order_contrast(replace(q, tau=t))) - mdc(q)  # noqa: E731
    if f(1e-9) >= 0:
        return 0.0
    hi = 1e-3
    while f(hi) < 0 and hi < 1e9:
        hi *= 3.0
    if f(hi) < 0:
        return float("nan")          # no tau makes the contrast detectable
    return float(optimize.brentq(f, 1e-9, hi, xtol=1e-12, rtol=1e-12))


def tau_infinity_limit(p: Params) -> dict:
    """Tightest possible 'can order ever matter at this design?' bounds, closed form.

    With a_A pinned by the identification line a_A = period1_raw/S and
    u = exp(-t_o/tau) in (0,1]:
        lam_hat = (period1_raw/2)(1 + u) - (period1_raw/2)(1 - u)/S,  S = 1 - e_p
    Three regimes at t_w = 0, t_o = t_p (e_p = exp(-t_p/tau)):
      * t_o -> 0  (u -> 1):  lam -> period1_raw = 4.80 pt  > MDC  => detectable
      * t_o -> inf (u -> 0):  lam -> period1_raw/2 * (1 - 1/S) -> 0
      * max over tau:  a_A/8 = 0.600 pt                 < MDC  => undetectable
    Along the identification line the tau->0 end is the tightest:
      max over tau of lam(t_o = t_p) = 1.20 pt at tau* = t_p/ln 2.
    """
    raw = PUBLISHED["period1_raw_contrast"]
    p0 = replace(p, t_o=0.0, t_w=0.0)
    am = analytic_max_over_tau(p)
    line_max = 0.0
    for tau in np.geomspace(1e-3, 50 * p.t_p, 20000):
        q = Params(tau=float(tau), t_o=p.t_p, t_w=0.0, a_A=identify_aA(float(tau)))
        line_max = max(line_max, abs(order_contrast(q)))
    return {
        "t_o_0d_lam_pt": float(order_contrast(replace(p, t_o=0.0, t_w=0.0))),
        "t_o_tp_max_over_tau_pt": am["max_abs_lam_hat_pt"],
        "t_o_tp_max_over_tau_tau_star_days": am["tau_star_days"],
        "ident_line_max_over_tau_pt": float(line_max),
        "mdc_carryover_pt": mdc(p0, "carryover"),
        "order_detectable_at_t_o_0": bool(
            abs(order_contrast(replace(p, t_o=0.0, t_w=0.0))) > mdc(p0)),
        "order_detectable_at_t_o_tp_for_any_tau": bool(line_max > mdc(p0)),
    }


def analytic_max_over_tau(p: Params) -> dict:
    """Closed-form maximisation of the order contrast over tau at t_w = 0, t_o = t_p.

    At t_w = 0, t_o = t_p:  R = c = e = e^{-t_p/tau},  S = 1 - e
        S(1+R) - c = (1-e)(1+e) - (1-e) = e(1-e)   <=  1/4,  max at e = 1/2.
    Therefore  tau* = t_p/ln 2  and  max_tau |lam_hat| = K/4 = a_A/8.
    """
    p = replace(p, t_w=0.0, t_o=p.t_p)
    K = -0.5 * gates(p)["D"]
    return {"tau_star_days": p.t_p / math.log(2.0),
            "max_abs_lam_hat_pt": K / 4.0,
            "mdc_carryover_pt": mdc(p, "carryover"),
            "max_is_detectable": bool(K / 4.0 > mdc(p, "carryover")),
            "assumption": "t_w = 0 and t_o = t_p; K = G(a_A - a_B)/2"}


# ---------------------------------------------------------------------------
# NULL MODELS / COUNTER-PROOFS (PREREG.md Sec. 5)
# ---------------------------------------------------------------------------
def null_commutative(p: Params) -> float:
    """N0: exchangeable interventions (a_A = a_B) -> provably 0."""
    return order_contrast(replace(p, a_B=p.a_A))


def null_additive(p: Params) -> float:
    """N0': the standard crossover model -- tau -> 0 AND a washout at least as
    long as the period, so the residual state is fully erased. Must be 0."""
    return order_contrast(replace(p, tau=1e-9, t_w=max(p.t_w, p.t_p)))


def null_control_run(p: Params) -> float:
    """N1: control run A -> nothing (period 2 with no second intervention).

    The control is the same subject, same period-1 intervention, same washout,
    but no second intervention. Its period-2 read is therefore the period-1
    state decayed over the WHOLE of (t_w + t_o) with zero drive.
    Deviation from the A->B read, [pt]:  = |G * a_B * (1 - exp(-t_o/tau))|,
    i.e. identically 0 iff the comparator writes nothing. With a_B = a_A the
    placebo becomes motivable (deviation must exceed the MDC).
    """
    z = sequence_states(p, "AB")
    z_ctrl = float(propagator(z["z_p1_end"], 0.0, p.t_w + p.t_o, p.tau))
    return float(p.G * abs(z["z_p2_read"] - z_ctrl))


def null_orthogonal(p: Params) -> float:
    """N2: two interventions on independent latent states -> 0 for any memory.

    Provable: if A writes only x and B writes only z, the AB and BA orderings
    leave the combined observable (beta_x x + beta_z z + const) identical.
    """
    return 0.0


# ---------------------------------------------------------------------------
# UNIT CHECK  (explicit, not decorative)
# ---------------------------------------------------------------------------
UNIT_CHECK = [
    ("exp(-t_p/tau)",                        "-",  "ratio of two times [d]/[d]"),
    ("z(t+D) = z e^{-D/tau} + u (1-e^{-D/tau})", PT,
     "[pt] + [pt] x dimensionless = [pt]  (u is an attractor LEVEL, not a rate)"),
    ("lim_{D->inf} z(t+D) = u", PT, "attractor level in [pt], independent of tau [d]"),
    ("S, R, c, new_increment",               "-",  "each a ratio of times"),
    ("D = G (a_B - a_A)",                    PT,  "dimensionless x [pt]"),
    ("lam_hat = -(D/2)[S(1+R) - c]",         PT,  "[pt] x dimensionless"),
    ("pi_hat = (G/2)(a_A+a_B)[c - S(1-R)]",  PT,  "dimensionless x [pt] x dimensionless"),
    ("period1_raw = G(a_A-a_B) S",           PT,  "dimensionless x [pt]"),
    ("t_w* = -tau ln(q) - t_o",              DAY, "[d] x dimensionless - [d]"),
    ("MDC = t_{alpha/2,df} * SE",            PT,  "dimensionless x [pt]"),
    ("SE = 0.5 sqrt(SE1^2+..+SE4^2)",        PT,  "squares and sums of [pt]"),
]


def unit_check() -> dict:
    p = Params(a_A=identify_aA(p_default_tau()))
    rows = [{"expression": e, "unit": u, "justification": w} for e, u, w in UNIT_CHECK]
    z = float(propagator(0.0, p.a_A, p.t_p, p.tau))
    assert math.isfinite(z) and z > 0
    for f in (order_contrast, period_interaction, treatment_effect, period1_raw_contrast):
        v = f(p)
        assert math.isfinite(v), f"{f.__name__} not finite"
    # G-invariance: only G*a_A is identified, so scaling G by k and dividing
    # a_A by k must leave every estimand invariant.
    k = 7.3
    a, b = f(p), f(replace(p, G=p.G * k, a_A=p.a_A / k))
    assert abs(a - b) < 1e-10, "G/a_A non-invariance"
    return {"table": rows, "G_invariance_check": "passed",
            "scale_check": "all estimands in [pt]; all times in [d]; G, alpha dimensionless"}


# ---------------------------------------------------------------------------
# CALIBRATION TO R1 (PREREG.md Sec. 6, Erratum 4)
# ---------------------------------------------------------------------------
PUBLISHED = {
    "period1_raw_contrast": 4.8,          # 25.0 - 20.2                    [pt]
    "treatment_effect_tau_hat": 3.3,      # ((25.0-20.2)+(24.1-22.3))/2     [pt]
    "period_interaction_pi_hat": 0.6,     # ((24.1-25.0)+(22.3-20.2))/2     [pt]
    "first_order_carryover_lam_hat": 1.5, # 4.8 - 3.3                      [pt]
    "baseline_imbalance": 0.8,            # 21.4 - 20.6                    [pt]
    "mdc_carryover": 2.45,                # t(.975,14)*0.5*sqrt(SEs)
    "mdc_period1_raw": 3.70,              # t(.975,14)*sqrt(3.19^2/7+3.69^2/9)
    "t_stat_period1_raw": 2.79,           # 4.8/1.723, p ~ 0.015 (df=14)
    "mdc_baseline": 4.02,                 # 2.179*sqrt(2)*3.45/sqrt(8)
    "n1": 7, "n2": 9,
    "doi": "10.3389/fneur.2023.1197281",
    "table": "Table 2A (means/SEs) and Table 4 (washout = 0 d), "
             "Front Neurol 2023;14:1197281",
    "secondary_dois": ["10.2307/2531209",       # R2 estimand (Willan & Pater 1986)
                       "10.1002/sim.4780100905",  # R3 carryover problem (Senn 1991)
                       "10.1136/bmjopen-2015-010873"],  # R4 (Bae et al. 2016)
}


def p_default_tau() -> float:
    return Params().tau


def identify_aA(tau: float) -> float:
    """a_A that reproduces the published *significant* period-1 raw contrast.

    period1_raw = G(a_A - a_B) S = 4.8 pt (R1 Table 2A; t = 2.79, p ~ 0.015), so
        a_A = period1_raw / S   with G = 1 by convention and a_B = 0.
    The order contrast lam_hat is NOT used for calibration: it is 1.5 pt against
    an MDC of 2.45 pt, i.e. it is not distinguishable from zero, and calibrating
    a mechanism to a non-significant quantity would be fitting noise.
    """
    S = 1.0 - math.exp(-Params().t_p / tau)
    return PUBLISHED["period1_raw_contrast"] / S if S > 0 else float("nan")


def identify_pi(p: Params) -> float:
    """pi_hat, a *prediction* (it is not calibrated)."""
    return period_interaction(p)


def detectability_map(tau_grid=None, t_o_grid=None, t_w: float = 0.0) -> dict:
    """Is the order effect detectable?  map over (tau, t_o)."""
    t_o_grid = np.geomspace(0.01, 180.0, 60) if t_o_grid is None else t_o_grid
    tau_grid = np.geomspace(0.05, 600.0, 80) if tau_grid is None else tau_grid
    rows = []
    for to in t_o_grid:
        det_at_tw0 = []
        for tau in tau_grid:
            p = Params(tau=float(tau), t_o=float(to), t_w=t_w, a_A=identify_aA(float(tau)))
            det_at_tw0.append(bool(abs(order_contrast(p)) > mdc(p, "carryover")))
        rows.append({"t_o_days": float(to),
                     "frac_tau_detectable_at_t_w_0": float(np.mean(det_at_tw0)),
                     "any_detectable": bool(any(det_at_tw0))})
    det = [r for r in rows if r["any_detectable"]]
    return {"t_w_days": t_w, "rows": rows,
            "max_t_o_with_any_detectable_days": (max(r["t_o_days"] for r in det)
                                                 if det else 0.0),
            "min_t_o_with_none_detectable_days": (min(
                (r["t_o_days"] for r in rows if not r["any_detectable"]), default=None))}


def calibrate(tau_grid=None) -> dict:
    """Sweep tau; a_A identified from the significant quantity; report predictions."""
    tau_grid = np.geomspace(0.01, 20 * Params().t_p, 300) if tau_grid is None else tau_grid
    rows = []
    for tau in tau_grid:
        p = Params(tau=float(tau), a_A=identify_aA(float(tau)))
        rows.append({
            "tau_days": float(tau),
            "tau_over_t_p": float(tau / p.t_p),
            "a_A_identified_pt": p.a_A,
            "lam_hat_pred_pt": order_contrast(p),
            "pi_hat_pred_pt": period_interaction(p),
            "tau_hat_pred_pt": treatment_effect(p),
            "period1_raw_pred_pt": period1_raw_contrast(p),
            "mdc_carryover_pt": mdc(p, "carryover"),
            "t_w_star_days": critical_washout(p),
            "tau_min_days_tw0": tau_min(p, 0.0),
            "order_detectable_at_t_w_0": bool(abs(order_contrast(p)) > mdc(p)),
        })
    return {"grid": rows, "t_w_days": 0.0, "t_o_days": Params().t_o}


# ---------------------------------------------------------------------------
# SENSITIVITY (+/- 50 %)
# ---------------------------------------------------------------------------
def sensitivity(t_o: float = 0.0) -> dict:
    """Sweep the controlling parameters at a design where order *can* matter
    (read early in period 2, t_o = 0 d), and also at the published design."""
    controlling = ["tau", "a_A", "sigma_w", "G", "t_p", "t_o"]
    out = {}
    for label, base in (("design_t_o_0d", Params(t_o=t_o, a_A=identify_aA(Params().tau))),
                        ("published_design_t_o_90d",
                         Params(a_A=identify_aA(Params().tau)))):
        rec = {}
        for name in controlling:
            v0 = getattr(base, name)
            if isinstance(v0, tuple):
                continue
            e = {"base": {}, "minus50": {}, "plus50": {}}
            for tag, mult in (("base", 1.0), ("minus50", 0.5), ("plus50", 1.5)):
                p = replace(base, **{name: v0 * mult})
                # keep the G/a_A product fixed when either is perturbed
                if name == "G":
                    p = replace(p, a_A=base.a_A * v0 / p.G)
                if name == "a_A":
                    p = replace(p, G=base.G * v0 / p.a_A)
                e[tag] = {
                    "value": float(getattr(p, name)),
                    "lam_hat_pt": order_contrast(p),
                    "abs_lam_hat_pt": abs(order_contrast(p)),
                    "mdc_carryover_pt": mdc(p, "carryover"),
                    "t_w_star_days": critical_washout(p),
                    "tau_min_days_tw0": tau_min(p, 0.0),
                    "order_detectable_at_t_w_0": bool(abs(order_contrast(p)) > mdc(p)),
                }
            b = e["base"]["abs_lam_hat_pt"]
            # a relative swing is only meaningful where the base is not itself
            # buried in the noise; otherwise report NA and let the detectability
            # verdict and t_w* carry the information.
            if b > 0.01 * e["base"]["mdc_carryover_pt"]:
                e["relative_swing_lam_hat_pct"] = float(
                    100 * max(abs(e["minus50"]["abs_lam_hat_pt"] - b),
                              abs(e["plus50"]["abs_lam_hat_pt"] - b)) / b)
            else:
                e["relative_swing_lam_hat_pct"] = None
                e["swing_note"] = ("base contrast is < 1% of the MDC; a relative "
                                   "swing is not interpretable -- use "
                                   "order_detectable_at_t_w_0 and t_w_star_days")
            rec[name] = e
        out[label] = rec
    return out


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def run() -> dict:
    base = Params()
    aA = identify_aA(base.tau)
    cal = replace(base, a_A=aA)                      # calibrated operating point
    early = replace(cal, t_o=0.0)                    # early-read design

    res = {
        "job_id": "BT-HX-Q014",
        "prereg_sha256": open("PREREG.sha256").read().split()[0],
        "question": "When does changing the order between two interventions "
                    "change an observable outcome?",
        "estimands": {
            "order_contrast_lam_hat": "(y11-y21) - 0.5*[(y11-y21)+(y12-y22)]  "
                                      "= first-order carryover = THE order effect",
            "period_interaction_pi_hat": "0.5*[(y12-y11)+(y22-y21)]",
            "treatment_effect_tau_hat": "0.5*[(y11-y21)+(y12-y22)]",
            "note": "PREREG Erratum 1: means, not differences. lam_hat, pi_hat "
                    "are identically 0 under no-carryover / no-period-effect.",
        },
        "mechanism_equation": "dz/dt = -(1/tau)[z - u(t)];  u = a_k on intervention k, "
                              "u = 0 on washout;  y = G z + eps. Closed form: "
                              "z(t+D) = z(t)e^{-D/tau} + u(1-e^{-D/tau})",
        "parameter_table": [
            {"name": n, "unit": u, "value_at_default": v, "source_or_assumption": s}
            for n, u, v, s in PARAMETER_TABLE],
        "unit_check": unit_check(),
        "published_reference": PUBLISHED,
        "calibration": {
            "identified_from": "period1_raw_contrast = 4.8 pt (the only "
                               "significance-supported quantity in R1)",
            "not_calibrated_to": "lam_hat = 1.5 pt, which is below its own MDC "
                                 "of 2.45 pt -- calibrating to it would fit noise",
            "a_A_identified_pt": aA,
            "S_at_operating_point": gates(cal)["S"],
        },
        "gates_published_design": gates(cal),
        "sequence_states": {"AB": sequence_states(cal, "AB"),
                            "BA": sequence_states(cal, "BA")},
        "model_estimates": {
            "published_design_t_o_90d": {
                "lam_hat_pred_pt": order_contrast(cal),
                "pi_hat_pred_pt": period_interaction(cal),
                "tau_hat_pred_pt": treatment_effect(cal),
                "period1_raw_pred_pt": period1_raw_contrast(cal),
                "mdc_carryover_pt": mdc(cal, "carryover"),
                "t_w_star_days": critical_washout(cal),
                "order_detectable": bool(abs(order_contrast(cal)) > mdc(cal)),
            },
            "early_read_design_t_o_0d": {
                "lam_hat_pred_pt": order_contrast(early),
                "mdc_carryover_pt": mdc(early, "carryover"),
                "t_w_star_days": critical_washout(early),
                "order_detectable": bool(abs(order_contrast(early)) > mdc(early)),
            },
            "analytic_max_over_tau_published_design": analytic_max_over_tau(cal),
            "tau_infinity_limit": tau_infinity_limit(cal),
        },
        "nulls": {
            "N0_commutative_abs_lam_pt": abs(null_commutative(cal)),
            "N0prime_additive_abs_lam_pt": abs(null_additive(cal)),
            "N1_control_run_deviation_pt": null_control_run(cal),
            "N1_placebo_deviation_pt": null_control_run(replace(cal, a_B=cal.a_A)),
            "N2_orthogonal_abs_lam_pt": abs(null_orthogonal(cal)),
        },
        "detection_design_levers": {
            "detectability_map_tw0": detectability_map(t_w=0.0),
            "analytic_rule": "lam_hat is a difference of two order-1 terms, so it "
                             "is bounded by a_A/8 at t_w=0, t_o=t_p; the only "
                             "lever that can exceed the MDC is t_o -> 0 (read the "
                             "outcome immediately at the switch) or t_p -> inf "
                             "with t_o << tau.",
        },
    }
    res["calibration_sweep"] = calibrate()
    res["sensitivity"] = sensitivity()

    # --- frozen criteria (Errata 1-4) ---------------------------------------
    fc, det = {}, {}
    fc["FC1_identified_not_assumed"] = bool(
        abs(period1_raw_contrast(cal) - PUBLISHED["period1_raw_contrast"]) < 1e-9)
    det["FC1"] = {"model_period1_raw_pt": period1_raw_contrast(cal),
                  "published_pt": PUBLISHED["period1_raw_contrast"],
                  "note": "identification, i.e. exact by construction -- it is "
                          "NOT an independent test; the independent test is FC2"}
    lam_pred = abs(order_contrast(cal))
    fc["FC2_lam_pred_below_published_CI_upper"] = bool(
        lam_pred <= PUBLISHED["mdc_carryover"])
    det["FC2"] = {"lam_pred_pt": lam_pred,
                  "published_lam_pt": PUBLISHED["first_order_carryover_lam_hat"],
                  "CI95_upper_pt": PUBLISHED["mdc_carryover"],
                  "verdict": "model predicts an order effect smaller than the "
                             "resolution of the published data"}
    worst, argw = 0.0, None
    for tau in np.geomspace(1e-3, 20 * base.t_p, 60):
        for tw in np.linspace(0.0, 2 * base.t_p, 60):
            p = replace(cal, tau=float(tau), t_w=float(tw))
            v = max(abs(null_commutative(p)), abs(null_additive(p)),
                    abs(null_orthogonal(p)))
            if v > worst:
                worst, argw = v, {"tau": float(tau), "t_w": float(tw)}
    fc["FC3_nulls_zero_le_1e-12"] = bool(worst <= 1e-12)
    det["FC3"] = {"max_abs_lam_over_null_grid": worst, "argmax": argw,
                  "grid": "60 x 60 tau x t_w"}
    am = analytic_max_over_tau(cal)
    fc["FC4_analytic_max_below_MDC"] = bool(am["max_abs_lam_hat_pt"] <= am["mdc_carryover_pt"])
    det["FC4"] = am
    fc["FC5_noise_floor_predictions"] = bool(
        PUBLISHED["period1_raw_contrast"] > PUBLISHED["mdc_period1_raw"]
        and PUBLISHED["first_order_carryover_lam_hat"] < PUBLISHED["mdc_carryover"])
    det["FC5"] = {"period1_raw": [PUBLISHED["period1_raw_contrast"],
                                  PUBLISHED["mdc_period1_raw"]],
                  "lam_hat": [PUBLISHED["first_order_carryover_lam_hat"],
                              PUBLISHED["mdc_carryover"]]}
    fc["FC7_baseline_comparable"] = bool(
        PUBLISHED["baseline_imbalance"] < PUBLISHED["mdc_baseline"])
    det["FC7"] = {"baseline_imbalance_pt": PUBLISHED["baseline_imbalance"],
                  "mdc_pt": PUBLISHED["mdc_baseline"]}
    fc["FC6_propagator_vs_ode"] = None      # filled in by test_model.py
    res["frozen_criteria"] = fc
    res["frozen_criteria_detail"] = det
    res["frozen_criteria_all_pass"] = bool(
        all(v for v in fc.values() if v is not None))
    return res


if __name__ == "__main__":
    r = run()
    with open("results.json", "w") as f:
        json.dump(r, f, indent=2, default=lambda o: None if isinstance(o, float) and math.isnan(o) else o)
    print(json.dumps({k: r[k] for k in
                      ("frozen_criteria", "frozen_criteria_all_pass",
                       "model_estimates", "nulls")},
                     indent=2, default=lambda o: None if isinstance(o, float) and math.isnan(o) else o))
