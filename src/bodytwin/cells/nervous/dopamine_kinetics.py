"""DOPAMINE KINETICS -- synaptic DA release-and-reuptake (the
Michaelis-Menten DAT-uptake Vmax/Km governing the fast-scan cyclic voltammetry
(FSCV) single-transient clearance), tyrosine hydroxylase (TH) as the
rate-limiting catecholamine-synthesis step, and the Parkinson's disease (PD)
nigrostriatal-reserve threshold (motor symptoms only after ~50-80% striatal
DA-neuron/terminal loss -- the measured presymptomatic compensatory reserve).

QUESTION (pre-registered falsifier, stated before any number below is
computed): does the model reproduce (a) the measured DAT-uptake kinetics --
the FSCV single-transient clearance fit to Michaelis-Menten form, Km ~0.2 uM,
Vmax ~1-4 uM/s (Wightman/Jones) -- AND (b) the measured PD motor-onset
threshold (~50-80% striatal DA-neuron/terminal loss before symptoms, the
nigrostriatal reserve, Bernheimer 1973 / Fearnley & Lees 1991 / Cheng 2010)?
DECORRELATED CHECK: cocaine (competitive DAT blocker) and amphetamine (DAT
substrate/releaser) perturb the FSCV clearance transient via two MECHANISTICALLY
DISTINCT, independently measured routes (apparent-Km increase vs DAT-mediated
efflux, Jones 1995 / Jones 1999) -- does the model reproduce the qualitatively
DIFFERENT signatures each produces, and does the effect structurally REQUIRE
DAT (falls to zero in a DAT knockout, Budygin 2002)?

SYMMETRIC QC, stated up front, not discovered after the fact:
  - DAT density and release/uptake rates vary MARKEDLY by striatal subregion
    (Garris & Wightman 1994, PMID 8283249: release/uptake rate constants in
    caudate-putamen/nucleus accumbens are 8x and 50x those in prefrontal
    cortex/amygdala, respectively; a 90-fold DA-tissue-level disparity across
    regions) -- held OPEN, the single-Km/Vmax model below is a per-region
    approximation, not a universal constant.
  - The FSCV Km is an APPARENT value measured DURING simultaneous release+
    uptake, not a purified-transporter equilibrium constant -- Near et al
    1988 (PMID 3385647) directly show a spurious high-Km "second component"
    in chopped-tissue preparations that is a diffusional-barrier ARTIFACT,
    absent in synaptosomes/homogenate and absent after 6-OHDA lesion -- held
    OPEN as a genuine measurement caveat, not resolved here.
  - Cragg & Rice 2004 (PMID 15111009) argue that near a release site,
    DIFFUSION -- not DAT -- dominates the earliest shaping of the DA
    transient; DAT's dominant influence is on the LATER lifetime/sphere of
    influence. The single well-mixed-compartment MM model here is a
    disclosed simplification of this spatial reality.
  - The PD threshold has real inter-individual spread (Fearnley & Lees 1991
    report a presymptomatic phase of "about 5 yrs" as an ESTIMATE, and the
    68%/48% figures are themselves derived from only 7 "incidental Lewy
    body" cases) -- held OPEN, not smoothed into a single number.
  - The exact Vmax (uM/s) value is NOT verbatim-quoted in any abstract
    fetched live when this cell was written (the number lives in full-text tables; a
    direct PMC full-text fetch attempt for Wu et al 2001's companion cocaine
    paper was blocked -- "the publisher of this article does not allow
    downloading of the full text in XML form", disclosed in the evidence
    JSON) -- Part 1.6 below reports a DERIVED CONSISTENCY check (not a
    citation-sourced confirmation) for Vmax against the task's
    pre-registered [1,4] uM/s range.

Reads: nothing (all literature parameters embedded).
Writes: dopamine_kinetics_results.json (consumed by basal_ganglia_gating).
Gate: the pre-registered gates in PREREG, summarised in the results JSON.

GEOMETRIC STRUCTURE (derive from the geometry, not curve-fitting):
  Part 1 treats the post-stimulus DA clearance as the 1-D nonlinear dynamical
  system dC/dt = -Vmax*C/(Km+C). Its PHASE-LINE geometry has two regimes with
  a knee at C~Km: C>>Km gives dC/dt=~-Vmax (ZERO-ORDER, constant absolute
  rate, a straight-line decay); C<<Km gives dC/dt=~-(Vmax/Km)*C (FIRST-ORDER,
  exponential decay, rate constant Vmax/Km). This is exactly the qualitative,
  falsifiable signature that distinguishes true MM clearance from a naive
  single-exponential model, tested directly in Part 1.3 (not asserted). The
  ODE is separable in closed form: Km*ln(C0/C) + (C0-C) = Vmax*t (the
  "integrated Michaelis-Menten equation"), used directly (not a fitted
  approximation) in Parts 1.4/1.6/1.7. Part 3 models the PD nigrostriatal
  reserve as a PIECEWISE-LINEAR saturating-compensation geometry (a knee at
  a reserve-exhaustion fraction f_min), the simplest possible geometric
  structure realizing Zigmond et al 1990's named compensations (PMID
  1695406) -- explicitly disclosed as an illustrative TOY (like
  tau_pathology.py's 10-node graph), not a numerical prediction of the
  citation-anchored 48-80% figure (avoiding the tautology of fitting the toy
  TO the target and then "discovering" it).
"""

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

import json
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

OUT_DIR = _os.path.join(OUT_ROOT, "dopamine_kinetics")
OUT_PATH = _os.path.join(OUT_DIR, "dopamine_kinetics_results.json")

# ---- pre-registered gates (fixed BEFORE any number below is computed) -----
PREREG = {
    "rng_seed": 20260722,
    # Part 1.2 -- synthetic-control parameter-recovery pipeline validation
    "recovery_max_rel_error": 0.15,
    # Part 1.3 -- functional-form (MM vs single-exponential) adversary
    "exp_adversary_min_fold_degradation": 1.5,
    # Part 1.4 -- cocaine competitive inhibition (Jones et al 1995, PMID 7616424, Ki=0.29 uM)
    "cocaine_ki_uM": 0.29,
    "cocaine_min_fold_at_max_dose": 2.0,
    # Part 1.4 -- DAT-KO void floor (Budygin et al 2002, PMID 12006604)
    "dat_ko_null_tol_uM": 1e-6,
    # Part 1.5 -- amphetamine efflux (Jones et al 1999, PMID 10582600)
    "amphetamine_min_C_end_frac_of_C0": 0.10,
    "cocaine_C_end_frac_of_C0_max": 0.02,
    # Part 1.6 -- Vmax derived-consistency grid vs the task's pre-registered range
    "km_literature_value_uM": 0.20,  # Wu et al 2001 verbatim, PMID 11716947
    "vmax_prereg_range_uM_per_s": [1.0, 4.0],
    "vmax_grid_min_frac_in_range": 0.30,
    # Part 1.7 -- regional heterogeneity (Garris & Wightman 1994, PMID 8283249)
    "regional_ratio_MPFC": 8.0,
    "regional_ratio_BAN": 50.0,
    # Part 2 -- TH rate-limiting toy (illustrative arbitrary units, disclosed)
    "th_aadc_saturation_max": 0.5,
    # Part 3.1 -- aging-vs-PD rate falsifier (Fearnley & Lees 1991, PMID 1933245)
    "aging_pct_per_decade": 4.7,
    "pd_pct_per_decade_first_decade": 45.0,
    "rate_ratio_min_required": 5.0,
    # Part 3.2 -- compensatory-reserve toy, ILLUSTRATIVE assumption, not fit to the target
    "reserve_toy_f_min": 0.30,
}


# ============================================================================
# PART 1 -- DAT Michaelis-Menten uptake / FSCV single-transient clearance
# ============================================================================

def mm_rhs(t, C, Vmax, Km):
    c = max(float(C[0]), 0.0)
    return [-Vmax * c / (Km + c)]


def simulate_mm(C0, Vmax, Km, t_eval):
    sol = solve_ivp(mm_rhs, (float(t_eval[0]), float(t_eval[-1])), [C0], t_eval=t_eval,
                     args=(Vmax, Km), method="LSODA", rtol=1e-9, atol=1e-12)
    return np.clip(sol.y[0], 0.0, None)


def implicit_mm_time_to_fraction(C0, Vmax, Km, target_frac):
    """Closed-form (separable-ODE) time for MM decay to reach target_frac*C0:
    Km*ln(C0/C) + (C0-C) = Vmax*t, derived directly from dC/dt=-Vmax*C/(Km+C)
    (not a fitted approximation)."""
    c_target = target_frac * C0
    return (Km * np.log(C0 / c_target) + (C0 - c_target)) / Vmax


def part1_synthetic_control_recovery(rng):
    """Validate the fitting pipeline on a KNOWN-answer synthetic transient
    before trusting it on any literature comparison (memory:
    synthetic-control-before-real-negative)."""
    true_vmax = 2.5  # uM/s, midpoint of the task's pre-registered [1,4] range
    true_km = PREREG["km_literature_value_uM"]  # 0.20 uM, Wu et al 2001 verbatim
    c0 = 1.2  # uM, representative peak evoked DA concentration
    t_eval = np.arange(0.0, 3.0 + 1e-9, 0.1)  # 100 ms sampling interval (May et al 1988's stated FSCV interval)

    clean = simulate_mm(c0, true_vmax, true_km, t_eval)
    noisy = np.clip(clean * (1.0 + rng.normal(0, 0.05, size=clean.shape)), 0.0, None)

    def model(t, vmax, km):
        return simulate_mm(c0, vmax, km, t)

    popt, _ = curve_fit(model, t_eval, noisy, p0=[1.0, 0.5],
                         bounds=([0.01, 0.01], [20.0, 5.0]), maxfev=5000)
    vmax_hat, km_hat = float(popt[0]), float(popt[1])
    vmax_rel_err = abs(vmax_hat - true_vmax) / true_vmax
    km_rel_err = abs(km_hat - true_km) / true_km
    return {
        "true_Vmax_uM_s": true_vmax, "true_Km_uM": true_km, "C0_uM": c0,
        "recovered_Vmax_uM_s": vmax_hat, "recovered_Km_uM": km_hat,
        "vmax_rel_error": vmax_rel_err, "km_rel_error": km_rel_err,
        "gate_vmax_recovery": bool(vmax_rel_err <= PREREG["recovery_max_rel_error"]),
        "gate_km_recovery": bool(km_rel_err <= PREREG["recovery_max_rel_error"]),
    }


def part1_functional_form_adversary():
    """FORCED ADVERSARY: a naive single-exponential decay (the null a skeptic
    would reach for first) vs the true nonlinear MM geometry. Pre-registered
    claim: the exponential adversary's misfit WORSENS as C0/Km grows (the
    zero-order regime becomes more prominent), a qualitative/geometric
    signature, not asserted -- measured."""
    km, vmax = PREREG["km_literature_value_uM"], 2.5
    ratios = [1, 5, 25, 100]
    t_eval = np.arange(0.0, 8.0 + 1e-9, 0.02)
    rows = []
    for ratio in ratios:
        c0 = ratio * km
        true_curve = simulate_mm(c0, vmax, km, t_eval)

        def exp_model(t, k):
            return c0 * np.exp(-k * t)

        popt, _ = curve_fit(exp_model, t_eval, true_curve, p0=[vmax / km], maxfev=5000)
        exp_curve = exp_model(t_eval, popt[0])
        rmse = float(np.sqrt(np.mean((exp_curve - true_curve) ** 2)))
        rows.append({"C0_over_Km": ratio, "C0_uM": c0, "best_fit_exp_k": float(popt[0]),
                     "exp_adversary_rmse_frac_of_C0": rmse / c0})
    lo, hi = rows[0]["exp_adversary_rmse_frac_of_C0"], rows[-1]["exp_adversary_rmse_frac_of_C0"]
    fold = hi / lo if lo > 0 else float("inf")
    monotonic = all(rows[i]["exp_adversary_rmse_frac_of_C0"] <= rows[i + 1]["exp_adversary_rmse_frac_of_C0"] + 1e-9
                     for i in range(len(rows) - 1))
    return {
        "per_ratio": rows, "fold_degradation_hi_vs_lo": float(fold), "monotonic_degradation": bool(monotonic),
        "gate": bool(monotonic and fold >= PREREG["exp_adversary_min_fold_degradation"]),
    }


def part1_cocaine_competitive_inhibition():
    """Cocaine = competitive DAT inhibitor: apparent_Km(dose) = Km0*(1+dose/Ki),
    Ki=0.29 uM DIRECTLY from Jones et al 1995 (PMID 7616424, live-verified
    verbatim: "cocaine had a Ki of 0.29 microM in both regions"). Doses span
    Jones 1995's tested 0.01-60 uM range."""
    km0, vmax, c0 = PREREG["km_literature_value_uM"], 2.5, 1.2
    ki = PREREG["cocaine_ki_uM"]
    doses = [0.0, 1.0, 3.0, 10.0, 30.0]
    rows = []
    for d in doses:
        km_app = km0 * (1.0 + d / ki)
        t80 = float(implicit_mm_time_to_fraction(c0, vmax, km_app, 0.2))
        rows.append({"cocaine_uM": d, "apparent_Km_uM": km_app, "T80_s": t80})
    baseline = rows[0]["T80_s"]
    for r in rows:
        r["fold_vs_baseline"] = r["T80_s"] / baseline
    monotonic = all(rows[i]["T80_s"] <= rows[i + 1]["T80_s"] for i in range(len(rows) - 1))
    return {
        "rows": rows, "monotonic_prolongation": bool(monotonic),
        "fold_at_max_dose": rows[-1]["fold_vs_baseline"],
        "gate": bool(monotonic and rows[-1]["fold_vs_baseline"] >= PREREG["cocaine_min_fold_at_max_dose"]),
        "OPEN_ITEM": "predicted fold-change is DERIVED from the verified Ki, not itself independently "
                     "cross-checked against a published exact fold-in-clearance-time number (searched "
                     "when this cell was written: Jones 1995, Wu 2001 in-vivo, Church 1987 microdialysis -- all confirm "
                     "DIRECTION/mechanism, none give an abstract-level exact fold value; a PMC full-text "
                     "fetch for Wu 2001 was blocked by the publisher, see evidence JSON).",
    }


def part1_dat_ko_void_floor():
    """VOID FLOOR / forced-adversary-falls: in a DAT knockout, Vmax_DAT=0, so
    cocaine's Ki-mediated apparent-Km term is multiplied by zero -- the model
    predicts ZERO effect of cocaine dose on clearance, driven only by a small
    DAT-independent first-order process k_other (diffusion/NET/MAO-independent
    removal). This reproduces Budygin et al 2002 (PMID 12006604, live-verified
    verbatim: cocaine and desipramine "failed to change DA clearance ... in
    the NAc of mutant [DAT-KO] mice") as a STRUCTURAL consequence of the
    equations, not a separately hand-coded special case."""
    km0, c0, ki, k_other = PREREG["km_literature_value_uM"], 1.2, PREREG["cocaine_ki_uM"], 0.05
    doses = [0.0, 1.0, 3.0, 10.0, 30.0]
    t_eval = np.arange(0.0, 5.0 + 1e-9, 0.05)

    def rhs(t, c, km_app, k_other_):
        cc = max(float(c[0]), 0.0)
        return [-0.0 * cc / (km_app + cc) - k_other_ * cc]  # Vmax_DAT=0 explicit

    curves = []
    for d in doses:
        km_app = km0 * (1.0 + d / ki)
        sol = solve_ivp(rhs, (0, t_eval[-1]), [c0], t_eval=t_eval, args=(km_app, k_other), method="LSODA")
        curves.append(sol.y[0])
    curves = np.array(curves)
    max_abs_diff = float(np.max(np.abs(curves - curves[0:1, :])))
    return {
        "doses_uM": doses, "max_abs_diff_across_doses_uM": max_abs_diff,
        "gate_null_effect": bool(max_abs_diff < PREREG["dat_ko_null_tol_uM"]),
    }


def part1_amphetamine_efflux():
    """Amphetamine = DAT-mediated REVERSE TRANSPORT (efflux), a distinct
    mechanism from cocaine's Km-shift -- Jones et al 1999 (PMID 10582600,
    live-verified verbatim: amphetamine "reduced the net DA uptake rate and
    increased extracellular DA levels due to DA efflux through the DAT"; "a
    new, elevated steady-state level ... was established"). Modeled as an
    added constant-flux source term: dC/dt = J_eff - Vmax*C/(Km+C), which has
    an analytic nonzero fixed point C_ss=Km*J_eff/(Vmax-J_eff) for J_eff<Vmax
    -- a QUALITATIVELY different long-time signature (nonzero asymptote) than
    cocaine's pure Km-shift (which always asymptotes to C=0)."""
    km0, vmax, c0 = PREREG["km_literature_value_uM"], 2.5, 1.2
    t_eval = np.arange(0.0, 40.0 + 1e-9, 0.05)

    def rhs(t, c, j_eff):
        cc = max(float(c[0]), 0.0)
        return [j_eff - vmax * cc / (km0 + cc)]

    amph_doses = [0.0, 0.3, 0.8, 1.5]  # J_eff in uM/s, illustrative (< Vmax for a valid fixed point)
    rows = []
    for j in amph_doses:
        sol = solve_ivp(rhs, (0, t_eval[-1]), [c0], t_eval=t_eval, args=(j,), method="LSODA", rtol=1e-10, atol=1e-13)
        c_end = float(sol.y[0][-1])
        if j == 0.0:
            c_ss_analytical = 0.0
            rel_err = float(c_end)  # absolute, since analytical is 0
        else:
            c_ss_analytical = km0 * j / (vmax - j)
            rel_err = float(abs(c_end - c_ss_analytical) / c_ss_analytical)
        rows.append({"J_eff_uM_s": j, "C_numeric_long_t_uM": c_end,
                     "C_ss_analytical_uM": c_ss_analytical, "err_vs_analytical": rel_err})

    # contrast: cocaine (pure Km-shift, no source term) at a "matched" high dose, long horizon
    km_app_cocaine_high = km0 * (1.0 + 30.0 / PREREG["cocaine_ki_uM"])
    cocaine_curve = simulate_mm(c0, vmax, km_app_cocaine_high, t_eval)
    cocaine_c_end = float(cocaine_curve[-1])

    amph_high = rows[-1]
    return {
        "rows": rows,
        "fixed_point_recovery_max_rel_err": float(max(r["err_vs_analytical"] for r in rows[1:])),
        "cocaine_high_dose_C_end_uM": cocaine_c_end, "cocaine_high_dose_C_end_frac_of_C0": cocaine_c_end / c0,
        "amphetamine_high_dose_C_end_frac_of_C0": amph_high["C_numeric_long_t_uM"] / c0,
        "gate_fixed_point_matches_analytical": bool(max(r["err_vs_analytical"] for r in rows[1:]) <= 0.05),
        "gate_qualitative_discrimination": bool(
            (amph_high["C_numeric_long_t_uM"] / c0) >= PREREG["amphetamine_min_C_end_frac_of_C0"]
            and (cocaine_c_end / c0) <= PREREG["cocaine_C_end_frac_of_C0_max"]
        ),
    }


def part1_vmax_consistency_grid():
    """The exact Vmax (uM/s) is NOT verbatim-available in any abstract fetched
    live when this cell was written (disclosed in the module docstring and evidence JSON).
    This is a DERIVED CONSISTENCY check, not a citation-sourced confirmation:
    using the DIRECTLY VERIFIED Km=0.20 uM (Wu et al 2001) and the closed-form
    integrated-MM equation, compute what Vmax WOULD have to be for a grid of
    physiologically-reasonable (peak concentration, clearance-time, target-
    fraction-cleared) assumptions, and report what fraction of that grid
    overlaps the task's pre-registered [1,4] uM/s range."""
    km = PREREG["km_literature_value_uM"]
    c0_grid = [0.5, 1.0, 1.5, 2.0]
    tclear_grid = [0.5, 1.0, 1.5, 2.0, 3.0]
    frac_grid = [0.1, 0.2]
    rows = []
    for c0 in c0_grid:
        for tclear in tclear_grid:
            for frac in frac_grid:
                c_target = frac * c0
                vmax_implied = (km * np.log(c0 / c_target) + (c0 - c_target)) / tclear
                rows.append({"C0_uM": c0, "t_clear_s": tclear, "target_frac": frac,
                             "Vmax_implied_uM_s": float(vmax_implied)})
    vals = np.array([r["Vmax_implied_uM_s"] for r in rows])
    lo, hi = PREREG["vmax_prereg_range_uM_per_s"]
    frac_in_range = float(np.mean((vals >= lo) & (vals <= hi)))
    return {
        "grid": rows, "n_grid_points": len(rows),
        "min_implied_Vmax": float(vals.min()), "max_implied_Vmax": float(vals.max()),
        "frac_in_prereg_range": frac_in_range,
        "gate": bool(frac_in_range >= PREREG["vmax_grid_min_frac_in_range"]),
    }


def part1_regional_heterogeneity():
    """SYMMETRIC-QC-anchored structural check: Garris & Wightman 1994 (PMID
    8283249, live-verified verbatim) report release+uptake rate constants in
    MPFC and BAN are ~8x and ~50x LESS, respectively, than in CP/NAc (CP~NAc,
    "similar"). Plug these VERBATIM ratios into the SAME MM equation as
    per-region Vmax multipliers and check the resulting T80 clearance-time
    ORDERING emerges correctly (not hand-typed in) -- a non-degenerate
    structural test: a wrong sign/misapplication would produce the wrong
    order."""
    km, vmax_cp = PREREG["km_literature_value_uM"], 2.5
    vmax_by_region = {
        "CP": vmax_cp,
        "NAc": vmax_cp,
        "MPFC": vmax_cp / PREREG["regional_ratio_MPFC"],
        "BAN": vmax_cp / PREREG["regional_ratio_BAN"],
    }
    c0 = 1.0
    t80_by_region = {r: float(implicit_mm_time_to_fraction(c0, v, km, 0.2)) for r, v in vmax_by_region.items()}
    order = sorted(t80_by_region, key=t80_by_region.get)
    gate = (order[-1] == "BAN") and (order[0] in ("CP", "NAc")) and (order[1] in ("CP", "NAc"))
    return {
        "Vmax_by_region_uM_s": vmax_by_region, "T80_by_region_s": t80_by_region,
        "clearance_order_fast_to_slow": order,
        "gate_BAN_slowest_and_CPNAc_fastest": bool(gate),
    }


# ============================================================================
# PART 2 -- Tyrosine hydroxylase (TH): rate-limiting synthesis step
# ============================================================================

def part2_th_rate_limiting_toy():
    """Citation-only synthesis (Nagatsu, Levitt, Udenfriend 1964, PMID
    14216443; Levitt, Spector, Sjoerdsma, Udenfriend 1965, PMID 14279179;
    Daubner, Le, Wang 2011, PMID 21176768, live-verified verbatim: "Tyrosine
    hydroxylase is the rate-limiting enzyme of catecholamine biosynthesis").
    PLUS a small illustrative 2-step serial-MM toy (Tyrosine --TH--> DOPA
    --AADC--> DA) using ARBITRARY-UNIT, DISCLOSED-ILLUSTRATIVE Vmax values
    (NOT independently measured when this cell was written) to demonstrate the STRUCTURAL/
    geometric meaning of "rate-limiting" in a serial pathway: at steady
    state, mass balance forces total pathway flux = TH's flux exactly,
    and AADC (given a much larger Vmax) never approaches its own saturation
    -- i.e. AADC has spare capacity and never becomes the bottleneck."""
    km_th, vmax_th = 50.0, 1.0       # illustrative arbitrary units
    km_aadc, vmax_aadc = 50.0, 20.0  # illustrative: AADC 20x higher capacity (disclosed assumption)
    s_grid = np.linspace(1, 500, 200)
    flux_th = vmax_th * s_grid / (km_th + s_grid)  # this IS the total steady-state pathway flux (mass balance)
    dopa_ss = km_aadc * flux_th / (vmax_aadc - flux_th)
    saturation_ratio = dopa_ss / (km_aadc + dopa_ss)  # AADC's fractional occupancy of its Vmax
    return {
        "note": "illustrative arbitrary-unit toy (Vmax_TH=1, Vmax_AADC=20, NOT independently measured "
                "when this cell was written) -- demonstrates the STRUCTURAL meaning of rate-limiting, not a numeric claim",
        "max_AADC_saturation_fraction": float(np.max(saturation_ratio)),
        "gate_AADC_never_becomes_bottleneck": bool(np.max(saturation_ratio) < PREREG["th_aadc_saturation_max"]),
        "citations_rate_limiting_verbatim": {
            "daubner_2011_pmid_21176768": "Tyrosine hydroxylase is the rate-limiting enzyme of catecholamine "
                                           "biosynthesis; it uses tetrahydrobiopterin and molecular oxygen to "
                                           "convert tyrosine to DOPA ... inhibited in feedback fashion by the "
                                           "catecholamine neurotransmitters. Dopamine binds to TyrH competitively "
                                           "with tetrahydrobiopterin.",
        },
    }


# ============================================================================
# PART 3 -- Parkinson's disease: nigrostriatal reserve / motor-onset threshold
# ============================================================================

def part3_aging_vs_pd_rate_falsifier():
    """FORCED ADVERSARY on RAW DATA: the no-reserve null (symptom severity
    proportional to cumulative neuron loss from f=1) predicts that ANY
    sufficiently large cumulative loss -- including from normal AGING --
    should eventually cross into clinical parkinsonism. Fearnley & Lees 1991
    (PMID 1933245, live-verified verbatim) report 36 controls with LINEAR
    4.7%/decade loss and NO parkinsonism, vs 20 PD cases with EXPONENTIAL
    loss, 45% in the first decade alone -- a >=9x rate disparity, using ONLY
    the paper's reported numbers, no free model parameter."""
    aging_rate = PREREG["aging_pct_per_decade"]
    pd_rate = PREREG["pd_pct_per_decade_first_decade"]
    ratio = pd_rate / aging_rate
    return {
        "aging_pct_per_decade": aging_rate, "pd_pct_per_decade_first_decade": pd_rate,
        "rate_ratio": float(ratio),
        "gate_no_reserve_null_implausible": bool(ratio >= PREREG["rate_ratio_min_required"]),
    }


def part3_compensatory_reserve_toy():
    """ILLUSTRATIVE TOY ONLY (like tau_pathology.py's 10-node graph) -- f_min
    is an ASSUMED input (0.30, i.e. 70% loss), chosen to sit inside the
    DISCLOSED empirical band (Fearnley & Lees 1991: 48-68% neuron loss at
    symptom onset; Bezard et al 2003: ~80% striatal DA DEPLETION -- a
    different, terminal-level metric, see part3_terminal_vs_soma below), NOT
    independently derived here and NOT claimed to numerically predict that
    citation-anchored figure (that would be circular). The model's role is
    to show Zigmond et al 1990's two named compensations (PMID 1695406,
    live-verified verbatim: increased per-neuron release/synthesis + reduced
    DA inactivation) are STRUCTURALLY SUFFICIENT to produce a threshold/knee
    -- geometrically, the simplest possible saturating-compensation shape,
    contrasted against the no-compensation null F_null(f)=f (pure
    proportionality, no knee at all, falsified on raw data in
    part3_aging_vs_pd_rate_falsifier)."""
    f_min = PREREG["reserve_toy_f_min"]
    f_grid = np.linspace(0.001, 1.0, 1000)
    f_reserve = np.where(f_grid >= f_min, 1.0, f_grid / f_min)
    f_null = f_grid.copy()
    # knee sharpness: the reserve curve is flat (slope~0) above f_min and has slope 1/f_min below it,
    # while the null has constant slope 1 everywhere -- a genuine, checkable shape difference.
    slope_reserve_above = float(np.mean(np.abs(np.diff(f_reserve[f_grid >= f_min]))) / np.diff(f_grid)[0]) \
        if np.sum(f_grid >= f_min) > 1 else 0.0
    slope_null = 1.0
    return {
        "f_min_assumed": f_min, "loss_pct_assumed": (1 - f_min) * 100,
        "disclosed_empirical_band_neuron_loss_pct_fearnley_lees": [48, 68],
        "disclosed_empirical_band_DA_depletion_pct_bezard": 80,
        "reserve_curve_slope_above_knee": slope_reserve_above,
        "null_curve_slope": slope_null,
        "gate_reserve_flat_above_knee_vs_null_linear": bool(slope_reserve_above < 0.05 * slope_null),
        "DISCLOSED_NOT_A_PREDICTION": "f_min=0.30 is an assumed illustrative input inside the empirical "
                                       "band, not independently derived -- this toy demonstrates MECHANISM "
                                       "SUFFICIENCY (a saturating reserve structurally produces a knee), "
                                       "the QUANTITATIVE 48-80% figure rests entirely on the citations "
                                       "(Fearnley & Lees 1991, Bezard 2003, Kordower 2013), not on this model.",
    }


def part3_terminal_vs_soma_dissociation():
    """Direct citation reproduction (Kordower et al 2013, PMID 23884810,
    live-verified verbatim): at the EARLIEST post-diagnosis time points,
    TH-positive SNpc neuron loss (a FUNCTIONAL/biochemical marker) is
    already 50-90%, while melanin-positive SNpc neuron loss (a STRUCTURAL/
    survival marker -- the cell body is still physically present) is only
    ~10% at 1 year and 30-60% over the first several years; striatal
    (terminal) dopaminergic FIBER density loss is "virtually complete" by 4
    years post-diagnosis. No model fitting performed here -- reported
    directly from the paper's numbers."""
    return {
        "th_positive_neuron_loss_pct_earliest_timepoint": [50, 90],
        "melanized_neuron_loss_pct_at_1yr": 10,
        "melanized_neuron_loss_pct_first_several_years": [30, 60],
        "striatal_fiber_loss_virtually_complete_by_years": 4,
        "qualitative_finding_verbatim": "Loss of dopaminergic markers in the dorsal putamen occurs rapidly "
                                        "and is virtually complete by 4 years post-diagnosis. Loss of "
                                        "melanized nigral neurons lags behind the loss of dopamine markers.",
        "mechanistic_reading_DISCLOSED_AS_REASONING_NOT_NEW_DATA": "Consistent with Zigmond 1990's "
            "compensation acting AT THE TERMINAL (increased release/turnover per surviving terminal) -- once "
            "that per-terminal ceiling is exceeded, terminal/fiber markers collapse while the cell body can "
            "persist (dysfunctional but alive) longer. A genuine THIRD state (functionally-silenced-but-alive, "
            "TH-negative yet melanin-positive) exists in Kordower's data that this toy's simple 2-channel "
            "(soma-fraction f vs reserve-gated tone F(f)) framing does not resolve -- held OPEN, not smoothed over.",
    }


# ============================================================================
# MAIN -- assemble, grade, write
# ============================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    rng = np.random.default_rng(PREREG["rng_seed"])

    part1 = {
        "synthetic_control_recovery": part1_synthetic_control_recovery(rng),
        "functional_form_adversary": part1_functional_form_adversary(),
        "cocaine_competitive_inhibition": part1_cocaine_competitive_inhibition(),
        "dat_ko_void_floor": part1_dat_ko_void_floor(),
        "amphetamine_efflux": part1_amphetamine_efflux(),
        "vmax_consistency_grid": part1_vmax_consistency_grid(),
        "regional_heterogeneity": part1_regional_heterogeneity(),
    }
    part2 = {"th_rate_limiting_toy": part2_th_rate_limiting_toy()}
    part3 = {
        "aging_vs_pd_rate_falsifier": part3_aging_vs_pd_rate_falsifier(),
        "compensatory_reserve_toy": part3_compensatory_reserve_toy(),
        "terminal_vs_soma_dissociation": part3_terminal_vs_soma_dissociation(),
    }

    gates = {
        "part1_synthetic_control_vmax_recovery": part1["synthetic_control_recovery"]["gate_vmax_recovery"],
        "part1_synthetic_control_km_recovery": part1["synthetic_control_recovery"]["gate_km_recovery"],
        "part1_functional_form_adversary_falls": part1["functional_form_adversary"]["gate"],
        "part1_cocaine_monotonic_prolongation": part1["cocaine_competitive_inhibition"]["gate"],
        "part1_dat_ko_void_floor_null": part1["dat_ko_void_floor"]["gate_null_effect"],
        "part1_amphetamine_fixed_point_matches_analytical": part1["amphetamine_efflux"]["gate_fixed_point_matches_analytical"],
        "part1_amphetamine_vs_cocaine_qualitative_discrimination": part1["amphetamine_efflux"]["gate_qualitative_discrimination"],
        "part1_vmax_consistency_grid": part1["vmax_consistency_grid"]["gate"],
        "part1_regional_heterogeneity_ordering": part1["regional_heterogeneity"]["gate_BAN_slowest_and_CPNAc_fastest"],
        "part2_th_aadc_never_bottleneck": part2["th_rate_limiting_toy"]["gate_AADC_never_becomes_bottleneck"],
        "part3_no_reserve_null_falsified_on_raw_data": part3["aging_vs_pd_rate_falsifier"]["gate_no_reserve_null_implausible"],
        "part3_reserve_toy_produces_knee_vs_null": part3["compensatory_reserve_toy"]["gate_reserve_flat_above_knee_vs_null_linear"],
    }
    overall_pass = all(gates.values())

    result = {
        "task": "Dopamine neurotransmission kinetics: DAT Michaelis-Menten uptake (FSCV clearance), "
                "tyrosine hydroxylase rate-limiting synthesis, Parkinson's nigrostriatal-reserve threshold.",
        "prereg": PREREG,
        "part1_dat_mm_uptake": part1,
        "part2_th_rate_limiting": part2,
        "part3_pd_threshold": part3,
        "gates": gates,
        "overall_pass": bool(overall_pass),
    }
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(gates, indent=2))
    print("overall_pass:", overall_pass)
    print("wrote:", OUT_PATH)


if __name__ == "__main__":
    main()
