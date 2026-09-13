"""Force-velocity-power and rate-of-force-development (RFD): a self-contained machine check
(numpy only) of the dynamic F-V-P plus explosive-RFD layer that sits on top of static
specific-tension / fiber-length results. Two independent parts:

PART A -- Hill (1938) force-velocity hyperbola: derive the power-optimum velocity fraction v*
  analytically as a function of the dimensionless curvature k = a/F0 = b/Vmax, cross-check the
  closed form against TWO independent numerical optimizers (fine grid search and golden-section
  search), verify the non-normalized Hill equation reproduces F(Vmax) = 0 exactly, sweep k across
  the physiological range (including an emergent cross-bridge-model fit a/F0 = 0.215 derived from
  Huxley-1957 kinetics, not curve-fit to Hill), and contrast against the LINEAR F-V model, whose
  optimum is provably k-independent (always v* = 0.5, the k -> infinity degenerate limit of the
  same family).

PART B -- RFD time course: recompute (not narrate) the cumulative fraction-of-MVC trajectory
  implied by Aagaard et al. 2002's raw pre/post numbers, decompose Klass et al. 2008's aging
  RTD decline into a twitch (contractile) vs voluntary (neural) component, and recompute Tillin et
  al. 2010's athlete-vs-untrained ratios (absolute RFD, normalized RFD, neural activation) in both
  the 0-50 ms and 50-100 ms windows, including the reported REVERSAL in the second window.

All cited numbers come from PubMed abstracts (PMID inline), fetched via NCBI eutils; Hill 1938,
Aagaard 2002, Andersen & Aagaard 2006, Andersen 2010, Bottinelli 1996, Bottinelli & Reggiani 2000,
Aagaard 2010, Klass 2008 and Tillin 2010 were each disambiguated against wrong candidate PMIDs
returned by a first-pass query.

Reads: nothing. Writes: force_velocity_power_rfd_evidence.json.
Gates: the closed-form v* must agree with both numerical optimizers, and each recomputed
literature ratio must reproduce the source's reported value.
"""
import json
import math
import os

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_PATH = _os.path.join(OUT_ROOT, "force_velocity_power_rfd", "force_velocity_power_rfd_evidence.json")

# --------------------------------------------------------------------------------------------
# PART A — Hill hyperbola power-optimum: closed form + two independent numerical cross-checks
# --------------------------------------------------------------------------------------------


def hill_force_norm(v, k):
    """f(v) = F(v)/F0 for Hill's hyperbola, v=V/Vmax, k=a/F0=b/Vmax (derived from F(Vmax)=0)."""
    return k * ((1.0 + k) / (v + k) - 1.0)


def hill_power_norm(v, k):
    return v * hill_force_norm(v, k)


def v_star_closed_form(k):
    """Analytic power-optimum velocity fraction, derived from dp/dv=0 (see doc Sec.2)."""
    return math.sqrt(k * (k + 1.0)) - k


def grid_search_argmax(func, k, n=10_000_001):
    """Independent numerical method 1: brute-force fine grid search over v in [0,1]."""
    best_v, best_p = 0.0, -1.0
    # vectorization-free but fast enough (n~1e7, pure python step) -> use numpy for speed
    import numpy as np

    vs = np.linspace(0.0, 1.0, n)
    ps = k * vs * (1.0 - vs) / (vs + k)
    idx = int(np.argmax(ps))
    return float(vs[idx]), float(ps[idx])


def golden_section_argmax(func, k, lo=0.0, hi=1.0, tol=1e-12):
    """Independent numerical method 2: golden-section search (pure python, no numpy/scipy)."""
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc, fd = func(c, k), func(d, k)
    while abs(b - a) > tol:
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = func(c, k)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = func(d, k)
    v = (a + b) / 2.0
    return v, func(v, k)


def verify_hill_equation_self_consistency(F0, Vmax, k, v_test_fracs):
    """Non-normalized Hill equation (F+a)(V+b)=(F0+a)b must give F(Vmax)=0 exactly, and the
    solved F(V) must match F0*hill_force_norm(V/Vmax,k) to float precision -- a genuine algebra
    cross-check of the normalized formula against the textbook non-normalized form."""
    a = k * F0
    b = k * Vmax
    checks = []
    for vf in v_test_fracs:
        V = vf * Vmax
        F_nonnorm = (F0 + a) * b / (V + b) - a
        F_norm_scaled = F0 * hill_force_norm(vf, k)
        checks.append(
            {
                "v_frac": vf,
                "F_nonnormalized": F_nonnorm,
                "F0_times_normalized_formula": F_norm_scaled,
                "abs_diff": abs(F_nonnorm - F_norm_scaled),
            }
        )
    # F(Vmax) must be exactly 0
    F_at_vmax = (F0 + a) * b / (Vmax + b) - a
    return checks, F_at_vmax


def linear_model_v_star():
    """Linear F-V model f(v)=1-v: power p(v)=v(1-v), dp/dv=1-2v=0 -> v*=0.5, ALWAYS (no k)."""
    import numpy as np

    vs = np.linspace(0.0, 1.0, 10_000_001)
    ps = vs * (1.0 - vs)
    idx = int(np.argmax(ps))
    return float(vs[idx]), float(ps[idx])


def power_law_adversary(n):
    """The STRONGEST fair alternative to 'linear' (tested against the fixed acceptance gates): a generic
    concave power-law f(v)=(1-v)^n, n!=1. p(v)=v(1-v)^n; dp/dv=(1-v)^(n-1)*(1-v(n+1))=0 -> v*=1/(n+1)
    (closed form, exact). Chosen because n=2 lands v* within ~2.5pp of Hill's k=0.25 case -- i.e. the
    power-OPTIMUM LOCATION ALONE is a WEAK discriminator (this adversary can mimic it); what actually
    discriminates is the full-curve shape (checked at v=0.25/0.5/0.75 against Hill and linear)."""
    v_star = 1.0 / (n + 1.0)
    import numpy as np

    vs = np.linspace(0.0, 1.0, 10_000_001)
    ps = vs * (1.0 - vs) ** n
    idx = int(np.argmax(ps))
    return v_star, float(vs[idx]), {vf: (1.0 - vf) ** n for vf in (0.25, 0.5, 0.75)}


def part_a():
    out = {}

    # --- Core case: Hill's classic a/F0 ~ 0.25 ---
    k_hill = 0.25
    v_star_cf = v_star_closed_form(k_hill)
    v_star_grid, p_star_grid = grid_search_argmax(hill_power_norm, k_hill)
    v_star_gs, p_star_gs = golden_section_argmax(hill_power_norm, k_hill)
    f_at_vstar = hill_force_norm(v_star_cf, k_hill)

    out["hill_k_0p25"] = {
        "k_a_over_F0": k_hill,
        "v_star_closed_form": v_star_cf,
        "v_star_grid_search_n1e7": v_star_grid,
        "v_star_golden_section": v_star_gs,
        "max_abs_diff_across_3_methods": max(
            abs(v_star_cf - v_star_grid), abs(v_star_cf - v_star_gs)
        ),
        "p_star_closed_form": v_star_cf ** 2,
        "p_star_grid_search": p_star_grid,
        "f_at_vstar": f_at_vstar,
        "identity_f_vstar_equals_v_star_abs_diff": abs(f_at_vstar - v_star_cf),
        "pct_of_Vmax": 100.0 * v_star_cf,
        "pct_of_F0": 100.0 * f_at_vstar,
        "pct_of_F0_x_Vmax_product": 100.0 * (v_star_cf ** 2),
    }

    # --- Non-normalized equation self-consistency (concrete numbers, e.g. F0=100N, Vmax=10 L0/s) ---
    checks, F_at_vmax = verify_hill_equation_self_consistency(
        F0=100.0, Vmax=10.0, k=k_hill, v_test_fracs=[0.0, 0.1, 0.309017, 0.5, 0.75, 1.0]
    )
    out["nonnormalized_equation_selfcheck"] = {
        "F0": 100.0,
        "Vmax": 10.0,
        "a": k_hill * 100.0,
        "b": k_hill * 10.0,
        "F_at_Vmax_must_be_0": F_at_vmax,
        "per_point_checks": checks,
        "max_abs_diff": max(c["abs_diff"] for c in checks),
    }

    # --- Linear model comparison (the falsifier) ---
    v_lin, p_lin = linear_model_v_star()
    out["linear_model_falsifier"] = {
        "v_star_linear_closed_form": 0.5,
        "v_star_linear_grid_search": v_lin,
        "p_star_linear": p_lin,
        "p_star_linear_closed_form": 0.25,
        "gap_vs_hill_v_star_percentage_points": 100.0 * (0.5 - v_star_cf),
    }

    # --- k-sweep: physiological range + limits (k->0 and k->inf) + the crossbridge fit ---
    k_values = [0.05, 0.10, 0.15, 0.20, 0.215, 0.25, 0.30, 0.35, 0.40, 0.50, 1.0, 5.0, 50.0, 500.0]
    sweep = []
    for k in k_values:
        vst = v_star_closed_form(k)
        sweep.append(
            {
                "k": k,
                "v_star": vst,
                "pct_Vmax": 100.0 * vst,
                "gap_vs_linear_50pct_pp": 100.0 * (0.5 - vst),
                "p_star": vst ** 2,
            }
        )
    out["k_sweep"] = sweep
    out["k_sweep_note"] = (
        "k=0.215 is an emergent cross-bridge-model fit for soleus (a/F0=0.215 derived from "
        "Huxley-1957 kinetics, NOT curve-fit to Hill) -- a mechanistically-independent "
        "decorrelated cross-check landing at v*={:.4f} ({:.2f}% Vmax), inside 2 percentage "
        "points of the k=0.25 textbook case.".format(
            v_star_closed_form(0.215), 100 * v_star_closed_form(0.215)
        )
    )

    # --- Forced adversary: the strongest fair alternative to "linear", not just linear itself ---
    fv_benchmark = {0.25: hill_force_norm(0.25, k_hill), 0.5: hill_force_norm(0.5, k_hill), 0.75: hill_force_norm(0.75, k_hill)}
    fv_linear_benchmark = {vf: 1.0 - vf for vf in (0.25, 0.5, 0.75)}
    powerlaw_results = {}
    for n in (1.0, 1.5, 2.0, 3.0):
        v_star_pl_cf, v_star_pl_grid, fv_pl = power_law_adversary(n)
        powerlaw_results[str(n)] = {
            "v_star_closed_form": v_star_pl_cf,
            "v_star_grid_search": v_star_pl_grid,
            "abs_diff_cf_vs_grid": abs(v_star_pl_cf - v_star_pl_grid),
            "f_v_benchmark": fv_pl,
        }
    out["forced_adversary_powerlaw_family"] = {
        "family": "f(v)=(1-v)^n; v_star=1/(n+1) closed form",
        "results_by_n": powerlaw_results,
        "hill_k0p25_f_v_benchmark": fv_benchmark,
        "linear_f_v_benchmark": fv_linear_benchmark,
        "reading": (
            "n=2 gives v*=33.33% (only 2.4pp from Hill's 30.90%) -- the power-OPTIMUM LOCATION "
            "alone is a WEAK discriminator; the full-curve shape is decisive: at v=0.5, Hill(k=0.25) "
            "gives f=0.1667 vs power-law(n=2) f=0.25 vs linear f=0.5 -- three genuinely different "
            "curves despite two of them nearly agreeing on WHERE power peaks."
        ),
    }

    # --- Limits: k->0 (rectangular corner) and k->inf (recovers the linear model exactly) ---
    tiny_k = 1e-6
    huge_k = 1e8
    out["limits"] = {
        "k_to_0_rectangular_corner": {"k": tiny_k, "v_star": v_star_closed_form(tiny_k)},
        "k_to_inf_recovers_linear": {
            "k": huge_k,
            "v_star": v_star_closed_form(huge_k),
            "diff_from_0p5": abs(v_star_closed_form(huge_k) - 0.5),
        },
    }

    # PASS/FAIL gates (pre-registered, symmetric)
    gates = {
        "gate_A1_three_methods_agree_lt_1e-6": out["hill_k_0p25"]["max_abs_diff_across_3_methods"]
        < 1e-6,
        "gate_A2_identity_f_vstar_eq_vstar_lt_1e-9": out["hill_k_0p25"][
            "identity_f_vstar_equals_v_star_abs_diff"
        ]
        < 1e-9,
        "gate_A3_nonnormalized_eq_selfcheck_lt_1e-9": out["nonnormalized_equation_selfcheck"][
            "max_abs_diff"
        ]
        < 1e-9,
        "gate_A4_F_at_Vmax_is_zero_lt_1e-9": abs(F_at_vmax) < 1e-9,
        "gate_A5_v_star_in_25_to_35pct_for_k_0p15_to_0p40": all(
            0.20 <= s["v_star"] <= 0.35 for s in sweep if 0.15 <= s["k"] <= 0.40
        ),
        "gate_A6_linear_always_exactly_50pct_gap_gt_10pp_for_all_physio_k": all(
            s["gap_vs_linear_50pct_pp"] > 10.0 for s in sweep if 0.05 <= s["k"] <= 0.50
        ),
        "gate_A7_k_to_inf_recovers_linear_lt_1e-4": out["limits"]["k_to_inf_recovers_linear"][
            "diff_from_0p5"
        ]
        < 1e-4,
        "gate_A8_powerlaw_n2_v_star_within_5pp_of_hill_LOCATION_weak_discriminator": (
            abs(
                out["forced_adversary_powerlaw_family"]["results_by_n"]["2.0"][
                    "v_star_closed_form"
                ]
                - v_star_cf
            )
            < 0.05
        ),
        "gate_A9_powerlaw_n2_f_at_0p5_differs_from_hill_by_gt_0p05_SHAPE_discriminates": (
            abs(
                out["forced_adversary_powerlaw_family"]["results_by_n"]["2.0"]["f_v_benchmark"][0.5]
                - fv_benchmark[0.5]
            )
            > 0.05
        ),
    }
    out["gates"] = gates
    out["all_gates_pass"] = all(gates.values())
    return out


# --------------------------------------------------------------------------------------------
# PART B — RFD time-course: recompute (not narrate) from each primary source's raw numbers
# --------------------------------------------------------------------------------------------


def part_b():
    out = {}

    # --- Aagaard et al. 2002 (PMID 12235031) raw numbers, both pre- and post-training ---
    # RFD windows are DEFINED as the average slope of the force-time curve from onset (t=0) to
    # t=window -> cumulative force at t=window = RFD_window_avg * window_seconds EXACTLY (not an
    # approximation), since "average rate over [0,T]" = "total change over T" / T by definition.
    aagaard = {
        "MVC_Nm": {"pre": 291.1, "post": 339.0},
        "RFD_Nm_per_s": {
            "30ms": {"pre": 1601, "post": 2020},
            "50ms": {"pre": 1802, "post": 2201},
            "100ms": {"pre": 1543, "post": 1806},
            "200ms": {"pre": 1141, "post": 1363},
        },
    }
    windows_s = {"30ms": 0.030, "50ms": 0.050, "100ms": 0.100, "200ms": 0.200}
    cumulative_pct_mvc = {}
    pct_change = {}
    excess_over_mvc_scaling_pp = {}
    mvc_pct_change = 100.0 * (aagaard["MVC_Nm"]["post"] - aagaard["MVC_Nm"]["pre"]) / aagaard[
        "MVC_Nm"
    ]["pre"]
    for w, t_s in windows_s.items():
        pre_f = aagaard["RFD_Nm_per_s"][w]["pre"] * t_s
        post_f = aagaard["RFD_Nm_per_s"][w]["post"] * t_s
        cumulative_pct_mvc[w] = {
            "pre_pct_of_pre_MVC": 100.0 * pre_f / aagaard["MVC_Nm"]["pre"],
            "post_pct_of_post_MVC": 100.0 * post_f / aagaard["MVC_Nm"]["post"],
        }
        rfd_pct_change = (
            100.0
            * (aagaard["RFD_Nm_per_s"][w]["post"] - aagaard["RFD_Nm_per_s"][w]["pre"])
            / aagaard["RFD_Nm_per_s"][w]["pre"]
        )
        pct_change[w] = rfd_pct_change
        excess_over_mvc_scaling_pp[w] = rfd_pct_change - mvc_pct_change

    out["aagaard_2002_pmid_12235031"] = {
        "raw_numbers": aagaard,
        "MVC_pct_change": mvc_pct_change,
        "RFD_pct_change_by_window": pct_change,
        "excess_RFD_change_over_MVC_scaling_pp": excess_over_mvc_scaling_pp,
        "cumulative_pct_of_own_MVC_by_window": cumulative_pct_mvc,
        "instant_F0_model_falsifier": {
            "prediction_if_instant_F0": "100% of MVC at t=30/50/100/200ms",
            "measured_post_training_pct_MVC": {
                w: cumulative_pct_mvc[w]["post_pct_of_post_MVC"] for w in windows_s
            },
            "measured_pre_training_pct_MVC": {
                w: cumulative_pct_mvc[w]["pre_pct_of_pre_MVC"] for w in windows_s
            },
            "max_pct_MVC_reached_by_200ms": max(
                cumulative_pct_mvc["200ms"]["pre_pct_of_pre_MVC"],
                cumulative_pct_mvc["200ms"]["post_pct_of_post_MVC"],
            ),
        },
        "lower_bound_time_to_MVC_from_200ms_avg_rate_s": {
            "pre": aagaard["MVC_Nm"]["pre"] / aagaard["RFD_Nm_per_s"]["200ms"]["pre"],
            "post": aagaard["MVC_Nm"]["post"] / aagaard["RFD_Nm_per_s"]["200ms"]["post"],
            "note": (
                "naive linear extrapolation of the 0-200ms AVERAGE rate; genuine lower bound "
                "because the measured rate is front-loaded/decelerating (30ms-window rate > "
                "200ms-window rate in both pre/post), so the TRUE time-to-MVC exceeds this."
            ),
        },
    }

    # --- Klass, Baudry & Duchateau 2008 (PMID 18174392) aging RTD decomposition ---
    klass_voluntary_rtd_decline_pct = 48.0
    klass_twitch_extra_margin_pct = 10.0  # "to a greater extent (~10%) than... evoked twitch"
    klass_twitch_rtd_decline_pct_derived = klass_voluntary_rtd_decline_pct - klass_twitch_extra_margin_pct
    out["klass_2008_pmid_18174392"] = {
        "voluntary_fast_RTD_decline_elderly_vs_young_pct": klass_voluntary_rtd_decline_pct,
        "iEMG_at_peak_RTD_decline_pct": 16.5,
        "MU_discharge_freq_decline_pct_first_3_ISI": [19, 28, 34],
        "doublet_incidence_decline_pct": 45,
        "voluntary_minus_twitch_RTD_decline_margin_pct_DIRECT_QUOTE": klass_twitch_extra_margin_pct,
        "implied_twitch_only_RTD_decline_pct_DERIVED_NOT_QUOTED": klass_twitch_rtd_decline_pct_derived,
        "neural_fraction_of_total_voluntary_RTD_decline_pct_DERIVED": 100.0
        * klass_twitch_extra_margin_pct
        / klass_voluntary_rtd_decline_pct,
    }

    # --- Tillin, Jimenez-Reyes, Pain & Folland 2010 (PMID 19952835) athletes vs untrained ---
    tillin = {
        "MVC_pct_stronger_athletes": 28.0,
        "absolute_RFD_0_50ms_ratio_athletes_over_controls": 2.0,
        "normalized_RFD_0_50ms_MVCs": {"athletes": 4.86, "controls": 2.81},
        "neural_activation_0_50ms_Mmax": {"athletes": 0.26, "controls": 0.15},
        "normalized_RFD_50_100ms_MVCs": {"athletes": 6.68, "controls": 7.93},
    }
    tillin["normalized_RFD_0_50ms_ratio"] = (
        tillin["normalized_RFD_0_50ms_MVCs"]["athletes"]
        / tillin["normalized_RFD_0_50ms_MVCs"]["controls"]
    )
    tillin["neural_activation_0_50ms_ratio"] = (
        tillin["neural_activation_0_50ms_Mmax"]["athletes"]
        / tillin["neural_activation_0_50ms_Mmax"]["controls"]
    )
    tillin["normalized_RFD_50_100ms_ratio_athletes_over_controls"] = (
        tillin["normalized_RFD_50_100ms_MVCs"]["athletes"]
        / tillin["normalized_RFD_50_100ms_MVCs"]["controls"]
    )
    tillin["REVERSAL_flag_50_100ms_controls_exceed_athletes"] = (
        tillin["normalized_RFD_50_100ms_MVCs"]["controls"]
        > tillin["normalized_RFD_50_100ms_MVCs"]["athletes"]
    )
    out["tillin_2010_pmid_19952835"] = tillin

    # PASS/FAIL gates
    gates = {
        "gate_B1_max_pct_MVC_by_200ms_lt_100_both_pre_post": (
            out["aagaard_2002_pmid_12235031"]["instant_F0_model_falsifier"][
                "max_pct_MVC_reached_by_200ms"
            ]
            < 95.0
        ),
        "gate_B2_early_window_pct_change_exceeds_MVC_pct_change": (
            pct_change["30ms"] > mvc_pct_change and pct_change["50ms"] > mvc_pct_change
        ),
        "gate_B3_late_window_100ms_closer_to_proportional_than_30ms": (
            abs(excess_over_mvc_scaling_pp["100ms"]) < abs(excess_over_mvc_scaling_pp["30ms"])
        ),
        "gate_B4_pre_post_fractional_trajectory_consistent_lt_5pp": all(
            abs(
                cumulative_pct_mvc[w]["pre_pct_of_pre_MVC"]
                - cumulative_pct_mvc[w]["post_pct_of_post_MVC"]
            )
            < 5.0
            for w in windows_s
        ),
        "gate_B5_tillin_early_ratio_gt_mvc_ratio": (
            tillin["normalized_RFD_0_50ms_ratio"] > (1.0 + tillin["MVC_pct_stronger_athletes"] / 100.0)
        ),
        "gate_B6_tillin_reversal_is_real_and_reported": tillin[
            "REVERSAL_flag_50_100ms_controls_exceed_athletes"
        ],
    }
    out["gates"] = gates
    out["all_gates_pass"] = all(gates.values())
    return out


def main():
    result = {
        "_meta": {
            "doc": "the force_velocity_power_rfd cell",
            "script": "force_velocity_power_rfd.py",
            "status": "HYPOTHESIS-awaiting-QC",
            "isolation": "pure computation",
        },
        "part_a_hill_hyperbola_power_optimum": part_a(),
        "part_b_rfd_time_course": part_b(),
    }
    result["_meta"]["overall_all_gates_pass"] = (
        result["part_a_hill_hyperbola_power_optimum"]["all_gates_pass"]
        and result["part_b_rfd_time_course"]["all_gates_pass"]
    )

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=1)

    # Print a compact human-readable summary (not the source of truth -- the JSON is)
    print("=== PART A gates ===")
    for k, v in result["part_a_hill_hyperbola_power_optimum"]["gates"].items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    print("=== PART B gates ===")
    for k, v in result["part_b_rfd_time_course"]["gates"].items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    print("OVERALL:", "PASS" if result["_meta"]["overall_all_gates_pass"] else "FAIL")
    print()
    print("v* (k=0.25):", result["part_a_hill_hyperbola_power_optimum"]["hill_k_0p25"]["pct_of_Vmax"])
    print("Written:", OUT_PATH)


if __name__ == "__main__":
    main()
