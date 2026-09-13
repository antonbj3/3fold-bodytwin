"""
TMJ lever model -- static-equilibrium bicondylar jaw model, machine-solved.

Every equilibrium is solved as an explicit linear system (numpy.linalg.solve),
never hand algebra. Three decorrelated sub-models:

  A. SAGITTAL single-fulcrum model (bilateral symmetric clench). Two unknowns
     (F_bite, R_condyle), two equations (moment about the condyle, vertical
     force). Derives the molar/incisor MAXIMAL BITE FORCE ratio -- and shows
     it depends ONLY on the ratio of bite-point distances from the condyle
     (d_incisor/d_molar), not on muscle-force magnitude, PCSA, or specific
     tension (those cancel algebraically -- verified numerically below, not
     asserted). Forces the "ignore-the-joint-reaction" adversary to commit to
     a number (R_condyle) and shows it is large and non-zero.

  B. CORONAL two-condyle model (unilateral molar clench). Two unknowns
     (R_L balancing, R_R working), two equations (vertical force, moment about
     the balancing condyle). Sweeps the working-side muscle-activation
     fraction alpha as the forced adversary: can a working-side EMG bias flip
     the balancing>working asymmetry Hylander (1975) reports from EMG data?
     Finds the exact crossover alpha, machine-computed, not assumed.

  C. KINEMATIC two-phase ROM model (rotation then translation). A rigid-body
     chain reproduces total interincisal opening from a rotation-phase chord
     length plus a real measured condylar-translation distance -- cross-
     checked against two independent ROM anchors. Chen (1998)'s finding that
     real motion is simultaneous (not a sharp two-phase switch) is folded in
     as an explicit, disclosed idealization gap, not hidden.

Reads: nothing. Writes: tmj_lever_model_results.json plus the sweep CSVs.

All literature numbers are NCBI-eutils-verified values carrying a DOI/PMID per
number (see ANCHOR below). Geometric nuisance
parameters NOT literature-pinned (r_muscle, intercondylar
geometry fractions, rotation angle) are swept over a disclosed plausible
range rather than fixed at one arbitrary point -- the falsifier is whether the
qualitative/ordinal claims survive across that whole swept instance-space.
"""
import csv
import json
from pathlib import Path

import numpy as np

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = Path(OUT_ROOT) / "tmj_lever_model"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. LITERATURE ANCHORS -- NCBI-eutils-verified (DOI/PMID per number)
# ---------------------------------------------------------------------------
ANCHOR = {
    # Waltimo & Kononen 1995, Acta Odontol Scand, PMID 7484109. n=129 (56M,73F).
    "bite_molar_men_N": 909.0, "bite_molar_men_sd": 177.0,
    "bite_incisor_men_N": 382.0, "bite_incisor_men_sd": 133.0,
    "bite_molar_women_N": 777.0, "bite_molar_women_sd": 168.0,
    "bite_incisor_women_N": 325.0, "bite_incisor_women_sd": 116.0,
    # Waltimo & Kononen 1993, Scand J Dent Res, PMID 8322012. n=30, 2nd device/cohort.
    "bite_molar_men_1993_N": 847.0, "bite_molar_women_1993_N": 597.0,
    # Lewis, Buschang & Throckmorton 2001, Am J Orthod Dentofacial Orthop, PMID 11552129. n=56.
    "incisor_opening_men_mm": 52.1, "incisor_opening_women_mm": 46.0,
    "condylar_translation_men_mm": (15.4, 17.6),
    "condylar_translation_women_mm": (12.4, 12.7),
    # Agrawal et al 2015, Indian J Dent Res, PMID 26481881. n=500. Independent ROM cross-check.
    "max_mouth_opening_men_mm": 50.3, "max_mouth_opening_women_mm": 49.9,
    # van Eijden, Koolstra & Brugman 1995, J Dent Res, PMID 7560404. n=8 cadavers.
    "pcsa_med_pterygoid_cm2": 2.47 + 3.53,  # anterior + posterior heads, direct sum
    # van Eijden, Koolstra & Brugman 1996, Anat Rec, PMID 8955797. n=8 cadavers, 6 AP portions.
    "pcsa_temporalis_portion_range_cm2": (1.82, 2.93),
    "n_temporalis_portions": 6,
    # Daboul et al 2018, J Nutr Health Aging, PMID 30080228, PMCID PMC12880510. n=747.
    "masseter_csa_men_cm2": 5.10, "masseter_csa_women_cm2": 3.90,
    # The already-verified specific-tension convention/range
    # (Arnold 2010 PMC2903973, Powell 1984 PMID 6511546,
    # Maganaris 2001) -- a cross-domain (limb-muscle) proxy, disclosed as such, NOT a
    # jaw-specific number verified.
    "specific_tension_range_Ncm2": (15.0, 100.0),
}

RNG = np.random.default_rng(20260722)  # seeded -- reproducible sweep


# ---------------------------------------------------------------------------
# 2. PART A -- sagittal single-fulcrum model (bilateral symmetric clench)
# ---------------------------------------------------------------------------
def sagittal_equilibrium(F_muscle_bilateral_N, r_muscle_mm, d_bite_mm):
    """Solve [moment about condyle; vertical force] for [F_bite, R_condyle].

    Unknowns x = [F_bite, R_condyle]. Equations:
      (moment about condyle)  d_bite * F_bite + 0 * R_condyle = F_muscle * r_muscle
      (vertical force)        1 * F_bite       + 1 * R_condyle = F_muscle
    """
    A = np.array([[d_bite_mm, 0.0],
                  [1.0, 1.0]])
    b = np.array([F_muscle_bilateral_N * r_muscle_mm, F_muscle_bilateral_N])
    F_bite, R_condyle = np.linalg.solve(A, b)
    return float(F_bite), float(R_condyle)


def total_elevator_pcsa_cm2(sex, temporalis_portion_mean_cm2):
    masseter = ANCHOR["masseter_csa_men_cm2"] if sex == "men" else ANCHOR["masseter_csa_women_cm2"]
    temporalis = temporalis_portion_mean_cm2 * ANCHOR["n_temporalis_portions"]
    med_pterygoid = ANCHOR["pcsa_med_pterygoid_cm2"]
    return masseter + temporalis + med_pterygoid, {
        "masseter_cm2": masseter, "temporalis_cm2": temporalis, "med_pterygoid_cm2": med_pterygoid,
    }


# ---------------------------------------------------------------------------
# 3. PART B -- coronal two-condyle model (unilateral molar clench)
# ---------------------------------------------------------------------------
def coronal_equilibrium(F_muscle_total_N, x_muscle_mm, F_bite_N, x_bite_mm, w_mm):
    """Solve [vertical force; moment about balancing condyle x=0] for [R_L, R_R].

    Left condyle (balancing side) at x=0, right condyle (working side) at x=w.
    Unknowns x = [R_L, R_R]. Equations:
      (vertical force)               R_L + R_R = F_muscle_total - F_bite
      (moment about x=0, balancing)  w * R_R    = F_bite*x_bite - F_muscle_total*x_muscle
    """
    A = np.array([[1.0, 1.0],
                  [0.0, w_mm]])
    b = np.array([F_muscle_total_N - F_bite_N,
                  F_bite_N * x_bite_mm - F_muscle_total_N * x_muscle_mm])
    R_L, R_R = np.linalg.solve(A, b)
    return float(R_L), float(R_R)


def coronal_equilibrium_vec(F_muscle_total_N, x_muscle_mm_arr, F_bite_N, x_bite_mm, w_mm):
    """Vectorized closed-form evaluation of the SAME 2x2 linear system solved by
    coronal_equilibrium() -- exact algebraic solution (R_R = b1/w; R_L = b0-R_R),
    evaluated element-wise over an array of x_muscle positions. Cross-checked
    against the per-point numpy.linalg.solve version in verify_vectorization()."""
    b0 = F_muscle_total_N - F_bite_N
    b1 = F_bite_N * x_bite_mm - F_muscle_total_N * x_muscle_mm_arr
    R_R = b1 / w_mm
    R_L = b0 - R_R
    return R_L, R_R


def verify_vectorization(n_checks=200):
    """Adversary check on the fast path itself: n_checks random points, compare
    the vectorized closed-form against the numpy.linalg.solve reference."""
    max_diff = 0.0
    for _ in range(n_checks):
        F_m = RNG.uniform(500, 3000)
        F_b = RNG.uniform(100, 1000)
        w = RNG.uniform(80, 130)
        x_bite = RNG.uniform(0.5, 1.0) * w
        x_musc = RNG.uniform(0.0, 1.0) * w
        R_L_ref, R_R_ref = coronal_equilibrium(F_m, x_musc, F_b, x_bite, w)
        R_L_vec, R_R_vec = coronal_equilibrium_vec(F_m, np.array([x_musc]), F_b, x_bite, w)
        max_diff = max(max_diff, abs(R_L_ref - R_L_vec[0]), abs(R_R_ref - R_R_vec[0]))
    return max_diff


def find_alpha_critical(F_muscle_total_N, F_bite_N, w_mm, x_bite_frac, xL_frac, xR_frac):
    """Vectorized dense grid search for the working-side force fraction alpha at
    which |R_L| == |R_R| -- the crossover point (machine-computed root, not
    hand-derived)."""
    alphas = np.linspace(0.5, 0.999, 20000)
    x_bite_mm = x_bite_frac * w_mm
    x_muscle_mm = ((1 - alphas) * xL_frac + alphas * xR_frac) * w_mm
    R_L, R_R = coronal_equilibrium_vec(F_muscle_total_N, x_muscle_mm, F_bite_N, x_bite_mm, w_mm)
    diffs = np.abs(R_L) - np.abs(R_R)
    sign_changes = np.where(np.diff(np.sign(diffs)) != 0)[0]
    if len(sign_changes) == 0:
        return None  # no crossover in [0.5, 0.999] -- balancing always wins
    idx = sign_changes[0]
    return float(alphas[idx])


# ---------------------------------------------------------------------------
# 4. PART C -- kinematic two-phase ROM model (rotation then translation)
# ---------------------------------------------------------------------------
def kinematic_opening_mm(L_condyle_to_incisor_mm, theta_rot_deg, d_translation_mm):
    """Rigid-body chord (pure rotation about a fixed condyle) plus a translation
    phase where the (now-translated) condyle center carries the ramus forward,
    adding its distance ~1:1 to interincisal opening (idealized, disclosed)."""
    theta = np.deg2rad(theta_rot_deg)
    rotation_chord_mm = 2.0 * L_condyle_to_incisor_mm * np.sin(theta / 2.0)
    total_mm = rotation_chord_mm + d_translation_mm
    return float(rotation_chord_mm), float(total_mm)


# ---------------------------------------------------------------------------
# 4b. PERTURBATION ADVERSARY -- muscle weakness (reuses Part A unchanged)
# ---------------------------------------------------------------------------
def muscle_weakness_adversary(central_params):
    """van Spronsen et al 1992 (PMID 1613176) measured long-face vs normal adults:
    masseter/medial-pterygoid/anterior-temporalis CSA 30%/22%/15% smaller, with
    correspondingly smaller max molar bite force (no inline % given). This model
    is LINEAR in each muscle's PCSA (F_muscle = Sigma PCSA_i * sigma), so it
    makes a sharp, falsifiable, per-muscle-weighted prediction with no new
    machinery -- reusing sagittal_equilibrium() exactly as in Part A."""
    p = dict(central_params)
    pcsa_total, breakdown = total_elevator_pcsa_cm2("men", p["temporalis_portion_mean_cm2"])
    weakened = dict(
        masseter_cm2=breakdown["masseter_cm2"] * (1 - 0.30),
        temporalis_cm2=breakdown["temporalis_cm2"] * (1 - 0.15),
        med_pterygoid_cm2=breakdown["med_pterygoid_cm2"] * (1 - 0.22),
    )
    pcsa_total_weak = sum(weakened.values())
    F_muscle_healthy = 2.0 * pcsa_total * p["specific_tension_Ncm2"]
    F_muscle_weak = 2.0 * pcsa_total_weak * p["specific_tension_Ncm2"]
    F_bite_molar_healthy, _ = sagittal_equilibrium(F_muscle_healthy, p["r_muscle_mm"], p["d_molar_mm"])
    F_bite_molar_weak, _ = sagittal_equilibrium(F_muscle_weak, p["r_muscle_mm"], p["d_molar_mm"])
    pct_drop = 100.0 * (1 - F_bite_molar_weak / F_bite_molar_healthy)
    pct_pcsa_drop = 100.0 * (1 - pcsa_total_weak / pcsa_total)
    return dict(pcsa_total_healthy_cm2=pcsa_total, pcsa_total_weakened_cm2=pcsa_total_weak,
                pct_pcsa_drop=pct_pcsa_drop,
                F_bite_molar_healthy_N=F_bite_molar_healthy, F_bite_molar_weakened_N=F_bite_molar_weak,
                pct_bite_force_drop=pct_drop,
                # identity check: this model is exactly linear in total PCSA (F_muscle is a
                # sum of independent per-muscle PCSA*sigma terms with a common sigma), so the
                # % bite-force drop must equal the % PCSA drop to machine precision.
                linear_identity_abs_diff=abs(pct_drop - pct_pcsa_drop))


# ---------------------------------------------------------------------------
# 4c. PERTURBATION ADVERSARY -- disc-displacement delayed translation onset
# ---------------------------------------------------------------------------
def disc_displacement_adversary(L_incisor_mm, theta_rot_deg, d_translation_mm, delay_frac):
    """Disc-displacement-with-reduction (Kalaykova et al 2010, PMID 21197509) is
    reported clinically as the moment-of-disc-reduction (MDR) shifting to a LATER
    point in opening as intermittent locking develops -- i.e. translation
    contributes LESS of the total opening at a given rotation, not a literal
    smooth delay (clicking is a discrete event; disclosed as a simplified proxy,
    not a claim to reproduce the click itself). delay_frac in [0,1) models the
    fraction of the normal translation distance that is withheld/delayed."""
    rot_chord, total_normal = kinematic_opening_mm(L_incisor_mm, theta_rot_deg, d_translation_mm)
    d_translation_effective = d_translation_mm * (1.0 - delay_frac)
    _, total_delayed = kinematic_opening_mm(L_incisor_mm, theta_rot_deg, d_translation_effective)
    return dict(delay_frac=delay_frac, total_normal_mm=total_normal, total_delayed_mm=total_delayed,
                opening_deficit_mm=total_normal - total_delayed)


# ---------------------------------------------------------------------------
# 5. SWEEP -- geometric nuisance parameters not literature-pinned
# ---------------------------------------------------------------------------
SWEEP_RANGES = {
    "r_muscle_mm": (15.0, 30.0),          # elevator-resultant moment arm about condyle
    "d_molar_mm": (40.0, 60.0),           # condyle-to-1st-molar sagittal distance
    "d_incisor_mm": (95.0, 120.0),        # condyle-to-central-incisor sagittal distance
    "specific_tension_Ncm2": ANCHOR["specific_tension_range_Ncm2"],
    "temporalis_portion_mean_cm2": ANCHOR["pcsa_temporalis_portion_range_cm2"],
    "w_intercondylar_mm": (90.0, 120.0),
    "x_bite_frac": (0.75, 0.95),          # working molar position, fraction of w
    "xL_muscle_frac": (0.05, 0.20),       # balancing-side muscle centroid, fraction of w
    "theta_rot_deg": (15.0, 30.0),
}
CENTRAL = {k: 0.5 * (lo + hi) for k, (lo, hi) in SWEEP_RANGES.items()}


def run_sagittal_sweep(n=4000):
    rows = []
    for _ in range(n):
        p = {k: RNG.uniform(lo, hi) for k, (lo, hi) in SWEEP_RANGES.items()}
        pcsa_total, pcsa_break = total_elevator_pcsa_cm2("men", p["temporalis_portion_mean_cm2"])
        F_muscle_one_side = pcsa_total * p["specific_tension_Ncm2"]
        F_muscle_bilateral = 2.0 * F_muscle_one_side
        F_bite_molar, R_molar = sagittal_equilibrium(F_muscle_bilateral, p["r_muscle_mm"], p["d_molar_mm"])
        F_bite_incisor, R_incisor = sagittal_equilibrium(F_muscle_bilateral, p["r_muscle_mm"], p["d_incisor_mm"])
        ratio_predicted = p["d_incisor_mm"] / p["d_molar_mm"]
        ratio_from_forces = F_bite_molar / F_bite_incisor
        rows.append(dict(
            F_muscle_bilateral=F_muscle_bilateral, F_bite_molar=F_bite_molar, F_bite_incisor=F_bite_incisor,
            R_molar=R_molar, R_incisor=R_incisor, ratio_predicted=ratio_predicted,
            ratio_from_forces=ratio_from_forces,
            eff_molar=F_bite_molar / F_muscle_bilateral, eff_incisor=F_bite_incisor / F_muscle_bilateral,
            load_per_bite_molar=R_molar / F_bite_molar, load_per_bite_incisor=R_incisor / F_bite_incisor,
            **p,
        ))
    return rows


def run_central_sagittal():
    p = CENTRAL
    pcsa_total, pcsa_break = total_elevator_pcsa_cm2("men", p["temporalis_portion_mean_cm2"])
    F_muscle_one_side = pcsa_total * p["specific_tension_Ncm2"]
    F_muscle_bilateral = 2.0 * F_muscle_one_side
    F_bite_molar, R_molar = sagittal_equilibrium(F_muscle_bilateral, p["r_muscle_mm"], p["d_molar_mm"])
    F_bite_incisor, R_incisor = sagittal_equilibrium(F_muscle_bilateral, p["r_muscle_mm"], p["d_incisor_mm"])
    return dict(
        pcsa_total_cm2=pcsa_total, pcsa_breakdown=pcsa_break,
        F_muscle_one_side_N=F_muscle_one_side, F_muscle_bilateral_N=F_muscle_bilateral,
        F_bite_molar_N=F_bite_molar, F_bite_incisor_N=F_bite_incisor,
        R_condyle_molar_N=R_molar, R_condyle_incisor_N=R_incisor,
        ratio_predicted=p["d_incisor_mm"] / p["d_molar_mm"],
        ratio_from_forces=F_bite_molar / F_bite_incisor,
        load_per_bite_molar=R_molar / F_bite_molar,
        load_per_bite_incisor=R_incisor / F_bite_incisor,
    )


def run_coronal_sweep(F_muscle_total_N, F_bite_N, n=1500):
    rows = []
    for _ in range(n):
        w = RNG.uniform(*SWEEP_RANGES["w_intercondylar_mm"])
        x_bite_frac = RNG.uniform(*SWEEP_RANGES["x_bite_frac"])
        xL_frac = RNG.uniform(*SWEEP_RANGES["xL_muscle_frac"])
        xR_frac = 1.0 - xL_frac  # mirror symmetry
        x_bite_mm = x_bite_frac * w
        # symmetric bilateral activation (alpha=0.5) -- the conservative baseline
        x_muscle_mm = 0.5 * (xL_frac + xR_frac) * w
        R_L, R_R = coronal_equilibrium(F_muscle_total_N, x_muscle_mm, F_bite_N, x_bite_mm, w)
        alpha_crit = find_alpha_critical(F_muscle_total_N, F_bite_N, w, x_bite_frac, xL_frac, xR_frac)
        rows.append(dict(w_mm=w, x_bite_frac=x_bite_frac, xL_frac=xL_frac, xR_frac=xR_frac,
                          R_L_balancing_N=R_L, R_R_working_N=R_R,
                          balancing_gt_working=abs(R_L) > abs(R_R),
                          working_side_tensile=R_R < 0,
                          alpha_critical=alpha_crit))
    return rows


def run_kinematic_check():
    L = CENTRAL["d_incisor_mm"]
    results = []
    for theta in np.linspace(*SWEEP_RANGES["theta_rot_deg"], 6):
        for d_trans in np.linspace(12.4, 20.7, 6):  # real measured range, Lewis 2001
            rot_chord, total = kinematic_opening_mm(L, theta, d_trans)
            results.append(dict(theta_rot_deg=float(theta), d_translation_mm=float(d_trans),
                                 rotation_chord_mm=rot_chord, total_opening_mm=total))
    return results


# ---------------------------------------------------------------------------
# 6. MACHINE GATES -- PASS/FAIL, pre-registered thresholds
# ---------------------------------------------------------------------------
def main():
    vec_check_max_diff = verify_vectorization()

    central = run_central_sagittal()
    sweep = run_sagittal_sweep()

    ratios_pred = np.array([r["ratio_predicted"] for r in sweep])
    ratios_from_forces = np.array([r["ratio_from_forces"] for r in sweep])
    # identity check: ratio_predicted (pure geometry) must equal ratio_from_forces
    # (full force solve) to numerical precision -- proves F_muscle/r_muscle truly
    # cancel, not just assumed algebraically.
    identity_max_abs_diff = float(np.max(np.abs(ratios_pred - ratios_from_forces)))

    measured_ratio_men = ANCHOR["bite_molar_men_N"] / ANCHOR["bite_incisor_men_N"]
    measured_ratio_women = ANCHOR["bite_molar_women_N"] / ANCHOR["bite_incisor_women_N"]

    gate_F1_ratio_band = dict(
        threshold="predicted ratio sweep range must CONTAIN measured molar/incisor ratio",
        measured_men=measured_ratio_men, measured_women=measured_ratio_women,
        predicted_min=float(ratios_pred.min()), predicted_max=float(ratios_pred.max()),
        verdict="PASS" if (ratios_pred.min() <= measured_ratio_men <= ratios_pred.max()
                           and ratios_pred.min() <= measured_ratio_women <= ratios_pred.max()) else "FAIL",
    )

    gate_F1b_identity = dict(
        threshold="ratio_predicted (geometry-only, d_incisor/d_molar) == ratio_from_forces "
                   "(full muscle-force solve) to <1e-9 -- proves F_muscle & r_muscle cancel",
        max_abs_diff=identity_max_abs_diff,
        verdict="PASS" if identity_max_abs_diff < 1e-9 else "FAIL",
    )

    # F2: absolute bite force regime sanity vs measured, band convention (0.3,3.0)x
    # a band convention also used for the ankle-force comparison.
    lo_ratio_band, hi_ratio_band = 0.30, 3.00
    molar_regime_ok = lo_ratio_band < central["F_bite_molar_N"] / ANCHOR["bite_molar_men_N"] < hi_ratio_band
    incisor_regime_ok = lo_ratio_band < central["F_bite_incisor_N"] / ANCHOR["bite_incisor_men_N"] < hi_ratio_band
    gate_F2_absolute = dict(
        threshold=f"central-estimate F_bite within ({lo_ratio_band},{hi_ratio_band})x measured (men)",
        F_bite_molar_N=central["F_bite_molar_N"], measured_molar_men_N=ANCHOR["bite_molar_men_N"],
        ratio_molar=central["F_bite_molar_N"] / ANCHOR["bite_molar_men_N"],
        F_bite_incisor_N=central["F_bite_incisor_N"], measured_incisor_men_N=ANCHOR["bite_incisor_men_N"],
        ratio_incisor=central["F_bite_incisor_N"] / ANCHOR["bite_incisor_men_N"],
        verdict="PASS" if (molar_regime_ok and incisor_regime_ok) else "FAIL",
    )

    # F3: joint reaction must be strictly positive (non-zero) across the ENTIRE sweep
    # -- forces the "ignore-the-joint-reaction" adversary to fail structurally.
    R_molar_all = np.array([r["R_molar"] for r in sweep])
    R_incisor_all = np.array([r["R_incisor"] for r in sweep])
    gate_F3_nonzero_reaction = dict(
        threshold="R_condyle > 0 N for both molar and incisor bite, across 100% of sweep",
        frac_molar_positive=float(np.mean(R_molar_all > 0)),
        frac_incisor_positive=float(np.mean(R_incisor_all > 0)),
        min_R_molar_N=float(R_molar_all.min()), min_R_incisor_N=float(R_incisor_all.min()),
        verdict="PASS" if (np.all(R_molar_all > 0) and np.all(R_incisor_all > 0)) else "FAIL",
    )

    # F4: joint-load-per-unit-bite-force higher at incisor than molar, across 100% of sweep
    load_molar_all = np.array([r["load_per_bite_molar"] for r in sweep])
    load_incisor_all = np.array([r["load_per_bite_incisor"] for r in sweep])
    gate_F4_efficiency_order = dict(
        threshold="R_condyle/F_bite (incisor) > R_condyle/F_bite (molar), across 100% of sweep "
                  "[molar bite is more joint-load-efficient than incisor bite]",
        frac_holds=float(np.mean(load_incisor_all > load_molar_all)),
        median_load_per_bite_molar=float(np.median(load_molar_all)),
        median_load_per_bite_incisor=float(np.median(load_incisor_all)),
        verdict="PASS" if np.all(load_incisor_all > load_molar_all) else "FAIL",
    )

    # --- Part B: coronal working/balancing asymmetry ---
    coronal = run_coronal_sweep(central["F_muscle_bilateral_N"], central["F_bite_molar_N"])
    balancing_wins = np.array([r["balancing_gt_working"] for r in coronal])
    working_tensile = np.array([r["working_side_tensile"] for r in coronal])
    alpha_crits = np.array([r["alpha_critical"] for r in coronal if r["alpha_critical"] is not None])
    gate_F5_balancing_gt_working = dict(
        threshold="|R_L| (balancing) > |R_R| (working) at alpha=0.5 (symmetric bilateral activation), "
                  "across 100% of geometric sweep -- reproduces Hylander 1975 PMID 1101706 EMG finding",
        frac_holds=float(np.mean(balancing_wins)),
        # Forced-adversary surfaced finding, NOT a bug: the unconstrained 2-force linear
        # solve gives R_R<0 (a tensile/distractive reaction) for a large fraction of the
        # sweep -- physically inadmissible for a compression-only (contact) joint. A real
        # TMJ cannot pull; the true constrained solution clamps R_R to 0 and redistributes
        # the remainder onto R_L, which makes the balancing>working asymmetry MORE extreme,
        # never less. Corroborated externally by Osborn & Baragar 1992 (PMID 1517273), whose
        # own model imposes the same compression-only/direction constraint on condylar
        # reaction. Disclosed explicitly, not hidden inside an absolute-value comparison.
        frac_working_side_tensile_unconstrained=float(np.mean(working_tensile)),
        note="where working_side_tensile=True, the reported |R_R| UNDERSTATES the true "
             "asymmetry under a proper compression-only contact constraint",
        verdict="PASS" if np.all(balancing_wins) else "FAIL",
    )
    gate_F6_crossover_margin = dict(
        threshold="crossover working-side-force-fraction alpha_critical (where balancing/working "
                  "reaction flips) must be > 0.80 (i.e. needs an extreme, EMG-inconsistent bias) "
                  "for at least 95% of the geometric sweep",
        n_with_crossover_in_range=int(len(alpha_crits)),
        n_total=int(len(coronal)),
        frac_no_crossover_or_extreme=float(np.mean(alpha_crits > 0.80)) if len(alpha_crits) else 1.0,
        median_alpha_critical=float(np.median(alpha_crits)) if len(alpha_crits) else None,
        min_alpha_critical=float(np.min(alpha_crits)) if len(alpha_crits) else None,
        verdict="PASS" if (len(alpha_crits) == 0 or float(np.mean(alpha_crits > 0.80)) >= 0.95) else "FAIL",
    )

    # --- Part C: kinematic ROM check ---
    kinematic = run_kinematic_check()
    totals = np.array([k["total_opening_mm"] for k in kinematic])
    measured_rom_lo = min(ANCHOR["incisor_opening_women_mm"], ANCHOR["max_mouth_opening_women_mm"])
    measured_rom_hi = max(ANCHOR["incisor_opening_men_mm"], ANCHOR["max_mouth_opening_men_mm"])
    frac_in_band = float(np.mean((totals >= measured_rom_lo) & (totals <= measured_rom_hi)))
    gate_F7_rom = dict(
        threshold=f"predicted total-opening sweep range must overlap measured band "
                  f"[{measured_rom_lo:.1f},{measured_rom_hi:.1f}]mm "
                  "(Lewis/Buschang/Throckmorton 2001 PMID 11552129; Agrawal 2015 PMID 26481881)",
        predicted_min_mm=float(totals.min()), predicted_max_mm=float(totals.max()),
        predicted_median_mm=float(np.median(totals)),
        measured_band_mm=[measured_rom_lo, measured_rom_hi],
        # Disclosed refinement beyond the pre-registered (weak) overlap test: what FRACTION
        # of the swept grid actually lands inside the measured band. The additive
        # rotation-chord + translation idealization is an upper-bound approximation (Chen
        # 1998, PMID 9590523, found real opening is SIMULTANEOUS rotation+translation, not
        # cleanly separable phases -- so summing a pure-rotation chord and a pure-translation
        # distance double-counts some shared displacement at the high end of the sweep).
        frac_grid_points_in_measured_band=frac_in_band,
        verdict="PASS" if (totals.min() <= measured_rom_hi and totals.max() >= measured_rom_lo) else "FAIL",
    )

    # --- Perturbation adversary: muscle weakness (van Spronsen 1992, PMID 1613176) ---
    weakness = muscle_weakness_adversary(CENTRAL)
    gate_F8_weakness_linearity = dict(
        threshold="% bite-force drop must equal % PCSA drop to <1e-9 (this model is exactly "
                  "linear in PCSA) -- a sharp, falsifiable structural prediction, not a fit",
        pct_pcsa_drop=weakness["pct_pcsa_drop"], pct_bite_force_drop=weakness["pct_bite_force_drop"],
        verdict="PASS" if weakness["linear_identity_abs_diff"] < 1e-9 else "FAIL",
    )

    # --- Perturbation adversary: disc-displacement delayed translation onset ---
    disc = [disc_displacement_adversary(CENTRAL["d_incisor_mm"], CENTRAL["theta_rot_deg"],
                                         16.0, delay) for delay in (0.0, 0.25, 0.5, 0.75)]
    gate_F9_disc_monotonic = dict(
        threshold="opening deficit must increase monotonically with delay_frac (mechanistic "
                  "sanity, not a literature-pinned magnitude -- disclosed as a simplified proxy)",
        deficits_mm=[d["opening_deficit_mm"] for d in disc],
        verdict="PASS" if all(disc[i]["opening_deficit_mm"] <= disc[i + 1]["opening_deficit_mm"]
                               for i in range(len(disc) - 1)) else "FAIL",
    )

    gates = dict(F1_ratio_band=gate_F1_ratio_band, F1b_identity=gate_F1b_identity,
                 F2_absolute_regime=gate_F2_absolute, F3_nonzero_reaction=gate_F3_nonzero_reaction,
                 F4_efficiency_order=gate_F4_efficiency_order,
                 F5_balancing_gt_working=gate_F5_balancing_gt_working,
                 F6_crossover_margin=gate_F6_crossover_margin, F7_rom=gate_F7_rom,
                 F8_weakness_linearity=gate_F8_weakness_linearity,
                 F9_disc_monotonic=gate_F9_disc_monotonic)

    all_pass = all(g["verdict"] == "PASS" for g in gates.values())

    evidence = dict(
        anchor=ANCHOR, sweep_ranges=SWEEP_RANGES, central_params=CENTRAL,
        vectorization_cross_check_max_abs_diff_N=vec_check_max_diff,
        central_sagittal=central,
        sagittal_sweep_summary=dict(
            n=len(sweep),
            ratio_predicted_min=float(ratios_pred.min()), ratio_predicted_max=float(ratios_pred.max()),
            ratio_predicted_median=float(np.median(ratios_pred)),
        ),
        coronal_sweep_summary=dict(
            n=len(coronal),
            frac_balancing_gt_working=float(np.mean(balancing_wins)),
            median_R_L_balancing_N=float(np.median([r["R_L_balancing_N"] for r in coronal])),
            median_R_R_working_N=float(np.median([r["R_R_working_N"] for r in coronal])),
        ),
        kinematic_sweep_summary=dict(
            n=len(kinematic), total_opening_min_mm=float(totals.min()), total_opening_max_mm=float(totals.max()),
            total_opening_median_mm=float(np.median(totals)),
        ),
        muscle_weakness_adversary=weakness,
        disc_displacement_adversary=disc,
        gates=gates,
        overall_verdict="PASS" if all_pass else "FAIL",
    )

    out_path = OUT_DIR / "tmj_lever_model_results.json"
    out_path.write_text(json.dumps(evidence, indent=1))

    # Raw per-sample audit trail (every swept point, not just summaries) --
    # deterministic given the fixed RNG seed, but persisted anyway for direct
    # inspection without re-running.
    for name, rows in (("sagittal_sweep_raw.csv", sweep),
                        ("coronal_sweep_raw.csv", coronal),
                        ("kinematic_sweep_raw.csv", kinematic)):
        if rows:
            with open(OUT_DIR / name, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)

    print(f"Wrote {out_path}")
    print(json.dumps({k: v["verdict"] for k, v in gates.items()}, indent=1))
    print("OVERALL:", evidence["overall_verdict"])
    print("\n--- central sagittal estimate ---")
    for k, v in central.items():
        print(f"  {k}: {v}")
    print("\n--- measured anchors ---")
    print(f"  molar/incisor ratio (men)  : {measured_ratio_men:.3f}")
    print(f"  molar/incisor ratio (women): {measured_ratio_women:.3f}")
    print("\n--- muscle weakness adversary (van Spronsen 1992 PCSA deltas) ---")
    print(f"  PCSA drop: {weakness['pct_pcsa_drop']:.2f}%  ->  bite-force drop: {weakness['pct_bite_force_drop']:.2f}%")
    print("\n--- disc-displacement delayed-translation adversary ---")
    for d in disc:
        print(f"  delay_frac={d['delay_frac']:.2f}: total={d['total_delayed_mm']:.1f}mm "
              f"(deficit {d['opening_deficit_mm']:.1f}mm)")


if __name__ == "__main__":
    main()
