#!/usr/bin/env python3
"""FINGERTIP TACTILE: the cutaneous mechanoreceptor / spatial-sampling layer.

Geometric model of the four glabrous-fingertip mechanoreceptor classes (Merkel SA-I, Meissner
FA-I, Ruffini SA-II, Pacinian FA-II): each afferent population is a 2D point lattice with a finite,
density-linked receptive footprint (sigma_r = spacing/2, density-derived, not fitted), and the cell
asks whether its own sampling geometry (density -> hexagonal nearest-neighbour spacing ->
footprint-limited resolution floor) can support the independently measured human psychophysical
grating-orientation acuity limit. Each receptor integrates the stimulus over its Gaussian footprint
before multiplicative response noise is added, so the decoder has a genuine, density-governed,
noise-limited resolution floor rather than the hyperacuity behaviour of an exact-frequency
point-sampling matched filter (that earlier variant never degraded, even 13x below receptor
spacing; the finite aperture is the fix and is disclosed in honest_gaps).

Literature numbers carry an explicit tier: PRIMARY (quoted from a fetched primary abstract),
WEAKER (citation confirmed and correctly attributed, but the number is a widely cited
textbook/tertiary figure), DERIVED (computed here from a PRIMARY number and a WEAKER ratio, and
re-flagged WEAKER).

Reads: nothing (all inputs embedded). Writes:
<BODYTWIN_OUT>/fingertip_tactile/fingertip_tactile_results.json
Gate: the grading block at the end; the cell prints its grading and exits 0.
"""

import json
import os

import numpy as np
from pathlib import Path

OUT_DIR = Path(os.environ.get("BODYTWIN_OUT", Path.cwd() / "outputs")) / "fingertip_tactile"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_JSON = OUT_DIR / "fingertip_tactile_results.json"

SIGMA_R_RATIO = 0.5          # receptor Gaussian-footprint radius = this x nearest-neighbor spacing
                              # (density-derived, NOT a free-fit parameter)
NOISE_STD_LEVELS = [0.1, 0.2, 0.3]   # relative response-noise sweep (robustness, not tuned)
NOISE_STD_PRIMARY = 0.2      # the pre-registered primary grading value (middle of the sweep)
CRITERION = 0.75             # standard 2AFC psychophysical threshold convention

# ============================================================================================
# 0. LITERATURE ANCHORS -- every PMID verified LIVE via NCBI eutils esearch/efetch
#    (never trusted from memory -- one recalled PMID (Corniani & Saal, guessed 32965153) was
#    CAUGHT WRONG live: esearch returned 32965159, the actual PMID.)
# ============================================================================================
LITERATURE_ANCHORS = {
    "johansson_vallbo_1979": {
        "citation": "Johansson RS, Vallbo AB. Tactile sensibility in the human hand: relative and "
                    "absolute densities of four types of mechanoreceptive units in glabrous skin. "
                    "J Physiol. 1979;286:283-300.",
        "pmid": "439026", "tier": "PRIMARY", "verified_live": True,
        "quoted": "mechanoreceptor density 'increase[s] in the proximo-distal direction', with "
                  "relative densities of 1, 1.6, 4.2 in palm, finger, fingertip respectively; "
                  "estimated ~241 units/cm^2 at the fingertip vs ~58 units/cm^2 in the palm "
                  "(ratio 241/58=4.16, matches the paper's stated ~4.2).",
        "used_for": "TOTAL (all 4 types pooled) density at fingertip (241/cm^2) and palm (58/cm^2) "
                    "-- the only per-region density numbers independently live-extracted.",
    },
    "vallbo_johansson_1984": {
        "citation": "Vallbo AB, Johansson RS. Properties of cutaneous mechanoreceptors in the human "
                    "hand related to touch sensation. Hum Neurobiol. 1984;3(1):3-14.",
        "pmid": "6330008", "tier": "PRIMARY (qualitative)", "verified_live": True,
        "quoted": "'The 17,000 tactile units in this skin area of the human hand are of four "
                  "different types: two fast adapting types, FA I and FA II (formerly RA and PC), "
                  "and two slowly adapting types, SA I and SA II. The receptive field "
                  "characteristics and the densities in the skin of the TYPE I UNITS (FA I and SA "
                  "I) indicate that these account for the detailed spatial resolution that is of "
                  "paramount importance for the motor skill and the explorative role of the hand.'",
        "used_for": "The qualitative claim under test: TYPE I (SA-I+FA-I combined), not TYPE II "
                    "(SA-II+FA-II), is the field's claimed substrate for fine spatial acuity. "
                    "Does NOT give a per-type numeric density table in machine-fetchable form "
                    "(no abstract table) -- per-type split is WEAKER/DERIVED below, disclosed.",
    },
    "johansson_1978_rf": {
        "citation": "Johansson RS. Tactile sensibility in the human hand: receptive field "
                    "characteristics of mechanoreceptive units in the glabrous skin area. J "
                    "Physiol. 1978;281:101-125.",
        "pmid": "702358", "tier": "PRIMARY (qualitative)", "verified_live": True,
        "quoted": "RA(FA-I) and SA-I units: 'receptive fields with several zones of maximal "
                  "sensitivity distributed over an approximately circular or oval area' spanning "
                  "'five to ten papillary ridges' with sensitivity 'diminish[ing] steeply' beyond "
                  "it. PC(FA-II) and SA-II units: 'receptive fields with a single zone of maximal "
                  "sensitivity and gentle continuous threshold increase outside this zone.'",
        "used_for": "CORRECTED this build's initial hypothesis (see honest_gaps): SA-I and "
                    "FA-I are BOTH multi-zone/ridge-linked (grouped together), while SA-II and "
                    "FA-II are BOTH single-zone/smooth (grouped together) -- the type-I/type-II "
                    "split, not an SA-I-alone-vs-everything split, is what the receptive-field "
                    "architecture itself actually supports. No mm-scale RF diameter given (spans "
                    "quoted in papillary-ridge counts; ridge-to-mm conversion NOT independently "
                    "verified live -- left as an honest gap, not fabricated).",
    },
    "van_boven_johnson_1994": {
        "citation": "Van Boven RW, Johnson KO. The limit of tactile spatial resolution in humans: "
                    "grating orientation discrimination at the lip, tongue, and finger. Neurology. "
                    "1994;44(12):2361-2366.",
        "pmid": "7991127", "tier": "PRIMARY", "verified_live": True,
        "quoted": "Grating-orientation groove-width discrimination thresholds: lip 0.51 mm, tongue "
                  "0.58 mm, finger 0.94 mm.",
        "used_for": "THE external psychophysical anchor. Finger threshold 0.94mm groove/bar-width "
                    "sits inside the task's pre-specified 0.8-1.5mm target range. NOTE units: this "
                    "is a groove/bar WIDTH -- period = 2x bar-width. Lip/tongue numbers reserved as "
                    "an honest gap: the cleanest possible same-paradigm cross-site test could not "
                    "be completed (no verified oral-mucosa type-I density).",
    },
    "johnson_phillips_1981": {
        "citation": "Johnson KO, Phillips JR. Tactile spatial resolution. I. Two-point "
                    "discrimination, gap detection, grating resolution, and letter recognition. J "
                    "Neurophysiol. 1981;46(6):1177-1192.",
        "pmid": "7320742",
        "tier": "WEAKER -- citation confirmed real live (title/journal/year match); no abstract "
                "indexed for live re-extraction (pre-indexing era, same situation as the "
                "Reilly&Burstein-1975 case). Cited "
                "for the task-defining methodology only; the operative NUMBER used is Van Boven & "
                "Johnson 1994's.",
        "verified_live": "citation only (PMID/title/journal confirmed; abstract not retrievable)",
        "used_for": "Establishes the grating-orientation paradigm this whole model targets.",
    },
    "corniani_saal_2020": {
        "citation": "Corniani G, Saal HP. Tactile innervation densities across the whole body. J "
                    "Neurophysiol. 2020;124(4):1229-1240.",
        "pmid": "32965159",
        "tier": "PRIMARY (qualitative) -- PMID CAUGHT WRONG once live: first-guess recall was "
                "32965153 (a different paper); esearch returned 32965159, used instead. Full "
                "numeric body-map table not retrieved (PMC full-text fetch blocked by a repeated "
                "elink EOF/connection error) -- qualitative abstract claim only.",
        "verified_live": True,
        "quoted": "'Innervation density correlates well with psychophysical spatial acuity across "
                  "different body regions, and, additionally, on hairy skin, with hair follicle "
                  "density.'",
        "used_for": "Independent modern (2020 reconciliation-of-all-sources) confirmation of this "
                    "model's core premise, at the qualitative/whole-body level.",
    },
    "pare_2003_ruffini_paucity": {
        "citation": "Pare M, Behets C, Cornu O. Paucity of presumptive Ruffini corpuscles in the "
                    "index finger pad of humans. J Comp Neurol. 2003.",
        "pmid": "12528190", "tier": "PRIMARY", "verified_live": True,
        "quoted": "'SAII fibers represent approximately 15% of the myelinated mechanosensitive "
                  "axons in the peripheral nerves innervating the volar surface of the hand ... "
                  "Only one presumptive Ruffini corpuscle was found ... very few SAII primary "
                  "afferents are likely to terminate as Ruffini corpuscles in human glabrous skin.'",
        "used_for": "IMPORTANT HONEST CAVEAT on SA-II/'Ruffini' depth claims below: SA-II as a "
                    "PHYSIOLOGICAL afferent class is solid (~15% of recorded hand mechanosensory "
                    "fibers) but 'Ruffini corpuscle' as its finger-PAD structural correlate is not "
                    "(only 1 found in an anatomical study) -- the weakest-anchored structural "
                    "identity of the four classes.",
    },
    "purves_neuroscience_bookshelf": {
        "citation": "Purves D et al. (eds). Neuroscience, 2nd ed. Ch. 'Mechanoreceptors Specialized "
                    "to Receive Tactile Information.' NCBI Bookshelf NBK10895.",
        "pmid": "n/a (textbook, NCBI Bookshelf)",
        "tier": "WEAKER (textbook secondary synthesis, live-fetched)",
        "verified_live": True,
        "quoted": "Population fractions among hand mechanoreceptors: Meissner/FA-I ~40%, "
                  "Merkel/SA-I ~25%, Ruffini/SA-II ~20%, Pacinian/FA-II ~10-15%. Depth: Meissner "
                  "'between the dermal papillae just beneath the epidermis'; Pacinian 'subcutaneous "
                  "tissue'; Merkel 'located in the epidermis'; Ruffini 'deep in the skin, as well "
                  "as in ligaments and tendons'.",
        "used_for": "DERIVED per-type density split (cross-check variant): applied to the PRIMARY "
                    "verified fingertip total (241/cm^2). Fraction is stated for the 'human hand' "
                    "generally, not fingertip-specific -- an approximation, disclosed.",
    },
    "task_given_split": {
        "citation": "Commonly-cited literature figures for fingertip density, as given in this "
                    "build's task brief (traceable to the Vallbo & Johansson 1984 / Johansson "
                    "& Vallbo 1979 lineage).",
        "pmid": "n/a",
        "tier": "WEAKER -- widely-cited figure, NOT independently re-extracted from a live-fetched "
                "primary-source numeric table (the two source papers' abstracts do "
                "not carry the per-type table in machine-fetchable form).",
        "verified_live": False,
        "quoted": "SA-I ~70/cm^2, FA-I ~140/cm^2 at the fingertip.",
        "used_for": "Primary working per-type split, run ALONGSIDE the independently-DERIVED split "
                    "above as a forced two-variant robustness check.",
    },
    "wikipedia_two_point_bates": {
        "citation": "Bates' Guide to Physical Examination (2007), via Wikipedia 'Two-point "
                    "discrimination' summary.",
        "pmid": "n/a", "tier": "WEAKER (tertiary clinical reference)", "verified_live": True,
        "quoted": "Fingertip 2-point discrimination 2-8mm; palm 8-12mm; lips 2-4mm; shins/back "
                  "30-40mm.",
        "used_for": "C4 cross-site bonus context ONLY (different paradigm than grating "
                    "orientation) -- never the falsifier.",
    },
}

# ============================================================================================
# 1. DENSITY TABLE (units/cm^2 at the fingertip unless noted) -- two independent split variants
# ============================================================================================
TOTAL_FINGERTIP_CM2 = 241.0   # PRIMARY, johansson_vallbo_1979
TOTAL_PALM_CM2 = 58.0         # PRIMARY, johansson_vallbo_1979

DENSITY_CITED = {"SA-I": 70.0, "FA-I": 140.0}
DENSITY_CITED["typeI_pooled"] = DENSITY_CITED["SA-I"] + DENSITY_CITED["FA-I"]           # 210
DENSITY_CITED["typeII_pooled"] = TOTAL_FINGERTIP_CM2 - DENSITY_CITED["typeI_pooled"]     # 31

DENSITY_DERIVED = {
    "SA-I": 0.25 * TOTAL_FINGERTIP_CM2, "FA-I": 0.40 * TOTAL_FINGERTIP_CM2,
    "SA-II": 0.20 * TOTAL_FINGERTIP_CM2,
    "FA-II_mid": 0.125 * TOTAL_FINGERTIP_CM2, "FA-II_lo": 0.10 * TOTAL_FINGERTIP_CM2,
    "FA-II_hi": 0.15 * TOTAL_FINGERTIP_CM2,
}
DENSITY_DERIVED["typeI_pooled"] = DENSITY_DERIVED["SA-I"] + DENSITY_DERIVED["FA-I"]
DENSITY_DERIVED["typeII_pooled_mid"] = DENSITY_DERIVED["SA-II"] + DENSITY_DERIVED["FA-II_mid"]
DENSITY_DERIVED["typeII_pooled_lo"] = DENSITY_DERIVED["SA-II"] + DENSITY_DERIVED["FA-II_lo"]
DENSITY_DERIVED["typeII_pooled_hi"] = DENSITY_DERIVED["SA-II"] + DENSITY_DERIVED["FA-II_hi"]
DENSITY_DERIVED["sum_of_four_check"] = (DENSITY_DERIVED["SA-I"] + DENSITY_DERIVED["FA-I"]
                                        + DENSITY_DERIVED["SA-II"] + DENSITY_DERIVED["FA-II_mid"])

# ============================================================================================
# 2. PRE-REGISTRATION -- frozen BEFORE the (corrected) simulation's authoritative numbers.
#    C1/C2 use a hard sufficiency/insufficiency (physical-impossibility) operator: a receptor
#    population's sampling floor exceeding measured human performance is a geometric
#    contradiction (cannot resolve finer than your own noise-limited aperture floor).
# ============================================================================================
PREREG = {
    "external_anchor_mm": 0.94, "external_anchor_source": "van_boven_johnson_1994",
    "psychophysical_criterion": CRITERION,
    "primary_noise_std": NOISE_STD_PRIMARY, "noise_std_sweep": NOISE_STD_LEVELS,
    "footprint_rule": "sigma_r = 0.5 x hex nearest-neighbor spacing (density-derived, not free-fit)",
    "C1_sufficiency_typeI": {
        "claim": "Pooled type-I (SA-I+FA-I) machine-measured sampling floor (bar-width) <= 0.94mm "
                 "at the PRIMARY noise level, for BOTH density-split variants (cited, derived).",
        "operator": "<=", "gate": "both_variants_at_primary_noise",
    },
    "C2_insufficiency_typeII": {
        "claim": "Pooled type-II (SA-II+FA-II) machine-measured sampling floor > 0.94mm at the "
                 "PRIMARY noise level, for BOTH density-split variants -- geometrically ruling out "
                 "type-II as sole substrate.",
        "operator": ">", "gate": "both_variants_at_primary_noise",
    },
    "C3_individual_channel_honesty_check": {
        "claim": "SA-I-alone and FA-I-alone reported WITHOUT a pre-committed pass/fail -- "
                 "exploratory, symmetric-QC probe of whether either channel ALONE is individually "
                 "sufficient. No cherry-picking: both variants x both channels reported plainly.",
        "gate": "exploratory_report_only",
    },
    "C4_cross_site_direction": {
        "claim": "Fingertip total-density floor (241/cm^2) < palm total-density floor (58/cm^2) "
                 "-- must match the documented DIRECTION of clinical two-point discrimination "
                 "(fingertip << palm). Magnitude (ratio) reported as bonus context only.",
        "gate": "direction_only",
    },
    "C5_lattice_irregularity_robustness": {
        "claim": "Jittering the lattice at realistic fractions of spacing changes the measured "
                 "floor by a bounded, reported factor -- exploratory sensitivity check, not gated.",
        "gate": "exploratory_report_only",
    },
    "unit_conversion_note": "Van Boven & Johnson report groove/bar WIDTH; grating PERIOD = 2x "
                            "bar-width. All internal quantities computed in PERIOD units, converted "
                            "to bar-width (/2) only when compared to the external anchor.",
}

# ============================================================================================
# 3. GEOMETRY -- closed-form hex-lattice nearest-neighbor spacing (derivation, not rote lookup)
# ============================================================================================
def hex_nn_spacing_mm(rho_per_mm2):
    """2D triangular (hexagonal) lattice, number density rho [pts/mm^2] -> nearest-neighbor
    spacing a [mm]. Primitive cell of a triangular lattice with side a is a rhombus of two
    equilateral triangles, area (sqrt(3)/2) a^2, containing exactly 1 point:
    rho = 1 / [(sqrt(3)/2) a^2]  =>  a = sqrt( 2 / (sqrt(3) * rho) )."""
    return np.sqrt(2.0 / (np.sqrt(3.0) * rho_per_mm2))


def closed_form_bar_width_mm(rho_per_mm2):
    """Analytic pure-Nyquist cross-check ONLY (a DIFFERENT, complementary mechanism to the
    footprint+noise simulation below -- expected to be same-ballpark, not identical; both scale
    with 1/sqrt(density), which is the shared geometric content). period_Nyquist = sqrt(3)*a
    (row spacing dy=(sqrt(3)/2)a; a resolvable grating needs >=1 row per half-period).
    bar-width = period/2 = (sqrt(3)/2)*a."""
    a = hex_nn_spacing_mm(rho_per_mm2)
    period = np.sqrt(3.0) * a
    return period / 2.0, period, a


def rho_cm2_to_mm2(rho_cm2):
    return rho_cm2 / 100.0


# ============================================================================================
# 4. SIMULATION -- finite-footprint, noisy point-sampling + ideal-observer (matched-filter DFT)
#    2AFC decoder. This is the MEASURED number (machine cross-check), not a formula taken on faith.
# ============================================================================================
def hex_lattice(spacing_mm, patch_mm, jitter_frac, rng):
    dy = spacing_mm * np.sqrt(3.0) / 2.0
    n_rows = int(patch_mm / dy) + 6
    n_cols = int(patch_mm / spacing_mm) + 6
    xs, ys = [], []
    for r in range(n_rows):
        y = r * dy
        x_off = (spacing_mm / 2.0) if (r % 2 == 1) else 0.0
        xs.append(np.arange(n_cols) * spacing_mm + x_off)
        ys.append(np.full(n_cols, y))
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    x = x - x.mean()
    y = y - y.mean()
    mask = (np.abs(x) <= patch_mm / 2.0) & (np.abs(y) <= patch_mm / 2.0)
    pts = np.stack([x[mask], y[mask]], axis=1)
    if jitter_frac > 0:
        pts = pts + rng.normal(0.0, jitter_frac * spacing_mm, size=pts.shape)
    return pts


def accuracy_at_period(points, period_mm, n_trials, rng, sigma_r_mm, noise_std):
    """Each receptor: Gaussian-footprint-attenuated sinusoidal grating response + additive noise.
    attenuation = exp(-2 pi^2 sigma_r^2 / period^2) is the standard Gaussian-aperture MTF applied
    to the grating's spatial frequency 1/period (same mechanism as a camera pixel's aperture
    rolloff) -- this is what makes the task genuinely harder as period shrinks, unlike the
    discarded pure point-sampling decoder (see module docstring)."""
    x = points[:, 0]
    y = points[:, 1]
    attenuation = np.exp(-2 * np.pi ** 2 * sigma_r_mm ** 2 / period_mm ** 2)
    true_labels = rng.integers(0, 2, size=n_trials)
    phases = rng.uniform(0, period_mm, size=n_trials)
    coord = np.where(true_labels[:, None] == 0, x[None, :], y[None, :])
    signal = attenuation * np.sin(2 * np.pi * (coord - phases[:, None]) / period_mm)
    noise = rng.normal(0.0, noise_std, size=signal.shape)
    vals = signal + noise
    exp_x = np.exp(-2j * np.pi * x / period_mm)
    exp_y = np.exp(-2j * np.pi * y / period_mm)
    Sx = np.abs(vals @ exp_x)
    Sy = np.abs(vals @ exp_y)
    pred = np.where(Sx > Sy, 0, 1)
    return float(np.mean(pred == true_labels))


def measure_threshold(rho_per_mm2, jitter_frac=0.0, patch_mm=20.0, n_trials=500,
                       criterion=CRITERION, seed=0, noise_std=NOISE_STD_PRIMARY,
                       period_grid=None):
    spacing = hex_nn_spacing_mm(rho_per_mm2)
    sigma_r = SIGMA_R_RATIO * spacing
    rng = np.random.default_rng(seed)
    points = hex_lattice(spacing, patch_mm, jitter_frac, rng)
    if period_grid is None:
        period_grid = np.linspace(0.2, 4.5, 87)
    rng2 = np.random.default_rng(seed + 10_000)
    accs = np.array([accuracy_at_period(points, P, n_trials, rng2, sigma_r, noise_std)
                      for P in period_grid])
    above = accs >= criterion
    if not above.any():
        return {"period_star_mm": None, "bar_width_star_mm": None, "n_points": len(points),
                "spacing_mm": float(spacing), "sigma_r_mm": float(sigma_r),
                "accuracy": accs.tolist(), "period_grid": period_grid.tolist(),
                "note": "never reached criterion within tested grid"}
    idx_above = np.where(above)[0]
    i_lo = idx_above.min()
    if i_lo == 0:
        period_star = float(period_grid[0])
    else:
        P_above, acc_above = period_grid[i_lo], accs[i_lo]
        P_below, acc_below = period_grid[i_lo - 1], accs[i_lo - 1]
        frac = 0.0 if acc_above == acc_below else (criterion - acc_below) / (acc_above - acc_below)
        period_star = float(P_below + frac * (P_above - P_below))
    return {"period_star_mm": period_star, "bar_width_star_mm": period_star / 2.0,
            "n_points": len(points), "spacing_mm": float(spacing), "sigma_r_mm": float(sigma_r),
            "accuracy": accs.tolist(), "period_grid": period_grid.tolist()}


# ============================================================================================
# 5. RUN
# ============================================================================================
def main():
    results = {"literature_anchors": LITERATURE_ANCHORS, "prereg": PREREG,
               "density_cited_cm2": DENSITY_CITED, "density_derived_cm2": DENSITY_DERIVED,
               "total_fingertip_cm2": TOTAL_FINGERTIP_CM2, "total_palm_cm2": TOTAL_PALM_CM2}

    cases_cm2 = {
        "SA-I_cited": DENSITY_CITED["SA-I"], "FA-I_cited": DENSITY_CITED["FA-I"],
        "typeI_pooled_cited": DENSITY_CITED["typeI_pooled"],
        "typeII_pooled_cited": DENSITY_CITED["typeII_pooled"],
        "SA-I_derived": DENSITY_DERIVED["SA-I"], "FA-I_derived": DENSITY_DERIVED["FA-I"],
        "typeI_pooled_derived": DENSITY_DERIVED["typeI_pooled"],
        "typeII_pooled_derived_mid": DENSITY_DERIVED["typeII_pooled_mid"],
        "typeII_pooled_derived_lo": DENSITY_DERIVED["typeII_pooled_lo"],
        "typeII_pooled_derived_hi": DENSITY_DERIVED["typeII_pooled_hi"],
        "all4_pooled_total_fingertip": TOTAL_FINGERTIP_CM2,
    }
    seeds = [0, 1, 2]

    # ---- 5a. Core cases at the PRIMARY noise level, 3-seed stability + noise-level sweep ----
    core = {}
    for name, rho_cm2 in cases_cm2.items():
        rho_mm2 = rho_cm2_to_mm2(rho_cm2)
        seed_runs = [measure_threshold(rho_mm2, seed=s, noise_std=NOISE_STD_PRIMARY) for s in seeds]
        bar_widths = [r["bar_width_star_mm"] for r in seed_runs if r["bar_width_star_mm"] is not None]
        noise_sweep = {}
        for ns in NOISE_STD_LEVELS:
            r = measure_threshold(rho_mm2, seed=0, noise_std=ns)
            noise_sweep[str(ns)] = r["bar_width_star_mm"]
        cf_bar, cf_period, cf_a = closed_form_bar_width_mm(rho_mm2)
        core[name] = {
            "rho_cm2": rho_cm2, "rho_mm2": rho_mm2, "hex_spacing_mm": cf_a,
            "closed_form_bar_width_mm": cf_bar,
            "sim_bar_width_mm_per_seed_primary_noise": bar_widths,
            "sim_bar_width_mm_mean": float(np.mean(bar_widths)) if bar_widths else None,
            "sim_bar_width_mm_std": float(np.std(bar_widths)) if bar_widths else None,
            "sim_bar_width_by_noise_level": noise_sweep,
            "closed_form_vs_sim_pct_diff": (
                100.0 * (np.mean(bar_widths) - cf_bar) / cf_bar if bar_widths else None
            ),
            "n_points_seed0": seed_runs[0]["n_points"],
        }
    results["core_fingertip_cases"] = core

    # ---- 5b. Grading against PREREG C1/C2 (mechanical PASS/FAIL @ primary noise, both variants) ----
    anchor = PREREG["external_anchor_mm"]

    def floor_of(name):
        return core[name]["sim_bar_width_mm_mean"]

    def floor_at_noise(name, ns):
        return core[name]["sim_bar_width_by_noise_level"][str(ns)]

    c1_cited = floor_of("typeI_pooled_cited") <= anchor
    c1_derived = floor_of("typeI_pooled_derived") <= anchor
    c2_cited = floor_of("typeII_pooled_cited") > anchor
    c2_derived_mid = floor_of("typeII_pooled_derived_mid") > anchor
    c2_derived_lo = floor_of("typeII_pooled_derived_lo") > anchor
    c2_derived_hi = floor_of("typeII_pooled_derived_hi") > anchor

    grading = {
        "C1_sufficiency_typeI": {
            "cited_variant_floor_mm": floor_of("typeI_pooled_cited"), "cited_variant_pass": bool(c1_cited),
            "derived_variant_floor_mm": floor_of("typeI_pooled_derived"), "derived_variant_pass": bool(c1_derived),
            "overall_pass_both_required": bool(c1_cited and c1_derived),
            "robust_across_noise_sweep": {
                name: {str(ns): bool(floor_at_noise(name, ns) <= anchor) for ns in NOISE_STD_LEVELS}
                for name in ["typeI_pooled_cited", "typeI_pooled_derived"]
            },
        },
        "C2_insufficiency_typeII": {
            "cited_variant_floor_mm": floor_of("typeII_pooled_cited"), "cited_variant_pass": bool(c2_cited),
            "derived_mid_floor_mm": floor_of("typeII_pooled_derived_mid"), "derived_mid_pass": bool(c2_derived_mid),
            "derived_lo_floor_mm": floor_of("typeII_pooled_derived_lo"), "derived_lo_pass": bool(c2_derived_lo),
            "derived_hi_floor_mm": floor_of("typeII_pooled_derived_hi"), "derived_hi_pass": bool(c2_derived_hi),
            "overall_pass_both_required": bool(c2_cited and c2_derived_mid),
            "robust_across_noise_sweep": {
                name: {str(ns): bool(floor_at_noise(name, ns) > anchor) for ns in NOISE_STD_LEVELS}
                for name in ["typeII_pooled_cited", "typeII_pooled_derived_mid"]
            },
        },
        "C3_individual_channel_honesty_check": {
            "SA-I_cited_floor_mm": floor_of("SA-I_cited"), "SA-I_cited_le_anchor": bool(floor_of("SA-I_cited") <= anchor),
            "FA-I_cited_floor_mm": floor_of("FA-I_cited"), "FA-I_cited_le_anchor": bool(floor_of("FA-I_cited") <= anchor),
            "SA-I_derived_floor_mm": floor_of("SA-I_derived"), "SA-I_derived_le_anchor": bool(floor_of("SA-I_derived") <= anchor),
            "FA-I_derived_floor_mm": floor_of("FA-I_derived"), "FA-I_derived_le_anchor": bool(floor_of("FA-I_derived") <= anchor),
            "note": "Exploratory, no pre-committed pass/fail gate (PREREG C3) -- reported plainly.",
        },
        "all4_pooled_floor_mm": floor_of("all4_pooled_total_fingertip"),
        "all4_pooled_le_anchor": bool(floor_of("all4_pooled_total_fingertip") <= anchor),
    }
    results["grading"] = grading

    # ---- 5c. Cross-site (C4): palm vs fingertip TOTAL density, direction + bonus magnitude ----
    fingertip_floor = measure_threshold(rho_cm2_to_mm2(TOTAL_FINGERTIP_CM2), seed=0)
    palm_floor = measure_threshold(rho_cm2_to_mm2(TOTAL_PALM_CM2), seed=0)
    direction_pass = fingertip_floor["bar_width_star_mm"] < palm_floor["bar_width_star_mm"]
    predicted_ratio = palm_floor["bar_width_star_mm"] / fingertip_floor["bar_width_star_mm"]
    two_point_fingertip_mm_range = (2.0, 8.0)
    two_point_palm_mm_range = (8.0, 12.0)
    two_pt_ratio_midpoint = np.mean(two_point_palm_mm_range) / np.mean(two_point_fingertip_mm_range)
    results["C4_cross_site"] = {
        "fingertip_total_cm2": TOTAL_FINGERTIP_CM2, "palm_total_cm2": TOTAL_PALM_CM2,
        "fingertip_bar_width_mm": fingertip_floor["bar_width_star_mm"],
        "palm_bar_width_mm": palm_floor["bar_width_star_mm"],
        "direction_pass_fingertip_finer_than_palm": bool(direction_pass),
        "predicted_palm_over_fingertip_ratio": predicted_ratio,
        "bonus_context_two_point_discrimination_mm": {
            "fingertip_range": two_point_fingertip_mm_range, "palm_range": two_point_palm_mm_range,
            "tier": "WEAKER (wikipedia_two_point_bates, different paradigm than grating "
                    "orientation -- illustrative only, NOT the falsifier)",
            "midpoint_ratio_palm_over_fingertip": float(two_pt_ratio_midpoint),
        },
    }

    # ---- 5d. Lattice-irregularity robustness (C5, exploratory, not gated) ----
    jitter_fracs = [0.0, 0.15, 0.30]
    jitter_sweep = {}
    for name in ["typeI_pooled_cited", "typeII_pooled_cited"]:
        rho_mm2 = rho_cm2_to_mm2(cases_cm2[name])
        row = []
        for jf in jitter_fracs:
            r = measure_threshold(rho_mm2, jitter_frac=jf, seed=0)
            row.append({"jitter_frac": jf, "bar_width_star_mm": r["bar_width_star_mm"], "n_points": r["n_points"]})
        jitter_sweep[name] = row
    results["C5_jitter_robustness"] = jitter_sweep

    # ---- 5e. Honest gaps ----
    results["honest_gaps"] = [
        "SELF-CAUGHT SIMULATION BUG (disclosed, fixed, re-run -- see module docstring): the first "
        "decoder version (pure point-sampling matched filter at exact known analytic positions "
        "and period) NEVER degraded with density -- it is mathematically a hyperacuity/exact-"
        "position computation, immune to sampling density, and does not test the intended "
        "phenomenon at all. Fixed by giving each receptor a finite Gaussian integration footprint "
        "(sigma_r = 0.5x its own density-derived spacing -- geometric, not a free parameter) plus "
        "response noise (swept at 3 levels, 0.1/0.2/0.3, for robustness -- 0.2 primary). This is "
        "disclosed as a live example of the OODA-forcing discipline (Observe degenerate result -> "
        "Orient: diagnosed as a hyperacuity/matched-filter exploit, not a real resolution limit -> "
        "Decide: add density-linked footprint+noise -> Act: re-measure), not papered over.",
        "Per-type (SA-I/FA-I/SA-II/FA-II) density split at the fingertip was NOT independently "
        "re-extracted from a live-fetched PRIMARY numeric table -- the two source "
        "papers' (1979, 1984) abstracts do not carry the itemized table in machine-fetchable "
        "form; full-text PDF retrieval was not attempted. Two independent WEAKER-tier variants "
        "(task-given commonly-cited figures; textbook-population-fraction-times-verified-total) "
        "are run side by side as a forced robustness pair -- never a single unchecked split.",
        "The single strongest possible over-determination test -- applying this SAME model to Van "
        "Boven & Johnson 1994's lip/tongue grating-orientation numbers (0.51mm, 0.58mm), a "
        "same-paradigm same-lab 3-point cross-site curve -- could NOT be completed: "
        "no verified oral-mucosa mechanoreceptor density-by-type number was found/verified live. "
        "Flagged as the single most valuable follow-up measurement, not silently dropped.",
        "SA-II/'Ruffini' depth and structural identity in the finger PAD specifically is the "
        "weakest-anchored of the four classes (Pare et al. 2003, PMID 12528190, PRIMARY, live-"
        "verified: only 1 presumptive Ruffini corpuscle found in an anatomical study of the index "
        "finger pad, despite SA-II being ~15% of recorded hand mechanosensory afferents).",
        "Receptive-field diameters in mm are NOT independently verified live for any of the 4 "
        "types -- Johansson 1978 (PMID 702358, PRIMARY) gives type-I RF spans in "
        "PAPILLARY RIDGE COUNTS (5-10 ridges), not mm; the ridge-to-mm conversion constant "
        "(commonly cited ~0.4-0.5mm) was searched for but NOT found in independently-quotable "
        "form -- left as an honest gap. Depth-in-mm figures (e.g. fingertip epidermis "
        "thickness) are similarly not independently verified live -- reported as qualitative layer "
        "location only, tiered WEAKER, no invented precision.",
        "The footprint+noise simulation is a LOWER-BOUND-STYLE geometric floor by construction "
        "(sigma_r tied to density; noise a plausible but not literature-pinned relative level) -- "
        "a floor exceeding the measured 0.94mm anchor is a genuine geometric/SNR contradiction; a "
        "floor at or below 0.94mm is a NECESSARY, not sufficient, condition for that population "
        "being the substrate. The closed-form pure-Nyquist number is reported alongside as an "
        "independent analytic cross-check (different mechanism, same 1/sqrt(density) scaling) -- "
        "the two are NOT expected to match exactly (see closed_form_vs_sim_pct_diff per case).",
        "Real Merkel-cell/Meissner-corpuscle mosaics are NOT perfect hexagonal lattices -- C5's "
        "jitter sweep is a bounded sensitivity check, not itself anchored to a measured real-mosaic "
        "irregularity number (none found/verified).",
        "Single geometric model -- no attempt to model actual measured Merkel/Meissner single-"
        "afferent spatial transfer functions (Phillips & Johnson 1981 II/III), no central/cortical "
        "stage, no adaptation/history-dependence, one body-region pair (fingertip vs palm) for the "
        "cross-site test.",
    ]

    OUT_JSON.write_text(json.dumps(results, indent=2))
    print(f"Wrote {OUT_JSON}")

    print("\n=== CORE FINGERTIP CASES (bar-width mm; external anchor = 0.94mm; primary noise=0.2) ===")
    for name, d in core.items():
        print(f"  {name:32s} rho={d['rho_cm2']:7.2f}/cm^2  spacing={d['hex_spacing_mm']:.3f}mm  "
              f"closed_form={d['closed_form_bar_width_mm']:.3f}mm  "
              f"sim_mean={d['sim_bar_width_mm_mean']:.3f}mm (std={d['sim_bar_width_mm_std']:.3f})  "
              f"noise_sweep={d['sim_bar_width_by_noise_level']}")

    print("\n=== GRADING ===")
    print(json.dumps(grading, indent=2))
    print("\n=== C4 CROSS-SITE ===")
    print(json.dumps(results["C4_cross_site"], indent=2))
    print("\n=== C5 JITTER ROBUSTNESS ===")
    print(json.dumps(jitter_sweep, indent=2))


if __name__ == "__main__":
    main()
