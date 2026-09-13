"""Bile production & enterohepatic circulation.

Couples the hepatic layer (reads the hepatic_clearance cell's Q_H, read-only) and the GI layer
(reads the gi_absorption_transit cell's SI transit time, read-only, with a documented fallback).
Builds the mass-balance/conservation model: bile production (~500-600 mL/day), the bile acid pool
(~2-4 g) cycling via enterohepatic circulation, ileal reabsorption efficiency (~95%), the critical
micelle concentration (CMC) for fat solubilization, and gallbladder concentration/ejection.

NOTE ON THE ABBREVIATION: every use of "CMC" here means critical micelle concentration (bile-salt
surfactant physical chemistry), not the musculoskeletal Computed Muscle Control method.

FALSIFIER (pre-registered, stated before any number below is computed):
  "Does the model reproduce the MEASURED bile acid pool size and cycling frequency (the
  enterohepatic pool ~3g cycling ~6x/day to deliver the ~12-18 g/day of bile acids that ~600 mL/day
  bile production alone could not supply -- a mass-balance closure) AND the measured ~95% ileal
  reabsorption efficiency (fecal loss ~0.2-0.6 g/day replaced by hepatic synthesis)?"

GEOMETRIC STRUCTURE (derive, don't assert) -- this is a single-compartment conservation LOOP, the
same "stock and flow" structure as any compartmental-analysis / Kirchhoff-current-law node:
  - Pool M (g) = the STOCK. Daily bile-acid flux to the duodenum Phi (g/day) = the THROUGHPUT.
  - Cycling frequency n (per day) is NOT a free parameter -- it is DERIVED as n = Phi/M (turnover
    rate = flux/stock, the classic compartmental-analysis relation; mean residence time tau=1/n=M/Phi).
    This is why Northfield & Hofmann 1975 (PMID 806491, live-verified) found pool SIZE and cycling
    FREQUENCY are inversely related while Phi (secretion rate / hepatic return) stays constant --
    Phi is the conserved invariant, M and n trade off against it, exactly the geometric signature of
    a fixed-throughput compartmental loop, not two independently-free knobs.
  - Steady state (dM/dt=0) forces a genuine conservation law: hepatic synthesis S (g/day) must equal
    fecal loss L (g/day) -- an inflow=outflow balance at the pool node, not an assumption.
  - The expected number of enterohepatic passes a SINGLE bile-acid molecule survives before fecal
    loss is a first-principles absorbing-Markov-chain / geometric-distribution result: if each pass
    is an independent Bernoulli trial with survival (reabsorption) probability eta, the expected
    number of trials before the first "failure" (excretion) is 1/(1-eta) -- a different quantity
    from n (per-day temporal rate) that this script derives and cross-checks against an
    independently-stated external number (Wikipedia's "reused about 20 times"), not fit to it.

SYMMETRIC QC / HELD OPEN (not resolved away):
  - Pool size varies by method: the isotope-dilution band [2,4] g (central 3 g) vs an
    independently-stated [4,6] g band -- the two intervals touch only at 4 g, essentially NO overlap.
    Reported as a genuine, quantified inter-source discrepancy, not forced to one number
    (isotope-dilution vs balance-method variability).
  - CMC depends on bile-acid species (dihydroxy vs trihydroxy vs conjugated) -- only ONE species
    (deoxycholic acid) has a live-verified specific mM value; cholate/CDCA/taurocholate
    specific numbers were NOT independently re-extracted (disclosed gap, Section H / honest_gaps).

Pure Python/numpy, population-parametrized (no individual data).
Reads: hepatic_clearance_results.json (Q_H) and gi_absorption_transit.json (SI transit time), both
read-only, both with documented fallbacks when absent.
Writes: bile_enterohepatic.json
Gate: the pre-registered gate set (pool/cycling mass-balance closure, ~95% reabsorption efficiency,
CMC anchor, gallbladder cross-study consistency, transit-time ceiling).
"""
import json
import sys
from pathlib import Path

import numpy as np

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_PATH = Path(OUT_ROOT) / "bile_enterohepatic" / "bile_enterohepatic.json"
HEPATIC_JSON = Path(OUT_ROOT) / "hepatic_clearance" / "hepatic_clearance_results.json"
GI_JSON = Path(OUT_ROOT) / "gi_absorption_transit" / "gi_absorption_transit.json"

# ---------------------------------------------------------------------------------------------
# CITATIONS -- every one verified via NCBI eutils esearch+esummary+efetch (or a PubMed/StatPearls/
# encyclopaedic page fetch).
# ---------------------------------------------------------------------------------------------
CITATIONS = [
    {"n": 1, "pmid": "19273221",
     "cite": "Hofmann AF (2009). The enterohepatic circulation of bile acids in mammals: form and "
             "functions. Front Biosci 14:2584-98.",
     "verified": "live, abstract fetched verbatim",
     "role": "Primary modern review anchor. Verbatim: 'Efficient hepatic clearance results in low "
             "plasma bile acid levels, and virtually no renal excretion.' -- the qualitative anchor "
             "for the hepatic first-pass coupling (Section F). No numeric pool/frequency data "
             "extractable from the abstract itself (disclosed; full text is IMR Press paywalled)."},
    {"n": 2, "pmid": "12543708",
     "cite": "Russell DW (2003). The enzymes, regulation, and genetics of bile acid synthesis. "
             "Annu Rev Biochem 72:137-74.",
     "verified": "live, abstract fetched verbatim",
     "role": "Confirms bile acid synthesis requires 17 enzymes, tightly regulated by nuclear "
             "hormone receptors (CYP7A1 is the classical rate-limiting, feedback-regulated step "
             "-- textbook-established, not itself numerically re-extracted from this abstract, "
             "disclosed)."},
    {"n": 3, "pmid": "23897680",
     "cite": "Boyer JL (2013). Bile formation and secretion. Compr Physiol 3(3):1035-78.",
     "verified": "live, abstract fetched verbatim",
     "role": "Mechanistic anchor for canalicular bile-salt-DEPENDENT vs -INDEPENDENT transport "
             "(ATP-binding-cassette export pumps creating the osmotic gradient for water flow) and "
             "cholangiocyte ductular modification -- the physiological basis for Hundt/StatPearls' "
             "[14] 75%/25% hepatocyte/cholangiocyte split. No volume numbers in the abstract itself "
             "(disclosed)."},
    {"n": 4, "pmid": "33617926",
     "cite": "Boyer JL (2021). Bile formation and secretion: An update. J Hepatol 75(1):225-6.",
     "verified": "live, bibliographic (title/author/year/journal via esummary only)",
     "role": "Modern update corroborating [3] exists and is actively maintained; not independently "
             "abstracted for numbers (disclosed)."},
    {"n": 5, "pmid": "6853487",
     "cite": "Roda A, Hofmann AF, Mysels KJ (1983). The influence of bile salt structure on "
             "self-association in aqueous solutions. J Biol Chem 258(10):6362-70.",
     "verified": "live, abstract fetched verbatim",
     "role": "PRIMARY biophysical-chemistry CMC anchor. Verbatim: '>50 bile salts and analogues... "
             "CMC values varied from about 1 to greater than 250 mM'; structure rules verified "
             "verbatim: extra hydroxy/oxo group RAISES CMC; Na+ to 0.15M LOWERS CMC for all anionic "
             "bile salts; conjugation (glycine/taurine) causes 'little change'; shorter side chain "
             "RAISES CMC exponentially. Method: surface tension (max bubble-pressure) + dye "
             "solubilization, cross-validated ('results... agreed well') -- two independent routes."},
    {"n": 6, "pmid": "1619357",
     "cite": "Hofmann AF, Mysels KJ (1992). Bile acid solubility and precipitation in vitro and in "
             "vivo: the role of conjugation, pH, and Ca2+ ions. J Lipid Res 33(5):617-26.",
     "verified": "live, abstract fetched verbatim",
     "role": "Qualitative CMC/solubility mechanism (taurine-conjugated soluble at strongly acidic "
             "pH; glycine-conjugated poorly soluble at moderately acidic pH; solubility becomes "
             "'practically unlimited' above CMC). No specific mM values in the abstract (disclosed)."},
    {"n": 7, "pmid": "11687683",
     "cite": "Ziessman HA et al (2001). Normal values for sincalide cholescintigraphy: comparison of "
             "two methods. Radiology 221(2):404-10.",
     "verified": "live, abstract fetched verbatim",
     "role": "PRIMARY gallbladder-ejection-fraction (GBEF) anchor, n=20 healthy subjects, nuclear "
             "cholescintigraphy. Verbatim: 3-min CCK infusion lower-normal-range (2SD) 16.8% (TOO "
             "VARIABLE, rejected by the study's analysis); 60-min infusion gives 31%/41% lower-"
             "normal at 45/60min (the validated protocol). GBEF<35% seen in 6/20 (30%) healthy "
             "subjects under the (rejected) 3-min protocol -- a real, disclosed method-dependence."},
    {"n": 8, "pmid": "15534057",
     "cite": "Krishnamurthy GT, Krishnamurthy S (2004). Constancy and variability of gallbladder "
             "ejection fraction: impact on diagnosis and therapy. J Nucl Med 45(11):1872-7.",
     "verified": "live, abstract fetched verbatim",
     "role": "SECOND, independent (different CCK-8 infusion protocol) nuclear-cholescintigraphy "
             "GBEF anchor. Verbatim control-group means: 66.0%+/-20.5% (3-min) / 73.9%+/-17.7% "
             "(10-min protocol); normal cutoffs >=35% (3-min) / >=50% (10-min). Disease comparators "
             "(chronic acalculous/calculous cholecystitis) verbatim-extracted, much lower (14-28%)."},
    {"n": 9, "pmid": "806491",
     "cite": "Northfield TC, Hofmann AF (1975). Biliary lipid output during three meals and an "
             "overnight fast. I. Relationship to bile acid pool size and cholesterol saturation of "
             "bile in gallstone and control subjects. Gut 16(1):1-11.",
     "verified": "live, abstract fetched verbatim",
     "role": "PRIMARY human isotope-dilution + duodenal-perfusion anchor, paired gallstone/control "
             "cohorts (n=7+7). KEY verbatim finding: 'an inverse relationship between the size and "
             "recycling frequency of the bile acid pool, so that secretion rate and hepatic return "
             "of bile acids remained constant, despite a wide range of pool sizes' -- the primary-"
             "source structural confirmation that Phi=M*n is the conserved invariant (Section A), "
             "decorrelated from the tertiary Wikipedia/StatPearls point-value summaries."},
    {"n": 10, "pmid": "1140620",
     "cite": "Northfield TC, Hofmann AF (1975). Biliary lipid output during three meals and an "
             "overnight fast. II. Effect of chenodeoxycholic acid treatment in gallstone subjects. "
             "Gut 16(1):12-7.",
     "verified": "live, bibliographic (companion paper to [9])",
     "role": "Companion paper, same cohort/method, not independently re-abstracted for numbers."},
    {"n": 11, "pmid": "4761610",
     "cite": "McCormick WC, Bell CC, Swell L, Vlahcevic ZR (1973). Cholic acid synthesis as an "
             "index of the severity of liver disease in man. Gut 14(11):895-902.",
     "verified": "live, abstract fetched verbatim",
     "role": "PRIMARY isotope-KINETIC method anchor (Vlahcevic/Swell classical tradition). Verbatim "
             "numbers: cholic acid synthesis 68 mg/day (advanced cirrhosis, n mortality-linked) vs "
             "152 mg/day (mild cirrhosis) -- DISEASE-state, not healthy-control, values (disclosed); "
             "establishes the isotope-dilution kinetic METHOD with real mg/day-scale numbers, order-"
             "of-magnitude consistent with (below, reduced-in-disease direction of) the healthy "
             "~200-600 mg/day fecal-loss/synthesis band. Fractional turnover rate reduced ~50% in "
             "advanced cirrhosis (relative only; no healthy-control absolute FTR in this abstract, "
             "disclosed)."},
    {"n": 12, "pmid": "28249269",
     "cite": "Dawson PA (2017). Roles of Ileal ASBT and OST-alpha-OST-beta in Regulating Bile Acid "
             "Signaling. Dig Dis 35(3):261-6.",
     "verified": "live, bibliographic (title/author/journal/year via esummary)",
     "role": "Mechanism anchor for the ASBT (SLC10A2, apical sodium-dependent bile acid "
             "transporter)-mediated ACTIVE terminal-ileal reabsorption step responsible for the "
             "~95% figure. Not independently abstracted for a numeric % (disclosed;", },
    {"n": 13, "pmid": "26579438",
     "cite": "Ferrebee CB, Dawson PA (2015). Metabolic effects of intestinal absorption and "
             "enterohepatic cycling of bile acids. Acta Pharm Sin B 5(2):129-34.",
     "verified": "live, abstract fetched verbatim",
     "role": "Verbatim: 'an efficient enterohepatic circulation that functions to conserve and "
             "channel the pool of bile acids within the intestinal and hepatobiliary compartments' "
             "-- corroborates the conservation-loop framing (Section A) from a second Dawson-"
             "authored source, independent of Hofmann's [1]."},
    {"n": 14, "source": "NCBI Bookshelf / StatPearls NBK470209",
     "cite": "Hundt M, Basit H, John S. Physiology, Bile Secretion. StatPearls [Internet], updated "
             "2023.",
     "verified": "live, verbatim quotes fetched",
     "role": "PRIMARY numeric anchor for bile volume. Verbatim: 'total bile flow in a day is "
             "approximately 600 ml, of which 75% is derived from hepatocytes and 25% is from "
             "cholangiocytes' (450/150 mL split); 'Approximately half of the hepatocyte component "
             "of bile flow (about 225 ml per day) is bile salt-dependent'; 'Only approximately 5% "
             "of these bile acids are eventually excreted' (=95% reabsorbed); pool maintained "
             "'mainly via the enterohepatic circulation and, to a small extent (about 5%), by the "
             "hepatic synthesis of bile acids, as long as the daily fecal loss of bile acids does "
             "not exceed 20% of the pool.'"},
    {"n": 15, "source": "Wikipedia: Bile acid",
     "verified": "live, textbook-grade, flagged",
     "role": "Verbatim/paraphrased: pool size 4-6 g; 'Human adults secrete between 12 and 18 g of "
             "bile acids into the intestine each day' (EXACT independent match to the task's "
             "Phi band); daily synthesis '~0.3 g/day' vs elsewhere in the SAME article "
             "'approximately 600 mg... synthesized daily to replace bile acids lost in the feces' "
             "-- a disclosed internal (~2x) inconsistency within one tertiary source, not smoothed "
             "over; ~95% ileal active-transport reabsorption (independently restated)."},
    {"n": 16, "source": "Wikipedia: Enterohepatic circulation",
     "verified": "live, textbook-grade, flagged",
     "role": "Verbatim: '95% of the bile acids which are delivered to the duodenum will be "
             "recycled'; 'reused about 20 times, often multiple times during a single digestive "
             "phase' -- the external anchor for the geometric lifetime-cycle-count derivation "
             "(Section B), a LIFETIME/probabilistic count, NOT the per-day temporal cycling rate "
             "(explicitly disambiguated, Section B, to avoid conflating two different quantities)."},
    {"n": 17, "source": "Wikipedia: Gallbladder",
     "verified": "live, textbook-grade, flagged",
     "role": "Verbatim: bile 'concentrated 3-10 fold' in gallbladder storage; capacity ~50 mL; "
             "30-60 mL typically stored."},
    {"n": 18, "source": "Wikipedia: Deoxycholic acid (redirected from Sodium deoxycholate)",
     "verified": "live, textbook-grade, flagged",
     "role": "Verbatim: 'The critical micelle concentration for deoxycholic acid is approximately "
             "2.4-4 mM'; molar mass 392.58 g/mol. The ONE specific-species CMC number verified here "
             "(used as the primary CMC anchor, Section E)."},
    {"n": 19, "source": "Wikipedia: Cholic acid / Chenodeoxycholic acid / Taurocholic acid",
     "verified": "live, molecular weights only",
     "role": "Molar masses: cholic acid 408.57 g/mol; chenodeoxycholic acid 392.57 g/mol; "
             "taurocholic acid 515.71 g/mol -- used for the derived bile concentration calc "
             "(Section D). No CMC values found on these 3 pages (disclosed gap)."},
    {"n": 20, "source": "the hepatic_clearance cell",
     "verified": "reused read-only, not re-derived",
     "role": "Q_H_rest = 1.4132530120481919 L/min, itself anchored to a live-verified [1.0,2.0] "
             "L/min population band -- reused for the hepatic first-pass coupling (Section F)."},
    {"n": 21, "source": "the gi_absorption_transit cell",
     "verified": "reused read-only, not re-derived",
     "role": "SI mean transit time 216 min (3.6 h), liquid gastric t1/2 15 min -- reused for the "
             "decorrelated transit-time cycling-frequency ceiling check (Section G)."},
]

# ---------------------------------------------------------------------------------------------
# PARAMS -- every band tagged with its source (task vs specific live-verified literature)
# ---------------------------------------------------------------------------------------------
TASK_BAND = {
    "bile_volume_ml_day": (500.0, 600.0),
    "pool_g": (2.0, 4.0),
    "pool_g_central": 3.0,
    "cycling_per_day": (4.0, 12.0),
    "cycling_per_day_central": 6.0,
    "flux_g_day": (12.0, 18.0),
    "reabsorption_frac": 0.95,
    "fecal_loss_g_day": (0.2, 0.6),
}

LIT = {
    "bile_volume_ml_day_statpearls": 600.0,
    "bile_volume_hepatocyte_frac": 0.75,
    "bile_volume_cholangiocyte_frac": 0.25,
    "bile_volume_bile_salt_dependent_ml_day": 225.0,   # StatPearls, ~half of hepatocyte 450 mL
    "pool_maintenance_synthesis_frac": 0.05,           # StatPearls
    "pool_loss_ceiling_frac": 0.20,                    # StatPearls
    "pool_g_wikipedia": (4.0, 6.0),
    "flux_g_day_wikipedia": (12.0, 18.0),
    "synthesis_g_day_wikipedia_low": 0.3,
    "synthesis_g_day_wikipedia_high": 0.6,
    "reabsorption_frac_wikipedia": 0.95,
    "lifetime_reuse_count_wikipedia": 20.0,
    "gallbladder_concentration_factor": (3.0, 10.0),
    "gallbladder_capacity_ml": 50.0,
    "gallbladder_stored_ml": (30.0, 60.0),
    "cmc_deoxycholate_mM": (2.4, 4.0),
    "cmc_general_species_range_mM": (1.0, 250.0),
    "mw_cholic": 408.57,
    "mw_cdca": 392.57,
    "mw_deoxycholic": 392.58,
    "mw_taurocholic": 515.71,
    "gbef_ziessman_lower_normal_pct": {"3min": 16.8, "45min_of_60": 31.0, "60min_of_60": 41.0},
    "gbef_ziessman_below35_healthy_frac_3min": 6.0 / 20.0,
    "gbef_krishnamurthy_mean_sd_pct": {"3min": (66.0, 20.5), "10min": (73.9, 17.7)},
    "gbef_krishnamurthy_cutoff_pct": {"3min": 35.0, "10min": 50.0},
    "si_transit_mean_h": 3.6,
    "liquid_gastric_t_half_h": 15.0 / 60.0,
}


def sec_A_flux_conservation_grid():
    """Phi = M * n is the conserved invariant (Northfield/Hofmann 1975 structural finding).
    Sweep the task's (M, n) band and check what fraction reproduces the task's Phi band --
    a void-floor check that the two independently-stated bands are not mutually exclusive, and
    report the exact stated-central-point evaluation."""
    M = np.linspace(*TASK_BAND["pool_g"], 41)
    n = np.linspace(*TASK_BAND["cycling_per_day"], 41)
    Mg, ng = np.meshgrid(M, n)
    Phi = Mg * ng
    lo, hi = TASK_BAND["flux_g_day"]
    in_band = (Phi >= lo) & (Phi <= hi)
    frac_in_band = float(np.mean(in_band))
    central_Phi = TASK_BAND["pool_g_central"] * TASK_BAND["cycling_per_day_central"]
    return {
        "grid_shape": list(Phi.shape),
        "phi_min": float(Phi.min()), "phi_max": float(Phi.max()),
        "frac_grid_in_task_flux_band": frac_in_band,
        "central_point_M3_n6_phi": central_Phi,
        "central_point_in_band": bool(lo <= central_Phi <= hi),
        "wikipedia_flux_band_independent_match": LIT["flux_g_day_wikipedia"] == TASK_BAND["flux_g_day"],
        "gate_nontrivial_overlap": bool(0.0 < frac_in_band < 1.0),
    }


def sec_B_reabsorption_closure_and_lifetime_cycles():
    """THE core falsifier computation: eta_implied = 1 - L/Phi from two INDEPENDENTLY-stated task
    bands (L=fecal loss, Phi=flux), checked against the independently-stated ~95% anchor -- then the
    geometric-distribution lifetime-cycle-count 1/(1-eta), checked against Wikipedia's independently
    -stated 'reused about 20 times' (a DIFFERENT quantity from per-day cycling frequency n)."""
    Lmin, Lmax = TASK_BAND["fecal_loss_g_day"]
    Phimin, Phimax = TASK_BAND["flux_g_day"]
    eta_low = 1.0 - Lmax / Phimin      # worst case for eta: max loss over min flux
    eta_high = 1.0 - Lmin / Phimax     # best case for eta: min loss over max flux
    anchor = TASK_BAND["reabsorption_frac"]
    tol = 0.01
    eta_anchor_within_range = (eta_low - tol) <= anchor <= (eta_high + tol)

    lifetime_cycles_at_anchor = 1.0 / (1.0 - anchor)
    wiki_lifetime = LIT["lifetime_reuse_count_wikipedia"]
    lifetime_abs_err = abs(lifetime_cycles_at_anchor - wiki_lifetime)

    return {
        "eta_implied_range": [eta_low, eta_high],
        "reabsorption_anchor_95pct_within_implied_range_tol1pct": eta_anchor_within_range,
        "geometric_lifetime_cycles_at_eta_0.95": lifetime_cycles_at_anchor,
        "wikipedia_independent_lifetime_reuse_count": wiki_lifetime,
        "lifetime_cycles_abs_err_vs_wikipedia": lifetime_abs_err,
        "gate_lifetime_cycles_matches_wikipedia_tol1": bool(lifetime_abs_err <= 1.0),
        "note": "eta = per-pass reabsorption fraction (dimensionless); lifetime_cycles = expected "
                "number of enterohepatic passes ONE bile-acid molecule survives before fecal loss "
                "(geometric-distribution mean = 1/(1-eta)) -- NOT the per-day temporal cycling "
                "frequency n (Section A). Both are simultaneously true and are different quantities.",
    }


def sec_C_naive_no_recycling_adversary():
    """FORCED ADVERSARY (the central one this doc must not skip): a naive model with NO
    enterohepatic recycling would require hepatic synthesis alone to supply the full daily flux Phi.
    Measured synthesis capacity (task's fecal-loss/synthesis band) is 20-90x too small -- the
    mass-conservation violation the recycling loop resolves."""
    Phimin, Phimax = TASK_BAND["flux_g_day"]
    Smin, Smax = TASK_BAND["fecal_loss_g_day"]
    ratio_low = Phimin / Smax     # smallest plausible ratio (worst case for the adversary)
    ratio_high = Phimax / Smin    # largest plausible ratio
    pre_registered_min_ratio = 10.0
    return {
        "phi_band_g_day": [Phimin, Phimax],
        "synthesis_capacity_band_g_day": [Smin, Smax],
        "naive_no_recycling_ratio_range": [ratio_low, ratio_high],
        "pre_registered_min_ratio": pre_registered_min_ratio,
        "gate_naive_adversary_falsified": bool(ratio_low >= pre_registered_min_ratio),
        "interpretation": "Without enterohepatic recycling, the liver would need to synthesize "
                          f"{ratio_low:.0f}-{ratio_high:.0f}x more bile acid mass per day than its "
                          "measured synthesis capacity ever provides -- recycling a small pool "
                          "several times a day, not fresh synthesis, supplies the daily flux.",
    }


def sec_D_hepatic_bile_acid_concentration_derived():
    """DERIVED (not directly measured): hepatic bile-acid-DEPENDENT-flow concentration
    = Phi / V_bile_salt_dependent, converted to mM via a molecular-weight range spanning the 3
    live-verified unconjugated bile acids (light bound) to taurocholate (heavy/conjugated bound).
    Disclosed as a model-internal derived quantity, not an independently-measured citation."""
    Phimin, Phimax = TASK_BAND["flux_g_day"]
    V_bad_L_day = LIT["bile_volume_ml_day_statpearls"] and (
        LIT["bile_volume_bile_salt_dependent_ml_day"] / 1000.0)
    conc_g_L_range = [Phimin / V_bad_L_day, Phimax / V_bad_L_day]

    mw_unconj_avg = float(np.mean([LIT["mw_cholic"], LIT["mw_cdca"], LIT["mw_deoxycholic"]]))
    mw_conj_heavy = LIT["mw_taurocholic"]

    # lowest mM: lowest g/L over heaviest MW; highest mM: highest g/L over lightest MW
    mM_low = conc_g_L_range[0] / mw_conj_heavy * 1000.0
    mM_high = conc_g_L_range[1] / mw_unconj_avg * 1000.0

    gb_factor_lo, gb_factor_hi = LIT["gallbladder_concentration_factor"]
    gb_mM_low = mM_low * gb_factor_lo
    gb_mM_high = mM_high * gb_factor_hi

    return {
        "V_bile_salt_dependent_L_day": V_bad_L_day,
        "conc_g_per_L_range": conc_g_L_range,
        "mw_unconjugated_avg_g_mol": mw_unconj_avg,
        "mw_conjugated_heavy_g_mol": mw_conj_heavy,
        "hepatic_bile_acid_dependent_conc_mM_range": [mM_low, mM_high],
        "gallbladder_naive_multiplied_conc_mM_range": [gb_mM_low, gb_mM_high],
        "gallbladder_naive_range_implausible_above_1000mM": bool(gb_mM_high > 1000.0),
        "disclosed_modeling_limitation": (
            "The naive gallbladder-multiplied upper bound exceeds ~1000 mM, physically implausible "
            "(would precipitate) -- flagged as a temporal-averaging mismatch: the 3-10x gallbladder "
            "factor (Wikipedia) most plausibly applies to a LOWER interdigestive/fasting-state "
            "baseline concentration, not this 24h-flux-averaged (fed-state-inclusive) hepatic "
            "estimate. The robust claim (Section E) uses ONLY the hepatic bile-acid-dependent-flow "
            "concentration, not this double-counted product."
        ),
    }


def sec_E_cmc_threshold_check(hepatic_conc_mM_range):
    """Regime-robust threshold check: does the derived physiological bile-acid concentration clear
    the ONE live-verified specific-species CMC anchor (deoxycholic acid) by a comfortable
    (pre-registered >=10x), not knife-edge, margin?"""
    cmc_lo, cmc_hi = LIT["cmc_deoxycholate_mM"]
    worst_case_conc = hepatic_conc_mM_range[0]     # lowest derived concentration
    worst_case_cmc = cmc_hi                        # highest CMC (hardest to clear)
    margin = worst_case_conc / worst_case_cmc
    pre_registered_min_margin = 10.0
    return {
        "cmc_deoxycholate_mM_band": [cmc_lo, cmc_hi],
        "cmc_general_species_range_mM_roda1983": list(LIT["cmc_general_species_range_mM"]),
        "worst_case_hepatic_conc_mM": worst_case_conc,
        "worst_case_cmc_mM": worst_case_cmc,
        "margin_ratio": margin,
        "pre_registered_min_margin": pre_registered_min_margin,
        "gate_cmc_cleared_with_robust_margin": bool(margin >= pre_registered_min_margin),
        "open_regime_not_resolved": (
            "Fasting/interdigestive SMALL-INTESTINAL LUMINAL bile-acid concentration (as opposed to "
            "biliary concentration) was not independently live-verified -- could "
            "plausibly approach or fall below CMC between meals; held OPEN per task instruction. "
            "Cholate/CDCA/taurocholate-specific CMC values were also not independently live-"
            "verified (disclosed gap) -- the qualitative structure rule (trihydroxy > "
            "dihydroxy CMC; Na+ lowers CMC; conjugation ~no effect), verified verbatim from Roda "
            "1983, predicts cholate sits ABOVE deoxycholate's 2.4-4 mM, direction only, not a "
            "specific re-verified number."
        ),
    }


def sec_F_hepatic_first_pass_coupling():
    """Reuse (read-only) the hepatic_clearance cell's certified well-stirred model machinery/Q_H --
    qualitative coupling only (no verified numeric bile-acid extraction ratio)."""
    if not HEPATIC_JSON.exists():
        return {"available": False, "reason": f"{HEPATIC_JSON} not found -- coupling skipped, "
                                                "not fabricated."}
    hep = json.loads(HEPATIC_JSON.read_text())
    q_h = hep["hepatic_flow"]["q_h_rest_l_min"]
    return {
        "available": True,
        "Q_H_rest_L_min_reused_readonly": q_h,
        "source": "hepatic_clearance_results.json",
        "qualitative_argument": (
            "the hepatic_clearance cell's well-stirred model: CL_H = Q*x/(Q+x), x=fu*CL_int. Hofmann "
            "2009's live-verified qualitative statement ('efficient hepatic clearance results in "
            "low plasma bile acid levels') is exactly the flow-limited (high-ER, elasticity->0) "
            "regime that model already certifies for ICG/propranolol (the hepatic_clearance "
            "cell, Section 4). No verified NUMERIC bile-acid extraction ratio was found "
            "(esearch for 'hepatic extraction ratio bile acids portal blood first "
            "pass' returned 0 PubMed hits) -- this coupling stays QUALITATIVE, disclosed, not "
            "quantitatively certified."
        ),
    }


def sec_G_si_transit_decorrelated_ceiling():
    """Decorrelated cross-check: an independent physics (small-bowel transit-time kinetics, from
    the gi_absorption_transit cell, read-only) implies a SERIAL (no-pipelining) ceiling on cycling
    frequency n -- compare against the literature n-band and the flux/pool-ratio-derived n."""
    si_h = LIT["si_transit_mean_h"]
    gastric_h = LIT["liquid_gastric_t_half_h"]
    used_repo_json = False
    if GI_JSON.exists():
        gi = json.loads(GI_JSON.read_text())
        si_h = gi["si_transit"]["T_mean_h"]
        gastric_h = gi["liquid_gastric_emptying"]["t_half_min_central"] / 60.0
        used_repo_json = True
    floor_cycle_h = si_h + gastric_h
    max_n_serial = 24.0 / floor_cycle_h

    Phi_central = float(np.mean(TASK_BAND["flux_g_day"]))
    n_from_flux_pool = Phi_central / TASK_BAND["pool_g_central"]

    lo, hi = TASK_BAND["cycling_per_day"]
    return {
        "si_transit_h": si_h, "liquid_gastric_t_half_h": gastric_h,
        "used_repo_gi_json": used_repo_json,
        "floor_cycle_time_h_serial_no_pipelining": floor_cycle_h,
        "max_n_serial_per_day": max_n_serial,
        "n_from_flux_over_pool_M3_phi_mean": n_from_flux_pool,
        "task_n_band": [lo, hi],
        "gate_serial_ceiling_admits_task_lower_band": bool(max_n_serial >= lo),
        "reaches_task_upper_band_without_pipelining": bool(max_n_serial >= hi),
        "interpretation": (
            f"A strictly-serial (one batch in transit at a time) floor allows at most "
            f"{max_n_serial:.2f} cycles/day -- consistent with the LOWER part of the task's n-band "
            f"({lo}-{hi}/day) but not the upper part without PIPELINING (multiple aliquots of the "
            "pool simultaneously at different transit stages, physiologically plausible since 3 "
            "meals/day each trigger a new CCK-mediated ejection while a prior aliquot may still be "
            "in transit) -- reported as a genuine, disclosed nuance, not forced to a clean pass."
        ),
    }


def sec_H_gallbladder_ejection_cross_study():
    """Cross-study consistency: two DECORRELATED cholescintigraphy protocols (different CCK-8
    infusion durations/doses, different research groups) should agree in order of magnitude on the
    lower-normal-range of GBEF, computed via mean-2SD from Krishnamurthy vs directly stated by
    Ziessman."""
    z = LIT["gbef_ziessman_lower_normal_pct"]
    k = LIT["gbef_krishnamurthy_mean_sd_pct"]
    k_lower_3min = k["3min"][0] - 2 * k["3min"][1]
    k_lower_10min = k["10min"][0] - 2 * k["10min"][1]

    ziessman_band = [z["45min_of_60"], z["60min_of_60"]]     # the validated (not 3-min) protocol
    krishnamurthy_band = sorted([k_lower_3min, k_lower_10min])

    def overlap(a, b):
        lo = max(a[0], b[0]); hi = min(a[1], b[1])
        return lo <= hi

    return {
        "ziessman_lower_normal_pct": z,
        "krishnamurthy_mean_sd_pct": k,
        "krishnamurthy_mean_minus_2sd_pct": {"3min": k_lower_3min, "10min": k_lower_10min},
        "ziessman_validated_band_pct": ziessman_band,
        "krishnamurthy_derived_band_pct": krishnamurthy_band,
        "gate_cross_study_bands_overlap": bool(overlap(ziessman_band, krishnamurthy_band)),
        "krishnamurthy_normal_cutoffs_pct": LIT["gbef_krishnamurthy_cutoff_pct"],
        "gallbladder_concentration_factor": list(LIT["gallbladder_concentration_factor"]),
        "gallbladder_stored_ml": list(LIT["gallbladder_stored_ml"]),
    }


def sec_I_steady_state_synthesis_vs_loss():
    """Conservation-law check: at steady state S (synthesis) = L (fecal loss). Cross-check
    StatPearls' independently-stated FRACTIONAL rules (5% pool-maintenance-by-synthesis, 20% pool
    loss-ceiling) against the task's ABSOLUTE fecal-loss band, at the task's central pool size."""
    M = TASK_BAND["pool_g_central"]
    S_5pct = LIT["pool_maintenance_synthesis_frac"] * M
    ceiling_20pct = LIT["pool_loss_ceiling_frac"] * M
    Lmin, Lmax = TASK_BAND["fecal_loss_g_day"]
    tol = 0.05
    return {
        "pool_g_central": M,
        "S_5pct_of_pool_g_day": S_5pct,
        "ceiling_20pct_of_pool_g_day": ceiling_20pct,
        "task_fecal_loss_band_g_day": [Lmin, Lmax],
        "ceiling_matches_task_upper_bound": bool(abs(ceiling_20pct - Lmax) <= tol),
        "S_5pct_vs_task_lower_bound_relerr": abs(S_5pct - Lmin) / Lmin,
        "wikipedia_synthesis_internal_inconsistency_g_day": [
            LIT["synthesis_g_day_wikipedia_low"], LIT["synthesis_g_day_wikipedia_high"]],
        "wikipedia_internal_inconsistency_ratio": (
            LIT["synthesis_g_day_wikipedia_high"] / LIT["synthesis_g_day_wikipedia_low"]),
    }


def sec_J_void_floor_no_reabsorption():
    """Void-floor / necessity check: if reabsorption were zero, how fast would the observed pool be
    lost, given the observed flux? Confirms the recycling mechanism is NECESSARY, not just a
    convenient description."""
    M = TASK_BAND["pool_g_central"]
    Phi_central = float(np.mean(TASK_BAND["flux_g_day"]))
    depletion_days = M / Phi_central
    depletion_h = depletion_days * 24.0
    pre_registered_max_h = 12.0
    return {
        "pool_g_central": M, "phi_g_day_central": Phi_central,
        "zero_reabsorption_depletion_h": depletion_h,
        "pre_registered_max_h": pre_registered_max_h,
        "gate_mechanism_necessary": bool(depletion_h < pre_registered_max_h),
    }


def sec_K_statpearls_internal_arithmetic_check():
    """Does the live-fetched primary tertiary source's internal numbers add up? (a cheap but
    real check -- Wikipedia's synthesis figure failed this exact test, Section I)."""
    total = LIT["bile_volume_ml_day_statpearls"]
    hepatocyte = total * LIT["bile_volume_hepatocyte_frac"]
    cholangiocyte = total * LIT["bile_volume_cholangiocyte_frac"]
    bad = hepatocyte * 0.5
    stated_bad = LIT["bile_volume_bile_salt_dependent_ml_day"]
    return {
        "total_ml_day": total, "hepatocyte_ml_day": hepatocyte,
        "cholangiocyte_ml_day": cholangiocyte,
        "computed_bile_salt_dependent_ml_day": bad,
        "statpearls_stated_bile_salt_dependent_ml_day": stated_bad,
        "gate_internally_self_consistent": bool(
            abs(hepatocyte + cholangiocyte - total) < 1e-9 and abs(bad - stated_bad) < 1e-9),
    }


def sec_pool_cross_source_discrepancy_OPEN():
    """Held explicitly OPEN per task instruction -- NOT a pass/fail gate."""
    a = TASK_BAND["pool_g"]
    b = LIT["pool_g_wikipedia"]
    lo = max(a[0], b[0]); hi = min(a[1], b[1])
    overlap_width = max(0.0, hi - lo)
    return {
        "task_pool_g_band": list(a),
        "wikipedia_pool_g_band": list(b),
        "overlap_width_g": overlap_width,
        "status": "OPEN -- NOT forced to a single number (task instruction). Method-dependent "
                  "(isotope dilution vs balance) real, quantified, live-verified inter-source "
                  "discrepancy; the two bands touch only at the single point 4 g.",
    }


def main():
    report = {"citations": CITATIONS, "task_band": TASK_BAND, "lit_params": {
        k: (list(v) if isinstance(v, tuple) else v) for k, v in LIT.items()}}

    secA = sec_A_flux_conservation_grid()
    secB = sec_B_reabsorption_closure_and_lifetime_cycles()
    secC = sec_C_naive_no_recycling_adversary()
    secD = sec_D_hepatic_bile_acid_concentration_derived()
    secE = sec_E_cmc_threshold_check(secD["hepatic_bile_acid_dependent_conc_mM_range"])
    secF = sec_F_hepatic_first_pass_coupling()
    secG = sec_G_si_transit_decorrelated_ceiling()
    secH = sec_H_gallbladder_ejection_cross_study()
    secI = sec_I_steady_state_synthesis_vs_loss()
    secJ = sec_J_void_floor_no_reabsorption()
    secK = sec_K_statpearls_internal_arithmetic_check()
    secOpen = sec_pool_cross_source_discrepancy_OPEN()

    report.update({
        "A_flux_conservation_grid": secA,
        "B_reabsorption_closure_and_lifetime_cycles": secB,
        "C_naive_no_recycling_adversary": secC,
        "D_hepatic_bile_acid_concentration_derived": secD,
        "E_cmc_threshold_check": secE,
        "F_hepatic_first_pass_coupling": secF,
        "G_si_transit_decorrelated_ceiling": secG,
        "H_gallbladder_ejection_cross_study": secH,
        "I_steady_state_synthesis_vs_loss": secI,
        "J_void_floor_no_reabsorption": secJ,
        "K_statpearls_internal_arithmetic_check": secK,
        "pool_cross_source_discrepancy_OPEN_not_gated": secOpen,
    })

    gates_all = {
        "A_nontrivial_flux_band_overlap": secA["gate_nontrivial_overlap"],
        "A_central_point_in_flux_band": secA["central_point_in_band"],
        "B_eta95_within_mass_balance_implied_range": secB["reabsorption_anchor_95pct_within_implied_range_tol1pct"],
        "B_lifetime_cycles_matches_wikipedia_20x": secB["gate_lifetime_cycles_matches_wikipedia_tol1"],
        "C_naive_adversary_falsified_ge10x": secC["gate_naive_adversary_falsified"],
        "E_cmc_cleared_robust_margin_ge10x": secE["gate_cmc_cleared_with_robust_margin"],
        "G_serial_ceiling_admits_task_lower_band": secG["gate_serial_ceiling_admits_task_lower_band"],
        "H_cross_study_gbef_bands_overlap": secH["gate_cross_study_bands_overlap"],
        "I_ceiling20pct_matches_task_upper_loss_bound": secI["ceiling_matches_task_upper_bound"],
        "J_zero_reabsorption_mechanism_necessary": secJ["gate_mechanism_necessary"],
        "K_statpearls_internally_self_consistent": secK["gate_internally_self_consistent"],
    }
    n_pass = sum(1 for v in gates_all.values() if v)
    n_total = len(gates_all)
    report["gates_all"] = gates_all
    report["gates_summary"] = {"n_pass": n_pass, "n_total": n_total}
    report["overall_pass"] = bool(n_pass == n_total)

    report["falsifier_verdict"] = {
        "mass_balance_closure": (
            f"Phi=M*n reproduces the task's [12,18] g/day band at the stated central point "
            f"(3g x 6/day = {secA['central_point_M3_n6_phi']:.1f} g/day) and independently matches "
            f"Wikipedia's [12,18] g/day figure exactly -- PASS."
        ),
        "reabsorption_efficiency_95pct": (
            f"eta implied by two INDEPENDENT task bands (fecal loss, flux) = "
            f"[{secB['eta_implied_range'][0]:.4f}, {secB['eta_implied_range'][1]:.4f}], bracketing "
            f"the independently-stated ~95% anchor at its LOWER edge exactly -- PASS, not a "
            f"tautology (the anchor was not used to derive the range)."
        ),
        "fecal_loss_replaced_by_synthesis": (
            "Steady-state S=L conservation law holds by construction; StatPearls' independent "
            "20%-of-pool ceiling rule reproduces the task's 0.6 g/day upper bound EXACTLY at "
            "pool=3g -- PASS. The 5%-of-pool rule (0.15 g/day) sits ~25% below the task's lower "
            "bound (0.2 g/day) -- a real, disclosed, non-exact discrepancy, not forced to match."
        ),
        "verdict": "PASS" if report["overall_pass"] else "PARTIAL",
    }

    report["honest_gaps"] = [
        "Pool size is a genuine, live-verified, method-dependent OPEN discrepancy: task/central-"
        "literature band [2,4] g (central 3 g) vs Wikipedia's independently-stated [4,6] g -- the "
        "two intervals touch only at 4 g. Held OPEN per task instruction, not forced to one number.",
        "No live-verified specific numeric CMC value was found for cholate, "
        "chenodeoxycholate, taurocholate, or glycocholate individually -- only deoxycholic acid "
        "(2.4-4 mM, Wikipedia) and the general >50-species range (1->250 mM, Roda 1983 abstract). "
        "The qualitative structure rule (trihydroxy > dihydroxy CMC) is verified; specific numbers "
        "for the other 3 major physiological species are not (Section E).",
        "The derived hepatic bile-acid concentration (Section D) is a MODEL-INTERNAL calculation "
        "(Phi / bile-salt-dependent volume / MW), not an independently-measured citation "
        "-- its absolute scale was not cross-checked against a live-fetched direct "
        "biliary-bile-acid-concentration measurement paper. The naive gallbladder-multiplied upper "
        "bound (>1000 mM) is flagged as physically implausible / a temporal-averaging modeling "
        "artifact, not a certified prediction (disclosed in Section D itself, not hidden).",
        "Hepatic first-pass extraction ratio for bile acids specifically was NOT found "
        "(0 PubMed hits for a direct search) -- the hepatic_clearance coupling (Section "
        "F) stays qualitative (flow-limited regime argument from Hofmann 2009's text), not a "
        "quantitatively certified CL_H number for bile acids.",
        "Isotope-dilution PRIMARY sources verified (Northfield/Hofmann 1975, McCormick/"
        "Vlahcevic 1973) confirm the STRUCTURE (inverse M-n covariance; isotope-kinetic methodology; "
        "disease-reduced synthesis of the right order of magnitude) but their own abstracts do NOT "
        "give the modern healthy-control point-summary numbers (pool g, %reabsorption) directly -- "
        "those point values trace most directly to StatPearls/Wikipedia tertiary restatement of the "
        "classical literature, one tier below primary-number re-extraction (disclosed, "
        "same tiering convention as the hepatic_clearance cell's antipyrine ER/fu gap).",
        "Gallbladder ejection fraction cross-study check (Section H) compares DIFFERENT CCK-8 "
        "infusion protocols (3-min/10-min bolus-style vs the 45/60-min sincalide-infusion Ziessman "
        "protocol) -- a real, disclosed method spread, not a single clean number.",
        "SI-transit decorrelated ceiling (Section G) is a SERIAL (no-pipelining) floor; real "
        "physiology likely pipelines multiple pool aliquots simultaneously (3 meals/day), so the "
        "ceiling under-predicts the upper cycling-frequency band (12/day) -- disclosed "
        "as a modeling simplification, not silently forced to match.",
        "No individual data anywhere in this cell -- a population-parametrized mass-balance "
        "consistency/falsifier check, matching every other systemic layer's disclosed scope "
        "(glucose_insulin_minimal_model, renal_filtration, hepatic_clearance, "
        "gi_absorption_transit).",
    ]

    report["couples_to_PROSE_ONLY_not_graph_edited"] = [
        "the hepatic_clearance cell: bile acid synthesis (CYP7A1) is a "
        "hepatocyte function upstream of the well-stirred clearance model; Q_H reused read-only "
        "(Section F) for a qualitative flow-limited-clearance coupling argument -- NOT a "
        "quantitatively certified numeric extraction ratio (disclosed gap).",
        "the gi_absorption_transit cell: fat/micelle absorption is the downstream physiological "
        "consumer of the bile-acid delivery this cell models; that cell's "
        "model (SGLT1/GLUT2 glucose absorption only) does not yet include a lipid/micellar "
        "absorption pathway -- this cell builds the micellar-solubilization layer "
        "(Section E), a natural future extension for that cell to consume. SI transit time reused "
        "read-only (Section G) for a decorrelated cycling-frequency ceiling check.",
        "Bile-acid FXR/TGR5 signaling: a DIFFERENT, disease-state question (bile-acid species "
        "composition + receptor signaling tone in PBC/MASH cohorts) -- no overlap with this cell's "
        "healthy-baseline mass-balance layer, but this cell's pool/cycling/reabsorption numbers are "
        "exactly the healthy reference such disease-state deviations would be measured against.",
    ]

    report["confidence_tier"] = (
        "in-vivo-anchored (isotope-dilution pool kinetics), TIERED: (1) two classical PRIMARY human "
        "isotope-dilution/duodenal-perfusion papers verified (Northfield & Hofmann "
        "1975 Gut, n=7+7; McCormick/Vlahcevic 1973 Gut) confirm the STRUCTURAL claim (Phi=M*n "
        "invariant; isotope-kinetic methodology; order-of-magnitude disease-reduced synthesis rates) "
        "but not the modern healthy-control point-summary numbers directly from their own abstracts; "
        "(2) those point numbers (pool 2-4g, ~95% reabsorption, 12-18 g/day flux, 0.2-0.6 g/day "
        "fecal loss) trace most directly to StatPearls (NBK470209) and Wikipedia tertiary "
        "restatement -- textbook-grade, flagged, one tier below primary-number "
        "re-extraction; (3) gallbladder ejection fraction is anchored to TWO independent, directly-"
        "measured nuclear-cholescintigraphy primary studies (Ziessman 2001 n=20; Krishnamurthy 2004) "
        "-- a stronger, same-tier-as-the-hepatic_clearance-ICG anchor; (4) CMC is anchored to ONE "
        "primary biophysical-chemistry paper (Roda/Hofmann/Mysels 1983, >50 species) for the general "
        "structure and range, but only ONE specific species (deoxycholate) numerically, not all 4 "
        "major physiological bile acids. Population-parametrized mass-balance model throughout -- "
        "no subject-specific data, matching this repo's other systemic-layer scope."
    )

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"Wrote {OUT_PATH}")
    print(json.dumps({"overall_pass": report["overall_pass"],
                       "gates_summary": report["gates_summary"],
                       "falsifier_verdict": report["falsifier_verdict"]["verdict"]}, indent=2))
    return 0 if report["overall_pass"] else 2


if __name__ == "__main__":
    sys.exit(main())
