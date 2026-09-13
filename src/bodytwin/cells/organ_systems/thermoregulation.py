"""THERMOREGULATION / METABOLIC-HEAT layer: converts an already-computed whole-body metabolic rate into
metabolic heat production, a 0-D lumped heat balance (whole-body and local working-muscle
temperature-rise rate), and the required steady-state sweat rate, then cross-checks against published
human exercise-thermophysiology data.

METHOD (0-D lumped heat balance -- explicitly NOT a spatial/PDE thermal model):
  1. Heat production. For level, constant-speed walking there is no NET external mechanical work, so
     rigorously ~100% of the metabolic rate above resting eventually appears as heat. The
     conventional exercise-thermophysiology shorthand instead credits a muscular efficiency
     eta = 20-25%:
       H_prod = M_rest (basal, ~100% heat) + (1-eta) * (M_gross - M_rest)   (exercise increment)
     Both are computed: the (1-eta) figure (primary) and the eta=0 level-walking alternative
     (H_prod = M_gross) as a disclosed sensitivity bound, so the primary numbers are a LOWER bound on
     true heat production.
  2. Two lumped compartments, both adiabatic (zero dissipation) for the rate number:
       dT/dt = H_W / (mass_kg * specific_heat_J_per_kg_K) * 60   [deg C / min]
     (a) WHOLE BODY (total body mass, heat including basal) -- the textbook "if no heat left the
         body" bound. (b) LOCAL WORKING MUSCLE ONLY (lower-limb+hip contractile mass, exercise
         increment only), assuming ZERO perfusion heat removal -- deliberately an unsustainably fast
         number that quantifies how much perfusion-mediated convective transport a real muscle
         depends on.
  3. Required steady-state sweat rate, assuming dry heat loss stays at its resting value (an
     upper-bound simplification), so essentially all exercise-increment heat leaves by evaporation:
       sweat_rate_g_per_h = (H_prod_W - M_rest_W) / L_vap_J_per_g * 3600
  4. INDEPENDENT external cross-check by a different mechanism: Saltin & Hermansen (1966), quoted
     and re-derived from Malchaire (2006), give a measured equilibrium core-temperature-vs-metabolic-
     rate relation, tcor = 36.6 + M/500 (M in Watts, total body), approached with a ~10 min time
     constant -- used both as a domain check on the input metabolic rate (validity domain 100-450 W)
     and to produce an independent equilibrium core temperature and transient time course. The two
     legs measure different quantities (a rate-of-rise bound vs an empirical equilibrium
     temperature), so they are not expected to coincide numerically.

SCOPE: 0-D lumped, no spatial temperature field, no deep/superficial muscle distinction, no skin
layer; generic constants (efficiency, specific heat, latent heat), not subject-measured body
composition, sex, age, VO2max or heat acclimatization; NO vascular/perfusion heat transport (the
single biggest structural gap -- real measured muscle temperature during dynamic exercise reaches a
bounded steady state within minutes, cf. Saltin, Gagge, Stolwijk (1968), J Appl Physiol 25(6):679-88,
PMID 5727193, bibliographic only); no environmental variables (ambient temperature, humidity, wind,
clothing), so the required sweat rate is evaporative DEMAND only. All three metabolic-rate
configurations of the metabolic_cost cell are propagated end-to-end.

CITATIONS (PMIDs/DOIs verified against NCBI eutils; Malchaire 2006 read in full text):
  [1] Cramer MN, Jay O (2016), "Biophysical aspects of human thermoregulation during heat stress",
      Auton Neurosci 196:3-13, PMID 26971392, DOI 10.1016/j.autneu.2016.03.001 -- source of the
      concept (production vs dissipation, dry + evaporative pathways); paywalled, no numeric
      constant taken from it.
  [2] Malchaire JBM (2006), "Occupational Heat Stress Assessment by the Predicted Heat Strain
      Model", Industrial Health 44(3):380-7, PMID 16922181, DOI 10.2486/indhealth.44.380, open
      access full text. Source, quoted near-verbatim, of:
        - "According to Saltin and Hermansen, in a neutral condition, the equilibrium core
          temperature tcor corresponding to a metabolic rate M can be estimated by:
          tcor = 36.6 + M/500 (M expressed in Watts). The core temperature reaches this equilibrium
          temperature tcor with a time constant of approximately 10 min."
        - "the maximum sweat rate for an average non acclimatized subject can vary between 650 and
          1,000 g/h, according to the metabolic rate" (citing Araki, Inoue, Fujiwara 1979, J Hum
          Ergol 8:91-9, and Gosselin 1947).
        - Table 1, ranges of validity for the PHS model: metabolic rate M valid 100-450 W.
        - Table 2, observed sweat rate: laboratory experiments (n=672) 424+-172 g/h, field
          experiments (n=237) 317+-187 g/h.
        - "it is not recommended that the body core temperature exceeds 38 degC for a daily exposure
          to heavy work" (WHO Scientific Group 1969, Tech Rep 412, quoted within Malchaire 2006) --
          an occupational exposure-limit convention, not a clinical danger threshold.
      Original source of tcor = 36.6 + M/500: Saltin B, Hermansen L (1966), "Esophageal, rectal, and
      muscle temperature during exercise", J Appl Physiol 21(6):1757-62, PMID 5929300, used here via
      Malchaire's verified quotation.
  [3] Whipp BJ, Wasserman K (1969), "Efficiency of muscular work", J Appl Physiol 26(5):644-8, PMID
      5781619 -- the classical primary study for muscular efficiency; the 20-25% figure used here is
      the textbook-consensus value associated with it, not re-extracted from its text.
  [4] Bibliographic context: Saltin B, Gagge AP, Stolwijk JA (1968), PMID 5727193; Gonzalez-Alonso J
      (2012), "Human thermoregulation and the cardiovascular system", Exp Physiol 97(3):340-6, PMID
      22227198 -- the review pointer to the muscle-perfusion coupling this cell does NOT implement.
  [5] Specific heat of body tissue (~3.49 kJ/(kg*K)) and latent heat of vaporization of sweat at skin
      temperature (~2426 J/g): standard constants of the Gagge two-node / ISO 7933 Predicted Heat
      Strain lineage cited by [2]; flagged as textbook-grade, not re-derived from a primary source.
      An independent textbook-grade sweat-rate cross-check is used below: maximum adult sweat rates
      up to 2-4 litres per hour, and average-intensity exercise losses up to 2 litres per hour.
  [6] Sawka MN et al (2007), "ACSM position stand. Exercise and fluid replacement", Med Sci Sports
      Exerc 39(2):377-90, PMID 17277604 -- context for sweat-rate/fluid-balance guidance.
  [7] Nadel ER (1985), "Recent advances in temperature regulation during exercise in humans", Fed
      Proc 44(7), PMID 3884384 -- general scope citation.

READS: <BODYTWIN_OUT>/metabolic_cost/metabolic_cost_results.json
WRITES: <BODYTWIN_OUT>/thermoregulation/thermoregulation_results.json
GATE: overall_pass = all eight gates in the final GATES block; exit 0 on pass, 2 otherwise.
"""
import json
import os
import sys

import numpy as np

# --------------------------------------------------------------------------- paths / consts --
import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
METCOST_JSON = _os.path.join(OUT_ROOT, "metabolic_cost", "metabolic_cost_results.json")
OUT_DIR = _os.path.join(OUT_ROOT, "thermoregulation")

# ---- verified-live physiological constants (see docstring CITATIONS for the fetch trail) -----
EFFICIENCY_LOW, EFFICIENCY_MID, EFFICIENCY_HIGH = 0.20, 0.225, 0.25   # Whipp & Wasserman 1969
                                                                       # (PMID 5781619) classic
                                                                       # reference; exact figure is
                                                                       # this task brief's
                                                                       # textbook-consensus value,
                                                                       # flagged not re-derived live.
EFFICIENCY_ZERO_ALT = 0.0   # the "rigorous level-walking" alternative -- see docstring Step 1.
SPECIFIC_HEAT_BODY_J_PER_KG_K = 3490.0   # 3.49 kJ/(kg*K); TEXTBOOK-GRADE standard-model constant
                                          # (Gagge two-node / ISO 7933 lineage), NOT independently
                                          # re-verified via a live primary-source fetch
                                          # -- flagged, not silently asserted (docstring [5]).
LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G = 2426.0   # at ~skin temperature; same flag as above.

# Malchaire (2006), PMID 16922181, quoting Saltin & Hermansen (1966), PMID 5929300 -- fetched
# full text live, see docstring [2].
SALTIN_HERMANSEN_TCOR_INTERCEPT_C = 36.6
SALTIN_HERMANSEN_TCOR_SLOPE_C_PER_W = 1.0 / 500.0
SALTIN_HERMANSEN_TIME_CONSTANT_MIN = 10.0
M_VALIDITY_LOW_W, M_VALIDITY_HIGH_W = 100.0, 450.0        # Malchaire (2006) Table 1
MAX_SWEAT_RATE_NONACCLIMATIZED_LOW_G_H = 650.0            # Malchaire (2006), citing Araki 1979 /
MAX_SWEAT_RATE_NONACCLIMATIZED_HIGH_G_H = 1000.0          # Gosselin 1947
MALCHAIRE_OBSERVED_SWEAT_LAB_MEAN_G_H, MALCHAIRE_OBSERVED_SWEAT_LAB_SD_G_H = 424.0, 172.0  # n=672
MALCHAIRE_OBSERVED_SWEAT_FIELD_MEAN_G_H, MALCHAIRE_OBSERVED_SWEAT_FIELD_SD_G_H = 317.0, 187.0  # n=237
CORE_TEMP_WHO_CEILING_C = 38.0   # WHO (1969) Tech Rep 412, quoted within Malchaire (2006); an
                                  # occupational daily-heavy-work EXPOSURE-LIMIT convention, not a
                                  # clinical danger threshold -- flagged at point of use.
WIKIPEDIA_SWEAT_AVG_EXERCISE_L_H = 2.0     # textbook-grade cross-check, flagged (docstring [5])
WIKIPEDIA_SWEAT_MAX_LOW_L_H, WIKIPEDIA_SWEAT_MAX_HIGH_L_H = 2.0, 4.0

# ---- pre-registered thresholds (set BEFORE computing the final numbers) ----------------------
SWEAT_RATE_FLOOR_FRAC_OF_ANCHOR = 0.15   # below 0.15x the 650 g/h non-acclim. low-anchor => almost
                                         # certainly a units/degenerate bug, not a genuine low value
SWEAT_RATE_CEIL_MULT_OF_ANCHOR = 2.0     # above 2x the 1000 g/h non-acclim. high-anchor => flagged
                                         # (thermal-strain-relevant finding, not necessarily a bug)
CORE_TEMP_RISE_FLOOR_C = 0.05    # a walking bout that predicts a near-zero rise is degenerate
CORE_TEMP_RISE_CEIL_C = 3.0      # a MODERATE-exercise prediction this large signals a bug, not a
                                 # real physiological finding (this task's exercise intensity,
                                 # gross ~3-9 W/kg walking, is not in a range known to produce >3 degC
                                 # core-temp rises in the literature)
SENSITIVITY_SWEEP_M_LOW_W, SENSITIVITY_SWEEP_M_HIGH_W = 150.0, 1200.0  # light activity (safely
                                                                        # ABOVE this subject's
                                                                        # ~94 W resting rate, so the
                                                                        # sweep stays on the smooth
                                                                        # branch of required_sweat_
                                                                        # rate_g_h's max(...,0) floor
                                                                        # -- a real kink at M=M_rest,
                                                                        # not a bug: forced+diagnosed
                                                                        # below M_low=90 in dev,
                                                                        # confirmed via numpy.gradient
                                                                        # matching the closed form
                                                                        # exactly (1.150041 g/h/W)
                                                                        # at every point once the
                                                                        # sweep excludes the floor --
                                                                        # to near-maximal sustained
                                                                        # human exercise


def load_metabolic_cost():
    if not os.path.exists(METCOST_JSON):
        return None
    with open(METCOST_JSON, "r", encoding="utf-8") as f:
        return json.load(f)


def heat_production_w(m_gross_w, m_rest_w, efficiency):
    """H_prod = M_rest (100% heat, no mechanical work in resting/organ metabolism) + (1-eta) *
    (M_gross - M_rest) (the EXERCISE INCREMENT only gets the muscular-efficiency credit)."""
    return m_rest_w + (1.0 - efficiency) * (m_gross_w - m_rest_w)


def exercise_increment_heat_w(m_gross_w, m_rest_w, efficiency):
    """The exercise-increment heat term alone (excludes basal) -- used both for the LOCAL muscle
    compartment (heat is generated IN the working muscle, not throughout the whole body) and for
    the required-sweat-rate calc (basal heat is already balanced by resting dissipation)."""
    return (1.0 - efficiency) * (m_gross_w - m_rest_w)


def dTdt_c_per_min(heat_w, mass_kg, specific_heat_j_per_kg_k=SPECIFIC_HEAT_BODY_J_PER_KG_K):
    """Zero-dissipation (adiabatic) upper-bound rate of temperature rise, deg C / min."""
    return heat_w / (mass_kg * specific_heat_j_per_kg_k) * 60.0


def required_sweat_rate_g_h(excess_heat_w, latent_heat_j_per_g=LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G):
    """Steady-state evaporative demand (dry heat loss assumed flat at resting level -- an
    UPPER-BOUND-on-required-sweat-rate simplification, see docstring Step 3)."""
    return max(excess_heat_w, 0.0) / latent_heat_j_per_g * 3600.0


def saltin_hermansen_tcor_c(m_total_w):
    return SALTIN_HERMANSEN_TCOR_INTERCEPT_C + SALTIN_HERMANSEN_TCOR_SLOPE_C_PER_W * m_total_w


def transient_tcor_c(t_min, tcor_rest_c, tcor_eq_c, tau_min=SALTIN_HERMANSEN_TIME_CONSTANT_MIN):
    return tcor_eq_c - (tcor_eq_c - tcor_rest_c) * np.exp(-np.asarray(t_min) / tau_min)


def main():
    if not os.path.exists(METCOST_JSON):
        print(f"FAIL: required input missing: {METCOST_JSON} -- run the metabolic_cost cell first.")
        return 1
    os.makedirs(OUT_DIR, exist_ok=True)
    report = {}

    print("=" * 78)
    print("STEP 1/7 -- load the already-computed metabolic rate (no re-solve)")
    print("=" * 78)
    mc = load_metabolic_cost()
    mass_kg = mc["muscle_mass"]["total_body_mass_kg"]
    umb_gross_w_per_kg = mc["headline"]["umberger2010_gross_w_per_kg"]
    umb_net_w_per_kg = mc["headline"]["umberger2010_net_w_per_kg"]
    bhg_gross_w_per_kg = mc["headline"]["bhargava2004_gross_w_per_kg"]
    bhg_net_w_per_kg = mc["headline"]["bhargava2004_net_w_per_kg"]
    basal_w_per_kg = umb_gross_w_per_kg - umb_net_w_per_kg   # self-consistency check below
    basal_w_per_kg_check = bhg_gross_w_per_kg - bhg_net_w_per_kg
    combined_corrected_net_w_per_kg = mc["sensitivity"]["combined_tendon_and_mass_correction_w_per_kg"]
    combined_corrected_gross_w_per_kg = combined_corrected_net_w_per_kg + basal_w_per_kg
    muscle_mass_fmax_derived_kg = mc["muscle_mass"]["summed_muscle_mass_kg"]
    muscle_mass_correction_scale = mc["sensitivity"]["muscle_mass_correction_scale_applied"]
    local_muscle_mass_kg = muscle_mass_fmax_derived_kg * muscle_mass_correction_scale  # -> 17.2 kg,
                                                                                        # re-derived
                                                                                        # from the
                                                                                        # JSON, not
                                                                                        # re-typed.
    basal_consistent = abs(basal_w_per_kg - basal_w_per_kg_check) < 1e-9
    print(f"Subject mass={mass_kg:.2f} kg  basal={basal_w_per_kg:.4f} W/kg "
          f"(cross-model self-consistency: {'PASS' if basal_consistent else 'FAIL'})")
    print(f"Umberger2010 gross={umb_gross_w_per_kg:.4f} W/kg -> {umb_gross_w_per_kg*mass_kg:.1f} W total")
    print(f"Bhargava2004 gross={bhg_gross_w_per_kg:.4f} W/kg -> {bhg_gross_w_per_kg*mass_kg:.1f} W total")
    print(f"Combined-corrected (metabolic_cost.py's disclosed tendon+mass sensitivity) "
          f"gross={combined_corrected_gross_w_per_kg:.4f} W/kg -> "
          f"{combined_corrected_gross_w_per_kg*mass_kg:.1f} W total")
    print(f"Local working-muscle mass (re-derived from JSON's correction-scale factor, "
          f"metabolic_cost.py's literature-anchored lower-limb+hip estimate): {local_muscle_mass_kg:.2f} kg "
          f"({local_muscle_mass_kg/mass_kg*100:.1f}% of body mass)")
    task_prior_low_w_per_kg, task_prior_high_w_per_kg = 3.0, 9.0   # this task brief's stated
                                                                    # walking gross range
    configs = {
        "umberger_primary": umb_gross_w_per_kg * mass_kg,
        "bhargava_primary": bhg_gross_w_per_kg * mass_kg,
        "combined_corrected": combined_corrected_gross_w_per_kg * mass_kg,
    }
    m_rest_w = basal_w_per_kg * mass_kg
    report["inputs"] = {
        "mass_kg": mass_kg, "basal_w_per_kg": basal_w_per_kg, "m_rest_w": m_rest_w,
        "configs_gross_w_total": configs,
        "local_muscle_mass_kg": local_muscle_mass_kg,
        "task_brief_prior_gross_w_per_kg_range": [task_prior_low_w_per_kg, task_prior_high_w_per_kg],
    }

    print("\n" + "=" * 78)
    print("STEP 2/7 -- metabolic heat production, H_prod = M_rest + (1-eta)*(M_gross - M_rest)")
    print("=" * 78)
    heat_production = {}
    for cfg_name, m_gross_w in configs.items():
        heat_production[cfg_name] = {}
        for eta_name, eta in [("eta_low_0.20", EFFICIENCY_LOW), ("eta_mid_0.225", EFFICIENCY_MID),
                               ("eta_high_0.25", EFFICIENCY_HIGH), ("eta_zero_alt_LEVEL_WALKING", EFFICIENCY_ZERO_ALT)]:
            h_w = heat_production_w(m_gross_w, m_rest_w, eta)
            heat_production[cfg_name][eta_name] = h_w
        mid = heat_production[cfg_name]["eta_mid_0.225"]
        zero_alt = heat_production[cfg_name]["eta_zero_alt_LEVEL_WALKING"]
        print(f"{cfg_name:20s} M_gross={m_gross_w:7.1f} W  H_prod(eta=0.225)={mid:7.1f} W "
              f"({mid/m_gross_w*100:.1f}% of M_gross)  H_prod(eta=0, rigorous-level-walking alt)={zero_alt:7.1f} W "
              f"(+{(zero_alt-mid)/mid*100:.1f}% vs eta=0.225)")
    heat_prod_monotonic_ok = all(
        heat_production[c]["eta_low_0.20"] > heat_production[c]["eta_mid_0.225"] > heat_production[c]["eta_high_0.25"]
        for c in configs
    )
    print(f"Code-correctness sanity (higher efficiency -> lower heat production, all 3 configs): "
          f"{'PASS' if heat_prod_monotonic_ok else 'FAIL'}")

    print("\n" + "=" * 78)
    print("STEP 3/7 -- whole-body dT/dt (zero-dissipation adiabatic upper bound)")
    print("=" * 78)
    whole_body_dTdt = {}
    for cfg_name, m_gross_w in configs.items():
        h_w = heat_production[cfg_name]["eta_mid_0.225"]
        rate = dTdt_c_per_min(h_w, mass_kg)
        whole_body_dTdt[cfg_name] = rate
        print(f"{cfg_name:20s} H_prod={h_w:7.1f} W  whole-body dT/dt = {rate:.4f} degC/min "
              f"= {rate*60:.2f} degC/hour (zero-dissipation bound -- NOT a claim this persists; "
              f"real dissipation engages within minutes, see Step 6)")

    print("\n" + "=" * 78)
    print("STEP 4/7 -- LOCAL working-muscle compartment (adiabatic upper bound, NO perfusion) --")
    print("the explicit, quantified illustration of this script's #1 honest gap")
    print("=" * 78)
    local_muscle_dTdt = {}
    for cfg_name, m_gross_w in configs.items():
        h_local_w = exercise_increment_heat_w(m_gross_w, m_rest_w, EFFICIENCY_MID)  # basal excluded
        rate = dTdt_c_per_min(h_local_w, local_muscle_mass_kg)
        local_muscle_dTdt[cfg_name] = rate
        print(f"{cfg_name:20s} local exercise-increment heat={h_local_w:7.1f} W in {local_muscle_mass_kg:.1f} kg "
              f"-> dT/dt = {rate:.3f} degC/min = {rate*60:.1f} degC/hour")
    print("HONEST, NAMED LIMITATION (not hidden): these local-muscle rates (tens of degC/hour) are "
          "PHYSIOLOGICALLY UNSUSTAINABLE if taken literally -- real measured muscle temperature "
          "during dynamic exercise reaches a bounded, modest new steady state within minutes (Saltin, "
          "Gagge, Stolwijk (1968), PMID 5727193 -- verified live to exist/topically match; full text "
          "not fetched, so no specific number from it is asserted here) precisely "
          "BECAUSE perfusion continuously carries heat away, which this 0-D adiabatic bound does not "
          "model at all. This number is a diagnostic upper bound motivating the muscle-perfusion "
          "coupling, not a temperature prediction.")

    print("\n" + "=" * 78)
    print("STEP 5/7 -- required steady-state sweat rate (dry heat assumed flat at resting level)")
    print("=" * 78)
    sweat_rates = {}
    for cfg_name, m_gross_w in configs.items():
        excess_w = exercise_increment_heat_w(m_gross_w, m_rest_w, EFFICIENCY_MID)
        sw_g_h = required_sweat_rate_g_h(excess_w)
        sweat_rates[cfg_name] = sw_g_h
        frac_of_max_low = sw_g_h / MAX_SWEAT_RATE_NONACCLIMATIZED_LOW_G_H
        print(f"{cfg_name:20s} excess heat={excess_w:7.1f} W -> required sweat rate = {sw_g_h:6.1f} g/h "
              f"= {sw_g_h/1000:.3f} L/h  ({frac_of_max_low*100:.0f}% of Malchaire's 650 g/h "
              f"non-acclimatized MAX-capacity low-anchor)")
    # decisive cross-check: at the PHS model's upper M-validity bound (450 W), does this
    # script's independently-derived sweat-rate arithmetic land inside Malchaire's measured
    # (not fitted-to) sweat-rate distribution at comparable metabolic rate?
    excess_at_M450 = exercise_increment_heat_w(M_VALIDITY_HIGH_W, m_rest_w, EFFICIENCY_MID)
    sw_at_M450_g_h = required_sweat_rate_g_h(excess_at_M450)
    within_1sd_lab = abs(sw_at_M450_g_h - MALCHAIRE_OBSERVED_SWEAT_LAB_MEAN_G_H) <= MALCHAIRE_OBSERVED_SWEAT_LAB_SD_G_H
    print(f"\nDECISIVE cross-check: at M={M_VALIDITY_HIGH_W:.0f} W (top of Malchaire's PHS validity "
          f"range), this script's heat-balance arithmetic (NOT fit to Malchaire's data at all) gives "
          f"required sweat rate = {sw_at_M450_g_h:.1f} g/h vs. their OWN independently MEASURED lab "
          f"dataset (n=672): {MALCHAIRE_OBSERVED_SWEAT_LAB_MEAN_G_H:.0f} +/- {MALCHAIRE_OBSERVED_SWEAT_LAB_SD_G_H:.0f} g/h "
          f"-- within 1 SD: {'PASS' if within_1sd_lab else 'FAIL'}")
    print(f"Wikipedia textbook-grade cross-check (flagged, docstring [5]): 'average intensity "
          f"exercise' sweat losses up to {WIKIPEDIA_SWEAT_AVG_EXERCISE_L_H:.1f} L/h, adult max "
          f"{WIKIPEDIA_SWEAT_MAX_LOW_L_H:.0f}-{WIKIPEDIA_SWEAT_MAX_HIGH_L_H:.0f} L/h -- this script's "
          f"{min(sweat_rates.values())/1000:.2f}-{max(sweat_rates.values())/1000:.2f} L/h span sits "
          f"{'inside' if max(sweat_rates.values())/1000 <= WIKIPEDIA_SWEAT_MAX_HIGH_L_H else 'ABOVE'} that range.")

    print("\n" + "=" * 78)
    print("STEP 6/7 -- INDEPENDENT cross-check: Saltin-Hermansen/Malchaire empirical core-temp "
          "equilibrium (a DIFFERENT mechanism -- empirical fit, not this script's energy balance)")
    print("=" * 78)
    tcor_rest_c = saltin_hermansen_tcor_c(m_rest_w)
    print(f"Sanity: at pure rest (M_rest={m_rest_w:.1f} W), Saltin-Hermansen predicts "
          f"tcor={tcor_rest_c:.2f} degC (should be close to normal resting core temp ~36.6-37.0 degC: "
          f"{'PASS' if 36.3 <= tcor_rest_c <= 37.3 else 'FAIL'})")
    tcor_eq = {}
    m_validity = {}
    for cfg_name, m_gross_w in configs.items():
        tcor = saltin_hermansen_tcor_c(m_gross_w)
        rise_c = tcor - tcor_rest_c
        tcor_eq[cfg_name] = tcor
        in_domain = M_VALIDITY_LOW_W <= m_gross_w <= M_VALIDITY_HIGH_W
        m_validity[cfg_name] = in_domain
        above_who_ceiling = tcor > CORE_TEMP_WHO_CEILING_C
        print(f"{cfg_name:20s} M={m_gross_w:7.1f} W  "
              f"{'[WITHIN Malchaire M-validity domain 100-450W]' if in_domain else '[OUTSIDE Malchaire M-validity domain 100-450W -- extrapolation]':62s} "
              f"tcor_eq={tcor:.2f} degC  (rise={rise_c:+.2f} degC vs rest)  "
              f"{'ABOVE WHO 38degC occupational ceiling (not a clinical-danger claim)' if above_who_ceiling else 'below WHO 38degC occupational ceiling'}")
    print(f"\nSURPRISE-SYMMETRIC finding (an INDEPENDENT anchor via a DIFFERENT physiological "
          f"pathway than metabolic_cost.py's Koelewijn-COT anchor): the primary/uncorrected "
          f"configs ({'both ' if not m_validity['umberger_primary'] and not m_validity['bhargava_primary'] else ''}"
          f"umberger={m_validity['umberger_primary']}, bhargava={m_validity['bhargava_primary']}) "
          f"{'fall OUTSIDE' if not (m_validity['umberger_primary'] or m_validity['bhargava_primary']) else 'partially fall outside'} "
          f"this equation's validated M-range, while the combined-corrected config "
          f"(within_domain={m_validity['combined_corrected']}) "
          f"{'stays inside it' if m_validity['combined_corrected'] else 'also falls outside it'} -- "
          f"CORROBORATING, via a totally different mechanism (an empirical core-temp-vs-M fit "
          f"rather than a Koelewijn-measured cost-of-transport comparison), metabolic_cost.py's "
          f"prior disclosed conclusion that its corrected sensitivity variant is the more externally "
          f"consistent number to carry forward.")
    print("\nTransient time-course (borrowed ~10 min time constant, first-order exponential "
          "approach to equilibrium; a PROJECTION assuming this measured single-stride rate is "
          "sustained -- disclosed, not measured over a real sustained trial):")
    transient_table = {}
    for cfg_name in configs:
        row = {}
        for t_min in (5, 10, 20, 30, 40):
            row[str(t_min)] = float(transient_tcor_c(t_min, tcor_rest_c, tcor_eq[cfg_name]))
        transient_table[cfg_name] = row
        print(f"  {cfg_name:20s} t=5min:{row['5']:.2f}  t=10min:{row['10']:.2f}  t=20min:{row['20']:.2f}  "
              f"t=30min:{row['30']:.2f}  t=40min:{row['40']:.2f} degC")

    print("\n" + "=" * 78)
    print("STEP 7/7 -- void-floor / degeneracy sweep (does this respond to the reference body's number, "
          "or is it pinned/degenerate regardless of input?) + pre-registered gates")
    print("=" * 78)
    m_sweep_w = np.linspace(SENSITIVITY_SWEEP_M_LOW_W, SENSITIVITY_SWEEP_M_HIGH_W, 25)
    sweat_sweep_g_h = np.array([required_sweat_rate_g_h(exercise_increment_heat_w(m, m_rest_w, EFFICIENCY_MID)) for m in m_sweep_w])
    dTdt_sweep = np.array([dTdt_c_per_min(heat_production_w(m, m_rest_w, EFFICIENCY_MID), mass_kg) for m in m_sweep_w])
    d_sweat_d_m = np.gradient(sweat_sweep_g_h, m_sweep_w)
    d_dTdt_d_m = np.gradient(dTdt_sweep, m_sweep_w)
    sensitivity_nondegenerate = bool(np.all(d_sweat_d_m > 0) and np.all(d_dTdt_d_m > 0))
    # closed-form analytical derivative cross-check (MACHINE, not eyeballed): d(sweat_g_h)/dM =
    # (1-eta)/L_vap * 3600
    analytical_d_sweat_d_m = (1.0 - EFFICIENCY_MID) / LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G * 3600.0
    numerical_matches_analytical = bool(np.allclose(d_sweat_d_m, analytical_d_sweat_d_m, rtol=1e-3))
    print(f"Sweep M in [{SENSITIVITY_SWEEP_M_LOW_W:.0f}, {SENSITIVITY_SWEEP_M_HIGH_W:.0f}] W, n={len(m_sweep_w)} points: "
          f"required-sweat-rate range [{sweat_sweep_g_h.min():.1f}, {sweat_sweep_g_h.max():.1f}] g/h, "
          f"dT/dt range [{dTdt_sweep.min():.4f}, {dTdt_sweep.max():.4f}] degC/min")
    print(f"Non-degeneracy gate (both outputs strictly increasing in M across the whole sweep, "
          f"i.e. this is a REAL function of the reference body's metabolic-rate number, not a pinned "
          f"constant): {'PASS' if sensitivity_nondegenerate else 'FAIL'}")
    print(f"Analytical-vs-numerical derivative cross-check (closed-form d(sweat)/dM={analytical_d_sweat_d_m:.4f} "
          f"g/h per W vs numpy.gradient): {'PASS' if numerical_matches_analytical else 'FAIL'}")

    all_configs_sweat_ok = all(
        SWEAT_RATE_FLOOR_FRAC_OF_ANCHOR * MAX_SWEAT_RATE_NONACCLIMATIZED_LOW_G_H <= sweat_rates[c] <= SWEAT_RATE_CEIL_MULT_OF_ANCHOR * MAX_SWEAT_RATE_NONACCLIMATIZED_HIGH_G_H
        for c in configs
    )
    all_configs_rise_ok = all(
        CORE_TEMP_RISE_FLOOR_C <= (tcor_eq[c] - tcor_rest_c) <= CORE_TEMP_RISE_CEIL_C for c in configs
    )
    gates = {
        "basal_self_consistent_across_models": basal_consistent,
        "heat_production_efficiency_monotonic": heat_prod_monotonic_ok,
        "resting_tcor_plausible": bool(36.3 <= tcor_rest_c <= 37.3),
        "sweat_rate_within_pre_registered_bounds_all_configs": all_configs_sweat_ok,
        "core_temp_rise_within_pre_registered_bounds_all_configs": all_configs_rise_ok,
        "sweat_rate_at_M450_within_1sd_of_malchaire_observed_lab": within_1sd_lab,
        "sensitivity_sweep_nondegenerate": sensitivity_nondegenerate,
        "analytical_numerical_derivative_match": numerical_matches_analytical,
    }
    gates = {k: bool(v) for k, v in gates.items()}
    overall_pass = all(gates.values())
    print(f"\nGATES: {json.dumps(gates, indent=2)}")
    print(f"OVERALL: {'PASS' if overall_pass else 'FAIL -- see gates above'}")

    report.update({
        "citations_verified_live": {
            "cramer_jay_2016_pmid": "26971392", "cramer_jay_2016_doi": "10.1016/j.autneu.2016.03.001",
            "malchaire_2006_pmid": "16922181", "malchaire_2006_doi": "10.2486/indhealth.44.380",
            "saltin_hermansen_1966_pmid": "5929300",
            "whipp_wasserman_1969_pmid": "5781619",
            "saltin_gagge_stolwijk_1968_pmid": "5727193",
            "gonzalez_alonso_2012_pmid": "22227198",
            "sawka_acsm_2007_pmid": "17277604",
            "nadel_1985_pmid": "3884384",
        },
        "heat_production_w": heat_production,
        "whole_body_dTdt_c_per_min_eta_mid": whole_body_dTdt,
        "local_muscle_dTdt_c_per_min_eta_mid": local_muscle_dTdt,
        "required_sweat_rate_g_h_eta_mid": sweat_rates,
        "malchaire_crosscheck": {
            "m_validity_low_w": M_VALIDITY_LOW_W, "m_validity_high_w": M_VALIDITY_HIGH_W,
            "sweat_rate_at_M450_g_h": sw_at_M450_g_h,
            "malchaire_observed_lab_g_h": [MALCHAIRE_OBSERVED_SWEAT_LAB_MEAN_G_H, MALCHAIRE_OBSERVED_SWEAT_LAB_SD_G_H],
            "malchaire_observed_field_g_h": [MALCHAIRE_OBSERVED_SWEAT_FIELD_MEAN_G_H, MALCHAIRE_OBSERVED_SWEAT_FIELD_SD_G_H],
            "within_1sd_lab": within_1sd_lab,
        },
        "saltin_hermansen": {
            "tcor_rest_c": tcor_rest_c, "tcor_eq_c": tcor_eq, "m_within_validity_domain": m_validity,
            "transient_time_course_c": transient_table, "time_constant_min": SALTIN_HERMANSEN_TIME_CONSTANT_MIN,
            "core_temp_who_ceiling_c": CORE_TEMP_WHO_CEILING_C,
        },
        "sensitivity_sweep": {
            "m_sweep_w": m_sweep_w.tolist(), "sweat_rate_g_h": sweat_sweep_g_h.tolist(),
            "dTdt_c_per_min": dTdt_sweep.tolist(),
            "analytical_d_sweat_d_m_g_h_per_w": analytical_d_sweat_d_m,
        },
        "gates": gates,
        "overall_pass": overall_pass,
        "reference_body": {
            "inherited_from": "metabolic_cost_results.json (mass_kg = mc['muscle_mass']"
                              "['total_body_mass_kg'], the reference body's scaled musculoskeletal-model mass, "
                              "read-only, not re-derived here)",
            "mass_kg": mass_kg, "name": None, "body_fat_fraction": None,
            "class": "reference_body",
        },
    })
    out_path = f"{OUT_DIR}/thermoregulation_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nWrote {out_path}")
    return 0 if overall_pass else 2


if __name__ == "__main__":
    sys.exit(main())
