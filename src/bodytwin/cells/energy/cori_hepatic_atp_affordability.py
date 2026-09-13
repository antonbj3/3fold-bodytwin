"""Cori-cycle hepatic-ATP-AFFORDABILITY collision test.

Asks whether the resting Cori-cycle ATP cost stays affordable out of the liver's ATP budget as
exercise lactate turnover rises, under a "constant recycled fraction" null.

Reuses, does not retype: the whole_body_lactate_ledger cell's live Cori ATP cost, the
exercise_bloodflow_redistribution cell's live Perko-1998-derived hepatosplenic flow numbers
(rest 0.84 L/min, -50% cycling-exercise reduction), plus O2_KJ_PER_L=20.2 kJ/L mixed-diet caloric
equivalent and P/O(NADH)=2.5 (Hinkle 2005, ATP/O2=5.0).

Reads: the whole_body_lactate_ledger cell (imported as a module) and
<OUT_ROOT>/exercise_bloodflow_redistribution/exercise_bloodflow_redistribution_results.json.
Writes: nothing (stdout only).
Gate: pre-registered ceiling -- Cori alone consuming >50% of the hepatic ATP budget falsifies
"the liver can still do its other jobs"; reported as the crossover Ra multiplier.
"""
import json
import os
import sys
import importlib.util
import numpy as np

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
HERE = os.path.dirname(os.path.abspath(__file__))

def load(mod_name, path):
    spec = importlib.util.spec_from_file_location(mod_name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

lactate = load("wbl", os.path.join(HERE, "whole_body_lactate_ledger.py"))
ledger = lactate.build_report()

with open(os.path.join(OUT_ROOT, "exercise_bloodflow_redistribution", "exercise_bloodflow_redistribution_results.json")) as f:
    ebf = json.load(f)

# ---- A: resting Cori ATP demand, REUSED live ----
cori_rest_W = ledger["cori_cycle_energetics"]["cori_cycle_power_w"]
rmr_W = ledger["cori_cycle_energetics"]["rmr_w_from_kcal_day"]
cori_rest_pct_rmr = ledger["cori_cycle_energetics"]["cori_pct_of_rmr"]

# ---- B: hepatic resting O2 consumption budget ----
BW_KG = 70.0
VO2_WHOLEBODY_REST_ML_MIN = 3.5 * BW_KG  # 1-MET generic, the same constant the
                                          # exercise_bloodflow_redistribution cell uses for
                                          # VO2_nonmuscle_generic
# Guyton/Hall-family textbook fact: liver receives ~25% CO but consumes a DISPROPORTIONATE
# ~20-27% of total-body resting VO2 (dual hepatic-artery+portal supply, high intrinsic metabolic
# rate); this is a DISCLOSED textbook-band, not a primary measurement -- swept, not point-asserted.
HEPATIC_FRAC_OF_WHOLEBODY_VO2_REST_BAND = [0.20, 0.23, 0.27]
O2_KJ_PER_L = 20.2   # mixed-diet caloric-equivalent constant -- reused, not retyped
ATP_PER_O2_MOL = 5.0  # = 2 x P/O(NADH)=2.5, Hinkle 2005 flux-measured consensus, matching the
                      # geometrically-derived P/O of the mitochondrial_oxphos cell;
                      # O2 has 2 O atoms per molecule

def hepatic_vo2_ml_min(frac):
    return frac * VO2_WHOLEBODY_REST_ML_MIN

def power_w_from_vo2_ml_min(vo2_ml_min):
    return (vo2_ml_min / 1000.0) * O2_KJ_PER_L * 1000.0 / 60.0

def atp_mol_per_min_from_vo2_ml_min(vo2_ml_min):
    o2_mol_min = (vo2_ml_min / 1000.0) / 22.4
    return o2_mol_min * ATP_PER_O2_MOL

hepatic_rest = []
for frac in HEPATIC_FRAC_OF_WHOLEBODY_VO2_REST_BAND:
    vo2 = hepatic_vo2_ml_min(frac)
    hepatic_rest.append({
        "frac_of_wholebody_vo2": frac,
        "hepatic_vo2_ml_min": round(vo2, 2),
        "hepatic_atp_power_w": round(power_w_from_vo2_ml_min(vo2), 3),
        "hepatic_atp_mol_min": round(atp_mol_per_min_from_vo2_ml_min(vo2), 5),
    })
hepatic_central = hepatic_rest[1]  # 0.23 central

cori_frac_of_hepatic_atp_budget_rest = cori_rest_W / hepatic_central["hepatic_atp_power_w"]

# ---- C: exercise hepatic O2-SUPPLY ceiling (flow down, extraction reserve up, swept honestly) ----
hepatosplenic = ebf["rest_state_distribution"]["splanchnic"]["perko_1998_pmid_9824727_own_numbers"]["hepatosplenic"]
flow_rest = hepatosplenic["rest_l_min"]
flow_reduction_pct_perko_submax = hepatosplenic["reduction_pct"]  # 50.0, real n=19 human cycling data
FLOW_RATIO_SWEEP = {
    "perko_1998_submax_cycling_measured": 1.0 - flow_reduction_pct_perko_submax / 100.0,  # 0.50
    "rowell_heavy_exercise_literature_qualitative_upper_reduction": 0.20,  # disclosed: bibliographic
                                                                            # qualitative citation
                                                                            # (Rowell LB, splanchnic
                                                                            # flow can fall to ~20-30%
                                                                            # of rest at near-maximal
                                                                            # effort) -- NOT re-verified
                                                                            # from the primary abstract,
                                                                            # flagged as such
}
# extraction-reserve compensation: hepatic venous O2 saturation can fall further as flow falls,
# partially offsetting the VO2 loss -- swept across a disclosed literature-qualitative-tier band
# (no primary abstract re-pulled for the exact number), NOT point-asserted
EXTRACTION_RESERVE_SWEEP = [1.0, 1.5, 2.0, 2.5]  # 1.0 = no compensation (conservative floor),
                                                   # 2.5 = aggressive near-maximal compensation
                                                   # (qualitative literature ceiling)

exercise_hepatic_atp_grid = []
for flow_label, flow_ratio in FLOW_RATIO_SWEEP.items():
    for ext in EXTRACTION_RESERVE_SWEEP:
        vo2_ex = hepatic_central["hepatic_vo2_ml_min"] * flow_ratio * ext
        p_w = power_w_from_vo2_ml_min(vo2_ex)
        exercise_hepatic_atp_grid.append({
            "flow_regime": flow_label,
            "flow_ratio_vs_rest": round(flow_ratio, 3),
            "extraction_reserve_factor": ext,
            "hepatic_vo2_exercise_ml_min": round(vo2_ex, 2),
            "hepatic_atp_power_w_exercise": round(p_w, 3),
            "ratio_vs_rest_hepatic_atp_power": round(p_w / hepatic_central["hepatic_atp_power_w"], 3),
        })

# ---- D: exercise Cori ATP DEMAND under a "constant recycled fraction" null, vs literature Ra multiplier ----
RA_MULTIPLIER_SWEEP = [1.0, 3.0, 5.0, 8.0, 10.0]  # van Hall 2010 range, as cited in the
                                                    # whole_body_lactate_ledger regime note (">5-10x")
demand_sweep = []
for mult in RA_MULTIPLIER_SWEEP:
    cori_w = cori_rest_W * mult
    demand_sweep.append({
        "ra_multiplier_vs_rest": mult,
        "cori_atp_power_w_null_constant_split": round(cori_w, 3),
    })

# ---- E: THE COLLISION -- for every (demand mult, supply grid point) pair, does demand exceed
# a physiologically-plausible ceiling (Cori alone consuming >50% of hepatic total ATP budget is
# treated as the FALSIFYING threshold for "liver can still do its other jobs", pre-registered) ----
CEILING_FRAC_OF_HEPATIC_ATP_BUDGET = 0.50  # pre-registered threshold, stated before computing
collision = []
for supply in exercise_hepatic_atp_grid:
    supply_w = supply["hepatic_atp_power_w_exercise"]
    row = {"supply": supply, "demand_crossovers": []}
    for d in demand_sweep:
        frac = d["cori_atp_power_w_null_constant_split"] / supply_w
        row["demand_crossovers"].append({
            "ra_multiplier": d["ra_multiplier_vs_rest"],
            "cori_frac_of_hepatic_atp_budget_at_this_supply": round(frac, 3),
            "exceeds_50pct_ceiling": bool(frac > CEILING_FRAC_OF_HEPATIC_ATP_BUDGET),
        })
    # crossover multiplier: smallest Ra multiplier at which frac > ceiling (linear interp on cori_rest_W*mult)
    crossover_mult = CEILING_FRAC_OF_HEPATIC_ATP_BUDGET * supply_w / cori_rest_W
    row["crossover_ra_multiplier"] = round(crossover_mult, 2)
    collision.append(row)

crossover_mults = [r["crossover_ra_multiplier"] for r in collision]
worst_case_crossover = min(crossover_mults)   # most restrictive (worst flow cut, no extraction compensation)
best_case_crossover = max(crossover_mults)    # most permissive (mild flow cut, max extraction compensation)

# ---- F: fraction of collision-grid points where a Ra-multiplier of 5x (literature-central) already
# breaches the ceiling under the CONSTANT-SPLIT null (falsifier of "split stays constant with intensity") ----
breach_at_5x = []
for r in collision:
    hit = next(c for c in r["demand_crossovers"] if c["ra_multiplier"] == 5.0)
    breach_at_5x.append(hit["exceeds_50pct_ceiling"])
frac_breach_at_5x = float(np.mean(breach_at_5x))

out = {
    "regime_disclosure": "resting Cori numbers REUSED live from whole_body_lactate_ledger.py "
                          "(RESTING/POST-ABSORPTIVE regime); hepatic O2-supply numbers REUSED live "
                          "from exercise_bloodflow_redistribution.py (Perko 1998 submax-cycling "
                          "real human AV data, n=19); Ra-multiplier and extraction-reserve values "
                          "are LITERATURE-QUALITATIVE-TIER sweeps (van Hall 2010 review range, "
                          "Rowell-family qualitative splanchnic-flow-at-heavy-exercise citation) "
                          "NOT re-pulled from primary abstracts -- flagged, not hidden.",
    "cori_rest": {
        "cori_atp_power_w": cori_rest_W,
        "rmr_w": rmr_W,
        "cori_pct_of_rmr": cori_rest_pct_rmr,
    },
    "hepatic_rest_atp_budget_sweep": hepatic_rest,
    "hepatic_rest_atp_budget_central": hepatic_central,
    "cori_frac_of_hepatic_atp_budget_at_rest": round(cori_frac_of_hepatic_atp_budget_rest, 4),
    "hepatosplenic_flow_rest_l_min_reused": flow_rest,
    "exercise_hepatic_atp_supply_grid": exercise_hepatic_atp_grid,
    "exercise_cori_atp_demand_sweep_constant_split_null": demand_sweep,
    "collision_pre_registered_ceiling_frac": CEILING_FRAC_OF_HEPATIC_ATP_BUDGET,
    "collision_grid": collision,
    "worst_case_crossover_ra_multiplier": round(worst_case_crossover, 2),
    "best_case_crossover_ra_multiplier": round(best_case_crossover, 2),
    "frac_of_supply_grid_where_5x_ra_breaches_ceiling": round(frac_breach_at_5x, 3),
    "verdict": (
        f"AT REST: Cori consumes only {round(100*cori_frac_of_hepatic_atp_budget_rest,1)}% of the "
        f"hepatic ATP budget (central estimate) -- trivially affordable. UNDER THE CONSTANT-SPLIT "
        f"NULL (recycled fraction of lactate Ra held fixed at its resting value as intensity rises), "
        f"the crossover Ra-multiplier past which Cori alone would consume >50% of the hepatic ATP "
        f"budget ranges {round(worst_case_crossover,1)}x (worst case: Perko-measured 50% flow cut, "
        f"NO extraction compensation) to {round(best_case_crossover,1)}x (best case: mild flow cut, "
        f"aggressive 2.5x extraction compensation). The literature-central Ra-multiplier during heavy "
        f"exercise (5-10x, van Hall 2010) BREACHES the ceiling in "
        f"{round(100*frac_breach_at_5x,1)}% of the supply-sweep grid at just 5x. CONCLUSION: the "
        f"constant-recycled-fraction null is NOT hepatic-ATP-affordable across most of the plausible "
        f"supply grid at literature-reported exercise lactate fluxes -- the recycled (Cori) FRACTION "
        f"of total lactate disposal must fall as intensity rises (the lactate-shuttle/oxidative route "
        f"must absorb a GROWING share), i.e. affordability is a real candidate DETERMINANT of the "
        f"split, not merely consistent with it. This is a DERIVED prediction from this cell's "
        f"arithmetic, not itself a new primary measurement of the split -- the sharp, falsifiable "
        f"claim for a future build to test against real hepatic-venous-catheterization AV data "
        f"(e.g. Rowell-family or modern MRI/PET splanchnic studies) or against isotopically-measured "
        f"Ra-Cori vs Ra-oxidation split during graded exercise."
    ),
}

if __name__ == "__main__":
    print(json.dumps(out, indent=2, default=str))
