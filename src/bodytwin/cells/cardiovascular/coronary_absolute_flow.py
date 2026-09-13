"""Absolute resting and maximal coronary blood flow via two decorrelated routes.

Route A: PET-measured per-gram myocardial blood flow x heart mass. Route B: Fick, using the
coronary arterio-venous O2 difference (CaO2 x extraction fraction) and two method-families of
myocardial O2 consumption (bottom-up mechanical / top-down allocation). Reports the overlap
between the routes, and feeds the result back into the resting and maximal cardiac-output
flow budgets.

Reads:  the blood_oxygen_transport, coronary_blood_flow and exercise_bloodflow_redistribution
        cell result JSONs.
Writes: nothing (prints its cross-checks).
Gate:   overlap width / point gaps between Route A and each Route B method family.
"""
import json
import os as _os

OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

# ---- pull real cross-cell numbers, read-only ----
bo = json.load(open(_os.path.join(OUT_ROOT, "blood_oxygen_transport",
                                  "blood_oxygen_transport_results.json")))
cao2_rest = bo["chemistry_rest"]["cao2_ml_dl"]  # mL/dL, the reference body's arterial O2 content

cbf = json.load(open(_os.path.join(OUT_ROOT, "coronary_blood_flow",
                                   "coronary_blood_flow_results.json")))
ext_pct = cbf["claim_B_o2_extraction_flow_reserve"]["lv_rest_o2_extraction_pct (midpoint of 70-80%)"]  # 75.0
ext_lo, ext_hi = 70.0, 80.0
flow_rise_factor = cbf["claim_B_o2_extraction_flow_reserve"]["lv_real_measured_flow_rise_x"]  # 5.0 (Duncker&Bache)

# ---- external anchors ----
HEART_MASS_G = 331.0   # Molina & DiMaio 2012, PMID 22182983, Am J Forensic Med Pathol 33(4):362-7,
                        # n=232 healthy men 18-35, mean 331+-56.7g, 95% ref range 233-383g
PET_MBF_LO, PET_MBF_HI = 0.6, 1.0   # mL/min/g resting whole-heart MBF, PET (13N-ammonia/Rb-82 family),
                                     # Murthy et al 2018 JNM 59(2):273-93 PMID 29242396 (SNMMI/ASNC joint position paper)

# ---- ROUTE A: direct (per-gram flow x heart mass) ----
routeA_lo = PET_MBF_LO * HEART_MASS_G
routeA_hi = PET_MBF_HI * HEART_MASS_G
routeA_mid = 0.5*(PET_MBF_LO+PET_MBF_HI) * HEART_MASS_G

# ---- ROUTE B: Fick, MVO2(mL/100g/min) x mass -> total O2 consumed, / coronary AVO2diff ----
def fick_flow(mvo2_per100g, ext_frac):
    avo2diff = cao2_rest * ext_frac          # mL O2 / dL blood
    mvo2_total = mvo2_per100g * HEART_MASS_G / 100.0   # mL O2/min
    return mvo2_total * 100.0 / avo2diff              # mL/min

# method-family MVO2 values (PVA/MVO2 and rate-pressure-product cross-check families)
BOTTOMUP_LO, BOTTOMUP_HI = 10.0, 20.0   # bottom-up mechanical (Suga PVA: LV=20.12,WH=14.93; coronaryFick=12.35)
TOPDOWN_LO, TOPDOWN_HI = 6.0, 10.0      # top-down allocation (wholebodyfrac=6.6-9.42; Gobel1978 direct=6.7+-1.3)
BOTTOMUP_PT = 14.93   # Suga-WH point (matches WH-scope mass used here)
TOPDOWN_PT_GOBEL = 6.7  # Gobel 1978 PMID624164 direct N2O-Fick measurement, n=27
TOPDOWN_PT_ALLOC = 8.01  # midpoint of wholebodyfrac 6.6-9.42

# full sensitivity (MVO2 range x extraction range) bounds
def fick_bounds(mvo2_lo, mvo2_hi):
    vals = [fick_flow(m, e/100.0) for m in (mvo2_lo, mvo2_hi) for e in (ext_lo, ext_hi)]
    return min(vals), max(vals)

bottomup_flow_lo, bottomup_flow_hi = fick_bounds(BOTTOMUP_LO, BOTTOMUP_HI)
topdown_flow_lo, topdown_flow_hi = fick_bounds(TOPDOWN_LO, TOPDOWN_HI)
bottomup_pt_flow = fick_flow(BOTTOMUP_PT, ext_pct/100.0)
topdown_pt_flow_gobel = fick_flow(TOPDOWN_PT_GOBEL, ext_pct/100.0)
topdown_pt_flow_alloc = fick_flow(TOPDOWN_PT_ALLOC, ext_pct/100.0)

def overlap(a_lo, a_hi, b_lo, b_hi):
    return max(0.0, min(a_hi, b_hi) - max(a_lo, b_lo))

ov_bottomup = overlap(routeA_lo, routeA_hi, bottomup_flow_lo, bottomup_flow_hi)
ov_topdown = overlap(routeA_lo, routeA_hi, topdown_flow_lo, topdown_flow_hi)

gap_bottomup_pt = abs(routeA_mid - bottomup_pt_flow)
gap_topdown_pt_gobel = abs(routeA_mid - topdown_pt_flow_gobel)
gap_topdown_pt_alloc = abs(routeA_mid - topdown_pt_flow_alloc)

print("=== ROUTE A (direct, PET mL/min/g x heart mass) ===")
print(f"heart_mass_g={HEART_MASS_G}, PET MBF {PET_MBF_LO}-{PET_MBF_HI} mL/min/g")
print(f"routeA rest flow mL/min: lo={routeA_lo:.1f} hi={routeA_hi:.1f} mid={routeA_mid:.1f}")

print("\n=== ROUTE B (Fick, coronary AVO2diff = CaO2*extraction) ===")
print(f"cao2_rest_ml_dl={cao2_rest:.4f}, extraction range {ext_lo}-{ext_hi}% (midpoint {ext_pct})")
print(f"bottom-up-mechanical family flow range mL/min: {bottomup_flow_lo:.1f}-{bottomup_flow_hi:.1f} (point Suga-WH={bottomup_pt_flow:.1f})")
print(f"top-down-allocation family flow range mL/min: {topdown_flow_lo:.1f}-{topdown_flow_hi:.1f} (point Gobel={topdown_pt_flow_gobel:.1f}, point alloc={topdown_pt_flow_alloc:.1f})")

print("\n=== CROSS-CHECK: which family does Route A (direct) agree with ===")
print(f"overlap width Route A vs bottom-up: {ov_bottomup:.1f} mL/min")
print(f"overlap width Route A vs top-down: {ov_topdown:.1f} mL/min")
print(f"point-gap Route A_mid vs bottomup_pt: {gap_bottomup_pt:.1f} mL/min ({100*gap_bottomup_pt/routeA_mid:.1f}%)")
print(f"point-gap Route A_mid vs topdown_pt(Gobel direct): {gap_topdown_pt_gobel:.1f} mL/min ({100*gap_topdown_pt_gobel/routeA_mid:.1f}%)")
print(f"point-gap Route A_mid vs topdown_pt(alloc): {gap_topdown_pt_alloc:.1f} mL/min ({100*gap_topdown_pt_alloc/routeA_mid:.1f}%)")

print("\n=== INTERNAL CONSISTENCY: extraction near-fixed -> flow should scale ~linearly with MVO2 ===")
demand_rise = 6.0
required_flow_rise_if_ext_maxes = demand_rise / (95.0/ext_pct)
print(f"measured flow_rise_x={flow_rise_factor}, demand_rise_x={demand_rise}, required-if-ext-maxes-at-95%={required_flow_rise_if_ext_maxes:.2f} (within 20% of measured, already gated PASS in the coronary_blood_flow cell)")

print("\n=== MAX EXERCISE coronary flow (x5, Duncker&Bache) ===")
maxA_lo, maxA_hi = routeA_lo*flow_rise_factor, routeA_hi*flow_rise_factor
maxA_mid = routeA_mid*flow_rise_factor
print(f"routeA max flow mL/min: lo={maxA_lo:.1f} hi={maxA_hi:.1f} mid={maxA_mid:.1f}")

# ---- FEED BACK: rest organ-flow sum ----
# Three resting-CO values are in circulation: 5.70 (Higginbotham et al. 1986, PMID 3948345, REAL
# direct right-heart catheterization, n=24), 5.5575 (the cardiac_output_geometric cell's
# population-geometric CMR-derived central estimate) and 5.653 (the cardiac_output cell's
# Fick-chain estimate). This cell deliberately anchors to the REAL measured value, a decorrelated
# cross-check against the two modeled routes -- collapsing it into either would destroy that
# decorrelation (max pairwise spread observed: 2.57%).
CO_REST_LMIN = 5.70   # Higginbotham 1986 anchor, already used by sibling cells
residual_pct = 27.17  # from the resting organ-flow-sum re-execution
residual_lmin = residual_pct/100.0 * CO_REST_LMIN
print("\n=== REST SUM FEEDBACK ===")
print(f"CO_rest={CO_REST_LMIN} L/min, residual={residual_pct}% = {residual_lmin:.3f} L/min")
for label, lo, hi, mid in [("routeA_direct", routeA_lo, routeA_hi, routeA_mid)]:
    pct_lo = lo/1000.0/CO_REST_LMIN*100
    pct_hi = hi/1000.0/CO_REST_LMIN*100
    pct_mid = mid/1000.0/CO_REST_LMIN*100
    frac_of_residual_lo = lo/1000.0/residual_lmin*100
    frac_of_residual_hi = hi/1000.0/residual_lmin*100
    frac_of_residual_mid = mid/1000.0/residual_lmin*100
    print(f"{label}: {lo:.1f}-{hi:.1f} mL/min = {pct_lo:.2f}-{pct_hi:.2f}% of CO_rest (mid {pct_mid:.2f}%)")
    print(f"  -> accounts for {frac_of_residual_lo:.1f}-{frac_of_residual_hi:.1f}% of the 27.17% residual (mid {frac_of_residual_mid:.1f}%)")

# ---- FEED BACK: max-exercise joint oversubscription (exercise + heat) ----
CO_MAX_LMIN = 23.6208
# The muscle flow fraction is READ from its producing cell rather than hardcoded: a hardcoded
# 91.647 came from a uniform-extraction tautology (systemic a-vO2diff applied to the muscle bed)
# and would go stale again.
import json as _json
_ebr = _json.load(open(_os.path.join(OUT_ROOT, "exercise_bloodflow_redistribution",
                                     "exercise_bloodflow_redistribution_results.json")))
muscle_max_pct = _ebr["forced_adversary_2_conservation_of_flow_budget"]["muscle_flow_frac_of_co_max_pct"]
skin_max_lo_pct = 6.0/CO_MAX_LMIN*100
skin_max_hi_pct = 8.0/CO_MAX_LMIN*100
existing_oversub_lo = muscle_max_pct + skin_max_lo_pct - 100.0
existing_oversub_hi = muscle_max_pct + skin_max_hi_pct - 100.0
print("\n=== MAX SUM FEEDBACK ===")
print(f"existing (muscle+skin) oversubscription: {existing_oversub_lo:.1f}-{existing_oversub_hi:.1f} pp over CO_max")
coronary_max_pct_lo = maxA_lo/1000.0/CO_MAX_LMIN*100
coronary_max_pct_hi = maxA_hi/1000.0/CO_MAX_LMIN*100
new_oversub_lo = existing_oversub_lo + coronary_max_pct_lo
new_oversub_hi = existing_oversub_hi + coronary_max_pct_hi
print(f"coronary adds {coronary_max_pct_lo:.2f}-{coronary_max_pct_hi:.2f}pp of CO_max")
print(f"new total oversubscription: {new_oversub_lo:.1f}-{new_oversub_hi:.1f} pp over CO_max (was {existing_oversub_lo:.1f}-{existing_oversub_hi:.1f})")
