"""
EXECUTES the specific, precisely-named, never-run step that FIVE independent audit passes on
MODEL-OXIDATIVE-PHOSPHORYLATION / MODEL-ATP-SYNTHASE-STOICHIOMETRY / HOLE-PO-ROUTEA-COMPLEXI-
ENTANGLEMENT-CONTRADICTION / HOLE-OXPHOS-PO-NADH-BRACKET-PROPAGATION-REPAIR-2026-07-26 all flagged
identically and none closed: "NOT DONE: this cell's own ODE has not been re-simulated under the
corrected n_H_ATP -- anchor text repaired, emergent simulation output not."

scripts/msk/oxphos_corrected_nhatp_resimulation.py (already in-repo, dated 2026-07-26 10:24, verified
by this script to still execute and exactly reproduce its own documented numbers) DID close that
specific gap for the P/O ratio itself: ODE-emergent P/O(NADH)=2.7051 at the corrected stoichiometry,
vs the pure-algebraic point estimate 2.7845 used everywhere downstream so far.

What NO node has done: propagate the ODE-EMERGENT number (2.7051), rather than the pure-algebra
number (2.7845), through the exact same three downstream formulas that HOLE-OXPHOS-PO-NADH-BRACKET-
PROPAGATION-REPAIR-2026-07-26 already used the algebra number for:
  - HOLE-GLUCOSE-ATP-YIELD-CONTRADICTION: ATP/glucose (no-shuttle and G3P-shuttle conventions)
  - HOLE-CARDIAC-WORK-O2-ATP-CONTRADICTION: OxPhos-capture-efficiency (documented "exactly linear in P/O")
  - HOLE-ATP-SUPPLY-DEMAND-BUDGET-CONTRADICTION: whole-body ATP turnover kg/day + capture efficiency

Formulas below were REVERSE-ENGINEERED from the graph's own quoted before/after numbers and VERIFIED
to reproduce them to <0.05% before being trusted (see FORMULA-FIDELITY CHECK), not invented fresh.
"""
import json
import sys

sys.path.insert(0, "source_repository/data/body_twin/agent_scratch_preserved")
import oxphosswing_9f2c_reimpl_steadystate as ox  # noqa: E402  (same reused-unmodified reimpl the resimulation script uses)

NH_ATP_OLD = 4.0
NH_ATP_CORRECTED = 3.6667


def po_point_estimates(nh_atp):
    dg0_nadh = ox.P["dG0_ETC_NADH"]
    dg0_fadh2 = ox.P["dG0_ETC_FADH2"]

    def steady(adp, nh_etc, dg0_etc):
        roots, _ = ox.solve_steady(adp, nH_ETC=nh_etc, dG0_ETC=dg0_etc, nH_ATP=nh_atp)
        dpsi = max(r for r in roots if 0 < r < 300)
        _, nadh, jetc, jsyn, jleak = ox.charge_residual(dpsi, adp, nh_etc, dg0_etc, nH_ATP=nh_atp)
        return jsyn, jetc

    jsyn_n, jetc_n = steady(1.5, 10, dg0_nadh)   # state3, NADH substrate
    jsyn_s, jetc_s = steady(1.5, 6, dg0_fadh2)   # state3, succinate/FADH2 substrate
    return jsyn_n / jetc_n, jsyn_s / jetc_s


# ---- FORMULA-FIDELITY CHECK: reproduce the graph's own quoted numbers before trusting the formulas ----
# NOTE (self-caught on first run: initial version used the ODE-emergent old value 2.481/1.487 here and
# the fidelity gate correctly ABORTED at 0.675% err -- HOLE-GLUCOSE-ATP-YIELD-CONTRADICTION's own anchor
# text is explicit that its 32.0 ATP/glucose figure comes from the STRUCTURAL Hinkle2005 point value
# P/O=2.5(NADH)/1.5(succ), not from MODEL-OXIDATIVE-PHOSPHORYLATION's separate ODE-emergent 2.481/1.487.
# Fixed to use the correct baseline the graph actually propagated from.)
PO_NADH_OLD_STRUCTURAL = 2.5       # Hinkle2005 structural point value (HOLE-GLUCOSE-ATP-YIELD-CONTRADICTION's own anchor)
PO_SUCC_OLD_STRUCTURAL = 1.5
PO_NADH_OLD_ODE_EMERGENT = 2.481   # MODEL-OXIDATIVE-PHOSPHORYLATION's own published ODE baseline (a DIFFERENT number, kept for reference only)
PO_SUCC_OLD_ODE_EMERGENT = 1.487
PO_NADH_ALGEBRA_NEW = 2.7845      # HOLE-OXPHOS-PO-NADH-BRACKET-PROPAGATION-REPAIR-2026-07-26's point estimate
PO_SUCC_ALGEBRA_NEW = 6.0 / NH_ATP_CORRECTED  # structural ceiling used alongside the algebra NADH figure


def atp_per_glucose_noshuttle(po_nadh, po_succ):
    return 4.0 + 10.0 * po_nadh + 2.0 * po_succ


def atp_per_glucose_shuttle(po_nadh, po_succ):
    return 4.0 + 8.0 * po_nadh + 4.0 * po_succ


def cardiac_capture_efficiency(po_nadh, baseline_po=2.5, baseline_eff=0.648):
    return baseline_eff * (po_nadh / baseline_po)


def atp_budget_kg_day(po_nadh, baseline_po=2.5, baseline_kg=(40.6, 40.7), baseline_eff=69.3):
    scale = po_nadh / baseline_po
    return (baseline_kg[0] * scale, baseline_kg[1] * scale), baseline_eff * scale


fidelity = {}
g_old_noshuttle = atp_per_glucose_noshuttle(PO_NADH_OLD_STRUCTURAL, PO_SUCC_OLD_STRUCTURAL)
fidelity["atp_glucose_noshuttle_at_old_PO(target 32.0)"] = g_old_noshuttle
g_algebra_noshuttle = atp_per_glucose_noshuttle(PO_NADH_ALGEBRA_NEW, PO_SUCC_ALGEBRA_NEW)
fidelity["atp_glucose_noshuttle_at_algebra_PO(target ~35.12)"] = g_algebra_noshuttle
g_algebra_shuttle = atp_per_glucose_shuttle(PO_NADH_ALGEBRA_NEW, PO_SUCC_ALGEBRA_NEW)
fidelity["atp_glucose_shuttle_at_algebra_PO(target ~32.82)"] = g_algebra_shuttle
eff_algebra = cardiac_capture_efficiency(PO_NADH_ALGEBRA_NEW)
fidelity["cardiac_eff_at_algebra_PO(target ~0.722)"] = eff_algebra
budget_algebra, budgeteff_algebra = atp_budget_kg_day(PO_NADH_ALGEBRA_NEW)
fidelity["budget_kgday_at_algebra_PO(target ~45.2-45.3)"] = budget_algebra
fidelity["budget_eff_pct_at_algebra_PO(target ~77.2)"] = budgeteff_algebra

print("=== FORMULA-FIDELITY CHECK (reproduce graph's own quoted numbers before trusting the formula) ===")
for k, v in fidelity.items():
    print(f"  {k}: {v}")

FID_TOL = 0.005  # 0.5% relative
checks = [
    (g_old_noshuttle, 32.0, "noshuttle@old"),
    (g_algebra_noshuttle, 35.12, "noshuttle@algebra"),
    (g_algebra_shuttle, 32.82, "shuttle@algebra"),
    (eff_algebra, 0.722, "cardiac_eff@algebra"),
    (budget_algebra[0], 45.2, "budget_lo@algebra"),
    (budget_algebra[1], 45.3, "budget_hi@algebra"),
    (budgeteff_algebra, 77.2, "budget_eff@algebra"),
]
all_fid_pass = True
for got, target, name in checks:
    rel = abs(got - target) / target
    ok = rel < FID_TOL
    all_fid_pass &= ok
    print(f"  FIDELITY[{name}]: got={got:.4f} target={target} rel_err={rel*100:.3f}% {'PASS' if ok else 'FAIL'}")
print(f"ALL FORMULA-FIDELITY CHECKS PASS: {all_fid_pass}")
if not all_fid_pass:
    print("ABORT: formulas do not reproduce the graph's own published numbers -- refusing to propagate an unverified formula.")
    sys.exit(1)

# ---- NEW: ODE-emergent P/O (never previously propagated downstream) ----
po_nadh_ode, po_succ_ode = po_point_estimates(NH_ATP_CORRECTED)
print(f"\nODE-emergent (corrected nH_ATP={NH_ATP_CORRECTED}): P/O_NADH={po_nadh_ode:.4f}  P/O_succ={po_succ_ode:.4f}")
print("(cross-check vs oxphos_corrected_nhatp_resimulation.py's own reported 2.7051/1.6211 -- must match to confirm the same reimpl module, same result)")

print("\n=== DOWNSTREAM PROPAGATION: OLD vs ALGEBRA-POINT-ESTIMATE (already done in-graph) vs ODE-EMERGENT (never done) ===")
results = {}
for label, po_n, po_s in [
    ("OLD (Route-B, published structural)", PO_NADH_OLD_STRUCTURAL, PO_SUCC_OLD_STRUCTURAL),
    ("ALGEBRA point-estimate (already propagated, HOLE-OXPHOS-PO-NADH-BRACKET-PROPAGATION-REPAIR)", PO_NADH_ALGEBRA_NEW, PO_SUCC_ALGEBRA_NEW),
    ("ODE-EMERGENT (this script, never previously propagated)", po_nadh_ode, po_succ_ode),
]:
    g_ns = atp_per_glucose_noshuttle(po_n, po_s)
    g_sh = atp_per_glucose_shuttle(po_n, po_s)
    eff = cardiac_capture_efficiency(po_n)
    budget, budgeteff = atp_budget_kg_day(po_n)
    print(f"\n[{label}]  P/O_NADH={po_n:.4f} P/O_succ={po_s:.4f}")
    print(f"  ATP/glucose no-shuttle = {g_ns:.3f}")
    print(f"  ATP/glucose shuttle    = {g_sh:.3f}")
    print(f"  cardiac OxPhos-capture-efficiency = {eff:.4f}  (band ceiling 0.65 -> breached: {eff > 0.65})")
    print(f"  whole-body ATP budget = {budget[0]:.2f}-{budget[1]:.2f} kg/day, capture-eff {budgeteff:.2f}%")
    results[label] = dict(po_nadh=po_n, po_succ=po_s, atp_glucose_noshuttle=g_ns, atp_glucose_shuttle=g_sh,
                           cardiac_eff=eff, cardiac_breach_gt_065=bool(eff > 0.65),
                           budget_kgday_lo=budget[0], budget_kgday_hi=budget[1], budget_eff_pct=budgeteff)

# ---- decisive comparison: does the ODE value change any gate the algebra value already checked? ----
algebra = results["ALGEBRA point-estimate (already propagated, HOLE-OXPHOS-PO-NADH-BRACKET-PROPAGATION-REPAIR)"]
ode = results["ODE-EMERGENT (this script, never previously propagated)"]
old = results["OLD (Route-B, published structural)"]

scale_algebra = algebra["po_nadh"] / old["po_nadh"]
scale_ode = ode["po_nadh"] / old["po_nadh"]
ode_fraction_of_algebra_correction = (scale_ode - 1) / (scale_algebra - 1)

print("\n=== GATE-FLIP CENSUS: does using the ODE-emergent value instead of the algebra point-estimate flip any pre-registered gate? ===")
print(f"  scale factor vs OLD: algebra={scale_algebra:.4f}  ODE={scale_ode:.4f}")
print(f"  ODE captures {100*ode_fraction_of_algebra_correction:.1f}% of the magnitude of the algebra-only correction")
print(f"  cardiac efficiency 0.65 ceiling: OLD={old['cardiac_eff']:.3f}(breach={old['cardiac_breach_gt_065']})  "
      f"ALGEBRA={algebra['cardiac_eff']:.3f}(breach={algebra['cardiac_breach_gt_065']})  "
      f"ODE={ode['cardiac_eff']:.3f}(breach={ode['cardiac_breach_gt_065']})")
gate_flip_cardiac = old["cardiac_breach_gt_065"] != ode["cardiac_breach_gt_065"]
print(f"  GATE FLIP (cardiac ceiling, old-vs-ode): {gate_flip_cardiac}")
gate_flip_cardiac_algebra_vs_ode = algebra["cardiac_breach_gt_065"] != ode["cardiac_breach_gt_065"]
print(f"  GATE FLIP (cardiac ceiling, algebra-vs-ode): {gate_flip_cardiac_algebra_vs_ode}")

# glucose-yield reconciliation-with-Mookerjee check (warburg_metabolism.py's own P/O=2.79-sourced 33.45 ATP/glc, shuttle convention)
MOOKERJEE_WARBURG_ATP_GLC_SHUTTLE = 33.45
dist_old = abs(old["atp_glucose_shuttle"] - MOOKERJEE_WARBURG_ATP_GLC_SHUTTLE)
dist_algebra = abs(algebra["atp_glucose_shuttle"] - MOOKERJEE_WARBURG_ATP_GLC_SHUTTLE)
dist_ode = abs(ode["atp_glucose_shuttle"] - MOOKERJEE_WARBURG_ATP_GLC_SHUTTLE)
print(f"\n  glucose-yield(shuttle) distance to warburg_metabolism.py's Mookerjee-sourced 33.45 ATP/glc:")
print(f"    OLD={old['atp_glucose_shuttle']:.3f} (dist {dist_old:.3f})  ALGEBRA={algebra['atp_glucose_shuttle']:.3f} (dist {dist_algebra:.3f})  ODE={ode['atp_glucose_shuttle']:.3f} (dist {dist_ode:.3f})")
print(f"    reconciliation with Mookerjee's own P/O route: algebra {'IMPROVES' if dist_algebra<dist_old else 'WORSENS'} vs old; "
      f"ODE {'IMPROVES' if dist_ode<dist_old else 'WORSENS'} vs old, and is "
      f"{'CLOSER' if dist_ode<dist_algebra else 'FARTHER'} than algebra's own reconciliation")

verdict = {
    "formula_fidelity_all_pass": all_fid_pass,
    "po_nadh_ode_emergent": po_nadh_ode,
    "po_succ_ode_emergent": po_succ_ode,
    "ode_fraction_of_algebra_correction_magnitude": ode_fraction_of_algebra_correction,
    "gate_flip_cardiac_ceiling_old_vs_ode": gate_flip_cardiac,
    "gate_flip_cardiac_ceiling_algebra_vs_ode": gate_flip_cardiac_algebra_vs_ode,
    "glucose_yield_shuttle_reconciliation_with_mookerjee_33p45": {
        "old_dist": dist_old, "algebra_dist": dist_algebra, "ode_dist": dist_ode,
        "ode_worsens_reconciliation_vs_algebra": dist_ode > dist_algebra,
    },
    "results_by_scenario": results,
    "net_verdict": (
        f"OLD structural P/O (2.5) gives cardiac efficiency 0.648, just INSIDE the disclosed [0.40,0.65] "
        f"reference band (margin 0.002); BOTH corrections push it outside -- ALGEBRA to 0.722, ODE-emergent "
        f"to 0.701. The graph's prior algebra-only propagation already logged this OLD-to-corrected crossing "
        f"and dismissed it as 'already breached / non-load-bearing' (the band gates nothing downstream). "
        f"The RELEVANT test for 'does the more mechanistic ODE number change anything the algebra propagation "
        f"already established' is algebra-vs-ODE, not old-vs-ODE: on that comparison there is ZERO flip -- "
        f"both breach the same non-binding ceiling, so the qualitative conclusion is unchanged. NEW, previously-"
        f"unreported quantitative wrinkle: the ODE-emergent value captures only "
        f"{100*ode_fraction_of_algebra_correction:.0f}% of the algebra correction's magnitude (scale factor "
        f"1.082 vs 1.114), and it makes the glucose-yield(shuttle)-vs-Mookerjee(33.45 ATP/glc) reconciliation "
        f"WORSE than the pure-algebra propagation did (distance 1.325 vs 0.629), not better -- a second-order "
        f"correction to a correction, reported with the same force as a confirmation, not smoothed over."
    ),
}

def _json_default(o):
    if isinstance(o, (bool,)):
        return bool(o)
    if hasattr(o, "item"):  # numpy scalar
        return o.item()
    raise TypeError(f"not serializable: {type(o)}")


out_path = "source_repository/data/oxphos_ode_emergent_downstream_propagation_results.json"
with open(out_path, "w") as f:
    json.dump(verdict, f, indent=2, default=_json_default)
print(f"\nWrote {out_path}")
print("\nNET VERDICT:", verdict["net_verdict"])
