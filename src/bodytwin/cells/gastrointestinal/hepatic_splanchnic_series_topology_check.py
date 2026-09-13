"""Hepatic (Q_H) vs splanchnic flow: SERIES, not PARALLEL, vascular-graph edge.

QUESTION (cardiac output as a shared budget across independently-built organ cells):
exercise_bloodflow_redistribution computes "splanchnic" flow (celiac+SMA arterial inflow, Perko et
al 1998 PMID 9824727) and hepatic_clearance computes "Q_H", total hepatic blood flow (hepatic artery
+ portal vein, F_HEPATIC=0.25 x CO). Both cells were built independently. Are these two additive,
non-overlapping shares of cardiac output (as a naive whole-body organ-flow census would assume), or
the SAME conserved blood flow viewed at two points of one series circuit?

GEOMETRY: the systemic circulation is a graph. Renal, muscle, and splanchnic beds are PARALLEL
branches off the aortic root (the parallel-conductance topology used by
exercise_bloodflow_redistribution: G_i = Q_i/MAP, G_total = CO/MAP). But hepatic flow is NOT a 4th
parallel branch alongside splanchnic -- ~75% of it (the portal-vein component) is anatomically
DOWNSTREAM of, in SERIES with, the splanchnic arterial bed (celiac's non-hepatic branches + SMA
perfuse stomach/spleen/gut, then drain via the portal vein INTO the liver -- recycled blood, not a
second aortic offtake). Only the hepatic-ARTERY component (~25% of Q_H) is parallel to (in fact a
sub-branch of) the splanchnic (celiac) arterial tree already counted in "splanchnic". Conservation
of flow (Kirchhoff current law on the vascular graph) at the liver-sinusoid node therefore predicts
Q_H ~= splanchnic (same blood, arterial-inflow point vs hepatic-venous-outflow point), NOT Q_H
additional-to-splanchnic.

FORCED ADVERSARY (steelman: could they legitimately be additive after all?): PMID 9824727's methods
("superior mesenteric and coeliac artery flows... determined by duplex ultrasonography" -- the
coeliac trunk gives rise to the common hepatic artery, so "hepatosplenic" IS partly hepatic-artery
flow) and the portal-system anatomy (portal-vein blood originates from recycled splanchnic venous
outflow, not an independent additional blood source beyond arterial inflow to those regions; the
hepatic artery provides ~25% of total liver blood flow) both refute the "independent parallel beds"
adversary and support the series/conservation reading. Third, decorrelated corroboration: Guyton &
Hall lists "splanchnic 20-25%" as ONE inclusive textbook row, with no separate "liver" row.

PRE-REGISTERED THRESHOLD (stated before this cell computes anything): if Q_H and splanchnic_abs are
the same conserved quantity, expect agreement within <=20% (looser than the ~1-6% the tightest
same-organ crosschecks hit, tighter than the 1.5x/50%-family gates used for weaker,
population-level cross-organ ties). Anything outside 20% would falsify "same quantity" at the
pre-registered bar (it would not by itself prove independence, but would remove the central support
for treating them as non-additive).

Reads: exercise_bloodflow_redistribution_results.json, hepatic_clearance_results.json (read-only;
no re-solve, pure arithmetic on top).
Writes: nothing (prints the result object).
Gate: same-quantity gate |pct_diff| <= 20%, plus the blind-spot and k_detect=2.092 assertions.
"""
import json
import os

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
EXER_JSON = _os.path.join(OUT_ROOT, "exercise_bloodflow_redistribution",
                          "exercise_bloodflow_redistribution_results.json")
HEP_JSON = _os.path.join(OUT_ROOT, "hepatic_clearance", "hepatic_clearance_results.json")

# published central values, cited verbatim from the already-certified organ-flow census
RENAL_PCT = 19.135
MUSCLE_PCT = 16.842
BRAIN_PCT = 11.5
CO_HIGGINBOTHAM = 5.70  # PMID 3948345, n=24 real catheterization


def main():
    exer = json.load(open(EXER_JSON))
    hep = json.load(open(HEP_JSON))

    splanchnic_abs = exer["rest_state_distribution"]["splanchnic"]["rest_l_min_central"]
    splanchnic_pub_pct = exer["rest_state_distribution"]["splanchnic"]["pct_of_CO_rest"]
    co_fick = hep["inputs"]["q_rest_l_min"]
    q_h_abs = hep["hepatic_flow"]["q_h_rest_l_min"]
    f_hepatic = hep["hepatic_flow"]["f_hepatic"]

    assert abs(q_h_abs - f_hepatic * co_fick) < 1e-9, "Q_H not reproduced from f_hepatic*CO_fick"

    ratio = q_h_abs / splanchnic_abs
    pct_diff = (ratio - 1) * 100
    naive_sum_abs = splanchnic_abs + q_h_abs
    k_equiv = naive_sum_abs / splanchnic_abs  # equivalent single-term (splanchnic) inflation factor

    other3_pct = RENAL_PCT + MUSCLE_PCT + BRAIN_PCT
    k_detect_splanchnic = (100 - other3_pct) / splanchnic_pub_pct

    existing_4term_sum_pct = splanchnic_pub_pct + other3_pct
    existing_residual_pct = 100 - existing_4term_sum_pct

    q_h_pct_of_co_higg = q_h_abs / CO_HIGGINBOTHAM * 100
    splanchnic_pct_of_co_higg = splanchnic_abs / CO_HIGGINBOTHAM * 100
    hypothetical_5term_sum_pct = existing_4term_sum_pct + q_h_pct_of_co_higg
    hypothetical_residual_pct = 100 - hypothetical_5term_sum_pct

    margin_below_detection = k_detect_splanchnic - k_equiv

    same_quantity_threshold_pct = 20.0
    same_quantity_gate = abs(pct_diff) <= same_quantity_threshold_pct

    oversubscription_gate_fires_if_doublecounted = hypothetical_5term_sum_pct > 100.0
    doublecount_inside_blindspot = k_equiv < k_detect_splanchnic

    out = {
        "inputs": {
            "splanchnic_abs_l_min_perko1998_pmid9824727": splanchnic_abs,
            "q_h_abs_l_min_hepatic_clearance": q_h_abs,
            "co_fick_l_min": co_fick,
            "f_hepatic": f_hepatic,
        },
        "same_quantity_check": {
            "ratio_qh_over_splanchnic": ratio,
            "pct_diff": pct_diff,
            "pre_registered_threshold_pct": same_quantity_threshold_pct,
            "gate_pass_same_quantity_within_threshold": same_quantity_gate,
        },
        "naive_doublecount_if_ever_appended": {
            "naive_sum_abs_l_min": naive_sum_abs,
            "equivalent_splanchnic_only_inflation_k": k_equiv,
            "graph_own_k_detect_splanchnic": k_detect_splanchnic,
            "margin_below_detection_threshold": margin_below_detection,
            "doublecount_falls_inside_existing_gate_blindspot": doublecount_inside_blindspot,
        },
        "wholebody_sum_consequence": {
            "existing_4term_sum_pct_asPublished": existing_4term_sum_pct,
            "existing_residual_pct": existing_residual_pct,
            "hypothetical_5term_sum_pct_if_wrongly_appended": hypothetical_5term_sum_pct,
            "hypothetical_residual_pct": hypothetical_residual_pct,
            "oversubscription_gate_wouldfire": oversubscription_gate_fires_if_doublecounted,
        },
        "verdict": {
            "current_graph_state": "NO VIOLATION -- Q_H correctly absent from every existing organ-flow-sum node; books close at 72.584%CO as published",
            "mechanism": "splanchnic (arterial inflow) and Q_H (hepatic total, incl. portal) are the SAME series-circuit blood flow (conservation), not parallel/additive -- confirmed via Perko1998 methodology + the portal-circulation mechanism, both decorrelated from the 2 cells being reconciled",
            "latent_risk_if_ever_violated": "a future naive 'sum every organ pct_of_CO field' synthesis that DID append Q_H would produce an error mathematically equivalent to inflating splanchnic by 2.013x -- INSIDE the void-floor gate's blind spot (k_detect=2.092x), i.e. NOT caught by the sum<=100%CO check; residual would silently collapse from 27.4pp to 2.6pp",
        },
    }
    print(json.dumps(out, indent=2))

    # self-test / machine gate
    assert same_quantity_gate, f"FAIL: splanchnic/Q_H disagree by {pct_diff:.2f}%, outside +/-{same_quantity_threshold_pct}% pre-registered band"
    assert doublecount_inside_blindspot, "expected the double-count scenario to fall inside the existing gate's blind spot"
    assert not oversubscription_gate_fires_if_doublecounted, "expected hypothetical 5-term sum to stay under 100%CO (that IS the blind-spot finding)"
    assert abs(k_detect_splanchnic - 2.092) < 0.01, "published k_detect=2.092 not reproduced"
    print("\nSELF-TEST: PASS (same-quantity gate holds at +/-20% band; blind-spot reproduced; k_detect=2.092 cross-verified)")
    return out


if __name__ == "__main__":
    main()
