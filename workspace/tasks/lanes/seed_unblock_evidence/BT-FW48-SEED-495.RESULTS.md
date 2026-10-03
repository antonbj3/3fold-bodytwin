# BT-FW48-SEED-495 — RESULTS (v6 definition/review; PENDING_INDEPENDENT_REVIEW)

## Status
`NOT_EXECUTED_DESIGN_REVIEW` → now plus **one small matched synthetic mechanism test**
(`m495_method_provenance.py`, seed 20260930, frozen in PREREG.md before the run).
No scientific admission. No physical validation. No universal certificate.

## What was done
The M495 new operation is: *Method provenance separates implementation, model assumptions and
empirical tests*. Observable = `fraction of correct reuse decisions` (correct method-reuse decisions /
planned queries, dimensionless). Consumer = recursive operator-/method-planner.
Discriminating contrast = **same methodname, different solver precision**, versioned executable
method registry vs informal doc retrieval.

A frozen 4-entry registry with **one shared methodname** `oct2_mate_renal_clearance` was built;
entries A–D differ only in version, declared `eps_solver`, cost and model-assumption tag.
12 frozen query cells (eps_y × context = two unseen task families F1 no_inhibitor, F2 inhibitor);
ground truth = min-cost entry matching context with `eps_solver <= eps_y`, else ABSTAIN.

## First actual result (N=100 planned reuse queries; matched rules)
| rule | correct fraction | wrong-answer frac | abstains (gt abstains=32) |
|------|-----------------|-------------------|---------------------------|
| R_candidate (implementation + assumptions provenance) | **1.00** | 0.00 | 32 (all correct) |
| R_doc_name_only (informal doc retrieval, canonical entry) | 0.17 | 0.83 | 0 |
| R_abl_no_impl_field (drop implementation field) | 0.43 | 0.57 | 0 |
| R_abl_no_assum_field (drop assumptions field) | 0.58 | 0.42 | 16 (16 wrong) |

Per-family (F1 no_inhibitor / F2 inhibitor), candidate = 1.00/1.00; doc = 0.34/0.00;
no-impl = 0.18/0.68; no-assum = 0.68/0.48.
N=1 and N=10 (planned query counts) agree in direction; N=10 candidate 1.00 vs doc 0.20.

Separation vs candidate at N=100: doc −0.83, no-impl −0.57, no-assum −0.42.
**No ablation tie**; candidate not beaten by any baseline → passes the separation/ablation
criterion of the M495 test **on this constructed adverse-case distribution**.

## What this does and does not establish
- ESTABLISHES (conditional, under PREREG assumptions): both the implementation field and the
  model-assumption field are *load-bearing*; dropping either destroys a large fraction of correct
  reuse decisions (0.57 / 0.42 loss). The two fields are not substitutes (ablation profiles differ
  by family), i.e. implementation ≠ assumptions.
- DOES NOT establish an comparative improvement: candidate's declared lookup cost is 4/query vs 1/query for
  doc retrieval. The frozen contract requires strictly lower total cost at equal eps_y **or**
  decision-loss reduction > delta_dec at equal total cost. Setup/tuning/discovery/validation cost
  is UNKNOWN, so neither branch is satisfied with a bound. **No cost-matched comparative improvement claim.**
- Separation is partly by construction (balanced discriminating cells), so this is a
  specification/ablation result, **not** an empirical measured gain.

## Missing inputs (UNKNOWN; do not invent)
- L3 UNKNOWN: `sigma` (noise) for `korrekta reusebeslut andel`; observation vs latent state.
- L4 UNKNOWN: transfer/reduction error to the real recursive operator-/method-planner; consumer bound.
- Real solver-precision ↔ runtime calibration on the frozen runtime replay (cold + N warm).
- The two *named* unseen task families from the supplied protocol (here F1/F2 are declared synthetic).
- Physical tau and sampling/followup window — bind to supplied protocol before any empirical run.
- Complete cost model (setup+tuning+discovery+validation+update) for the candidate registry.

## Inherited/failed criteria preserved
- Native CERTIFIED/seed 7×, 0.938, ARI 0.177, DEFF 2.70×, robot gain are inherited claims, not new
  measurements or a general gate (BRIEF §60). Not re-measured here.
- Q005 is an assumption-based mechanistic scenario; nothing Q005 was relabeled as a physical result.
- Historical parent results are not v6 tests. Rate/network interruption is not scientific refutation.

## Next test (frozen proposal)
Bind a complete cost ledger on a frozen runtime replay and require: at equal eps_y and equal total
cost, candidate decision-loss reduction > delta_dec = 0.05, with F1/F2 held out. Reject on tie,
coverage loss, or invalid cost. See FOLLOWUPS.json.
