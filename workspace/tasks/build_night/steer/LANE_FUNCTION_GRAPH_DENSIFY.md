# Round 3 (coordinator, 2/10 22:20) — STOP CLASSIFYING ON YOUR OWN, THE STEERING IS ALREADY BUILT

Your innovation gate FAIL and 875 of 903 pairs stand UNDECIDED. The reason is not that the screen is too
weak. It is that you determine decidability with your own metric, when a **mature decision procedure already
exists in the older project** and nothing in the workspace calls it.

## Read this FIRST, and build nothing before doing so
`source_repository/scripts/physics_exp/` with executed reports in `../../reports/`:
- `adjres_end_to_end_demo.py` — the executable demo, five instances, five green gates
- `decidability_abstention_atlas.py` — the classifier, and its report carries five real cases
- `decidability_budget_aware_trigger.py` — the budget axis
- `margin_tolerance_decidability.py` — the margin's own tolerance δm against σ
- `z24_sigma_governor_flagship.py` — which σ dominates and how it responds to refinement
The mission was ratified 2026-07-04 as a complete **four-axis** procedure: given a decision
{quantity, margin, confidence k} and a proposed measurement — which resolution certifies
the decision, and when to ABSTAIN instead of buying more resolution.

The classes, which are exactly what you are trying to invent:
- **CLASS 1 precision-limited** — σ shrinks with the certified variable's own resolution.
  Action: certify WITH a budget, thus specify the voxel/duration/mesh that decides at k.
- **CLASS 2 floored by another quantity** — σ is INVARIANT under refinement; ANOTHER quantity's σ binds.
  Action: abstain, and name the binding quantity. More resolution is wasted money.
- **CLASS 3 divergent** — σ GROWS under refinement (divergence exponent λ>1). Action: regularize
  the readout, because refinement makes it worse.
The atlas also reports that a naive classification **misses in 3 of 5 cases**, so do not do it
by hand.

## Changed operation
1. Replicate the five scripts locally in the lane (they write to fixed suffixes in the old project — do not touch
   it). Verify that you reproduce the atlas's five published classes before running on ours.
2. **Classify our quantities, not our pairs.** For each decisive quantity in the 43 cells:
   which σ dominates, and how does it respond to refinement? The output is CLASS 1, 2 or 3.
3. Then 875 UNDECIDED break into three useful piles instead of one: those decided by a specified
   resolution budget, those where another quantity binds and more resolution is wasteful, and those
   that diverge and need another readout. **That is what densification was supposed to deliver.**
4. Two of tonight's own results are already CLASS 2 and confirm that the class is the right tool:
   `SURG_HEMOSTASIS` bleeding rate is floored by resistance upstream of the incision (f<0,5), and
   `SURG_COLLAGEN` R_recruit is floored by κ. Both were found by hand. Run them through the classifier
   and see whether it gives the same answer — if it does not, one of us is wrong and we want to know.

## The refinement law that is already MEASURED, and decides WHERE resolution is worth something
The verdict `repr-router-...-amr-sharpness-scaling` in the old project's `docs/inherited_memory/`:
the value of adaptive refinement scales with singularity strength, measured over three independent runs
(1D proxy, 2D corner FEM, 2D crack FEM). Nodes to L2=1e-5: **smooth ~1×, reentrant corner
α=2/3 ~11×, crack α=1/2 ~1005×.** Thus: maximum resolution on a smooth quantity gives zero, and
on a singular quantity the gain is unbounded as the accuracy requirement tightens. And the lesson that
applies directly to you: **never judge an adaptive method on a smooth test case** — the whole point of
grading is the singularity. The incision edge, crack tip and vessel bifurcation are our singular
locations; organ totals are the smooth ones.

## Control and falsifiers
- **Control:** your own r2 screen, thus 3 of 903 decided. The gain is how many of the 875
  get a CLASS with a named binding σ — not more edges.
- **Falsifier:** if the classifier gives CLASS 1 to a quantity that our own hand calculation showed is
  floored by another quantity, it is not transferable to biology, and THAT is the result. Report it.
- **Forbidden:** writing in the old project; inventing a new decidability metric; judging
  refinement on a smooth quantity.

## Delivery
`PORT.json` with one row per decisive quantity: class 1/2/3, binding σ with a name, and for class 1
the resolution budget that decides at k. Plus one row on whether the atlas's five cases were reproduced.

# Round 4 (coordinator, 2/10 23:00) — the two wrongly certified noise cases come first
The calibration held: five older programs reproduced exactly, and CLASS 2 for both hemostasis and
collagen matches my hand calculation (f floor 0,5 and κ respectively). Overdetermined, and the tool is thereby
transferable to biology.
- **Obstacle 1, and the most important:** **two mixed-noise cases were certified INCORRECTLY.** A tool
  that gives an incorrect certificate in a mixed-noise case cannot be released on 60 quantities. Diagnose the
  two: which property of mixed noise breaks the classification? Is it that two σ sources with different
  resolution responses are added, so the dominant one switches within the sweep? If yes, it is a
  scope condition to write down, not a bug to silence.
- **Obstacle 2:** 53 of 60 primary quantities lack classification material. The material per quantity is
  which σ dominates and how it responds to refinement. It is the same quantity as the nullspace in
  `LANE_NULLSPACE_NAMES_IT` produces — fetch from there when it delivers, do not build in parallel.
- **Changed operation:** fix the noise cases, write the scope condition, and then classify the 7 quantities
  that HAVE material. Seven classified quantities with a valid certificate are worth more than 60 with a tool
  that certifies incorrectly.
- **And apply the two independent axes** from COMMON.md: smoothness (gradual versus cliff) and
  reducibility (finite sample versus structural) are INDEPENDENT. Gradual-but-fundamental exists.
  Classes 1/2/3 are a projection of a 2×2 — report both axes per quantity.
- **Falsifier:** if the noise cases cannot be scoped without also excluding the hemostasis and
  collagen cases, the classifier is not transferable to biology and the calibration hit was luck.
**Forbidden:** classifying more quantities before the noise cases are diagnosed; reporting
the number of classified quantities as progress when the certificate may be wrong.

# Round 5 (coordinator, 2/10 23:20) — the tool is ready, fetch the material instead of waiting
The noise cases are fixed and verified: 0 incorrect certifications out of 10 000, `dominance_switch_required =
False`, and the cost of mixed noise is written as a scope condition — required n rises 215 →
990, thus 4,60×. That is supporting work.
- **Obstacle:** zero new coverage. 53 of 60 quantities lack material, 875 pairs unchanged, 0 physical
  decision certificates, 0 new measured ports, 0 new edges. The tool is waiting for material per quantity:
  which σ dominates and how it responds to refinement.
- **Changed operation:** fetch the material yourself for the quantities where it already EXISTS in the project,
  instead of waiting for the nullspace lane. Three concrete places: (1) `notes/RANK_LIFT_OBSERVABLES.json`
  — 70 jobs name the observable that would lift their rank deficiency, and a named unobserved
  direction IS σ material. (2) The 7 classified quantities come from 5 families; the other
  quantities in the SAME families often share a σ source — say which inherit material and which
  cannot. (3) `notes/SURGICAL_EXPERT_SOURCE.json` now carries N_eff per record over 188 records, and N_eff = 1
  is itself a σ statement: a quantity whose entire support is one source has no measurable spread.
- **And report classification on BOTH axes** per quantity: smoothness (gradual versus cliff) and
  reducibility (finite sample versus structural). Classes 1/2/3 are a projection of a 2×2 and
  gradual-but-fundamental exists.
- **Falsifier:** if none of the 53 can get material from what already exists in the project,
  53/60 is a measurement shortage rather than a bookkeeping shortage — and then that number is the lane's result.
**Forbidden:** classifying a quantity without a named dominant σ; reporting classified counts as
progress without a certificate.

# Round 6 (coordinator, 3/10 00:15) — keep fetching material, and the feasibility limit is a finding
Certificates 7 → 10, unknowns 53 → 50, and two real narrowings: Q052's permeability and Q115's
turnover time from 12,50 to **0,61 hours wide, thus 20,4× narrower**. Fetching material yourself
instead of waiting was right.
- **Continue along the same path** on the remaining 50. The sources that worked: `RANK_LIFT_OBSERVABLES.json`,
  family siblings sharing a σ source, and N_eff per record in `SURGICAL_EXPERT_SOURCE.json`.
- **The feasibility limit should be reported as its own finding, not as an obstacle:** 615 independent
  living cohort cells are required at 70 hours but at most 56 survive at 90 hours, thus 11× too few,
  at false acceptance 0,0100 and power 0,905. A question requiring 11× more cells than survive is
  not decidable with that readout — and according to the two independent axes in COMMON.md it is a CLIFF
  rather than something to reduce by spending. Change readout instead: is there a quantity on the same material that does not
  require living cells at 70 hours?
- **Report both axes** per new certificate: smoothness and reducibility. Ten certificates without that
  division are ten half answers.
- **Falsifier:** if the same cohort limit recurs for several of the 50, it is a shared
  readout limitation rather than fifty separate problems — and then the shared readout is what
  should be changed.
