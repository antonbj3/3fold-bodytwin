BT-FW48-SEED-493

# RESULTS — Autocomposition seeks a type-/proof-/cost-valid operator DAG

Status: **PENDING_INDEPENDENT_REVIEW**. This is a *definition/review* job. The scientific
claim (does autocomposition improve a real molecular-to-organ error budget) is
**NOT_EXECUTED_DESIGN_REVIEW**. What *was* executed is the planner/instrument layer as the
brief's own design standard prescribes for software/planner cases: a synthetic mechanism test
on a frozen replay round with **declared** model scales, N in {1,10,100} as *planned* query
counts, on the supplied but **unaudited** Q005 model.

## 1. Actual result (first measurement)

Four deterministic questions over three families, exact reference targets
(`frozen reference procedure`, PREREG §1), on the bit-exact Q005 reproduction
(`audit/audit_findings.json` all-PASS):

| question | family | observable / unit | exact reference target (`p*`) |
|---|---|---|---|
| q1 | R1 | `CL_Cr` / mL min^-1 | `gfr_mL_min` |
| q2 | R1 | `CL_Met` / L h^-1 | `K_OCT2_ref_uM` |
| q3 | B1 | joint normalized error, budgeted **set** | `{Jmax, K_OCT2}` (cost 3) |
| q4 | D1 | transient peak cell conc. / microM | `P_MATE_cm_s` |

Interface-validity predicate (kernel of the composed map, derived, not fitted):
`R1_Cr = {P_MATE_cm_s}`, `R1_Met = {P_MATE_cm_s}`, `D1 = {gfr_mL_min}`.

Valid-solve count / 4 (abstentions recorded separately; primitive output = the count):

| arm | N=1 | N=10 | N=100 | operator invocations @N=100 | comm. bytes @N=100 |
|---|---|---|---|---|---|
| `reuse_only` (parent 2 alone) | 1 | 1 | 1 | 0 (replay) | 0 |
| `manual_strong` (cold additive pipeline) | 2 | 3 | 3 | 60 | 9.60e6 |
| `greedy_control` (additive, no predicate/joint) | 2 | 3 | 3 | 60 | 9.60e6 |
| `autocompose_no_joint` (predicate ON, joint OFF) | 2 | 3 | 3 | 44 | 7.04e6 |
| `joint_control` (**strongest conventional**: exact joint, no predicate) | 3 | **4** | **4** | 50 | 8.00e6 |
| `autocompose` (predicate ON, joint ON) | 3 | **4** | **4** | **38** | **6.08e6** |

## 2. Decision under the frozen rejection criteria

- **R-J1 no separation — FAIL.** `autocompose` (4/4) is **not strictly greater** than the
  strongest comparator `joint_control` (4/4) at equal budget. There is **no accuracy
  separation** against an equally informed conventional planner that evaluates the joint
  effect directly.
- **R-J2 ablation — PASS.** `autocompose` 4 > `autocompose_no_joint` 3: the *joint term* is
  the operative mechanism; removing it costs the B1 question.
- **R-J3 stronger baseline — not triggered** (no comparator strictly exceeds 4/4).
- **R-J4 coverage — PASS** (every arm names a target or records ABSTAIN; coverage 100 %).
- **R-J5 risk — PASS** (0 wrongly named targets for the candidate on any family).

**Outcome: the candidate is NOT an accuracy comparative improvement.** Its only measured advantage is
**strictly lower operator cost at equal accuracy**: 38 vs 50 operator invocations and
6.08e6 vs 8.00e6 communication bytes at N=100, because the interface-validity predicate
prunes the kernel parameter `P_MATE` (and `gfr` on D1) *by proof, without paying to evaluate
it*. Under the brief's second branch of the comparison contract — "at the same total cost the
decision loss must fall, / at equal eps_y total cost must be strictly lower" — this is a
**cost-only, accuracy-matched capability**, not a stronger answer.

This is a small, honest **negative finding on the headline claim** (auto beats manual on
validly solved questions) plus a **positive but narrow capability** (proof-based kernel
pruning saves evaluated queries at identical accuracy).

## 3. Why additive search failed — and a corrected direction of the frozen hypothesis

B1 exact set scores (full pool). The frozen PREREG §5 predicted **sub**-additivity and
"an additive scorer systematically over-values multi-parameter sets". The measurement shows
the opposite direction:

- `{Jmax, K_OCT2}`: exact joint reduction **0.25078**, additive sum of marginals **0.02377**,
  `sub_additive = false` — **super-additive (complementary)**, because K and J are
  log-correlated at 0.70 and refining K alone barely moves the joint error
  (`K` marginal = **-0.01928**, i.e. *negative*).
- `gfr`: additive marginal **0.16040** (passes the 0.10 gate) but exact joint reduction
  **0.05874** (fails the gate) — additive scoring inflates it.

So the additive scorer fails for **two** reasons (under-valuing a complementary pair and
over-valuing a marginal), and it lands on ABSTAIN on q3; the direction of the non-additivity
is the *opposite* of the frozen expectation. The **conclusion** (additive marginal scoring is
insufficient; an explicit joint term is required) survives; the **direction claim is
falsified**.

## 4. Amendments to the frozen PREREG (listed, not silent)

1. **Arm 4 degenerate as written.** PREREG §6 defines `autocompose_no_predicate` as
   "predicate and joint term removed" — identical to `greedy_control`. It was replaced by
   `autocompose_no_joint` (predicate ON, joint OFF) so the *operative* mechanism (the joint
   term) is actually ablated. Original wording preserved in `PREREG.md`.
2. **Added `joint_control`** (predicate OFF, joint ON, exact). R-J3 requires the strongest
   fair conventional planner; the frozen arm set did not include a joint-evaluating control
   and would otherwise have measured a strawman.
3. **PREREG §5 direction of non-additivity corrected** per §3 above.

## 5. What remains incomplete / missing inputs (unchanged, NOT simulated into existence)

- **Physical `tau` and the real Q005 sampling/follow-up window are UNKNOWN.** The design
  window `[0, 10 tau]` is a *declared model scale*, never a measured constant. Empirical
  execution stays blocked.
- **`BT-IM1-OBSERVABILITY`** — the observation instrument for "count of validly solved
  questions" is not bound; the count here is a planner-layer primitive, not a physical
  readout.
- **B1 is constructed by this job**, not supplied by Q005. Q005 defines only single-target
  one-at-a-time scoring; the budgeted set-selection family with the joint term is this job's
  extension and its status as an admissible "question" is a modeling choice.
- **Q005 is unaudited** and reports `NOT_ADMITTED`; all targets inherit that status.
- **No EXTERNALLY_MEASURED leaf created.** No physical, biological or clinical validation.

## 6. Unit definitions preserved verbatim

`CL_Cr` mL/min; `CL_Met` L/h; secretion mL/min; transient peak cell concentration microM;
declared planning cost units (K=2, J=1, f_u=1, P_MATE=2, gfr=3); valid-solve count is an
integer number of questions. Cost columns (cpu_s, peak bytes, comm bytes, operator
invocations) are reported separately and **never rolled up**.

## 7. Failed criteria and negative evidence preserved

- Original Q005 failures retained unchanged: no confirmed winner; `f_u` **misses** the 10 %
  gate at 9.397 %; `K_OCT2` negative; `P_MATE` unidentifiable at steady state;
  `corr(CL_Cr, CL_Met)=0.585363`.
- Frozen headline (auto beats manual on validly solved questions) is **not supported** by
  this executed instrument. The rate interruption of the prior attempt is an infrastructure
  event, **not** a scientific refutation.
