# BT-MAT-OBS0 — what material measurements actually differentiate the model branches?

Status: **mission prepared, not started**. Sol6 medium, no own sub-agents. Output:
`results/BT-MAT-OBS0/` ; sources and older results are read only. Ordered
as part of Anton's continued innovation hunt; this is a model study.

## Read before the trial plan is locked

- Gemensamma `AGENTS.md`, BodyTwin `START.md`.
- `results/C4/{README.md,RESULTS.md,c4_results.json,c4_budget.py}`.
- `results/X3c/{PREREG.md,README.md,RESULTS.md,model_x3c.py,run.py,sens.py,sensitivity.json}`.
- `results/BT-C1b/` for the contract around correlation, units and descent.
- `results/BT-IM1/RESULTS.md` and `BT-IM2/RESULTS.md`: more functional data did not automatically provide more geometry information.

## Starting point from our work

C4's known damage variance is dominated by the choice of density–modulus relation, but
unknown boundary conditions, calibration errors, and model errors remain. The X3c's
H1 character toggles within the sensitivity box. The previous simple marker metrics
separated the branches weakly. Calibrated q values ​​are partly outside the box that
was sampled. It justifies an identifiability test before a major optimization.

Note the concrete source: `X3c/sens.py:eval_row` builds marker features like
`log(marker[A+B]/marker[A])` for the same parameter set. These are paired counterfactual
model outcomes. They are not automatically two observable measurements in
the same person. Group contrasts or time sequences require separate contracts
for individual differences, time, treatment, and residual effects.

## Hypothesis and decomposition from first principles

Hypothesis: a combination of mechanical response and time-resolved markers
can distinguish some material/repair branches where a single marker cannot.
It may also prove to be impossible under the information available.

Start in a pronounced simplified elastic beam model,
not in a new clinical statement: `E = a (rho/rho_ref)^n`
, `epsilon = F l c / (E I)` . All terms must have
declared units. Logarithmic differentiation gives

`d log epsilon = d log F + d log(l c/I) - d log a
                 - log(rho/rho_ref) d n - n d log rho`.

At the same density and the same load shape, a, n, load scale and geometry
factor are partially mixed. Multiple repetitions of the same measurement operator
do not remove its exact null space. A load ratio for linear response can
also eliminate the very absolute stiffness we wanted to learn. Derive and
test this before an information matrix is ​​given an interpretation.

Mechanism → material law/repair → calculation operator → observable signal
→ representation/units → explicit physical and statistical prerequisite.
Have each branch end in verified algebra, a runnable test, or UNKNOWN.

## Limited delivery

1. **Observation contract.** List each proposed measurement with what
   actually observed, when, in which arm and on which sample/individ. Separation
   absolute measurement, repeated measurement, group contrast, and counterfactual oracle.
   Do not assign independent noise or precision where there is no basis.
2. **Rank and nullspace.** Derive low-dimensional counterexamples and check
   Jacobianrang/SVD with scaled variables. View a couple of different parameters like
   gives the same allowed signal but different targets, if any. Handle discreet
   material laws as alternative branches, not as independent extra measurements.
3. **Lite frozen model sample.** Determine before initial evaluation whether existing
   X3c data is enough. The 640 Saltelli rows are not independent training cases:
   keep all A/B/AB-families together at each split. They saved
   marker ratios are only good as an explicit oracle test. An executable real
   observation view requires separate arm output; do not fabricate these out of quota.
4. **Baselines and countersamples.** Compare prior without measurement, best individual
   allowed measurement and at most one small preselected pair. Kalibrering/val and
   evaluation is kept apart. A duplicate signal with fully shared error is not allowed
   counts as new information. Perfect noise-free oracle is an upper reference.
   Unknown correlation shall yield refusal or a declared sensitivity range.
5. **Decision.** Report whether a measurement breaks a degeneracy or just
   reduces noise within the same branch. Remaining model equivalence must be shown
   explicitly. Any improvement is within the simulator; it does not validate
   X3c's biological law or any treatment.

## Stops and resources

Lock `PREREG.md` and hash before new numerical evaluations. This document
is a brief, not a complete pre-registration. Primary subtask is contract
and identifiability; no new FE runs, downloads or experiments on people
needed. CPU2, nice19, under 2 GB, max 100 MB new output. If new ODE runs
are needed: a timed pilot case, then report budget before a larger
matrix. Do not reuse old parameter distribution as a calibrated posterior.
No clinical recommendation, graph upgrade or publication.

Supply `RESULTS.md` , `observations.json` , `nullspace_checks.json`
, any small run, and the exact next experiment. Separate
agent reviews before the results are posted as holding.
