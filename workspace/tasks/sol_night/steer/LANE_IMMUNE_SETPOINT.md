# Steering LANE_IMMUNE_SETPOINT — round 1 (coordinator, 2/10 20:00)

New lane, and the only one of Anton's seeded directions that lacked a deriving lane. It takes the slot from LANE_DEJ_BOND_INVENTORY, whose premise fell in the first round.

- **Obstacle:** underactivation and overactivation are treated as two sides of a scale with a point between them, and nobody has tested whether it is instead a region where both failures are possible simultaneously.
- **Changed operation:** derive the set's SHAPE from a mechanism with published parameters, instead of looking for a better threshold.
- **Strongest control:** the usual threshold description calibrated on the same data.
- **Falsifier:** if the set where both failures are met is empty for every published parameter combination, the point picture is correct, and that should be reported as clearly as the opposite.

## Start with the move that has already worked for us
A swarm job today used CRP's half-life invariance, 18,8 ± 3,9 h constant across inflammatory states, as a conservation constraint and rejected **four of seven arms** in an acute-phase model — all four were clearance-channel arms. An invariance thus rejects an entire mechanism class without measuring the mechanism. Look for such an invariance first; it is cheaper and sharper than a new parameter sweep.

## Three error classes that cost us today, all worth avoiding here
1. **The wrong regime gives correct arithmetic on the wrong question.** One of today's biggest claims fell because a requirement was calculated in a saturated regime and compared against a value in a linear one. Before comparing two numbers: say which regime and denominator each is in.
2. **Pooled is not paired.** Another lane showed that a gap of 2,9–9,1× to 2,17–7,76× was an artifact of amount and activity coming from different individuals. If you compare two quantities, say whether they were measured in the same individuals.
3. **Structural absence is not independence.** A third lane found that 25 of 27 anchors do not reverse sign under a physiological modifier — because the modifier is not even a port in the model. Mark such cases as untested, never as robust.

## If a path is missing where you run
`~` is not available on cloud hosts. Missing: write `missing_prerequisite` of the first subparagraph by: RESULTS.md and work from what is in the lane letter, labeled as quoted. Three jobs today drew the wrong conclusion by quietly reconstructing missing input, and one of them cost an incorrect headline.

Allt PENDING_INDEPENDENT_REVIEW. Inga interna data, inga patientdata.

# Round 2 (coordinator, 2/10 20:50) — make the comparison INDEPENDENT of dose

## Review of round 1: you answered the structural question and refused to call it a finding
`persistent_joint_failure_area = 3,7365` model units/h, thus a nonempty region where underactivation and overactivation persist simultaneously, with BOTH as a realizable class. The theorem identity held in 100 000 counterexample tests without a violation. And you still set `innovation_gate = TIE` and `biological_shape = UNKNOWN`, on the grounds that coexistence in a model is already known. That is correct: you could have used the headline "setpoint is a region" and did not.

## The assignment: 2000× in dose is no comparison, so build one that does not depend on dose
Your own blocker is named exactly: `external_assay_MOI_ratio = 2000`. A model and an assay differing 2000× in multiplicity of infection cannot be compared in level. But they can be compared in **shape**.

1. **Find a dimensionless combination that is invariant under dose.** A ratio between two outcomes in the same experiment, an exponent, one timescale divided by another timescale — something whose value does not move when MOI changes 2000×. Test several and report which are invariant and which are not; showing that a candidate ratio is NOT dose-independent is equally useful.
2. **Test the shape, not the level.** Does the region where both failures persist remain in the dimensionless description? If yes, the statement is independent of the dose mismatch that blocked round 1. If no, the region was an artifact of dose choice, and that is a stronger result than a TIE.
3. **The invariance move again, now as a tool rather than an example.** CRP's half-life invariance rejected four of seven arms without measuring the mechanism. Look for an invariance in immune response dynamics that rejects a class of setpoint descriptions — a quantity that should be constant across dose or state if a particular model class is correct.
4. **Declare the MOI regime for every number you compare**, just as the regime must be declared for concentration versus K_m. That error, saturated versus linear regime, cost us today's biggest claim.

## Control, falsifiers, prohibitions
- **Control:** round 1's comparison of levels, which was blocked by the 2000× mismatch. The gain should be a statement that survives the dose difference.
- **Falsifier:** if no dimensionless combination is invariant under MOI over the tested interval, the model and assay are not comparable at all, and then the deliverable is the assay design that would make them comparable. That is a fully adequate outcome.
- **Forbidden:** comparing levels across 2000× in dose; calling model coexistence a biological finding; calibrating against the outcome metric; treating a clinical threshold as a measurement.

# Round 3 (coordinator, 2/10 23:00) — the source conflict is the finding, not the obstacle

## Review of round 2: you made it dose-independent, and found an error in the publication
`dimensionless_local_joint_rate_area = 1,3194`, thus the region where both failures persist remains in the dimensionless description — it survives the dose mismatch. Terminally: BOTH 27, NEITHER 34 of 183 trajectories. The identity held in 100 000 instances and the maximum corrected assay error is 6,1e-13.

And: **`source_caption_Methods_conflict = True`** with `actual_external_MOI_ratio = None`. The publication's caption and Methods contradict each other about multiplicity of infection, so the number I asked you to compare against is not unambiguous in the source. That is a finding about the source, not an obstacle on your part.

## The assignment
1. **Document the source conflict so it can be cited.** Which section says what, verbatim, and which number is compatible with the article's own other numbers? A publication that contradicts itself about its dose is external ground truth that cannot be used directly, and we should be able to show that.
2. **Make the statement independent of the conflict.** You already have the dimensionless form. Does the region hold for BOTH interpretations of MOI? If yes, the conflict is irrelevant to the conclusion, and that is the strongest possible outcome.
3. **Search the index.** 1309 completed cells are in `notes/OLD_CELL_INDEX.json`, of which 431 REFUTED. Run `--quantity immun` and `--quantity uptake` — an already rejected question should not be rerun.

**Forbidden:** choosing the MOI interpretation that gives the desired answer; calling model coexistence a biological finding.

# Round 4 (coordinator, 2/10 22:45) — which initial condition is physiological?
The source conflict is documented: Methods 0,005 versus caption 10, factor 2000. Good. And the region exists for some initial conditions and not others — three class changes between `quasi_steady`, `dynamic_M0_zero` and `dynamic_M0_2_or_2_5`.
- **Obstacle:** the question has moved from dose to initial state, and we do not know which is physiological.
- **Changed operation:** determine which initial condition real tissue has. M0 = 0 means no macrophages present before the stimulus, which is an assumption rather than a fact — tissue-resident macrophages exist before any stimulus. Look for published resident density with a unit, and say which of your three cases it corresponds to.
- **Control:** the three cases as they stand, thus without designating any of them as physiological.
- **Falsifier:** if published resident density spans all three cases, the region's existence is not decidable from the initial state either, and then that is the missing measurement.
Forbidden: choosing the initial condition that gives the region; using either MOI interpretation without reporting both.

# Round 5 (coordinator, 2/10 22:40) — convert density to M0, it is the only gap left
You have the strongest outcome: `M0_001_02_05_both_doses = BOTH`, thus the region exists under BOTH
MOI interpretations for M0 at 0,01 / 0,2 / 0,5, and the publication's 2000× contradiction is thereby
irrelevant to the conclusion. And you fetched an external anchor: 1,73e6 ± 0,11e6 cells in 52,14 ± 4,02 µL
= 33,18 million per mL.
- **Obstacle:** `physiological_M0 = None` and `density_uncertainty = None`. The density has been fetched but not
  converted to the model's M0, so the anchor carries no weight yet.
- **Changed operation:** make the conversion explicit. Which volume or area normalizes the model's M0,
  and what is 33,18e6/mL in that unit? Propagate BOTH SEMs — 0,11/1,73 and 4,02/52,14 — to
  an interval in M0, not a point value. Then say whether that interval lies within 0,01–0,5, where
  the region exists under both doses. **If yes, the conclusion is anchored in real tissue rather than in a
  chosen initial condition, and that is the lane's whole point.**
- **And explain the 5902:** `scalar_uptake_false_clearance = 5902` are cases where scalar uptake is
  read as clearance. Is it the same error class as the kidney certificate's flow limitation? If yes, it is
  a shared failure mechanism across two lanes and worth its own line.
- **Falsifier:** if the propagated M0 interval extends outside 0,01–0,5, the region's
  existence is not anchored in the measurement, and then that is the missing measurement.
**Forbidden:** choosing an M0 within the interval that gives the region; reporting a point where two
SEMs exist.

# Round 6 (coordinator, 2/10 23:05) — the model's UNIT is the obstacle, not the data
You did the right thing and got a structural answer: `density_to_M0_calibration_gate =
FAIL_MISSING_FUNCTIONAL_NORMALIZATION`, `source_M_unit = PHENOMENOLOGICAL M-unit; no physical
count/volume normalizer`, and thereby `functional_M_identified_set = [0, +infinity)`. An external
anchor that cannot reach the model's unit carries no weight. The measurement exists — the model cannot receive
it.
- **This is tonight's system finding in your subsystem.** Measured over 5 984 swarm reports, 13,3 % make an
  explicit claim of structural unidentifiability versus 0,8 % complaining of parameter uncertainty.
  Your M unit is the same thing: a direction the observables do not see.
- **Changed operation, in this order:** (1) search ALL FOUR pools (see COMMON.md) for a
  published normalization of the same model family's M unit — somebody else has compared the same model to
  cell counts. (2) If none exists: rewrite the requirement as a unit specification. Which quantity with which
  unit would the model need to declare for 33,18 million cells/mL to become a valid M0? You have
  the interval 2,08e-10 to 7,85e-09 mL per cell and hour already — say what that coefficient IS
  physically, not just what it must lie between.
- **Keep what you have gained:** the region exists under BOTH MOI interpretations for M0 at 0,01–0,5, and
  the publication's 2000× contradiction is thereby irrelevant to the conclusion. That statement stands independently of
  the unit problem and should not be withdrawn.
- **The 60 of 10 000:** 0,6 % counterexamples against instantaneous clearance notification have the same form as the plug job's
  finding tonight — a scalar is almost always enough and fails for a measurable fraction. Say whether the 60 lie
  in a recognizable region of parameter space or are scattered. If they cluster, the cluster is
  the statement.
- **Falsifier:** if no published normalization exists and no unit specification can be
  formulated without choosing a free parameter, the model's M unit is unusable for external anchoring, and
  that should be said plainly as a requirement for a reformulated model.
**Forbidden:** choosing a normalization constant that makes the anchor fit.
