# LANE_NULLSPACE_NAMES_IT — turn 181 silent rank deficiencies into named observables

Resultatmapp `results/LANE_NULLSPACE_NAMES_IT/`.

## The finding the lane builds on, measured and not guessed
Across **5 984** swarm reports, **794 (13,3 %)** make an explicit claim of structural
non-identifiability, against **49 (0,8 %)** complaining of parameter uncertainty — ratio 16 to 1. The forms:
392 non-identifiability, 261 rank deficiency, 134 "not identifiable from", 84 structurally
unidentifiable, 31 unbounded above, 19 proved symbolically.

Thus: the project's binding limitation across subsystems is **observability**, not computation and
not coverage. Every seventh job says that the quantity it was asked to compute cannot be determined from the
observables it received. (Strict matching on purpose — the word "structurally" alone gave 38,1 %, an
overestimate by 2,9×.)

Of the 261 rank-deficient ones: **70 name the observable that would lift the rank**, 9 only exclude
insufficient routes, and **181 are SILENT**. Harvested in `notes/RANK_LIFT_OBSERVABLES.md` and `.json`
by `tasks/build_night/rank_lift.py`.

## Why the silence is the valuable part
A rank deficiency is not a lack of data. It is a property of the observation set, and every
rank-deficient matrix has a **NULLSPACE**. The nullspace is precisely the direction in parameter space that
the observables do not see. A job that has found a rank deficiency therefore already has everything needed to
name what is missing — and 181 of them left it unsaid. That is 181 dead negatives that can become
181 named observables without a single new measurement.

And it reinterprets seed 1: the value of an external biology model is not that it is another model,
but that it delivers an observable we do not have. An external model that reproduces a quantity we already
compute lifts no rank. **The nullspace direction is the requirements specification for which external model
is worth connecting.**

## Do this
1. Read `notes/RANK_LIFT_OBSERVABLES.json`. The 70 that name a lift are your CALIBRATION: for three of
   them, compute the nullspace from the report's own quantities and check that your direction
   corresponds to the observable the job itself proposed. If it does not, your method is wrong, and you
   must discover that on the 70 and not on the 181.
2. Then take the 181 silent ones, in the order `rank_lift.py` ranks them. For each: set up
   the observation matrix from the report's own declared quantities and observables, compute
   the nullspace, and **translate the direction into a physical quantity with a unit**. A direction in an
   unnamed parameter vector is not an answer; "the combination of k_on and K_off that a
   steady-state measurement does not distinguish" is an answer.
3. **State the order of the nullspace.** Dimension 1 means a single added independent observable
   suffices. Dimension 3 means three. That number is what an acquisition list needs, and it exists today
   for one (1) job out of 261.
4. **Distinguish liftable from structural.** A direction that no physically measurable observable can see is
   not an acquisition entry but a reformulation of the question. Use the two independent axes from
   `tasks/build_night/COMMON.md`: smoothness versus reducibility by spending. Gradual ⇒ pay down. Cliff ⇒ change
   readout.
5. Run the already graded tools instead of inventing: `decidability_abstention_atlas.py`
   and its siblings in `source_repository/scripts/physics_exp/` (replicate within the lane, do
   not write in the old project). CLASS 2 there is the same as a rank deficiency whose nullspace no
   finer measurement of the same quantity can see.

## Strongest control and falsifiers
- **Control:** the 70 jobs that named their own rank-lifting observable. Your method must reproduce
  their answers before it is applied to the silent ones. Report hits on 3 of 3, or which ones missed.
- **Falsifier:** if the nullspace for the silent jobs systematically CANNOT be translated into a
  physical quantity with a unit, the rank deficiencies are artifacts of how the jobs set up their matrices
  and not statements about the body. That would be a substantial finding about our own job generation and must
  be reported directly.
- **Forbidden:** ordering a measurement that is not derived from a nullspace; merging several
  jobs' deficiencies into a common statement without showing that they share a direction; writing in the old
  project.

## Delivery
`PORT.json`: one row per job with the nullspace dimension, the direction translated into a physical quantity with
a unit, classification as liftable or structural, and for the liftable ones the smallest set of independent
observables that suffices. Plus the calibration row: 3 of 3 against the self-named ones, or which ones failed.

Inga interna data. Allt PENDING_INDEPENDENT_REVIEW.
