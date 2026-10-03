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

# Round 2 (koordinatorn, 2/10 23:35) — scale from 20 till alla 181, the method is calibrated
Kalibrering 3 av 3 against the self-named, and a real lift: Q168 from rank 7 till 14 med sju
extra-independent shellers; 3 000 exact counter-examples without a ranking identity counter-example.
Holds.
- **Hinder, ur din egen utdata:** 20 av 181 rows have an estimated zero room dimension; 157 Wait, 156
  is not reconstructed; and `measured_application_gain = None`It's... 11 % after a game.
- **Changed operation:** automate reconstruction instead of doing it per job. 156
  which could not be reconstructed — Group them according to WHY. Missing observation matrix completely,
  missing units, or is the parameter vector nameless? Each group has its own action, and the groups'
  size tells you where the work is. Report the distribution before continuing per job.
- **And deliver `measured_application_gain` pr MINST ett fall.** A lift from 7 till 14 is a
  ranking statement; winnings are what becomes decisive. For Q168: which quantity could not be determined at
  rang 7 and can be at rank 14? Name it with unity. Without that number, the lift is a linear algebra exercise.
- **De 24 counter-examples of population source headlines** This item shall be reported separately by which number they:
  felling — a counter-example against a headline number is a find of its own and should not be buried in a bill.
- **Falsifierare:** om de 156 cannot be reconstructed for a common cause that is in OURS
  job generation and not in jobs, are the rank deficiencies in part artifacts of how we ask questions,
  And that's the heaviest possible find here.

# New round (coordinator, 3/10 00:25)
De 94 which requires source binding of an executable operator is: EN sorts arbete — automate that binding and take all 94 In one round, instead of job for job. And deliver the application value for three FLER cases you did for P (256× → ±3,76 %): that number is what an acquisition list needs. Falsifier: if the source binding cannot be automated for more than half of the 94, lies the obstacle in job generation and not in jobs.

# New round (coordinator, 3/10 00:55)
Two things in your output weigh heavier than the progress figures.

HILL-LAGEN FIELD I ALLA TRE MOTFALL. Hill_law_adverse_cases = 3 and Hill_law_counterexamples = 3.
A law that falls in any case against which it is tried is not a law in our context, and it is used
in several places in the directory. It goes first: what condition in the Hill form is broken, and in how
many of our cells is the condition assumed without being tested? Measure the distribution over
the cells before it becomes a system statement. If the Hill adoption is unfounded in a significant proportion, it is
a heavier find than the entire rank work.

NOLL STARK FYSISK AXEL BUNDEN, trots 66 The difference is the whole point: tying up a job
operator to its source is not tying it to a physical axis with unity, and only the latter does
a zero room translated into a measurement. Changed operation: take the 50 interdependent and say for each
what is missing for the axis to become physical. If the unit is missing, missing the name of the quantity, or missing
the link to a measurable observable? The size of the groups tells you where the work is, just like yours
egen gruppering av de 156 gjorde.

And qualify the win honestly: dissociation constant goes from 16× indeterminate to a conditional
intervall BARA if the number of binding places is independent known. It is a new requirement, not a free
win. Say if that requirement is met anywhere in our four pools, or if it is a new
acquisition item.

Falsifierare: om ingen av de 50 is capable of receiving a physical axis without a new measurement, the zero-space method is not one of:
road to determinability without a path to a measurement order list. It is still valuable, but
It should be said cleanly and it changes what the lane is.

# New Round (koordinatorn, 3/10 01:55)
You corrected me and the correction is better than my statement: it was a FIXERAD n = 1 who fell, not
The Hill family, and with free n 1 000 fall noll motexempel. Jag har skrivit om min loggrad.

Changed operation, and it follows from your own distribution: Hill mentioned in 365 av 6 082 snapshots med 271
The question is now how many of them 271 som antar n = 1 WITHOUT to have tried it.
The cost is measured — error Hill shape distorts semi-saturation constant +208 % in a case and
−37 procent i ett annat — So each hit is a number that can be wrong with a factor of two or three.
Count the number, it is a lexical search plus a check of whether n is free in the code.

And build on the one thing that succeeded: a physical axis bound without remeasurement.
technology that worked there, and try the same technology on them 16 OPERATOR_INSTANCE-axlarna. De 18
However, the interference and assay axes shall be explicitly depreciated: — they are not physical quantities; and
should not be counted as failures.

Falsifierare: om ingen av de 16 the operator's control shafts can be tied by the same technique, each
the hit is a special case and the method is not peeling.

# New Round (koordinatorn, 3/10 05:50) — makes the sign-teaching general, it's bigger than your case
Du fann att signerade designer passerar precisionskravet i 25 av 30 mot 18 av 30 unsigned; and
det ogynnsamma intervallet smalnar 1,673 administration — pr 3 mikromol vittnesmaterial mot 627 standard.
Sign information is therefore almost free and is normally discarded.

Changed operation: test whether the lesson applies to the ANDRA the zero rooms and not just this. 55 rader
with calculated zero-room dimension. For how many of them would a signed reading reduce the dimension
or narrow the interval? It's a cheap check per line and the answer is a distribution, not a
speech. If character information helps in a large proportion, it is the single most cost-effective
the recommendation night produced, and it belongs in the acquisition list.

And connect to notes/ACQUISITION_TARGETS.json, as I wrote last night: nine acquisition entries in node form
with quantity, drive and keywords. Tell each of them if a signed reading is possible and what it
would buy. Several of them request magnitude where the sign may be the one carrying the information —
especially the paired GFR- and the secretarial post, where it is the correlation TECKEN That's the decision.

Falsifier: if signed reading only helps for binding constants and not for flows
or arenas, the lesson is narrow and should be reported as such rather than as a design principle.
