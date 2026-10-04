# Direction LANE_ROBOT_PORT_CONNECT — round 1 (coordinator, 2/10 18:30)

New lane. It takes over the slot from LANE_CONSTANT_TO_MECHANISM, which closed in round 5 with its main gap resolved.

- **Obstacle:** the robot branch validates cleanly but only 7 of 202 ports connect to a quantity we compute, and all 7 are on the tissue side. Every robot-side port hangs loose, and 80 port rows lack a value entirely.
- **Changed operation:** connect instead of expanding. Every robot-side port should either carry a computed quantity with unit and origin, or be declared `REQUIRES_NEW_COMPUTATION` with the missing computation named.
- **Strongest control:** the branch as it stands, 41 validating nodes with 7 connected ports. The gain is measured in the number of robot-side ports carrying a computed number afterward.
- **Falsifier:** if no robot-side port can be connected, the branch is premature, and the delivery becomes the list of what we must compute first. That is a fully valid outcome.

## Start with the lever, not the list
Registration time is the binding term: 150 s per registration is 99,998 % of the model update chain, and the acquisition rate is 0,002 %. Start there. A connection concerning registration time is worth more than twenty concerning tremor, because tremor contributes 1,0 % against the target budget 1,08 mm while brain shift contributes 195 %.

## Two traps that already cost us today
1. **Bench-only requirements are not robot gaps.** Our measurement spec is an ex vivo bench experiment with matched adjacent specimens and before/after calibration of the same tool. A machine can neither fulfill nor miss these. Do not count them as gaps — there was an entire class of apparent gaps in the material.
2. **The wrong regime gives correct arithmetic for the wrong question.** Today one of our largest claims failed because a requirement was calculated in the saturated regime and compared against a measured value in the linear one. Before comparing two numbers: state which regime and quantity each belongs to, and whether the denominator is the same.

## The file you need is not on cloud hosts
The branch is in `external_research_path`. If `local_path` does not exist where you run: state this as `missing_prerequisite` in the first paragraph of RESULTS.md and work from what is actually in the lane brief, with every number marked as quoted from the brief and not from the branch. Three jobs today drew the wrong conclusion by silently reconstructing missing inputs, and one of them cost an incorrect headline.

Allt PENDING_INDEPENDENT_REVIEW. Inga interna data.

# Round 2 (coordinator, 2/10 19:00) — a single port, the one that binds

## Review of round 1: the falsifier fired, and that was the right answer
`robot_computed_bindings = 0`, `requires_new_computation = 167`. The branch is a blueprint and this is now quantified. The labeling is valuable: 152 ports are the robot’s own PRECISION, 40 moved to the intelligence engine, and 15 were marked CONTEXT_OR_OBLIGATION_NOT_ROBOT_GAP — i.e. bench requirements that a machine can neither fulfill nor miss. Separating them removes an entire class of apparent gaps. Good also that you refused to sell numerical exactness as TRE calibration.

**167 computations are not a Sol task and must not become one.** They go to the swarm. Your round does one thing.

## Task: compute the decomposition of registration time, and only that
It is the binding term and the only one whose improvement moves the capability: at 150 s per registration, acquisition is 0,002 % of the chain, so 312 → 6417 Hz buys 3 ms out of 150 s while 150 s → 1 s is 149,5×. Round 1 did not compute it.

1. **Break down the 150 s into terms with a source each:** preprocessing, similarity evaluation or feature extraction, optimization iterations and their count, any deformable field solution, and data transfer. State each term’s provenance class, and that a vendor number is not a measurement.
2. **Separate the ALGORITHMIC lower bound from the implementation.** For the dominant term: derive the operation count from the problem size, i.e. number of voxels or features times iteration count, and show the arithmetic. A time that drops 150× with better implementation is a different claim from one requiring another algorithm, and only the second changes what we must build. Say which it is.
3. **Check whether 150 s is per registration or per procedure**, and whether any human step is included. A work cycle sold as a frequency has already been found once in this material.

## Control, falsifier, prohibitions
- **Control:** the quoted acquisition improvement 312 → 6417 Hz, i.e. 20,6× on a nonbinding term. The gain should be a decomposition of the term that does bind.
- **Falsifier:** if the dominant term is already at its algorithmic lower bound for the stated problem size, 150 s → 1 s requires another problem formulation and not faster code. Then name the formulation change — that is a stronger result than a decomposition.
- **Prohibited:** quoting a vendor speed as measured; comparing against another anatomy or problem size without saying so; treating a GPU port as an algorithmic bound; touching the other 166 ports this round.

# Round 3 (coordinator, 2/10 19:30) — what gets past 6,9×?

## Review of round 2: you corrected the ambition level, and that is the finding
124 s solve + 21 s resampling = 145 s, solve fraction 0,8552, and `solve_only_max_total_gain = 6,9048` — I checked 145/21 and got exactly the same number. **The target 150 s → 1 s requires 149,5×, and the dominant term can give at most 6,9×.** It is a hard bound of Amdahl’s form and says that faster solve is insufficient, however good it gets. Good also that you confirmed per registration and refused to guess the algorithmic floor when N, K, optimizer and timers are unreported.

## Task: attack resampling, because it is now the binding term
After your decomposition, 21 s resampling sets the ceiling. The question is not how it is optimized but what makes it **unnecessary**.

1. **What is resampling for?** Name the quantity it produces and why it must be recalculated per registration. If it recalculates the same image onto a new grid every time, the question is whether it can be represented once and transformed instead — i.e. a formulation change, not an optimization.
2. **Calculate Amdahl for three formulations separately:** (a) solve → 0, resampling remains: 6,9×; (b) resampling → 0, solve remains: calculate it; (c) both reduced by a factor each: which combination gives 149,5×, and is any required factor physically unreasonable? The last is the important part — if (c) requires both terms to drop 20× each, say so, because then 1 s is not a target but a wish.
3. **Test whether 1 s is the right target.** What is actually required of the update rate to track brain shift, which contributes +195 % against the target budget 1,08 mm? If brain shift develops on a timescale of minutes, 10 s may suffice, and then the entire 149,5× requirement is self-imposed. **Calculate it from brain shift’s own timescale, not from a round number.** It may be the single most valuable computation in the entire robot branch.

## Control, falsifier, prohibitions
- **Control:** the 6,9× ceiling from round 2. The gain should be either a formulation that gets past it, or a demonstration that 1 s is not the required target.
- **Falsifier:** if brain shift’s timescale requires subsecond updates, 149,5× remains a requirement and no formulation we know reaches it. Say so directly — that would be a hard bound and an important result.
- **Prohibited:** quoting a vendor speed; assuming an unreported problem size; touching the other 166 ports; treating a GPU port as a formulation change.

# Round 4 (coordinator, 2/10 20:00) — change the metric from seconds to error frequency

## Review of round 3: you decided the question and rightly refused the one I actually asked
Both terms 20× faster gives 7,25 s, not 1 s — my check (124+21)/20 gave exactly the same. For 1 s, 145× on both is required, or 248× and 42×. **1 s is thus not an ambitious target but an unattainable one**, and that is the finding. Setting `clinical_1s_requirement_gate = UNKNOWN` instead of reading a speed from a figure was right: the anchor compares pre-first with post-second in staged procedures, not a second needle pass, and the inter-image times are unpublished. `physical_target_velocity_upper_bound_mm_s = None` is the honest answer.

## Task: there is a metric that does NOT require a speed
You extracted it yourself: the fraction of cases with error over 3 mm drops from **20,3 % to 4,1 % upon rescanning**. It connects updating to the clinical quantity without going through any speed, and that is the route that is open.

1. **Compute the value of an update in error frequency, not seconds.** If a rescan takes 20,3 % to 4,1 %, what do two give, and what are the diminishing returns? State the assumption form you use to extrapolate and why, because a single measurement point does not determine a curve.
2. **Reverse the question:** at 7,25 s per update, how many updates fit in a procedure of published duration, and what error frequency does that correspond to? That makes 7,25 s a clinical number instead of a technical bound.
3. **Say whether 7,25 s suffices.** Against the target budget 1,08 mm and the measured error frequency: is the achievable update rate sufficient, insufficient or undetermined? All three are valid answers, and undetermined requires you to name what is missing.
4. **Keep 20,3 and 4,1 as held numbers.** They are an incidence in staged procedures, not a speed, and must not become anything else along the way.

## Control, falsifier, prohibitions
- **Control:** the seconds metric, i.e. 150 s versus 1 s. The gain should be that the metric changes to something clinically decidable.
- **Falsifier:** if error frequency cannot be linked to update rate without an unpublished quantity, say which and stop there. Two rounds have already shown that a speed cannot be obtained.
- **Prohibited:** extrapolating a curve from one point without declaring the form; turning an incidence into a speed; touching the other 166 ports.

# Round 5 (coordinator, 2/10 20:30) — FINAL: write the port, because the question is answered

## Review of round 4: the metric change gave the answer, and it reverses the brainstorm’s picture
At 7,25 s per update, a procedure of published mean duration fits **1192 updates** — I recalculated 144,1 × 60 / 7,25 and got 1192,5, i.e. your number. Computation time is thus nowhere near binding at clinical rates. It is the opposite of the brainstorm’s 1,0e8 gap in acquisition rate, and it is now computed and not asserted.

Equally important is that you left extrapolation undecided: `repeat_update_identification_gate = FAIL` with unbounded `[0, 1]` for the second update, the two assumption forms 5× apart, and both marked SENSITIVITY with `hypothetical_sensitivity_fits_use_held_outputs = True` printed. One point determines no curve, and you said so instead of fitting one.

**Four rounds suffice. This is the lane’s last.** Everything remaining is unreported clinical quantities, and they are not obtained with Sol time.

## Deliver two things and close
1. **PORT.json for the robot branch**, with what is now decided and what is not:
   - 202 ports, 0 connected to a computed quantity, 167 require new computation.
   - 152 ports are the robot’s own precision; 40 moved to the intelligence engine; 15 are bench requirements and not robot gaps.
   - 150 s = 124 s solve + 21 s resampling; maximum gain on solve alone 6,90×, on resampling alone 1,17×; both 20× gives 7,25 s; 1 s requires 145× on both.
   - At 7,25 s a procedure fits 1192 updates, so computation does not bind at clinical rates.
   - The binding term is the actual update rate, which is unreported.
   - Every number with provenance class and type of spread. No vendor number as a measurement.
2. **The acquisition list, short and ranked.** At the top: actual clinical update rate and number of updates per procedure, with the publication or registry source that would give them. Then the inter-image times that would make a target speed computable. Formulate each row so a swarm job can take it, with external_referent filled in. **No template repeated per port** — another lane today produced 24 identical follow-up proposals and only one of 31 carried external ground truth.

## Prohibited
Fitting a curve to one point. Reporting 1192 as a clinical recommendation — it is a computation ceiling, not a rate anyone should use. Opening any of the other 166 ports.
