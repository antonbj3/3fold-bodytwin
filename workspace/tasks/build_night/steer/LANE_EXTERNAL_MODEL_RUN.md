# Control LANE_EXTERNAL_MODEL_RUN — round 1 (the coordinator, 2/10 21:00)

New lane, and it is a shift in what Sol time is used for. Anton: "I kind of don't really know what you do if you work purposefully on goals and a greatly improved bodytwin". He's right — the evening has been spent verifying numbers and closing lanes, and five lanes were closed, four of which were because the question was wrong. It's useful but it's cleaning. This lane builds.

- **Hinders:** interoperability is untested. The only true solver hookup was rejected by our own gate at 0 % valid frames, and FEBio is runnable but never executed here.
- **Changed operation:** take the already FAILED connection and make it a run that passes, instead of demoing a new model.
- **Strongest control:** our own reimplementation of the same mechanism. The benefit is that a model we didn't write drives through our gates, not that the speech gets better.
- **Falsifier:** if no external model can run end-to-end without manual steps, "drop-in" is not a feature we have — then supply the exact list of defects in the port contract.

## The order is important: the gate first, the model later
0 % valid frames is a gate outcome. **Determine if the gate was correct before fixing the model.** A gate that rejects a correct model is a bug in the gate; one who rejects a wrong did his job. Bypassing it without deciding which is the only way to make this round worthless.

## Internal data: Sol can read it, the swarm can't
`results/CX-SOLVER`, `CX-SOLVER2` and `CX-SOLVER-ROBUST` contain internal identifiers. You drive locally and get to read them. Do not quote them in any file that may leave the machine, and do not write internal names in RESULTS.md. This is precisely why this task is a Sol task and not a swarm job.

## Three error classes that cost us today
1. **Reimplementation is not plugging in.** It has already been done twice in the project and counted as interop. If you write the model yourself, you have answered another question.
2. **Missing file, silently reconstructed.** Four jobs today silently reimplemented a missing module or reconstructed missing constants, and one of them gave an incorrect header that I amplified further. Something is missing: `missing_prerequisite` in the FIRST paragraph of RESULTS.md.
3. **Prediction read as measurement.** A field with *conditional* in the name was treated as measured, by me. Label each number MEASURED or HELD.

Allt PENDING_INDEPENDENT_REVIEW.

# Round 2 (the coordinator, 2/10 22:15) — determines whether 41 % is the comparison or the model
A model ran end-to-end against a measurement: 1722,10 N against 1222,26 N, peak error 0,4089, RMSE 531,75. 39 of 39 frames converged. It's a real drive.
- **Obstacle:** your own obstacle says that the observable doesn't match the held index and that an own mechanism check is missing. Then 41 % can be a comparison miss.
- **Changed operation:** determines the comparison FIRST. Same size, same frame, same body part, same normalization? Do not touch the model until it is settled.
- **Control:** our own reimplementation of the same mechanism. The gain is that a model we didn't write drove through the gates, not that the number gets better.
- **Falsifier:** if the observable matches and the error persists, it is the model, and then 41 % is the result — report it as such.
Prohibited: scaling the curve to fit; to change the anchor when the current one does not fit.

# Round 3 (the coordinator, 2/10 22:40) — split the 40,9 %, three explanations are already excluded
You have the project's first external model run against the hold measurement, and the error is **40,89 %** in peak power
(1722,0983 N vs 1222,2613 N, RMSE 531,75 N). That's the result — not that the drive went.

**You've already ruled out three things yourself, and that's why the number will do:**
- geometry: CAD proxy moves the fault by 6e-5
- sensor direction: a 5° cone takes 40,9 % to 39,5 %
- timestamps: the sector holds at most 0,881 N, so 0,18 % of the missing 499,84 N
Do not open any of them again.

**Changed operation:** constitutive law and boundary conditions remain. Varies EN parameter at a time over
its **published** range and report how much of 499,84 N it closes. Out comes one
table parameter → closed share. If no single parameter within published range closes
the error, it is the **form of the model** that is wrong and not its parameters, and that is a much stronger one
results than a calibration. In that case, say it outright.

**Prohibited:** adjusting a parameter outside its published range to close the bug; to
report a calibrated fit as validation; to count the drive as the find.
**Falsifier:** if each parameter individually closes under 10 % and no combination within
published range exceeds 50 %, the model shape is collapsed for our case.
**Control:** the own spread of the measurements held. Report 40,9 % against it, not against zero.

# New round (the coordinator, 3/10 02:20)
The condition in my falsifier is met by your numbers: best single parameter closes 18,6 percent of
gap and best combination 29,2, both under 50. You still set the model shape gate to UNKNOWN i
instead of asserting the fall, and it was right — a sweep of tried combinations is not one
exhaustive search.

Changed operation: make the search exhaustive enough for the gate to be set. You have 20 fasting
native cases and 25 attempted solutions of which 16 converged. The question is not more combinations but
if there IS any point in the published parameter space that closes over half the gap. It is one
optimization with a bound, not a sweep: maximize closed proportion over the rectangle and report
maximum. If the maximum is below 50 percent, the shape is folded within published ranges, and then can
the gate is set.

And keep the proportion: the baseline error is 1,88 hall standard deviations, so below two sigma.
The model is wrong at its peak value but not absurd against the spread of the measurement, and that should be in the gate
so no one reads 40,9 percent as the model is unusable.

The nine solutions that did not converge must be reported separately: are they outside the physical area or
is it a numerical question? A non-converged point can hide the very one that closes the gap.

Falsifier: if the maximum over the rectangle exceeds 50 percent, the shape is not collapsed but
underparameterized, and then it is calibration that is missing — that is completely different and cheaper
problem.
