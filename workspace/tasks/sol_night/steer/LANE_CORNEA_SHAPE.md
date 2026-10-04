# Styrning LANE_CORNEA_SHAPE — efter r21

## The first thing to fix: two numbers with the same name and different meanings
Your `same89` block carries `recomputed: False` and three toric numbers: `toric_field_magnitude_MAE_D`
0,284333, `toric_same_information_collapsed_MAE_D` 0,282503 and `toric_population_MAE_D` 0,587433,
with 50/69 to the measured posterior. The assembly chain in `tasks/assembly/toric_decision.py` gives on
the same 69 eyes and with the same 50/69 split 0,2808 and **0,6209** for the population arm. I checked
our own history: before tonight's correction of `from_vector` it was at 0,2817/0,6264, after it at
0,2808/0,6209 — so your numbers are not an old copy of ours, but another estimation step under a
name read as the same quantity. The population arm differs by 5,4 %.

**Operation:** either rerun the block in this workspace, or rename each number so that the predictor appears
in the name (`field_model_*` versus `thin_element_vector_*`) and print which quantity is compared. A number
read as a check of another number, but computed with another predictor, is the most dangerous
kind of entry in the network: three corpus searches have already failed on precisely that.

## Hindret i sak
`scope`: "HeLa native camera functional, not original phantom/cornea. No new measured index/clinical
endpoint." The lane has thus acquired a native reference/sample field but not on the specimen the question
concerns, and `R17_residual_D = 0,867062644` remains without being bound to any measured reference answer.

## Changed operation
Bind 0,867 D to a quantity or delete it. The direction the lane owns is still the one r17
measured: two corneas with identical mean thickness, index field, hydration AND surface allocation differ by
0,3056 D in 4 of 4 pairs, and the rotationally symmetric control case gives −1,36e-08 D. The next construction
is to make the non-symmetric structure an INPUT to the decision, not a residual entry: which measurable
surface quantity carries the 0,3056 D, and what happens to the residual when it is introduced?

## The strongest control
The same 69 eyes, the same metric, but the predictor without the new surface quantity. The gain counts only against it, and
only when both arms are computed in this workspace in the same run.

## The falsifier
If the new surface quantity does not reduce the residual below the acquired population arm on the same eyes,
the non-symmetric structure does not support the cylinder decision, and that is a result that must
be printed in plain text.


## Addition after r22 (this takes precedence over the naming question)
The diagnostic fails on 150 of 150 images. I tested the explanation that had been close all night —
that disagreement between reference choices is a branch-cut artefact, hence phase wrapping at π — and
**it does not hold**: of the 150 differences, 17 lie near zero, 2 near π and 131 in between, with
median 1,0409 rad and mean 1,2082 rad. The maximum 3,1348706 rad is merely the extreme value in a
spread-out distribution bounded by π, not a wrap. The spread is real and continuous.

The decisive number is therefore this: the lane's own `conditional_radius_rad` is 0,0026 while
the reference-choice ambiguity has median 1,0409 rad — **400 times larger**. No tightening of the
conditional interval means anything while the reference is unfixed, and reporting a radius that
is 400× smaller than the dominant ambiguity invites reading the precision as
resolution. Your own gate says where the path goes: `same_time_reference_ray_support:
NOT_ACQUIRED`. Acquire it, or stop reporting the conditional radius without the ambiguity beside it.
