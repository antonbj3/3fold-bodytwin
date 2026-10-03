# LANE_REGISTRATION_SCOPE — the binding term may have the wrong scope in four lanes

Resultatmapp `results/LANE_REGISTRATION_SCOPE/`.

## Why the lane exists, and why it goes first
The bleeding lane found something in round 2 that weighs more than its own question:
**`published_timing_excludes_spatial_registration = True`.** The whole robot branch rests on a
registration taking 150 s, of which 124 s solve and 21 s resampling, and on 1 192 updates fitting
in a procedure of published mean duration 144,1 min at 7,25 s per update. Four lanes have built on
this being the binding term. If the published timing excludes spatial registration,
the number has the wrong scope, and then everything resting on it is wrongly posed — including the conclusion that computation
does not bind.

## Do this
1. Find and quote the primary source for the 150 s, with locator. State exactly what the timing
   **includes** and what it **excludes**. The bleeding lane’s new locator is a place to start:
   DOI 10.1038/s41598-026-71616-w.
2. If spatial registration is excluded: find what it costs separately, and recalculate both 150 s and
   1 192. Report the old and new number side by side.
3. **Trace the consequence upward.** Which conclusions in `results/LANE_ROBOT_PORT_CONNECT/`,
   `results/LANE_BLEEDING_VISIBLE/` and `notes/ACQUISITION_LIST.md` point 4 change if the number changes?
   List them with the direction they move. The conclusion "1 192 updates fit, so computation
   does not bind" is what is at stake.
4. If the timing actually INCLUDES spatial registration, say so and close the question. A
   confirmed number is a fully valid outcome and cheaper than all four lanes continuing in uncertainty.

## Control and falsifier
- **Control:** the 150 s as they stand today, with their current scope.
- **Falsifier:** if recalculation changes the number of updates per procedure by less than 10 %,
  the scope question is irrelevant to the conclusion and this should be stated plainly.
- **Prohibited:** estimating registration time from your own run; this is a source question.

## Delivery
`PORT.json`: the 150 s with scope and exclusions quoted, recalculated number if it changes, and the list
of conclusions that move. Everything PENDING_INDEPENDENT_REVIEW.
