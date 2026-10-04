# LANE_REGISTRATION_SCOPE — the binding term may have the wrong scope in four lanes

Result directory `results/LANE_REGISTRATION_SCOPE/`.

## Why the lane exists, and why it goes first
The bleeding lane found something in round 2 that outweighs its own question:
**`published_timing_excludes_spatial_registration = True`.** The entire robot branch rests on a
registration taking 150 s, of which 124 s solving and 21 s resampling, and on 1 192 updates fitting
in a procedure of published mean duration 144,1 min at 7,25 s per update. Four lanes have built on
that being the binding term. If the published timing excludes spatial registration,
the number has the wrong scope, and then everything resting on it is ill-posed — including the conclusion that computation
does not bind.

## Do this
1. Find and cite the primary source for the 150 s, with locator. Say exactly what the timing
   **includes** and what it **excludes**. The bleeding lane’s new locator is a place to start:
   DOI 10.1038/s41598-026-71616-w.
2. If spatial registration is excluded: find its cost separately, and recompute both 150 s and
   1 192. Report the old and new number side by side.
3. **Trace the consequence upwards.** Which conclusions in `results/LANE_ROBOT_PORT_CONNECT/`,
   `results/LANE_BLEEDING_VISIBLE/` and `notes/ACQUISITION_LIST.md` item 4 change if the number changes?
   List them with the direction they move. The conclusion "1 192 updates fit, so computation
   does not bind" is what is at stake.
4. If the timing actually INCLUDES spatial registration, say so and close the question. A
   confirmed number is a fully adequate outcome and cheaper than all four lanes continuing in uncertainty.

## Control and falsifier
- **Control:** the 150 s as they stand today, with their current scope.
- **Falsifier:** if recomputation changes the number of updates per procedure by less than 10 %,
  the scope question is irrelevant to the conclusion and that must be said plainly.
- **Forbidden:** estimating registration time from your own run; this is a source question.

## Delivery
`PORT.json`: the 150 s with scope and exclusions cited, recomputed number if it changes, and the list
of conclusions that move. Everything PENDING_INDEPENDENT_REVIEW.

# Round 2 (coordinator, 3/10 00:15) — reformulate the number in every lane that rests on it
You established it: 145 s is phantom core, 150 s cadaver approximation, and **the one-time overhead is
strictly greater than 866,75 s, hence 5,98× the core**, plus at least 0,81 s per repetition. And spatial
registration time is absent from the literature, so `full_cycle_capacity` is unknown. The laser lane came
independently to the same conclusion, so three lanes now point in the same direction.
- **Changed operation:** write the corrected statement in a form other lanes can consume. "1 192
  updates" is CORE CAPACITY. State three numbers instead: core capacity, a lower bound on
  cycle capacity that includes the 866,75 s, and an explicit UNKNOWN for the rest.
- **Trace the consequence:** `notes/ACQUISITION_LIST.md` item 4 says that computation does not bind at
  clinical pace, and that rests on 1 192. Recompute with the one-time overhead included and say whether
  the conclusion stands. If it reverses, it is one of tonight’s most consequential findings, because four lanes built on it.
- **Acquisition item:** spatial registration time is now a named missing measurement. Formulate it
  as a reporting request and not an experiment — someone has run it, they have simply not printed
  it. Say exactly which number to request, with unit.
- **Falsifier:** if cycle capacity with the one-time overhead included still exceeds the number of
  updates a procedure needs, computation still does not bind and the scope question was
  irrelevant to the conclusion. Say so plainly — it is a fully adequate outcome.

# New round (coordinator, 3/10 00:25)
Recomputation gave −9,98 %, hence just below my own 10 % threshold, and the conclusion that computation does not bind STANDS. Do not finish on that: the overhead 866,75 s has measured866_75_overhead = False, hence hypothetical and unmeasured. Changed operation: acquire a MEASURED number for the one-time overhead and spatial registration, formulated as a reporting request and not an experiment. And deliver the consumable triple: core capacity 1192, lower bound 1072, and UNKNOWN for the rest, so other lanes can cite the right number. Falsifier: if no publication reports spatial registration time separately, it is an acquisition item and the lane closes with that.
