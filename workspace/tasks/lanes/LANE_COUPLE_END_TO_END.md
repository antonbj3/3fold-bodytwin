# LANE_COUPLE_END_TO_END

Get ONE chain all the way from lower-level inputs to a held-out measurement, without a human moving anything. Results directory `results/LANE_COUPLE_END_TO_END/`.

## Why this lane

The goal is to outperform. Measured tonight, what is missing is not correctness and not computing power — it is that **almost nothing reaches anything else**:

- 0 of 202 robot ports carry a calculated quantity.
- 27 proved signs have 0 consumers, and 25 of the anchors have no decisions at all.
- External models are executable but have never been executed here.
- A typed `memory` field is not read by the claim layer.

The parts exist in thousands. The connections barely exist. A twin gets its value from the connections, not from the number of parts, and **no lane has worked on that axis**.

Packaging is not the obstacle: the controller copies `sources/<key>` for each key in `source_keys`, and measured over 112 new jobs, 47 % name two or more families and get them. Two families with artifacts on disk can therefore be coupled today.

## Do this

1. **Choose two families that both have artifacts and where one produces what the other consumes.** Candidates to test, not assume: Q115 (crypt/villus renewal) against Q168 (intestinal lumen loss, which declares `k_pre` as a synthetic assumption); Q012 against Q052 (metabolite concentration as a boundary condition); Q005 against Q012. Say why the chosen pair has a real producer-consumer relationship and not just two names in the same list.
2. **Build the chain so it RUNS.** One script, one command, no manual intermediate steps. Input through the typed ports, units converted explicitly and once. If a step requires someone to paste a number, the chain is not coupled.
3. **End in a held-out measurement.** Choose the anchor from `results/LANE_EXTERNAL_FACIT_HUNT/FACIT_INDEX_v2.json`. The anchor is held-out data and must never be input — that is the entire difference between a chain and calibration.
4. **Report the coupling's cost:** number of glue lines, number of unit conversions, number of assumptions that must be added that were not in either family. That number says what "frictionless" is worth for us today.
5. **Leave the chain reusable.** The next lane should be able to replace one end. That means a named port contract, not a script with two hardcoded families.

## Strongest control and falsifiers

- **Control:** the two families run separately, each against its own anchor. The gain should be that the coupled chain says something neither says alone — a quantity that only arises in the coupling.
- **Falsifier:** if the chain runs but says nothing beyond the two parts, the coupling is technical rather than scientific. Say it plainly; having an executable chain is still progress, but claim no more than it gives.
- **Forbidden:** reimplementing one family to make them fit (done twice in the project and counted as interop); letting the anchor become input; calling a chain coupled if a step requires a human; hardcoding unit conversions without naming them.

## Delivery

One executable command, one number, one comparison against a held-out measurement, and the coupling's cost in lines and assumptions. `PORT.json` with the port contract so the next lane can replace an end. If it cannot be done: the exact list of deficiencies in the contract.

Everything PENDING_INDEPENDENT_REVIEW. No internal data.
