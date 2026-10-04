# LANE_HISTORY_REPRESENTATION

Build what the claim layer cannot express: time, order and state. Result directory `results/LANE_HISTORY_REPRESENTATION/`. Form A.

## Why this lane

Anton's second seed was about compressing a hundred surgical experts into a dense expert. The crossing of fourteen agents showed that the search algorithm already exists in the graph engine — all twelve proposed operators point to existing functions. **The only thing genuinely missing was that the representation cannot express time, order or state:** "third passage through the same tissue", "after 40 minutes of retraction".

And that is not a theoretical objection. We measured it in our own material: in `results/LANE_SURGICAL_BINDINGS/PORTS_R4_FINAL.json`, two events, heat and opening, give work **670,58 J/m² in one order and 7559,25 J/m² in the reverse** — a factor 11,27 — with **identical final state**, max_opening 0,004 m in both. Same events, same endpoint, an order of magnitude apart.

A `Box` in the claim layer is `dict[str, tuple[float, float]]`, thus a parameter rectangle. Two states with identical parameter values and different history are the same point to it. The concept also already exists one layer down: `typed_throw_ops/signature.py` has a `memory` field with the values `memoryless` and `path_dependent`, and `margin_net` never reads it.

## Do this

**Do not change the graph engine** — the graph lane (anton-4d) owns it, and they added four verdicts for shared sign boxes today. Build in your own directory, against pinned commit `352c6d3`, branch `research/typed-throws-20261001`, under `/opt/agents/graph_engine` or the cited path.

1. **Formulate the criterion sharply.** The box form suffices exactly when a history compresses into a declared scalar: accumulated dose, number of previous passages. It does not suffice when two paths share every scalar total and still differ — which the heat-opening case does. **Decide which of the two our own case is**, and show it with the numbers above.
2. **Build the smallest representation that carries the difference.** A path identity, an ordered event list, or a state variable with an update rule — choose and justify. The requirement is that two paths with identical scalar sums get different identities, and that two paths that really are equivalent get the same one. The second requirement is the hard one: a representation that distinguishes everything is as useless as one that distinguishes nothing.
3. **Show that it is consumed.** A verdict nobody reads is worthless — we measured that four times today. Write at least one test that FAILS if a consumer ignores path identity, just as the graph lane's mutant test fails on 28 and 25 assertions respectively.
4. **Test against the surgical case.** "Third passage through the same tissue" must be expressible, and it must differ from "first passage" even when all declared scalars are equal.

## Strongest control and falsifier

- **Control:** a box with history lifted into an extra axis, thus "number of previous passages" as a parameter. It is cheap and often sufficient. The gain must be a case where it demonstrably does not suffice, and the heat-opening case is the candidate.
- **Falsifier:** if every path dependence we actually have compresses into a scalar, no new representation is needed and a declared axis suffices. Report that — it would be a strong negative that saves us a whole build.
- **Forbidden:** changing the graph engine; building a representation without a test that fails when it is ignored; treating 11,27× as a physical fact (it is our own synthetic thermo-cohesive model and the final states' equality may be a property of its construction — test that first).

## Deliverable

The representation, the test that fails, and the verdict on whether our own path dependence requires it. `PORT.json` to the graph lane if the result justifies a change in the engine — as a proposal, not as a change.

Everything PENDING_INDEPENDENT_REVIEW. No internal data.
