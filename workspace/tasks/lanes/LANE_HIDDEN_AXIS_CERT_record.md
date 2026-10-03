# LANE_HIDDEN_AXIS_CERT

Build the chain nobody runs: hidden axis → contradiction federation → cert with margin. Results folder `results/LANE_HIDDEN_AXIS_CERT/`.

## Why this lane

Anton 1/10: hidden variables in the data anchors are the hardest thing about biology, and they must be **vibrated out and certified**, not described in prose. The machinery exists in the graph engine, but the graph lane (anton-4d) has measured that **nobody runs the chain from beginning to end**: `hidden_variable.from_mechanism` is consumed only by `disagreement_field`, `claim_types` and its own test. The sequence must be built once, reusable, and run on a real BodyTwin anchor.

This is not infrastructure for its own sake. A biological anchor is a projection of a state we do not see, and its spread is the hidden variables' fingerprint. Matching the mean may mean being wrong about every individual. Without this chain, we cannot distinguish "the model is correct" from "the model is fitted at population level".

## The tools, read them before you build

Path: `~/projects/graph_workspace/graph-typed-throws-20261001/src/graph_engine/` (109 moduler), **pinnad commit `352c6d3`** on the branch `research/typed-throws-20261001` (suite green: 3401 passerade, 68 deserted, 0 fallerade)Pinnad instead of just a worktree path, because a moved file now twice got agents to report it as missing.

CORRECTED PATH 2/10 16:35: the brief previously pointed to `graph-coordinator-oed-batching-20260923`, which has 97 modules and does NOT contain `hidden_axis_chain.py`. Two agents reported it missing, and they were right about that path — the error was mine. The chain with `instrument_can_decide` (threshold-free, returns sigma_required through bisection) and `coverage_audit` is in typed-throws. Warning: there is also a `coverage_audit` in `scene_eyes_coverage_audit.py`, which is a DIFFERENT function for scene/eyes — do not confuse them.

- `hidden_variable.from_mechanism(model, x_range, params, n_x=41, n_p=21, fixed=None)` — the mechanism side, before any data. Sweeps every parameter over its interval and measures over what fraction of (x, θ) the **sign** of dy/dx differs from the midpoint's, plus the threshold where reversal first appears. Score 0 means the parameter never reverses the sign, so it can be omitted from the validity box without cost to the sign question.
- `hidden_variable.candidates(reports, min_side=2, min_known=0.7, n_perm=200)` — the data side. **Two thresholds that silence:** empty list when n < 2·min_side, thus **four** reports are required, and the attribute must be present on at least 70 % of them. An empty list is therefore NOT absence of a hidden variable. Reports: `[{"margin": float, "sigma": float, "attributes": {...}}]`, or `"sign": ±1` instead of margin.
- `claim_federation` — **the entry point** for comparing contradictory studies. It traces every source to its root sources and calculates N_eff, so four papers citing the same original support the claim as ONE, and it distinguishes REGIME-BOUNDARY (disjoint validity boxes, thus different regimes) from CONTRADICTION (overlapping boxes, thus an undeclared axis). **Do not touch the file** — the graph lane connected a typed gate there today.
- `decision_cert.certify(decision, state)` with `flip_attribution`, and `margin_net` for margin and alerts.
- `disagreement_field` is downstream of `hidden_variable`, not an alternative.

## Do this

1. **Build the sequence as a reusable module in your own directory** (import the graph engine, do not change it): mechanism sweep → federation of reports → hidden-axis candidates → cert with margin and reversing quantity. Write it so that a swarm job can call it with a model and a list of reports.
2. **Run it on a real anchor.** Choose one from `results/LANE_EXTERNAL_FACIT_HUNT/FACIT_INDEX_v2.json` (43 records, errata-checked by our own lane, thus PENDING_INDEPENDENT_REVIEW) where spread is reported and where there are at least four independent reports to acquire. Nasal nicotine, systemic fraction 0,53 ± 0,16, is a candidate with CYP2A6 genotype as the presumed dominant axis — but check attribute coverage before interpreting an empty answer.
3. **Report three things per anchor:** which parameters can reverse the sign within their own intervals and at what threshold; whether the model can produce the anchor's **spread** by vibrating the identified axes within published intervals (only the mean ⇒ label POPULATION_FITTED, resolution not reached); and the decision as a cert with margin plus the quantity that would reverse it.
4. **Name the mechanism that would give the same mean for other reasons**, and which part of the spread distinguishes them. Two mechanisms with the same mean are not separated by the anchor.

## Strongest control and falsifiers

- **Control:** an ordinary mean test against the anchor, thus what we have done so far. The gain must be of another kind: a statement about which axis carries the spread and a cert that says when the decision reverses. That is a guarantee a mean test cannot give at any budget.
- **Falsifier:** if the mechanism sweep gives zero sign reversals for every parameter, the anchor's box is sufficient and the chain adds nothing for that anchor — report it as an honest negative and test the next anchor in the same round. If the data side gives an empty list: show that it is n or attribute coverage, not absence.
- **Forbidden:** writing prose about hidden variables instead of running the sweep. Interpreting an empty list as absence. Changing the graph engine.

## Leverans

The module, one run per tested anchor, and a `PORT.json` that says which of our 43 anchors are POPULATION_FITTED and which carry an identified axis. It is a map of where resolution is actually reached.

Inga interna data. Allt PENDING_INDEPENDENT_REVIEW.


## PIN UPDATED (the coordinator 2/10 21:05)
Pinned commit moved from `cf9c1e9` to `352c6d3`, three commits later, seven files changed. The reason is that the three contain precisely the corrections we have pursued tonight: *Abstain when the downstream cert says nothing, instead of reading silence as a pass*, *Stop a validity filter from passing a point it never checked*, and *Say whether the declaration could separate two discordant reports, and count unconsumed verdicts*. The older pin would have let silence pass as approved.
The engine is now also on both cloud hosts under `/opt/agents/graph_engine`, with the cited path symlinked there, so a job no longer needs to reimplement it. Verified: import through the cited path works on both.
