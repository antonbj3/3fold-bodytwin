# MECHANISM — SCALING DOCTRINE (how the structure grows 100×)

> **Status: forward-looking architecture, NOT current implementation.** The live system today is a single
> coordinator + weak subagents (see `MECHANISM_HARDENED_CONVENTIONS.md` for what actually runs). This doc is
> the target we reason toward when scope + dataflow grow ~100×. Do not mistake it for present state; the
> hardened-conventions doc + the hardening check remain the authority on what IS.
>
> Developed 2026-07-19 in a coordinator↔operator design conversation. Refine it; it is a thinking tool.

## 1. The problem
At 100× scope the "coordinator = designer" model breaks: one strong-model main loop cannot design at that
throughput. But design is where quality lives and it resists the weakest models. So: **which class of model
sits where, and where does the system move over time?** The answer is two orthogonal scaling axes that
compose, plus one principle that keeps the expensive axis bounded.

## 2. The core inversion — you don't scale the designer, you scale the FUNNEL
"Design needs a strong model" is true for *one* design done well — but the quality lives in a **judge**, not
the generator, and verifying is cheaper than generating. So at scale you replace *one strong design* with
**many weak proposals + a strong FILTER**, and the filter is mostly **machine** (code), not model. The
verdict comes from the gate, never from the generator's say-so (CS doctrine: "make QC a lookup, not a
judgement"). Decompose design by *judgement difficulty*:

| Substep | Model class | Why |
|---|---|---|
| Enumerate holes (which occluded states are missing) | graph + weak | the graph surfaces its own gaps (honest_gaps, missing edges, frontier) |
| Propose legs + anchor (candidate decorrelation) | weak, massively parallel + retrieval | generation is cheap; many scouts per hole |
| **Admission gate** (is the anchor truly decorrelated / held-out, not a tautology) | **MACHINE** (`coupling_admission_precheck`, anchor-independence) | this is where quality is *forced*; as code it scales for free |
| Measure the decorrelation (run the legs, check they meet on the hidden state vs the anchor) | weak + compute | the cert is objective once designed — verdict is the number, not judgement |
| Adjudicate the residual (gate-ambiguous cases) | **strong, rarely** | only what the machine + consensus can't decide (= CS's ESCALATE-LLM) |

**Consequence:** the strong model's leverage migrates *up a level* — from designing *cells* to building the
reusable **objective gates** (the judges). Build `anchor_independence` once; it then judges millions of weak
proposals for free. (This is "systematize a repeated fix into infra" at architecture scale.) Per-cell strong-
model cost → 0.

## 3. Axis A (horizontal) — a CLUSTER of coordinator/designer main-sessions
You cannot thin the coordinator to ~1 if the *novelty* is wide, so parallelise the strong layer too:
**N coordinator main-sessions, each owning a REGION of the anchor-graph** (a cluster/subtree), each doing the
high-judgement design + prioritisation for its region. This is already what the CS fleet is (10+ lanes + a
hub). mechanism is currently ONE lane; at 100× it becomes a *fleet of coordinator-lanes*.

Two things keep a coordinator-cluster coherent (and we already have the pieces):
- **Coordinators coordinate THROUGH the graph, not via chatter.** Each *claims* a region (lock/claim system),
  designs + folds there; the graph's `depends_on` / `couples_to` edges surface cross-region couplings. Fibre
  doctrine: glass in between, LLM only at the boundaries. **Integrator cells are the boundary objects** where
  two coordinators' regions meet — that is why they matter (they are the seams, pre-built).
- **The shared graph is the single substrate** (one `MECHANISM_ANCHOR_GRAPH.json`, canonical schema, machine
  gates). No coordinator holds private truth; everything folds back into the shared graph.

## 4. Axis B (vertical) — weak-flood + machine gates under each coordinator
Within each coordinator's region: the massively-parallel weak-subagent flood (the four roles — acquire /
propose / measure / verify) + the machine gates carrying the verdict. This is what runs today, just replicated
per region. Volume scales here, cheaply.

## 5. The principle that bounds the expensive axis — the MOVING NOVELTY-FRONTIER
Coordinators (strong, expensive) are a **moving frontier, not a fixed assignment**. A strong model is needed
only where judgement is not yet mechanised into a gate:
- **Mature region** (failure-modes mechanised, gates cover it) → no coordinator, just weak-flood + machine.
- **Frontier region** (new domain, new cert structure, un-mechanised judgement) → one coordinator, doing what
  the gates can't yet.
- As a frontier region **matures** (its failure-modes get mechanised into gates) it **transitions** from
  coordinator-owned to weak-flood. **Coordinators advance forward.**

Therefore: **volume scales with weak+machine (cheap); the number of coordinators tracks the *width of the
un-mechanised novelty front*, not the total scale.** The architect's job is to shrink that front from behind
by building gates, so N does not grow linearly with scale. Design-selection ("which holes are worth
cornering") *partially* resists mechanisation — but only at the frontier, and the frontier is finite-width.

## 6. Where the strong-model leverage actually sits (two places)
- **(a) Being a frontier coordinator** — Axis A, design + judgement at the edge.
- **(b) Building the gates that let the front move** — Axis B's machine layer. The optimal system maximises
  (b) so that (a) stays a thin, moving front instead of an army.

## 6b. Federating stages 2 & 3 down to the weakest model — the substep → model-class map
Both DESIGN and MEASURE/CERTIFY federate across many small models run **sequentially through gates**. The rule:
**a substep runs on the weakest model whose output is caught by an existing machine gate** (Haiku, or weaker
still where the gate is strong); it needs a strong model only where the *coordinator is still the gate*.

**MEASURE/CERTIFY (stage 3) — mostly weak/compute TODAY:**
- pick + parametrise a forward-model → weak (Haiku) **if** a library of reusable forward-models exists (CS's
  "reusable judge" pattern); building a NOVEL forward-model → strong (the only residual here).
- run it → pure compute. compute agreement + verdict → pure code (`fold_gate` §5). No model.
- Binding constraint: forward-model **library coverage**.

**DESIGN (stage 2) — ~5 propose+gate substeps:**
| substep | model | gate that makes it safe |
|---|---|---|
| 2a enumerate candidate occluded states | weak + graph | graph: real occluded state, not already-measured / trivial |
| 2b propose legs | weak + retrieval, parallel | leg independently measures the state |
| 2c propose held-out anchor | weak | **anchor-independence** (different source AND family from every leg) |
| 2d check leg↔leg decorrelation | MACHINE | common-mode probe / `coupling_admission_precheck` |
| 2e falsifier | medium | triviality gate (tight inequality, not vacuous) — `fold_gate_v2` already does this |

**Residual that still needs a strong model** (the coordinator's irreducible job): (i) VALUE-prioritising which
enumerated candidate states are worth cornering, and (ii) genuinely NOVEL cert framings. Small — it is
SELECTION over an enumerated list, not generation from scratch.

**What runs on Haiku (or weaker) NOW vs still-strong:**
- NOW (already machine-gated): acquire-extraction (MT-JSON sentinel + `extract_cell_json` + refuse), schema,
  fold verdict, triviality, isolation. This session ran on a strong coordinator but these parts did not need it.
- STILL STRONG (the coordinator IS the gate): design QUALITY — is the anchor truly decorrelated, is the hidden
  state worth cornering, is the falsifier sharp. The RICHNESS of this session's cells (honest_gaps, the
  decorrelation reasoning) was enforced by the coordinator, not a machine, so it can't drop to Haiku until
  those judgements are mechanised.

**Gates are the quality↔quantity conversion.** A gate turns a one-time strong-model investment (build the
gate) into unbounded weak-model throughput under it. So the optimisation is not "how many strong models" but
"**how many judgements have we cast into gates**." Every judgement mechanised moves a substep strong → Haiku →
(weaker still). And the next gates are found by DOING (hit a failure-mode → mechanise it), not by more upfront
theory — past this point, upstream decomposition has diminishing returns vs. starting execution.

## 7. Trajectory + where we are now
- **Now:** one coordinator-designer (me) + weak subagents + `fold_gate_v2` + a first anchor-independence
  check. The 158-node graph is a *design frontier*, not yet certified (decorrelation measured for 0 cells).
- **Next concrete steps toward the doctrine:**
  1. ★**#1 highest-leverage:** harden `anchor_independence` + common-mode / `coupling_admission_precheck` into
     the objective **design-admission gate**, so a weak proposal's anchor is machine-judged for tautology/leak
     (not coordinator-judged). This is substeps 2c+2d — hardening it drops the whole hypothesis wave to Haiku.
     (My bt_memory flags anchor-independence as the wave-1 design weakness: subagents attach tautological
     anchors — so this gate is the known-weak one and the biggest single federation unlock.)
  2. Let weak models flood proposals; machine-gate them; reserve the coordinator for adjudication +
     gate-building.
  3. Stand up the **execution** stage (build forward-models, measure the decorrelation) so cells transition
     `OPEN → measured` — this is the currently-missing certification step; it is objective/compute, scales.
  4. Only then split into regions + a coordinator-cluster (Axis A), coordinating through the graph, with the
     integrator cells as the region boundaries.
- **Where we move:** the binding constraint shifts from *design capacity* to **gate coverage** — how much
  judgement is mechanised. Certification (decorrelation measurement) is already objective, so it scales for
  free once designs + data exist; the real frontier is mechanising the *design-quality* judgement
  (tautology detection, leg-independence) into gates.

## 8. Open questions (unresolved — the honest residual; a brainstorm-resume list)
- ★ **THE interesting one — "what is genuinely hard?"** The irreducible residual is (i) VALUE-prioritising
  which enumerated holes are worth cornering and (ii) genuinely NOVEL cert framings. Is even (i) eventually
  mechanisable — a **value/priority model** that scores candidate holes (reach × novelty × decorrelation-
  headroom) so selection also drops off the strong model — or is it the permanent strong-model floor? This
  sets the floor on coordinator count. (ii) novel framing feels harder to mechanise than (i) selection.
- Could the *falsifier* substep (2e) go below "medium" — is sharpness fully capturable by the triviality gate,
  or does it need judgement the gate can't express?
- Region partition + rebalancing across the coordinator-cluster as the frontier moves (load-balancing a
  moving front) — the fleet-hub problem, not yet solved for mechanism.
- Cross-region integrator (boundary) cells: who owns a coupling that spans two coordinators' regions?
- Region partition + rebalancing across the coordinator-cluster as the frontier moves (load-balancing a
  moving front) — the fleet-hub problem, not yet solved for mechanism.
- Cross-region integrator (boundary) cells: who owns a coupling that spans two coordinators' regions?
