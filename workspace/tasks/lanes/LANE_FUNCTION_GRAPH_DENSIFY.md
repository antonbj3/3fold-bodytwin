# LANE_FUNCTION_GRAPH_DENSIFY

Push connections: densify the function graph with edges, not nodes. Results directory `results/LANE_FUNCTION_GRAPH_DENSIFY/`.

## Why this lane

Anton 2/10: "push connections and densify the graph, perhaps we can call it the function graph". We have nodes; **we have almost no edges**, and this is measured:

- 0 of 202 robot ports carry a computed quantity.
- 27 proven signs have 0 consumers.
- **900 of 903 directory pairs are undecidable**, and only 3 decidable.
- 43 cells exist in the workspace, while the old project `~/projects/bodytwin` has **156 physiology domains with data**.

Three are large and entirely unused in the workspace:

| domain | data files | what it concerns |
|---|---|---|
| `data/glioma_fisher_v1` | 2266 | tumor growth, Fisher-KPP type |
| `data/bonemet_v1` | 1621 | bone metastasis |
| `data/hfpef_diastolic_v1` | 1500 | diastolic cardiac function |

You run locally and have read permission on them. **Read them without changing anything** — the old project is read-only and this workspace is an adapter.

## Uppgiften: leverera KANTER

A node without an edge is a blueprint. The lane is judged by the number of edges carrying a named shared quantity, not the number of domains it has touched.

1. **For each domain, determine which FUNCTION it computes.** Which quantity in, which out, in which units, and where the data state it. A directory name is not a function.
2. **Find the shared quantity against an existing cell.** Candidates to test, not assume: `hfpef_diastolic_v1` against Q021/Q022 (pulmonary flow, heart chamber) and against Q140 (renal plasma flow, 600 mL/min declared) via cardiac output; `bonemet_v1` against Q090 or Q127 (load-bearing bone, fiber/satellite-cell state) via mechanical load; `glioma_fisher_v1` against Q044 (radial spread in the head) and Q084 (transport) via diffusivity. **An edge requires the same quantity in the same unit on both sides** — otherwise it is an association and not an edge.
3. **Write the edge so it can run.** Producer, consumer, quantity, unit, and the conversion needed, once and named. An edge requiring a human to move a number is no edge.
4. **Say what each edge enables.** An edge is worth something if it makes a question askable that neither side can answer alone. Name the question.
5. **Report decidable pairs before and after.** The number is 3 of 903 today. It is the lane’s only primary metric.

## Control and falsifier

Classify according to COMMON.md: this is `information_link`. The control is what the cells can answer **without** the edge, not an equally informed counterpart — giving the control the same connection guarantees a tie.

- **Falsifier:** if none of the three domains shares a quantity in the same unit with any existing cell, they cannot be connected today and the delivery is which quantity is missing per domain. It is a fully valid outcome and better than an invented edge.
- **Prohibited:** creating an edge between two quantities sharing only a name (the unit decides, not the word); changing anything in `~/projects/bodytwin`; counting a domain as connected before the edge has run; building a new cell before the edges are mapped — Anton asked for connections, not more nodes.

## Leverans

`PORT.json` with the edges: producer, consumer, quantity, unit, conversion, and the question the edge makes askable. Plus the number of decidable pairs before and after. Narrow follow-ups in FOLLOWUPS.json — no template repeated per domain.

No internal data leave the machine. Everything PENDING_INDEPENDENT_REVIEW.
