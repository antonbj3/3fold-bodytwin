# Claims-graph excerpt

A small, filtered excerpt of a physiology claims graph, included so the framework and the graph
tools can be exercised on real (not synthetic) claim nodes.

- `ANCHOR_GRAPH.json` — 26 claim nodes, one per physiological claim.
- `FOLD_LEDGER.jsonl` — 17 evidence records, one per node that has at least one evidence artifact.

Each node is a claim plus the shape of its certification: a status (`OPEN` = designed but not
measured, `ASSUMED` = measured with a partial result, `REFUTED` = the claim failed its own test,
`DEFERRED`), a regime note stating where the claim is meant to hold, and the number of evidence
artifacts behind it. The vocabulary (`leg`, `anchor`, `n_eff`) is defined in `_meta.vocabulary`.

Validate it with the graph engine's tools:

    ANCHOR_GRAPH=ANCHOR_GRAPH.json FOLD_LEDGER=FOLD_LEDGER.jsonl \
      python anchor_graph_tools.py validate
    ANCHOR_GRAPH=ANCHOR_GRAPH.json FOLD_LEDGER=FOLD_LEDGER.jsonl \
      python anchor_graph_tools.py next-actions

Expected: 0 errors, 1 warning (`goal_weights_missing` — the excerpt carries no goal node).

## What this excerpt is not

It is not the full graph and not a biological atlas. Nodes were kept only if their cluster is a
physiology subsystem and their full text contains nothing about any individual: every node
mentioning a person, a personal measurement, an individual body model, a device or implant, or an
internal file path was dropped in full rather than edited. Evidence file paths were replaced by
counts, and internal annotations, dates and process labels were removed from the surviving text.
Dependency edges are empty here because every dependency target of a kept node was itself dropped
by those rules.
