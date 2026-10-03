# Map: which graphs are there, where, and how they are connected

Inventoryed from disk 2026-09-17 (script + raw data: `private_runs/graph_inventory.json`). Private file, not included in any repo.

## One machine, two file types

| file type | what it is | tool |
|---|---|---|
| `ANCHOR_GRAPH.json` | **The graph.** Nodes = targets and what the targets are on. Fields: `claim`, `type`, `status`, `depends_on`, `evidence`, `risk`, `cost`. This is README's "unlock graph". | `anchor_graph_tools.py` (`validate`, `next-actions`, `priority`) |
| `FOLD_LEDGER*.jsonl` | **The measurement log.** One row = a completed work with results. Points to the graph via `consumers`; the graph's nodes point back via `evidence`. `supersedes`/`refines` = corrections. | `fold_ledger_tools.py`, `fold_gate*.py` |

## ★ Antons egen bild av grafen (hans prompts jun–jul, minnesfilen `project_graph_as_anton_conceives_it`)

"Variables that have relations to each other" + unlock graph + densification. NOT accounting. The README
was written by an agent from the code and mostly describes the accounting layer. The layers, from top:

| lager | vad | fil i CS-I `data/` |
|---|---|---|
| Unlock graph | target → load-bearing subnodes | `ANCHOR_GRAPH.json` (97 nodes, 14 targets) |
| **Variables in relation (the core)** | quantities linked by quantified relations; edge status TIGHT / OPEN / UNKNOWN; evidence as list {file, key, value} | `CONSTRAINT_NET_VEHICLE.json` (20 variables, 34 edges: 32 TIGHT), `CONSTRAINT_NET_PROJECTOR.json` , `CONSTRAINT_NET_FLEET_FACTORY.json` ; stress map `CONSTRAINT_STRESS_MAP_V1.json` |
| Fuel | each measurement tightens an edge | `FOLD_LEDGER*.jsonl` |
| Densification (casts / conf / boundary) | done by agents reading full text; conf is determined with a discriminator cell | document (`CONFLICT_MAP_A/B.md`, `COMBINATORIAL_LIST.md`), no data structure |

Volvo: node `VEH-REDBLOCK-BASE` (root under `VEH-CONSTRAINT-NET` → `VEHICLE-GOAL` ), 46 folds, base in `~/research/volvo_redblock/` . No
migration needed: everything above is correct. What is missing is a READ VIEW that shows the layers together (reads, never writes).

## cad-to-simulation (ett git-repo, 16 worktrees = en gren per lane)

| vad | var | storlek |
|---|---|---|
| The graph | **only in CS-I** (`feat/breakthrough-hyperreal`, the hub) `data/ANCHOR_GRAPH.json` | 97 nodes (14 targets), 81 dependencies, 660 evidence pointers. Created 07-16 with 11 nodes, 364 commits |
| Measurement logs, one per lane | CS-I 976 · CS-L 704 · CS-F 626 · CS-H 319 · CS-G 124 · CS-B 8 · CS-C 8 · CS-J 7 | ≈ 2 770 folds total |
| Older generation (Lane D, 07-07…07-13) | `docs/demand/kernel_gen_field/`: `MASTER_UNLOCK_GRAPH.md` (3 202 lines of text), `GRAPH_STORE.json` 457 nodes (finding/law/hole/instrument), `FEDERATED_GRAPH.json` 430 nodes (descent seed→code→commit) | replaced by the anchor graph 07-16, remains |
| Not knowledge graphs | `seed_graph.json` (log seeds→code output), `substrate_graph.jsonl` (machine telemetry) | — |

The graph is small in number of nodes and heavy in evidence: 17 PROVEN, 78 of 97 nodes have evidence.

## BodyTwin / Mechanism (ett repo, gren `bodytwin`; `mechanism/` = samma repo vid 07-23, `bodytwin/` vid 08-14)

| vad | storlek |
|---|---|
| Inherited copy of the CS graph (mode 07-19) `data/ANCHOR_GRAPH.json` | 66 nodes — all 66 ids remain in CS today |
| Inherited measurement log `data/FOLD_LEDGER.jsonl` | 451 rows |
| **Own graph** `data/MECHANISM_ANCHOR_GRAPH.json`, same schema + `cluster`, `cert_design` | 1 666 nodes (07-23) → **3 950 nodes** (08-14) |

Mirror of CS, as intended. The difference in size: 3 873 of 3 950 nodes are OPEN (designed, not measured), 1 PROVEN; 1 876 comes from the cluster "autonomous-expansion". CS: few nodes, very saturated. BodyTwin: many designed nodes, a bit saturated.

## Staging / publikt (`3fold_staging/`)

The engine, not the graphs. `3fold-graph-engine` : all tools above + an example with 4 nodes and 4 folds. `3fold-bodytwin` : extract
45, 26 and 2 nodes with 35 resp. 17 folds. The README says the same: nodes and overlays were created privately and are not included.

## Not investigated
- "Product nodes": not looked up. Candidates: `data/DATA_GRAPH.json`, `capability_index.json`, `cell_registry_v1.json`, `fabriksspec_v*.json`, `koncept-V1` (29 folds point there).
- The lane loggers (L, F, H, G) have no anchor graph of their own: their `consumers` point to names that do not exist as nodes in the CS-I's graph. Put together, it becomes 3 700 names of which 24 % are related — free text names, not missing graph.
