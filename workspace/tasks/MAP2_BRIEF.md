# MAP2 — overlap map: the day's results against existing BodyTwin (step 1 tonight, Anton 23:20 "Start with 1")

Purpose: for every result from the day 22–23/9 decide, with evidence, which existing graph node and existing code it belongs to, and whether it agrees with, contradicts, extends or is new. No new experiments. This is mapping; correctness comes before coverage.

## Sources (read-only)
- Source graph: `~/projects/bodytwin/data/MECHANISM_ANCHOR_GRAPH.json` (3 950 nodes), also `data/ANCHOR_GRAPH.json`, `data/MECHANISM_DATA_GRAPH.json`, `data/FOLD_LEDGER*.jsonl`, `data/ARROW_LEDGER.jsonl`.
- Code: `~/projects/bodytwin/scripts/{msk,tissuetwin,physics_exp,cad,kernel,video_index,mechanism,local,hand_leg}` and `docs/`. Machine inventory with docstrings: `~/projects/3fold-workspaces/bodytwin/results/MAP/private_inventory.tsv`; domain nodes: `results/MAP/graph_domain_nodes.tsv`.
- Project memory: `~/projects/bodytwin/bt_memory/` (grep, read relevant files).
- The day's results: `~/projects/3fold-workspaces/bodytwin/notes/RESULTS_INDEX.md` (A/B rows), `notes/DAY1_REVIEW.md`, `results/<id>/`.

## Method (mandatory)
1. For every result-id in your area: read the row in RESULTS_INDEX and the result folder's main file (summary/JSON), so you know what was actually measured.
2. Search with several word choices — Swedish and English, synonyms, dataset and author names, method names (e.g. "scaling", "morph", "statistical shape", "landmark", "moment arm", "attachment", "insertion", "origin", "GRF", "ground reaction", "contact force", "OpenCap", "Grand Challenge", "eTibia", "Fregly"). Search all of the nodes' text fields (id, claim, notes, evidence, scripts), not only id. Search the code (file names + docstrings + grep in content).
3. For every hit: read the node/script enough to decide the relation. A keyword hit is not enough.
4. "New" may only be stated if searches are listed and empty. State exact search terms in the table.
5. If the graph contains a number for the same quantity: compare it with the day's number (same data? same conditions?).

## Output (only under `~/projects/3fold-workspaces/bodytwin/results/MAP2/<ditt-id>/`)
- `overlap.tsv`: result_id | what the day measured (1 row, with number) | node IDs (status) | script (path) | relation {same, agrees, contradicts, extends, new, unclear} | evidence (node's claim excerpt / script's docstring) | search terms used | proposal {record as trial on node X, fold proposal, new node justified, remove as duplicate}.
- `SUMMARY.md` (short, ≤ 1 page): per area what already existed, what the day added, contradictions, and the 3–5 most promising "build on" points (existing node/script + the day's result → concrete next trial). No adjectives, numbers with sources.
- Change nothing in `~/projects/bodytwin`, in `notes/` or in generated graphs.

## Resources
Load is high (load ~32/20). Only reading and grep/python-json; no computation jobs, no downloads, thread ceiling 2, `nice -n 19`. Never `pgrep -f`, never `nvidia-smi -q`. No own subagents. No emails, no publication, no lane_runner.
