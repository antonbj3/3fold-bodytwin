# CX-QGRAPH — bind the 168 question models to BodyTwin's graph, and let the GRAPH drive the The swarm refill (instead of neighbouring Q numbers)

Anton: the pairing/priority shall come from the graph. Today tasks/refill_v2.py pairs DATX/CPL by neighbouring Q numbers and selects INT by word search. Crude.

## Tasks
1. **Mapping Q → graph nodes:**
   - For each BT-HX-Q001..Q168 (the question in results/BT-HX-Q*/inputs/QUESTION.md, the model and RESULTS.md), find 1–3 matching nodes in BodyTwin's canonical graph: GRAPH.json (3,950 nodes), MERGED_GRAPH.json (depends_on is the only unlock edge), CONSTRAINT_NETS.json, COUPLING_ANNOTATIONS.json. Use `./graph working rank --query '<words>'` and your own text matching.
   - Give a score and the reason for each match. Where no node fits, mark it "new node proposed" with its parents.
   - Output: `results/CX-QGRAPH/Q_GRAPH_MAP.json`.
2. **Graph-driven pairs:**
   - DATX/CPL candidates = pairs of Qs whose nodes are connected by depends_on, a coupling annotation, or a shared quantity in a constraint net (distance ≤ 2).
   - Priority = the node's value in PRIORITY.json/NEXT_ACTIONS.json (heuristic).
   - Output: `results/CX-QGRAPH/PAIRS.json`, with the pair, the edge/shared quantity, and the priority.
3. **Rewrite tasks/refill_v2.py** so that it uses Q_GRAPH_MAP.json + PAIRS.json:
   - DATX = graph pairs where both BT-DAT are finished;
   - CPL = graph pairs where both models met their criterion;
   - DAT order = by graph priority;
   - INT = the robust models, from the synthesis packets' lists, and not by word search;
   - add a BT-AUD-<HHMM> review packet per refill: a random sample of 10 finished DAT/DATX/INT (does the sample load, are the units right, does the reference say what is claimed).
   Keep the rest of the script's contract: real files, ALLOW_WEB, ≤45 min, the queue format with A/B only.
   Test it: `python3 tasks/refill_v2.py 0` must run without errors and print which packets it WOULD create. Do not queue anything yourself.
4. **Graph binding (feedback):** for the finished INT/DATX results, write proposed feedback records (the format per notes/GRAPH_WORKFLOW.md: target ID, result path/hash, quantity, unit, uncertainty, gate, outcome, review_state PENDING_INDEPENDENT_REVIEW) to `results/CX-QGRAPH/GRAPH_FEEDBACK_PROPOSED.json`. Do NOT run `./graph feedback` and do not change the graph files.

Deliver RESULTS.md starting with `# CX-QGRAPH` (the number mapped, how many pairs, examples) + pytest for the new refill logic. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only; GRAPH.json/current/generations must NOT be edited.
