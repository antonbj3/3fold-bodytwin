# CX-REFILL3 — the queue is starving: refill_v2.py now builds only ~1 packet (everything already exists). Extend it with NEW, well-reasoned packet types so that The swarm stays full until 21:50

Read: tasks/refill_v2.py (the current logic + templates), results/CX-QGRAPH (Q_GRAPH_MAP.json, PAIRS.json), results/CX-SWARMPLANNER (QUALITY_LOG.md: which errors The swarm makes; the templates are already tightened), NEXT_ACTIONS.json/PRIORITY.json, and the strategy (datasets + combining datasets first, then integration and coupling, with the graph driving).

## Tasks
1. Add at least 4 new packet families to refill_v2.py. Together they must give ≥ 300 potential packets. Suggestions, choose and improve:
   - **BT-DATX2:** combinations across domains along the graph's edges in PAIRS.json (distance 2, not just 1), where both DAT are finished;
   - **BT-DATG:** a data search for GRAPH NODES in NEXT_ACTIONS/PRIORITY that lack data, rather than for Q questions;
   - **BT-VAL:** validate an integrated module (BT-INT) against a SECOND independent dataset (from BT-DAT) that it was not built on;
   - **BT-DATQ:** data quality: re-run a DAT that the quality log flagged, with the stricter template;
   - **BT-MULTI:** couplings of 3 models along a graph path;
   - **BT-DENT2:** continue the dental coupling (results/BT-DENT-*: LOADIF/DATX) with the next step that the results propose.
2. Order within each refill: roughly 60 % data/combination, then validation/integration, then coupling. Use the graph priority.
3. Test with `python3 tasks/refill_v2.py 0`, which must list ≥ 80 packets it WOULD create, and pytest. Keep a backup copy of the current script (refill_v2.py.bak_<time>). Do not queue anything yourself: the timer bt-swarm-keepalive runs refill every 15 min when there are fewer than 40 runnable.
Deliver RESULTS.md starting with `# CX-REFILL3`: the new families, how many potential packets each, and examples. lane runner has full permissions in the workspace; graph files must not be edited; internal data (the collaborator/GC/restricted model data) must never go into ALLOW_WEB packets.
