# CX-SWARMPLANNER — Sol runs the The swarm machinery WITHOUT the coordinator (coordinator's limit is exhausted) until the OVH lease ends (22:24)

Anton: "it has to be properly filled; put a Sol agent on it right away."
Infrastructure (already running):
- the timer `bt-swarm-keepalive` runs `tasks/swarm_keepalive.sh` every 15 min: it restarts the drivers, stops hung OVH agents, and runs `tasks/refill_v2.py 80` when the queue has fewer than 40 runnable jobs;
- the queue is tasks/lanes/bt_queue.txt; the OVH driver is tasks/lanes/ovh_agents/bt_queue_ovh.sh (cap in tasks/lanes/ovh_agents/cap); the log is tasks/lanes/swarm_keepalive.log and tasks/lanes/ovh_agents/queue_ovh.log;
- the strategy (Anton): DATASETS and COMBINING datasets first, then integration into bodytwin_core, then couplings. The graph drives pairs and priority (results/CX-QGRAPH, if it is finished).

## Your job (loop over rounds; each round ≈ 20–30 min)
1. **Check the flow:**
   - the success rate of the last 40 `fetched` lines (results=yes/no);
   - the running count against the cap; free memory on OVH (ssh as in the keepalive script).
   Adjust the cap in tasks/lanes/ovh_agents/cap:
   - high success and memory > 8 GB → raise it by 5 (max 40);
   - many `results=no` or The swarm server errors in agent.log → lower it.
   Accounts: A and B are primary. Test C with a trivial job via `/opt/agents/run_profile_ovh.py C swarm` on OVH; if it works, add C back in refill_v2.py's round-robin.
2. **Quality before quantity:**
   - Read the BT-PLAN-*/NEXT_PACKETS.md that the refill will consume. Remove low-value or meaningless packets BEFORE they are built (edit NEXT_PACKETS.md; mark the reason with a comment).
   - Read a random sample of 5 finished BT-DAT/BT-DATX/BT-INT/BT-AUD per round. Where The swarm's work is poor (the sample was not loaded, the units are wrong, the reference does not say what is claimed), write it into `results/CX-SWARMPLANNER/QUALITY_LOG.md` and adjust the templates in refill_v2.py (the brief texts) so that the error does not repeat.
3. **Keep the queue properly full:**
   - If refill_v2.py produces too few packets (e.g. everything already exists), extend it with new, well-reasoned packet types that follow the strategy. Examples: combinations across domains according to the graph; data packets for nodes in NEXT_ACTIONS that lack data; validation against measured values.
   - Test with `python3 tasks/refill_v2.py 0` before the timer uses it. Never break the script: keep a working backup copy (refill_v2.py.bak_<time>).
4. **At ~21:30:** write `results/CX-SWARMPLANNER/RESULTS.md` starting with `# CX-SWARMPLANNER`. Report:
   - the number of packets per type that ran and succeeded;
   - the success rate per account;
   - the top 10 results in value to BodyTwin (with paths);
   - the quality problems and what you changed.
   After 21:50: queue nothing new.
Never edit GRAPH.json/current/generations. No emails, pushes or credentials beyond what already exists. Internal data (the collaborator, GC, restricted model data) must NEVER go into packets with ALLOW_WEB.
