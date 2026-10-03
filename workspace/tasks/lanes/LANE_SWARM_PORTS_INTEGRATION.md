# LANE_SWARM_PORTS_INTEGRATION

Turn tonight's reviewed swarm couplings into tested, graph-bound BodyTwin modules. Results directory `results/LANE_SWARM_PORTS_INTEGRATION/`.

## Underlag

`results/LANE_SWARM_FINDINGS_REVIEW/PORTS_FROM_SWARM.json` — 14 conditional model ports from 25 reviewed swarm results (3 CONFIRMED, 20 PARTIAL, 2 WRONG), with validity conditions and correction status. The source families in `tasks/free48/sources/` were corrected tonight (15 families; see the catalogue context and `results/LANE_AUDIT_REPAIR_REVIEW/PATCH_PLAN_R*.json`); Q160 uncorrected, Q084/BIORESP undetermined.

## Leverans

1. Package `ports/` in the results directory: one module per port, with declared input/output quantities (unit), validity domain and status (CONFIRMED/CONDITIONAL), built on the **corrected** source models. The port refuses (raise) outside its validity domain.
2. Tests per port: dimensional check, conservation where relevant, the source job's decisive number reproduced (or a discrepancy explained by the correction), and edge cases.
3. Graph binding: write GRAPH_COORDINATOR_FEEDBACK_<port>.json against the respective BT-CTX targets (kind review/define) as SURGICAL_SYNTHESIS did; the coordinator imports.
4. RESULTS.md for readers: which ports are now usable building blocks in BodyTwin and under which conditions.

Strongest control: the source family's own model without the port. No internal data.
