# CX-ADJUDICATE — settle the 7 deviating The swarm reproductions + rerun the 21 KKT gaps in the what-if curves

Source: `results/CX-BOOKKEEP3/RESULTS.md` ("The swarm reproductions that DEVIATE", "Three recommended next steps").

## 1. Adjudication (model hierarchy: a The swarm finding against swarm_worker/lane_runner does not count until it has been verified by a stronger model)
For each pair below:
- read the original and the reproduction, including the code and data in the packets;
- recompute the disputed number with your OWN minimal code;
- give a verdict: ORIGINAL HOLDS / ORIGINAL WRONG (correction with number) / A MATTER OF DEFINITION (both reasonable; state which one should apply) / UNKNOWN.

| Reproduction | Row | Issue |
|---|---|---|
| BT-Q-BT-D018 | A1148 | the density factor M3 deviates 3.64 % |
| BT-S-BT-D018 | A1236 | the smallest rank-changing density factor 1.0000587 vs 1.0378289; the original took the largest? |
| BT-S-BT-AN-G68 | A1171 | A577 at a 16 % threshold |
| BT-S-BT-B381 | A1182 | rank change 2/4 vs the 3/4 requirement |
| BT-S-BT-D017 | A1235 | inertia normalisation +20 % gives 4/20 ≥ 1 % |
| BT-S-BT-D024 | A1242 | a +20 % null-model threshold flips the right-hip decision |
| BT-S-BT-D028 | A1245 | strict source-unit normalisation gives 0/3 PASS |

If ORIGINAL WRONG: add an inline correction to the original's row in `notes/RESULTS_INDEX.md` in the form `[CORRECTED 24/9 by CX-ADJUDICATE: …]`, as with A170. This is an exact string replacement in that row only; do not rewrite anything else.

## 2. Rerun the 21 KKT gaps
The gaps are in the CX-WHATIF2 curves (A369; points with 140/141 KKT; list in `results/CX-WHATIF2/RESULTS.md`).
- Rerun on OVH (`tasks/cloud_run.sh`, BodyTwin ≤ 12 vCPU, finish by 23:59; the bundle is at /opt/bt/shared/cxw_bundle on OVH, or use the cloud_shard.py pattern) with a tighter tolerance or more iterations.
- Fill in the curves and report whether the derivative intervals change.

## Deliverables
`results/CX-ADJUDICATE/RESULTS.md` starting with `# CX-ADJUDICATE`: a verdict table plus the KKT gaps before/after, and results.json.

lane_runner has full permissions in the workspace.
