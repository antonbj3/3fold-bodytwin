# Koordinatorhandoff: fyra reviewbindningar

Four dispatch/feedback-par are executed **locally**, with frozen canonical receipt code. No import to workspace tasks/graph_runs has been run because it is outside the lane write limit. The targets are definition-only; no scientific admission. The following concrete bindings can be imported by the coordinator within its mandate. The result's SHA and originalpath are already verified in the feedback files. No new calculation is needed.

```bash
./graph dispatch --id BT-CTX-SURG-INCISION --lane LANE_SURGICAL_SYNTHESIS_R1_INCISION_IMPORT --kind review --reason "Import synthesisR1 lane-local review binding; negative results preserved, no scientific admission" --task-file tasks/lanes/LANE_SURGICAL_SYNTHESIS.md
./graph feedback --lane LANE_SURGICAL_SYNTHESIS_R1_INCISION_IMPORT --record results/LANE_SURGICAL_SYNTHESIS/GRAPH_COORDINATOR_FEEDBACK_R1_INCISION.json
./graph dispatch --id BT-CTX-SURG-COLLAGEN --lane LANE_SURGICAL_SYNTHESIS_R1_COLLAGEN_IMPORT --kind review --reason "Import synthesisR1 lane-local review binding; negative results preserved, no scientific admission" --task-file tasks/lanes/LANE_SURGICAL_SYNTHESIS.md
./graph feedback --lane LANE_SURGICAL_SYNTHESIS_R1_COLLAGEN_IMPORT --record results/LANE_SURGICAL_SYNTHESIS/GRAPH_COORDINATOR_FEEDBACK_R1_COLLAGEN.json
./graph dispatch --id BT-CTX-SURG-HEMOSTASIS --lane LANE_SURGICAL_SYNTHESIS_R1_HEMOSTASIS_IMPORT --kind review --reason "Import synthesisR1 lane-local review binding; negative results preserved, no scientific admission" --task-file tasks/lanes/LANE_SURGICAL_SYNTHESIS.md
./graph feedback --lane LANE_SURGICAL_SYNTHESIS_R1_HEMOSTASIS_IMPORT --record results/LANE_SURGICAL_SYNTHESIS/GRAPH_COORDINATOR_FEEDBACK_R1_HEMOSTASIS.json
./graph dispatch --id BT-CTX-SURG-HEALING --lane LANE_SURGICAL_SYNTHESIS_R1_HEALING_IMPORT --kind review --reason "Import synthesisR1 lane-local review binding; negative results preserved, no scientific admission" --task-file tasks/lanes/LANE_SURGICAL_SYNTHESIS.md
./graph feedback --lane LANE_SURGICAL_SYNTHESIS_R1_HEALING_IMPORT --record results/LANE_SURGICAL_SYNTHESIS/GRAPH_COORDINATOR_FEEDBACK_R1_HEALING.json
```
