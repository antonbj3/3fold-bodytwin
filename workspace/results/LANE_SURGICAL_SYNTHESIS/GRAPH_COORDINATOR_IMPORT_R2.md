# R2 coordinator handoff: four review bindings

R1 import is verified. R2 is locally dispatched and bound to feedback, no source graph mutation/science admission. The coordinator's write authority is needed for canonical tasks/graph_runs; the following concrete records are ready:

```bash
./graph dispatch --id BT-CTX-SURG-INCISION --lane LANE_SURGICAL_SYNTHESIS_R2_INCISION_IMPORT --kind review --reason "Import synthesis R2 diagnosis and conditional acquisition planning; proxy and information failures retained" --task-file tasks/build_night/steer/LANE_SURGICAL_SYNTHESIS.md
./graph feedback --lane LANE_SURGICAL_SYNTHESIS_R2_INCISION_IMPORT --record results/LANE_SURGICAL_SYNTHESIS/GRAPH_COORDINATOR_FEEDBACK_R2_INCISION.json
./graph dispatch --id BT-CTX-SURG-COLLAGEN --lane LANE_SURGICAL_SYNTHESIS_R2_COLLAGEN_IMPORT --kind review --reason "Import synthesis R2 diagnosis and conditional acquisition planning; proxy and information failures retained" --task-file tasks/build_night/steer/LANE_SURGICAL_SYNTHESIS.md
./graph feedback --lane LANE_SURGICAL_SYNTHESIS_R2_COLLAGEN_IMPORT --record results/LANE_SURGICAL_SYNTHESIS/GRAPH_COORDINATOR_FEEDBACK_R2_COLLAGEN.json
./graph dispatch --id BT-CTX-SURG-HEMOSTASIS --lane LANE_SURGICAL_SYNTHESIS_R2_HEMOSTASIS_IMPORT --kind review --reason "Import synthesis R2 diagnosis and conditional acquisition planning; proxy and information failures retained" --task-file tasks/build_night/steer/LANE_SURGICAL_SYNTHESIS.md
./graph feedback --lane LANE_SURGICAL_SYNTHESIS_R2_HEMOSTASIS_IMPORT --record results/LANE_SURGICAL_SYNTHESIS/GRAPH_COORDINATOR_FEEDBACK_R2_HEMOSTASIS.json
./graph dispatch --id BT-CTX-SURG-HEALING --lane LANE_SURGICAL_SYNTHESIS_R2_HEALING_IMPORT --kind review --reason "Import synthesis R2 diagnosis and conditional acquisition planning; proxy and information failures retained" --task-file tasks/build_night/steer/LANE_SURGICAL_SYNTHESIS.md
./graph feedback --lane LANE_SURGICAL_SYNTHESIS_R2_HEALING_IMPORT --record results/LANE_SURGICAL_SYNTHESIS/GRAPH_COORDINATOR_FEEDBACK_R2_HEALING.json
```

No import into source_repository's tasks/graph_runs has been executed by this lane inR2.
