# BodyTwin: use the graph in research

Use BodyTwin's own `./graph working` view. `john_geometry` and `general_bodytwin` are separate goals. The initial nine packets bind geometry, observability, joint fitting, material/FE uncertainty and arithmetic certificate results. They are a selected research view, not full coverage of the 3,950-node native graph.

```sh
./graph working status
./graph working rank --query 'the collaborator geometry' --limit 5
./graph working packet --id BT-IM2-JOINTFIT
./graph dispatch --id BT-IM2-JOINTFIT --lane <NEW_LANE> --kind review --reason '<why this test now>' --task-file tasks/lanes/<NEW_LANE>.md
./graph feedback --lane <NEW_LANE> --record results/<NEW_LANE>/GRAPH_FEEDBACK.json
./graph history --id BT-IM2-JOINTFIT
```

`dispatch` saves a packet, brief hash, selection reason and prior feedback in `tasks/graph_runs/`; it does not launch a process. Read `next_action_advice` and `previous_feedback` before repeating a branch. Choose `define`, `review` or `experiment` according to the actual task. Anatomical and physiological validity remain distinct from arithmetic guarantees.

Feedback fields: target ID, result path/hash, measured quantity, units, uncertainty, population/regime, frozen gate, baseline, outcome and `negative_result` boolean. Start with `review_state: PENDING_INDEPENDENT_REVIEW`. Reviewed task guidance requires a scoped JSON review, hashes, and coordinator registry admission. Self-asserted REVIEWED is insufficient. This mechanism updates task context and does not promote native scientific status or perform numerical fusion.

The failed BT-IM2 held-out demonstration gate is narrowly reviewed in `results/GRAPH_WORKFLOW_REVIEW_20260923/`. Persistent feedback changes the next proposed task to diagnosing the failure. It does not rule out every motion-to-anatomy method. The followup receipt is not an executed experiment.

Collect source-citation errors, repeated failed branches, strongest-baseline use, falsifiable next tests, scope errors, context bytes, observed token counts and elapsed time. Frozen paired inputs are in `results/GRAPH_WORKING_VIEW_20260923/trial/`; the current routing proxy is not an agent-performance result. Unobserved tokens remain unknown.

Current limits: curated packet coverage and heuristic priorities. New results require explicit binding/review. `./graph working build` regenerates only this view after source review. Native graph generations and the pinned engine remain separate.
