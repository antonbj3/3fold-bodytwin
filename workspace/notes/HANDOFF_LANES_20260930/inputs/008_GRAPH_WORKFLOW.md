# Dental: use the graph in research

The native BodyTwin/manufacturing import and the dental research view have separate status. `./graph working` now reads dental drafts, proposed goal paths and bound recent results. Two conflicting draft IDs are exposed as conflicts. Definitions and candidate experiments are listed separately.

```sh
./graph working status
./graph working rank --kind experiment_ready --limit 5
./graph working packet --id DENT-DESIGN-J3-CREST-STRESS-SURROGATE
./graph dispatch --id DENT-DESIGN-J3-CREST-STRESS-SURROGATE --lane <NEW_LANE> --kind experiment --reason '<why this test now>' --task-file tasks/lanes/<NEW_LANE>.md
./graph feedback --lane <NEW_LANE> --record results/<NEW_LANE>/GRAPH_FEEDBACK.json
./graph history --id DENT-DESIGN-J3-CREST-STRESS-SURROGATE
```

`dispatch` saves an immutable packet, task hash, selection reason and prior feedback under `tasks/graph_runs/`; it does not launch a process. Read its `next_action_advice` and `previous_feedback` when writing the next brief. Use `--kind define` for missing variables/relations, and `review` for evidence checks. Experiment readiness means a test is specified; it is not a validated design or accepted physical model. Normal lane resource/approval rules still apply.

`GRAPH_FEEDBACK.json` uses the packet's result-binding fields: target ID, result path/hash, measured quantity, units, uncertainty, population/regime, frozen gate, baseline, outcome and `negative_result` boolean. Initially use `review_state: PENDING_INDEPENDENT_REVIEW`. Reviewed task guidance additionally requires a scoped review JSON, exact hashes, and the coordinator review registry. A lane cannot promote itself by writing REVIEWED. Source statuses and statistical fusion are untouched. Hash drift blocks reuse until reviewed.

The narrow 173/180 GEOM coverage failure is independently checked and recorded in `results/GRAPH_WORKFLOW_REVIEW_20260923/`; its graph receipts demonstrate persistent feedback into the next proposed task. The followup receipt does not mean a followup experiment ran.

Measure graph usefulness with source-citation accuracy, repeated failed branches, baseline parity, falsifiable next tests, scope errors, context bytes, observed tokens and elapsed time. Unknown token measurements stay unknown. Frozen comparison inputs are in `results/GRAPH_WORKING_VIEW_20260923/trial/`. Retrieval diagnostics alone do not show improved agent performance.

Current limits: incomplete dental mathematical coverage, manually curated result bindings and heuristic priorities. New result files require explicit binding/review; an output directory alone does not refresh the view. `./graph working build` rebuilds only this research view after reviewing changed inputs. Native `./graph refresh` is a separate import operation.
