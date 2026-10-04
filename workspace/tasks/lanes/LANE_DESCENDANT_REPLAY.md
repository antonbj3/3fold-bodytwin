# LANE_DESCENDANT_REPLAY

Swarm jobs created before the corrections (06:15, 06:45, 07:45, 09:00) ran source models with 32 known errors in 21 families. The source review's rerun of 7 conclusions reversed 1. This lane corrects the archive for the most important ones. Results directory `results/LANE_DESCENDANT_REPLAY/`.

## Material

`results/LANE_SOURCE_CONSERVATION_AUDIT/`: EXPOSURE_JOB_LISTS_R2/, DESCENDANT_FIELD_EXPOSURE_R2.jsonl, CONCLUSION_REASSESSMENT_R3.json, AUDIT_TABLE_R4.json. Corrected sources in `tasks/free48/sources/` (backups in `tasks/free48/sources_backup_20261001_audit*/`). Swarm jobs: `external_mount<ID>/`.

## Task

1. Choose the exposed jobs with the highest research value (research_value.assess) and whose conclusion directly relies on a corrected field — up to ~30.
2. Rerun each job's decisive calculation with the corrected source (copy in the results directory; never change the job directory). Report: original number → corrected number → does the decision change?
3. Write `CORRECTIONS.jsonl` (job ID, field, old/new number, decision changed yes/no, file:line) and a summary: fraction of reversed conclusions, per error class.
4. For jobs whose decision reverses: write GRAPH_COORDINATOR_FEEDBACK_<ID>.json (review, negative_result where applicable) so the graph history is corrected; the coordinator imports it.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
