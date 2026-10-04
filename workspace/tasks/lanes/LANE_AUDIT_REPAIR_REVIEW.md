# LANE_AUDIT_REPAIR_REVIEW

Independent review of the source audit's findings and proposed corrections, so the coordinator can apply the corrections to the swarm's source families with a clear conscience. The reviewer must not have written the proposals. Results directory `results/LANE_AUDIT_REPAIR_REVIEW/`.

## Inputs

`results/LANE_SOURCE_CONSERVATION_AUDIT/`: AUDIT_TABLE.json (11 errors in 8 families: MITOSTRESS, Q005, Q012, Q052, Q090, Q115, Q146, Q168), proposals_r1/ (tested proposed corrections), witness_<KEY>.py, RAW_*_SOURCE_FINAL_R1.json. Sources: `tasks/free48/sources/<KEY>/` (read-only for you). The SOURCE_CONSERVATION_AUDIT lane is running round 2 at the same time (reviewable diffs + downstream exposure) — read its new files if they appear, but make your own assessment.

## For each finding

1. **Is the error real?** Derive it from the code itself (dimensional analysis, flow balance), without trusting the witness. The coordinator has already confirmed Q146 mg/mL (C·MW in mol/m³·g/mol = mg/L, labelled mg/mL).
2. **Is the correction right and minimal?** Does it change only what is wrong? Does it introduce new errors? Does it change dynamics or just reporting?
3. **Consequence**: which output fields change and by how much at the model's default parameters.
4. Verdict: APPLY (safe), APPLY WITH CHANGE (provide the amended diff), REJECT (incorrect finding or wrong correction), UNDETERMINED.

## Deliverables

REPAIR_REVIEW.json (finding → verdict → reason → final diff file), a directory `approved_patches/<KEY>.patch` for what should be applied, and a script `apply_and_test.sh` that applies the changes to a copy of the sources and runs the witnesses (red before, green after). Do not touch the actual sources — the coordinator applies the changes.
