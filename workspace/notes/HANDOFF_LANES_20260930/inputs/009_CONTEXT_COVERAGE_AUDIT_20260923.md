# Context gap before new BodyTwin experiments — 2026-09-23

Anton has pointed out that new knee work started without a full comparison against the private implementation. No new coordinator lanes before a finished plan. This note changes no scientific status.

## Checked on disk

`MERGED_GRAPH.json` is a typed index. Its anchor records refer with `detail_file: GRAPH.json` to the imported detail graph. A script not being mentioned in the index file does not mean that the reference is missing from the import.

| Exact path | Source graph | Workspace GRAPH.json | MERGED_GRAPH.json |
|---|---:|---:|---:|
| scripts/msk/static_opt_knee.py | 19 | 19 | 0 |
| scripts/msk/cartilage_contact.py | 5 | 5 | 0 |
| scripts/msk/contact_waveform_fmax_corrected.py | 1 | 1 | 0 |
| scripts/msk/jam_contact_decorr.py | 9 | 9 | 0 |
| scripts/msk/moment_arm_validation.py | 5 | 5 | 0 |

For `KNEE-CELL`, `MSK-KNEE-6DOF-JAM-ACL` and `HOLE-MOMENT-ARM-SENSITIVITY-UNWIRED-AT-KNEE-CELL` is `claim` and `cert_design` equal in source graph and import. This is a delimited check, not a proof of complete code index or scientific correctness. `~/research/sol6_recovery_20260923/BODY_CONTEXT_COVERAGE_20260923/coverage_probe.json`.

`results/GRAPH_WORKING_VIEW_20260923/working_graph.py` has nine manually selected records. `notes/GRAPH_WORKFLOW.md` already states that they do not cover the entire private graph. These records are not sufficient as the sole search for existing implementations and the strongest comparison methods.

The other session's ongoing inventory is in `results/MAP/private_inventory.tsv`, `graph_domain_nodes.tsv` and `overlap_today_vs_graph.tsv`. Reuse it, and distinguish thematic overlap from verified reusable implementation.

## Required material before new implementation

1. Search both the imported detail graph and private `source_repository/{scripts,docs,data,reports,evidence}` for the relevant measure, model and load case. Use the existing MAP inventory as an entry point. Absence of a hit in a selected working view is not evidence that work is missing.
2. Bind relevant nodes to exact scripts, results, data sources and previous corrections. Record file hash, input, output/unit, load case, limitations and what has actually been run. Distinguish literature seed, implementation, synthetic test, data comparison and independent review.
3. Write a short baseline table before a new experiment: what is reused, which previous attempt is tested and which still untested mechanism distinguishes the candidate. Preserve negative results and older corrections.
4. Compare on the same physical question, data, observable input and cost. Example: `cartilage_contact.py` describes passive flexion; X1b's Test 4 concerns contact force during walking. The files' existence proves neither equivalence nor that X1b is worse. That requires a concrete comparison.
5. Register missing code/result connections as coverage gaps in the working view. Do not upgrade hypotheses or result status through indexing. Generated graphs and source graphs remain read-only; changes go through the regular readers and review.

What is needed next is therefore a traceable check of prior work as part of selection, before further innovation jobs are defined.
