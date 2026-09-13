# Coupled activation source replay

## Source mechanism before execution

| Stage | Required input | Output | Source acceptance |
| --- | --- | --- | --- |
| cardiac_output_geometric | None | cardiac_output_geometric_results.json | gates and overall_pass |
| renal_filtration | Complete cardiac output record | renal_filtration_results.json | embedded gates.overall_pass |
| calcium_pth_vitd | Complete renal record | calcium_pth_vitd_results.json | gates and overall_pass |
| vitamin_d_activation | Complete calcium record | vitamin_d_activation_results.json | gates and overall_pass_strict_all |

No manually fabricated sibling is used. Each of two isolated chains executes all
four unchanged scripts in order, with one numerical thread and a 200 second limit
per child. This first measurement records the exact emitted inventories before
any integration adapter is designed. No clinical or performance acceptance.

Frozen measurement gates before execution: all eight children complete; every
complete stage record is byte-identical between chains; all emitted gate values
are literal booleans and true; all numeric floats are finite; downstream stages
leave earlier full artifacts unchanged; the final coupling has exactly the three
declared non-null comparisons, all true. Any failure is retained numerically;
no projection or timestamp exemption, no source or tolerance changes.

## Measured table before contract design

| Stage | Literal passing gate entries | Bytes | Complete repeat |
| --- | --- | --- | --- |
| Cardiac output | 16/16 | 15319 | exact |
| Renal filtration | 18/18, including two summaries | 14296 | exact |
| Calcium regulation | 15/15 | 40087 | exact |
| Vitamin activation | 30/30 | 41134 | exact |

Six measurement gates pass; all eight children exit0. Total per-chain full output
110836 bytes. Every source remains unchanged; three final comparisons are present
and true. The earlier empty/partial sibling false positives remain failures.

## Frozen integration contract before execution

Pin the complete measurement report and source files. Replay the chain, require
the exact two-leg/four-stage inventory and literal gate-key inventory per stage,
true literal overall values, and all six literal measurement gates. Compare the
entire fresh report against the frozen report, which contains complete artifact
hashes. Refuse missing/duplicate legs, renamed/extra gates, integer booleans,
false overall, false source gate, and missing artifact bytes. A selected recipe
must repeat this observer twice independently. This proves the supplied complete
chain only, not malformed-input safety or biological validity.

## Verification result

The strict observer passes4/4, including8malformed-record refusals. The selected
recipe runs it twice independently and obtains identical complete report hashes.
Two further isolated Modal CPU invocations match the local complete observer and
measurement reports exactly. Each invocation rebuilds both four-stage chains;
no cached sibling artifact replaces a source execution.

Complete observer SHA256:
`f0b2bf701d7bd86ef02c3b7e7fa7c064af0c004d3dffd0018dcae27caa5d0d79`.
Complete measurement SHA256:
`fc2085a3ecc242b1e27bd0f4972d8a4eaffb8fd0abe9a3c1b2ba2fdb2f636c95`.

Coverage is32of239fresh declarations in15explicit recipes;207remain unmapped.
The vitamin source remains failed for incomplete sibling inputs. This successful
complete-input chain neither overrides that negative nor claims universal CPU
identity, clinical validation, or throughput.
