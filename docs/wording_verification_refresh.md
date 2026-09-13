# Verification bindings after source wording changes

Before repair, nine recipes refused the edited source hashes before starting numerical children.
The numerical code and literal gate keys are unchanged; full serialized descriptions
can change, so old complete-output receipts remain archived rather than overwritten.

| Recipe | Prior passing observer gates | Binding before repair |
| --- | ---: | --- |
| aging_cells | 4/4 | Source mismatch |
| explicit_all_gates | 4/4 | Source mismatch |
| pipeline_cells | 4/4 | Source mismatch |
| strict_endocrine | 4/4 | Source mismatch |
| declared_pipeline | 4/4 | Source mismatch |
| additional_cells | 4/4 | Source mismatch |
| embedded_overall | 5/5 | Source mismatch |
| optional_coupling_measurement | 5/5 | Source mismatch |
| coupled_activation | 4/4 | Source mismatch |

Before execution: preserve every existing numerical gate, refusal control and timeout;
require two independent complete probe reports to match byte for byte; compare all
reported gate values with the archived receipt; retain old reports at their original
paths. Only source hash bindings may change in the existing observer code. The coupled
measurement is rebuilt before its exact-report contract is rebound and rerun.
The binding repair remains unaccepted until all source/input/receipt bindings validate.
Full-scope coverage is a separate requirement.
No new solver, verifier suite or tolerance change is introduced.


## Rebuilt evidence

All nine existing observers pass their unchanged gates twice, and the dependent
coupled-chain contract passes 4/4 twice after the complete measurement is frozen.
This is ten existing commands and twenty complete reports; every paired report is
byte-identical. All reported gate values and nonmetadata fields equal the archived
observations. Changed fields are source hashes, complete-output hashes and byte
lengths affected by prose. This does not claim that the complete serialized output
is identical across source revisions.

| Recipe | Rebuilt gates | Pair differs | Source binding |
| --- | ---: | ---: | --- |
| aging_cells | 4/4 | 0 bytes | Valid |
| explicit_all_gates | 4/4 | 0 bytes | Valid |
| pipeline_cells | 4/4 | 0 bytes | Valid |
| strict_endocrine | 4/4 | 0 bytes | Valid |
| declared_pipeline | 4/4 | 0 bytes | Valid |
| additional_cells | 4/4 | 0 bytes | Valid |
| embedded_overall | 5/5 | 0 bytes | Valid |
| optional_coupling_measurement | 5/5 | 0 bytes | Valid |
| coupled_activation | 4/4 | 0 bytes | Valid |

The underlying coupled measurement additionally passes 6/6 in both runs. Full
receipts and field-by-field comparisons are in `reports/wording_refresh/receipts.json`
and `reports/wording_refresh/comparison.json`. Original reports remain at their
original paths. Only exact source pins, the coupled reference-report path/hash,
and recipe metadata changed; numerical code, gate inventories and timeouts did not.

All fifteen recipe bindings now validate, including the six unaffected recipes.
The existing selected aging verifier runs twice locally with identical complete
reports matching the rebuilt cloud receipt; its eight malformed-binding controls
also pass. Full verification still exits 2 with zero numerical children because
only 32 of 242 current declarations have recipes: **210 remain unmapped**. This
repair restores existing coverage rather than adding new claims. The earlier
three cloud pytest timeouts remain recorded in `reports/wording_tests.json`.

The complete shared-origin batch arrays from the previous CPU measurement are also
retained in `reports/tetra_ray_batch_arrays.npz`; all eight array hashes match the
already committed report. No new timing measurement was taken during this repair.
