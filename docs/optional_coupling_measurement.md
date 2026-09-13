# Optional coupling mechanism measurement

Before execution: replay the unchanged vitamin activation source twice in empty,
separate output directories. Each child has a 200 second limit, one numerical
thread, and no visible CUDA device. No sibling result is supplied. This is a
source-contract diagnostic, not a biological or performance claim.

| Mechanism | Source inspection | Measurement to collect |
| --- | --- | --- |
| Sibling comparison | Filters null comparisons before all() | Available comparisons and emitted gate |
| Optional self-test | Null when sibling absent; omitted from final gates | Null count and final inventory |
| Overall result | all() over non-null final gates | Literal booleans and overall result |
| Repetition | Whole JSON output | Bytes and SHA256 for both complete records |

Frozen gates, written before execution: children complete; complete result bytes
repeat; final gates are literal booleans and overall equals their conjunction;
all numerical floats are finite; coupling PASS requires at least one non-null
comparison and an available sibling. The last gate must not interpret all(empty)
as evidence. A failure is retained; no baseline edit or acceptance exception.

## Measured absent-input baseline

| Children | Whole record per child | Numerical gates | Optional nulls | Comparisons | Coupling | Overall |
| --- | --- | --- | --- | --- | --- | --- |
| 2, exit0 | 39684 bytes, identical SHA256 12ac3fd12f8f7b7c05268f5ff4bd46e9d53704b3e8a15bcecf90538a12b19808 | 28/29 | 1 | 0 | false | false |

All five measurement gates pass. The source's availability guard prevents a
vacuous pass for an absent file. This does not establish safe handling of an
existing file with an empty or incomplete data block.

## Preregistered synthetic schema controls

Before execution, test absent, empty-object, partial-half-life, complete-matching,
and complete-mismatching sibling records. Each arm runs the unchanged source twice
with the same 200 second child limit. Values are synthetic schema controls only.
The required evidence inventory is the three literal source half-life keys.
Frozen diagnostic gates: all ten children complete; complete bytes repeat within
each arm; records finite; source gate/overall values are literal booleans and
consistent; every emitted coupling PASS has all three non-null comparisons;
complete matching passes and complete mismatch fails. No missing key, null or
omitted self-test is accepted as evidence. This is a separate observer beside the
absent-input baseline, which remains unchanged.

## Measured schema table

All ten children exit0, are finite, and repeat complete bytes within each arm.

| Synthetic arm | Bytes | Source gates | Source overall | Available comparisons |
| --- | --- | --- | --- | --- |
| Absent | 39684 | 28/29 | false | 0/3 |
| Empty object | 39765 | 30/30 | true | 0/3 |
| Partial block | 40184 | 30/30 | true | 1/3 |
| Complete match | 40240 | 30/30 | true | 3/3 |
| Complete mismatch | 40246 | 28/30 | false | 3/3 |

The schema observer fails exactly one of six frozen gates: two distinct arms
emit PASS without complete evidence. The source remains unchanged and unpromoted.

## Separate evidence contract design and frozen gates

The measured empty/partial acceptance motivates a new explicit evidence contract,
not an exception to a required numerical gate. Declare a nonempty unique tuple of
comparison names. A missing/null declared comparison gives ABSTAIN; a complete
literal-boolean vector gives PASS only when all values are true, otherwise FAIL.
Unknown keys, non-boolean/non-null values, or malformed declarations raise.
ABSTAIN is never accepted. The caller must retain the original source result.

Before execution, require: complete positive and negative controls; empty/partial/
null abstention controls; malformed input refusal; nonmutation; full deterministic
verdict repetition; agreement on all ten frozen schema-observer legs. No source
row is upgraded by these synthetic contract controls.
