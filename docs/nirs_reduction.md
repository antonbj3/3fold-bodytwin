# NIRS reduction mechanism

Before design, the existing chain calls signal1320 times per leg with1052 distinct
optical input tuples. Complete observed and unobserved leg records repeat exactly.

| Flat4119-path fixture stage | First | Second |
| --- | ---: | ---: |
| Unobserved full leg | 294.737 ms | 297.621 ms |
| Path matrix product | 19.256 ms | 19.329 ms |
| Exponential | 63.758 ms | 64.169 ms |
| Two-column sum and normalization | 110.155 ms | 109.917 ms |
| Weighted moments and derivative | 31.947 ms | 31.873 ms |

The two-column strided reduction is the largest measured substage. A separate
chain sibling will use the last row of np.add.accumulate(weights,axis=0), preserving
ascending row addition order, then keep every remaining expression unchanged.
This adds one transient array with the shape of weights; it is not a memory saving.

Frozen acceptance before execution: all eight existing chain gates on flat and
curved fixtures; exact complete leg dictionaries and propagation arrays against the
original; exact repeats; unchanged inputs; both full-chain candidate timings below
both baseline timings for each fixture. No floating tolerance, reassociation, new
extinction model or clinical calibration is introduced. If exactness fails, retain
the mismatch count and do not promote the candidate.

Results: all five comparison gates and all eight original chain gates pass on
both fixtures. Fourteen reduction controls pass. All eight full propagation arrays
and complete leg dictionaries match exactly; source hashes match the measured run.

| Fixture | Baseline seconds, two legs | Candidate seconds, two legs |
| --- | --- | --- |
| Flat, 4119 paths | 0.262834 / 0.261987 | 0.210012 / 0.208811 |
| Curved, 3590 paths | 0.238953 / 0.239906 | 0.192976 / 0.193998 |

Modal CPU measurements use one numerical thread and reverse arm order for the
second leg. Exactness is scoped to the pinned runtime and these fixtures and
controls. The additional transient allocation is 65904 or 57440 bytes; no memory
improvement is claimed. Evidence: reports/nirs_reduction.json and its complete
NPZ capture; mechanism receipts: reports/nirs_next_phase.json and
reports/nirs_reuse_phase.json. Clinical calibration remains outside this result.

Post-change measurement: the flat fixture uses63.48/64.85ms for exponential,
49.88/50.11ms for accumulation and38.73/38.91ms for absorption. Of1320 signal
calls,1052 optical tuples are unique. A second per-propagation observer isolates
255 duplicates in the256-row zero-band control:48.123/47.934ms flat and
44.304/59.637ms curved. All four256-row physical variation batches have256
unique inputs. Observed/unobserved complete records and repeats match.

No duplicate-result candidate is implemented: the dominant repeated-input cost
belongs to a correctness control, so bypassing it would not establish a useful
physical-chain improvement. Retain the original control and its independent
calculations. Next useful target is the exponential stage on unique input batches.
Receipts: reports/nirs_reduced_phase.json and reports/nirs_duplicate_phase.json.
