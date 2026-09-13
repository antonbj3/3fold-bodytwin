# Photon query-distance tail

Before design, the complete captured query set contains187882queries. Median
limit0.0546264 and99th percentile0.591913 contrast with maximum650.712.

| Window threshold | Queries exceeding threshold |
| ---: | ---: |
| 1 | 774 |
| 2 | 721 |
| 4 | 701 |
| 8 | 670 |

The launch64 API costs321–323ms on L4 and216–233ms on H100 at1024photons.
The current fixed8 search window may admit excess triangle candidates in the
small long-query tail. This is a hypothesis, not measured candidate counts.
A temporary observer compares windows8/4/2/1 with reverse second-leg order on
both devices,1024and8192photons, preserving128window attempts and the exact
full-query fallback. All triangle intersection arithmetic and tie-breaking stay
unchanged. Full output bytes are compared, not only aggregate photon counts.

Frozen gates for any subsequent sibling: all original five transport gates,
full five-array parity, exact repeats, exact energy, zero caps/leaks, unchanged
caller inputs and both candidate times below both baseline times for all four
count/device cases. Existing boundary/tie/fallback controls must also pass on
CPU and GPU. Failed cases remain explicit; no window/event limit loosening.

Measured API milliseconds, before sibling design:

| Device | Photons | Window | First | Second |
| --- | ---: | ---: | ---: | ---: |
| L4 | 1024 | 8 | 327.070 | 330.101 |
| L4 | 1024 | 4 | 285.291 | 283.831 |
| L4 | 1024 | 2 | 279.771 | 269.724 |
| L4 | 1024 | 1 | 273.123 | 275.607 |
| L4 | 8192 | 8 | 379.493 | 381.805 |
| L4 | 8192 | 4 | 334.632 | 333.999 |
| L4 | 8192 | 2 | 321.846 | 319.549 |
| L4 | 8192 | 1 | 321.427 | 321.698 |
| H100 | 1024 | 8 | 221.610 | 216.682 |
| H100 | 1024 | 4 | 186.686 | 224.281 |
| H100 | 1024 | 2 | 177.061 | 177.701 |
| H100 | 1024 | 1 | 178.599 | 179.875 |
| H100 | 8192 | 8 | 220.688 | 219.910 |
| H100 | 8192 | 4 | 190.289 | 233.790 |
| H100 | 8192 | 2 | 182.945 | 177.915 |
| H100 | 8192 | 1 | 177.980 | 179.022 |

All32observations have zero output-byte mismatches. Window4 fails the strict
H100 timing criterion in both sizes:224.281ms versus216.682ms and233.790ms
versus219.910ms. Window2 passes all four observed strict comparisons and is
selected for a separate sibling. The128attempt bound remains unchanged; the
full-query fallback occurs earlier in distance, preserving exact intersection
and face-tie arithmetic. Fresh complete transport comparison is still required.

Initial fresh comparison upload omitted the prepared-launch baseline module:
ImportError before the first timed transport sample, zero samples. Boundary
controls ran, but no transport comparison is claimed from that attempt. The
public source manifest was corrected; numerical source and gates unchanged.

Fresh sibling result: all5existing gates pass in all4cases. The10existing
boundary/tie/fallback cases pass twice on CPU and each GPU, with exact arrays
across CPU/GPU. All80downloaded transport arrays,919996928bytes, match their
receipt hashes. The transport kernel AST is identical to the frozen segmented
kernel; only its separate query-window function changes.

| Device | Photons | Baseline ms, two legs | Sibling ms, two legs |
| --- | ---: | --- | --- |
| L4 | 1024 | 319.578 / 318.315 | 268.563 / 267.903 |
| L4 | 8192 | 374.043 / 374.105 | 311.544 / 315.047 |
| H100 | 1024 | 218.059 / 234.086 | 164.843 / 181.271 |
| H100 | 8192 | 225.505 / 224.018 | 184.077 / 183.920 |

Energy1099511627776/8796093022208 and zero leaks/caps hold in every arm.
Receipts: reports/photon_window2_* JSON; mechanism reports/photon_window_phase_*.
Measurements are bounded to these two counts and the fixed fixture/seed. No
million-photon, universal throughput or clinical claim. Baseline files unchanged.
