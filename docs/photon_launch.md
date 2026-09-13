# Photon launch mechanism

Measured prepared transport residual before launch investigation:

| Stage | L4 milliseconds | H100 milliseconds |
| --- | ---: | ---: |
| Frozen transport kernel, 1024 photons | 382.569 / 382.667 | 210.689 / 210.700 |
| Output readback | 24.752 / 24.742 | 15.246 / 15.320 |

Source inspection finds default block_dim=256. At 1024 photons this produces
four blocks; at 8192 it produces 32. These are launch-shape counts, not measured
occupancy. A temporary Modal observer varies only the launch block size through
256, 128, 64 and 32, reversing order on the second leg. The frozen source files
and per-thread arithmetic stay unchanged. A launch sibling is not designed until
this table has measured results.

Acceptance fixed before candidate execution: all five existing full-transport
gates, complete five-array parity against the prepared baseline, exact repeats,
zero leaks and caps, exact packet energy, unchanged caller inputs, and both
candidate times below both baseline times for each tested count and device.
No tolerance, seed, transport window or event cap changes. Counts are 1024 and
8192 on each device. Any failing case remains recorded and unpromoted.

Observed API milliseconds before sibling design:

| Device | Photons | Block | First | Second |
| --- | ---: | ---: | ---: | ---: |
| L4 | 1024 | 256 | 413.535 | 410.917 |
| L4 | 1024 | 128 | 342.744 | 346.009 |
| L4 | 1024 | 64 | 318.805 | 316.024 |
| L4 | 1024 | 32 | 304.271 | 306.284 |
| L4 | 8192 | 256 | 420.327 | 418.013 |
| L4 | 8192 | 128 | 432.731 | 435.500 |
| L4 | 8192 | 64 | 395.067 | 372.364 |
| L4 | 8192 | 32 | 371.510 | 370.948 |
| H100 | 1024 | 256 | 236.401 | 272.232 |
| H100 | 1024 | 128 | 220.188 | 222.366 |
| H100 | 1024 | 64 | 217.537 | 217.700 |
| H100 | 1024 | 32 | 235.573 | 236.681 |
| H100 | 8192 | 256 | 240.833 | 243.243 |
| H100 | 8192 | 128 | 230.341 | 261.632 |
| H100 | 8192 | 64 | 224.779 | 222.111 |
| H100 | 8192 | 32 | 223.429 | 222.749 |

All full output arrays match for all observed launches. Block32 fails the strict
speed criterion on H100 at1024:236.681ms exceeds baseline236.401ms. Block128
fails on L4 at8192:432.731/435.500ms versus420.327/418.013ms.
The separate sibling will fix block_dim=64, which beats both default timings in
all four observed cases. Selection uses these observations; a fresh paired run
must pass the preregistered gates. No universal or large-count speed claim.

Fresh sibling comparison: all five existing gates pass in each of four cases.
All80 full output arrays were fetched and checked against their recorded hashes.

| Device | Photons | Default256 milliseconds | Sibling64 milliseconds |
| --- | ---: | --- | --- |
| L4 | 1024 | 412.805 / 411.262 | 322.675 / 321.052 |
| L4 | 8192 | 417.088 / 416.776 | 376.222 / 377.347 |
| H100 | 1024 | 235.466 / 252.957 | 233.408 / 216.259 |
| H100 | 8192 | 239.644 / 334.878 | 224.120 / 224.843 |

H100 baseline timing varies substantially; at1024 the strict worst-candidate
margin is only2.058ms. These are bounded two-leg observations, not a stable
universal throughput claim. Energy totals are1099511627776 and8796093022208,
with zero leaks/caps in every arm. Evidence: reports/photon_launch64_* JSON
receipts and reports/photon_launch_phase_* mechanism receipts. Baseline and
transport kernel remain unchanged; caller must explicitly select this sibling.
