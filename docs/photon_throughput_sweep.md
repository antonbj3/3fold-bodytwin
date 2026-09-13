# Bounded photon throughput sweep

The final directive permits one sweep, then stop. No engine code changes.
Existing certified-cell engine and frozen window2 reference, same public voxel
fixture, seed20260913, complete five-array outputs and transport window unchanged.
One GPU thread per photon, fixed64-thread blocks; no per-photon host loop.
Counts fixed before execution:10000,100000,1000000 on L4 and H100.

Warm compilation uses1024photons. Construction/compilation are reported separately.
At each count run the reference once, then the unchanged candidate twice. Report
synchronized full API seconds and photons/s; launch timing includes synchronization
and dispatch. Compare complete candidate arrays against the reference and repeat,
exact packet energy, zero leaks/caps and immutable caller launch inputs. Hashing,
comparison and scalar summaries occur outside the timed API. First failed physical,
parity or launch-contract gate stops that worker; no retry or replacement target.

Report throughput ratios across the two decades. A saturation point is only claimed
if supported by observed flattening; no extrapolated large-count speed claim.
Original numerical tolerances and event/window bounds remain unchanged. This is a
full-output transport workload, not an equivalence claim to other photon engines.
The existing broad final pre-push validation budget is reserved; no push authorized.

Results: all18measured calls pass physics, exact full-array reference parity,
repeat and launch-contract checks. Each candidate has two bit-identical full
outputs per count; comparisons occurred in the worker before hashes were returned.
No full-array transfer to the local workspace was needed for this measurement.

| Device | Photons | API seconds, two legs | Photons/s, two legs | Reference seconds, one leg |
| --- | ---: | --- | --- | ---: |
| L4 | 10000 | 0.265819 / 0.261240 | 37619.5 / 38279.0 | 0.322320 |
| L4 | 100000 | 1.765083 / 1.776934 | 56654.6 / 56276.7 | 1.947885 |
| L4 | 1000000 | 16.548926 / 16.618574 | 60426.9 / 60173.6 | 18.054821 |
| H100 | 10000 | 0.129444 / 0.127474 | 77253.8 / 78447.2 | 0.191599 |
| H100 | 100000 | 0.299685 / 0.298468 | 333683.7 / 335044.6 | 0.508426 |
| H100 | 1000000 | 2.453203 / 2.446641 | 407630.3 / 408723.6 | 4.327816 |

| Device | Throughput gain10000 to100000 | Gain100000 to1000000 |
| --- | ---: | ---: |
| L4 | 1.488x | 1.068x |
| H100 | 4.295x | 1.221x |

L4 approaches a plateau around100000–1000000photons: the last decade adds
about6.8% throughput. H100 still gains about22.1%; the sweep does not establish
its final saturation point. The best measured rate is408724photons/s, so this
workload does not reach100million photons/s. No further attempt is authorized
by the bounded directive. At1000000photons,1368967976events are processed;
energy1073741824000000, zero leaks/caps, and192873120output bytes per call.
Full API timing includes allocation and output copies; candidate launch times at
1000000are16.459/16.526s on L4 and2.3683/2.3675s on H100.

Initial scene construction baseline/candidate: L42.332/1.330s, H1002.183/0.859s.
The first construction includes device initialization, so these are lifecycle
observations, not a constructor speed comparison. Compilation plus1024-photon
warmup baseline/candidate: L45.464/5.213s, H1003.797/3.512s. Neither cost is
hidden in a cold-speed claim. The reference has one timing leg per size: its
comparison is for full-output parity, not a repeated speedup claim.

Source-bound receipts: reports/photon_throughput_l4.json and
reports/photon_throughput_h100.json. No engine code changed. Sweep complete;
stop. The final pre-push pytest/full-tree review remains reserved, not claimed
as completed here. No push.
