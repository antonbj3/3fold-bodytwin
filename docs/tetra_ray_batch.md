# Shared-origin CPU traversal

The frozen `tetra_ray_walk_v1.py` locates the launch point over all cells on every
`trace` call. The existing CUDA comparison supplies 512 directions with one origin.
First measure two complete baseline batches and separately the identical point-location
expression repeated 512 times. No candidate exists at this measurement stage.

Before candidate execution, the gates are fixed: all complete float64 interval arrays
and cell sequences must equal the frozen baseline bit for bit, both complete repetitions
must agree, original singular/invalid refusals must remain, input arrays must remain
unchanged, and both paired batch timings must improve. The existing independent oracle
parameter bound remains 1e-12. No change to geometric arithmetic or acceptance tolerances.

## Measured mechanism before design

| Quantity, 1296 cells / 512 rays | First | Second |
| --- | ---: | ---: |
| Complete baseline batch | 83.204 ms | 71.051 ms |
| Repeated point location alone | 50.647 ms | 50.023 ms |
| Complete interval records | 2735 | 2735 |
| Different output bytes between repeats | 0 | 0 |

The repeated search costs about 61–70% of the observed batch, so eliminating it
cannot imply a 512-fold gain for the complete traversal. Raw values and the frozen
source hash are in `reports/tetra_ray_batch_baseline.json`.

## Candidate design after measurement

A separate `tetra_ray_batch_v1.py` exposes `TetraRayBatch.trace_batch(origin,directions)`.
It prepares the same complete barycentric array once per call, then uses the unchanged
per-ray normalization and neighbor traversal. It returns a list of complete float64
record arrays, in input order. Inherited scalar `trace` remains the baseline. A fresh
origin is located on each batch call; no mutable cross-call location cache is retained.
The method performs no artifact writes. Empty batches return an empty list after shape
and finiteness validation. Any invalid ray or singular traversal raises `ValueError`.

## Results

| Paired batch, 512 rays / 2735 interval records | First | Second |
| --- | ---: | ---: |
| Frozen scalar calls | 70.991 ms | 69.761 ms |
| Shared-origin batch | 16.597 ms | 16.526 ms |
| Complete batch gain | 4.277x | 4.221x |
| Different record/offset bytes | 0 | 0 |

All six preregistered gates pass. The 32 original/inverted, alternate-origin,
singular/invalid and empty-batch controls pass. The independent exhaustive interval
oracle differs by at most 4.440892098500626e-16, below the unchanged 1e-12 bound.
Both full record arrays and offsets match baseline and repeats bit for bit.
Timing order is baseline/candidate then candidate/baseline after one warmup each;
local CPU work used the shared exclusive hardware lock and one numerical thread.
This is a two-pair synthetic CPU observation, not a general throughput guarantee.

Run `PYTHONPATH=src python probes/optics/tetra_ray_batch_probe.py`.
The probe explicitly writes `reports/tetra_ray_batch.json` and the complete
`reports/tetra_ray_batch_arrays.npz`; the engine itself writes nothing.
Existing scalar and CUDA callers remain unchanged. Callers with shared-origin rays
can adopt the new batch API; varying-origin and photon workloads remain unmeasured.
