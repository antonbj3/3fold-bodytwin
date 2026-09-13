# Photon RT query experiment

## Preregistered mechanism measurement

The existing photon path implementation calls a Warp AABB traversal and performs
its own double-precision triangle test for every transport event. There is no
OptiX photon implementation in this repository. SER requires that missing native
backend; published SER speedups are not predictions for this implementation.

Before candidate design, capture the first256events from each of128photons in the
unchanged public atlas transport, seed20260912. Extra observations must leave all
five complete transport output arrays identical to the frozen source. Capture
twice; preserve every captured start, direction, limit and selected triangle.
The finite prefix is a query microbenchmark, not an unbiased photon workload or
full-transport speed claim. Misses and endpoint ties remain included.

Frozen measurement gates: original/capture full transport identity; two complete
captures identical; standalone original query returns every captured face;
nonempty finite query inventory. Two timing legs follow warmup, with explicit
synchronization and no competing job in that worker. Write measured table before
implementing the candidate. No physical tolerance or transport baseline changes.

## External prerequisites

OpenOximetry1.1.1 is a restricted-access resource requiring registration and a
signed data-use agreement. No dataset files or matched calibration observations
were obtained. Its public metadata is not an optical calibration dataset.
Source: https://physionet.org/content/openox-repo/1.1.1/ (accessed2026-09-13).

OptiX reordering requires an OptiX ray-generation program. The current Warp
transport cannot call optixReorder directly. Reordering benefits depend on the
workload; the native query boundary must be measured before a full transport port.
Source: https://developer.nvidia.com/rtx/ray-tracing/optix

## Prefix result and next measurement

| Queries | Surface hits | Mesh triangles | Baseline leg0 | Baseline leg1 |
| --- | --- | --- | --- | --- |
| 22175 | 130 | 2754304 | 0.03724375 ms | 0.03446445 ms |

All4measurement gates passed. This prefix alone does not justify an RT port.
Before deciding, repeat the same capture with capacity50000per photon, covering
every event allowed by the unchanged transport. Keep128photons, seed and optical
properties fixed, the same4gates, and report the full workload separately. No
candidate design or acceptance threshold is changed by this added measurement.

A complementary baseline measures complete API cost, synchronized mesh creation,
and synchronized transport launch, two legs after warmup at the same128photons.
Full five-output hashes must repeat. This distinguishes setup from traversal
before choosing between persistent geometry and a native transport port.

## Complete workload and API baseline

| Scope | Count | Leg0 | Leg1 |
| --- | --- | --- | --- |
| Full queries | 187882,2985hits | 2383.75063045 ms | 2367.37153455 ms |
| Complete API | 128photons | 28454.687156 ms | 28446.377466 ms |
| Transport kernel | 187882events | 27376.717381 ms | 27374.963972 ms |
| Mesh creation | 2754304triangles | 12.961541 ms | 17.646344 ms |

All measurement identities pass, caps0. Query limits exceed8in670cases,
32in546cases,128in237cases; maximum650.71213022. The prefix was not representative.
Mesh caching cannot address the dominant kernel cost on this fixture.

## Candidate design and frozen gates

Create a sibling transport with segmented AABB candidate enumeration only.
Segments have length8; each retains the original global ray origin and the exact
original double-precision triangle arithmetic and face-index tie break. Short
queries use the original enumeration. Stop at the first interval with a hit;
include both interval endpoints and conservative existing AABB padding. At most
128intervals, then fall back to the original whole query, so long/missing rays
cannot introduce an unbounded loop. No physical offset, property or RNG changes.

Before execution: require every captured face identical, two full face arrays
identical, no input mutation, and both query timing legs faster than the original
complete workload. Then require all five complete transport arrays identical to
the frozen source in two independent runs, exact energy/leak/cap counts, and both
complete API timing legs faster. Preserve any failed gate; do not loosen equality.
This is an alternative to a large OptiX port, not SER or RT acceleration.

After the query candidate passes, a separately declared H100 run uses1024photons
and independent seed20260913, two full transports per variant after warmup, with
the same complete-array/energy/leak/cap/input and strict speed gates. Keep the
original128photon/L4 result separately; no throughput extrapolation or automatic
cross-seed array comparison. Ten analytic boundary/tie/fallback cases run twice
on CPU and GPU before either full-transport worker.

## Full transport result

| Worker and fixed sample | Original API, two legs | Segmented API, two legs | Own gates |
| --- | --- | --- | --- |
| L4,128photons,seed20260912 | 28387.162689 / 28403.689517 ms | 1781.297489 / 1439.786640 ms | 5/5 |
| H100,1024photons,seed20260913 | 30625.705040 / 30600.244214 ms | 892.559194 / 888.359482 ms | 5/5 |

All five complete output arrays are identical to the original within each sample
and between repeats: absorbed field, terminal ledger, counters, paths and exits.
The L4sample has187882events; the independent H100sample has1410250events.
Both have zero leaks/caps and exact integer energy balance. Ten analytic
boundary/tie/fallback cases repeat on CPU and each GPU with exact arrays.
The query-only candidate passes4/4:187882face indices identical,1.620288 and
1.664016ms versus2352.589950 and2395.178370ms in its paired worker.

L4 used the same numerical candidate before its CLI was redirected to its own
transport probe. Reconstructing that CLI-only difference matches the recorded
source hash; all imports/functions are AST-identical. The independent H100 run
uses the final complete source hash recorded in its report.

These are synchronized API wall times after warmup; geometry construction and
output copies are included. They are not steady-state kernel throughput, a
million-photon result, a same-input cross-GPU certificate, or clinical acceptance.
The known atlas topology and missing calibration limitations remain. No SER or
OptiX implementation was added: smaller candidate boxes removed the measured
bottleneck without changing the transport physics.

## Reproduction

Use the pinned public `colin27_v3.mat` atlas with SHA256
`8f6c71ced645db796de1736ee77eb7514cdf6e7c6315407ecb2bebf81ee01672`,
from https://github.com/fangq/mcx/blob/master/mcxlab/examples/colin27_v3.mat.
Read `colin27` using `scipy.io.loadmat`, pass its labels to
`bodytwin.geometry.label_interfaces_v1.label_interfaces`, and save all returned
arrays plus `shape=np.array(labels.shape)` in `atlas_fixture.npz`.

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python probes/optics/photon_segmented_transport_probe.py --fixture atlas_fixture.npz --photons 1024 --seed 20260913 --output transport.json
CUDA_VISIBLE_DEVICES='' PYTHONPATH=src python -m pytest -q tests/test_photon_segmented.py
```

Recorded dependencies: Warp1.13.0, NumPy2.4.1, SciPy1.17.0, Python3.12.
A supplied fixture is read without mutation. The probe writes only the explicit
report; the engine writes no files and returns all five arrays. Large photon
counts and further anatomies still require their own measured physical gates.

The local isolated unit-test command (`pytest --noconftest` on the named test)
passes1test containing all10analytic cases in3.18s. It is not a claim that the
entire repository test suite passed. Downloaded full numerical captures from the
two transport workers total455611648bytes; every array was rechecked locally
against its recorded SHA256 and corresponding baseline/repeat array.
