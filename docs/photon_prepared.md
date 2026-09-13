# Prepared photon geometry: mechanism measurement

| Measured stage | L4 / CPU first–second | H100 first–second | Decision |
| --- | --- | --- | --- |
| Photon float normals | 392.683 / 391.730 ms | 244.463 / 242.947 ms | Prepare immutable geometry once |
| Photon double planes | 528.496 / 526.141 ms | 292.714 / 288.400 ms | Same preparation target |
| Photon transport kernel | 382.569 / 382.667 ms | 210.689 / 210.700 ms | Preserve arithmetic; profile next |
| Photon input upload excluding kernel | 44.942 / 45.054 ms | 15.624 / 15.638 ms | Keep prepared geometry on device |
| NIRS attenuation and reduction | 187.467 / 187.829 ms flat;164.758 / 163.199 ms curved | CPU observation only | Next after photon preparation |

Full photon API: L4 1418.419 / 1411.005 ms; H100 803.831 / 790.626 ms.
1024 photons,2754304 faces,57012384 complete output bytes; observed/unobserved and
both repeats match exactly on each device. All five hashes also agree across these
two device runs. Full histories, energy exact,zero leaks/caps. Prefix results were
not substituted for full histories. NIRS1320 signal calls per leg; complete leg
records/hash dictionaries unchanged by instrumentation and exact twice.

The next implementation is a separate prepared photon scene. Baseline files stay
unchanged. Warm reuse and cold construction+first-call costs must be reported separately.
No new verifier suite, README change, or numerical tolerance change.


Before candidate execution, require full five-array bit parity with the segmented
baseline and two exact repeats, exact packet-energy conservation, zero leaks/caps,
unchanged caller inputs, and both warm candidate times below both baseline times.
Record preparation and cold total separately; no cold-speedup claim is required.
Retained outputs must survive subsequent calls. The prepared scene owns uploaded
geometry and optical properties. Geometry mutation requires a new scene. Launch
source, direction, photon count and seed remain per-call arguments.

## Complete prepared-scene results

| Device | Baseline API, two runs | Reused scene, two runs | Preparation, two runs | Preparation plus first call |
| --- | --- | --- | --- | --- |
| L4 | 1434.129 ms / 1444.608 ms | 410.189 ms / 410.189 ms | 996.505 ms / 1007.567 ms | 1407.530 ms / 1418.702 ms |
| H100 | 941.978 ms / 958.444 ms | 229.856 ms / 229.957 ms | 695.240 ms / 808.035 ms | 925.207 ms / 1038.225 ms |

All five original full-transport gates pass on both devices. The complete five
arrays (57012384 bytes per result) match baseline and repeats exactly, including
1410250 events and energy1099511627776 with zero leaks/caps. All five hashes also
match across the two devices for this fixed1024-photon seed. Ten query boundary
controls pass on CPU and each GPU. Prepared ownership controls pass on CPU and each
GPU: mutate caller geometry/properties after preparation, retain earlier output,
change the next launch count/seed, and reject five invalid launch configurations.
The scene is reusable but not a mutable-geometry update service.

Warm gains are about3.50x on L4 and4.10–4.17x on H100. The second H100 cold total
1038.225ms exceeds the baseline941.978–958.444ms: **no general cold-speedup claim**.
No million-photon or clinical/anatomy claim follows from this fixed sample.

Use the existing transport probe with `--prepared --fixture <fixture.npz>`;
`PreparedPhotonScene(...).simulate(...)` takes the original per-launch arguments
except geometry, properties and device, which are supplied at construction.
The original segmented module and kernel are unchanged.
