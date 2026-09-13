# Certified open-cell query rejection

Measured before design:

| Item | L4 | H100 |
| --- | ---: | ---: |
| Window2 kernel1024, two legs ms | 240.504 / 240.492 | 156.941 / 156.702 |
| Window2 kernel8192, two legs ms | 286.043 / 285.992 | 158.403 / 158.407 |
| Remaining host/readback1024 ms | 26.63 / 29.40 | 23.15 / 23.03 |

Complete187882-query inventory has19452 triangle visits; median and90th
percentile are zero, maximum48. Precomputing all triangle edges would allocate
132206592additional bytes for little repeated intersection arithmetic, so that
candidate is not selected.160389queries have their entire conservatively padded
box strictly inside one integer cell; all160389 have expected face -1.

Design: a separate prepared scene accepts only meshes where every triangle lies
on at least one integer-coordinate axis plane, checked after the same float32
vertex conversion as the baseline mesh. Its query function first forms the same
padded float32 box as the baseline. If that entire box lies strictly inside one
unit cell in every axis, no certified triangle can intersect it. Return no-hit;
otherwise use the unchanged window2 query. No labels or approximate distance
model, no relaxed boundary comparison, and no default routing change.

Gates fixed before execution: all five original transport gates at1024/8192 on
L4/H100, full five-array parity against window2, exact repeats, zero leaks/caps,
exact energy and unchanged inputs; both candidate times below both baseline times
in every case. Existing10 boundary/tie/fallback controls must pass CPU/GPU.
Noninteger-plane geometry must be refused before simulation. Construction cost
and geometry-certificate overhead are outside warm timing and must be reported.
Evidence before design: reports/photon_window2_phases_* and
reports/photon_candidate_inventory.json with full repeated inventory arrays.

Fresh result: all five gates pass in all four cases. All80 complete transport
arrays919996928bytes were downloaded and checked against their recorded hashes.
The10 existing boundary/tie/fallback cases pass twice on CPU and each GPU with
exact arrays; the explicit noninteger-plane refusal passes before device setup.
The transport kernel AST remains identical to window2.

| Device | Photons | Warm baseline ms | Warm sibling ms | Construction baseline / sibling ms |
| --- | ---: | --- | --- | --- |
| L4 | 1024 | 251.698 / 250.170 | 183.527 / 185.334 | 692.011 / 807.995 |
| L4 | 8192 | 298.359 / 299.624 | 239.262 / 241.491 | 632.016 / 873.613 |
| H100 | 1024 | 174.572 / 174.430 | 127.669 / 117.086 | 708.349 / 772.686 |
| H100 | 8192 | 191.743 / 176.021 | 121.820 / 115.402 | 634.976 / 808.356 |

Construction is slower in all four measured cases: added115.985/241.597ms on
L4 and64.337/173.380ms on H100. No cold-call speedup is claimed. Geometry
validation uses transient host arrays and does not add persistent device arrays.
Exact packet energy1099511627776/8796093022208 and zero leaks/caps hold.
Receipts: reports/photon_cell_* JSON. The geometry restriction is mandatory;
this sibling is not a replacement for arbitrary mesh transport. No large-count
or clinical claim. Both numerical baselines remain unchanged.

API: PreparedPhotonScene(vertices,faces,properties,face_front,face_back,device=...)
snapshots the validated float32 geometry through the existing prepared owner.
Each simulate call returns independent full absorption, terminal, counters, paths
and exits arrays. No implicit artifact writes. The module CLI routes to the
existing transport probe with --cell-certified and explicit fixture/output paths.
