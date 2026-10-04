# CX-FLUOROLINK — Fluoroscopy as a verifiable CT↔marker link

Builds on `results/CX-CT2MARKER` (A368: only two femur points, 12,8 mm best JW residual, zero verified frames), `results/CX-JWGEOM`, and Grand Challenge JW4/PS5's inventoried synchronized fluoroscopy/marker files. No new arbitrary CT↔marker rotation may be called observed.

## Task
1. Reuse the archive catalog; identify actual numerical time series, calibration matrices, time base and implant/CT coordinate definition. Distinguish relative femur–tibia kinematics from absolute lab transform.
2. Build an identifiability analysis: which 6-DOF can be determined, with rank/Jacobian and measurement uncertainty? Choose at least one JW or PS sequence with known phase pairs. If absolute CT↔lab cannot be identified, deliver an exact minimal extra measurement data contract and abstain downstream.
3. If identifiable: robust synchronization/rigid fitting, hold out 20% frames, build transform for femur and tibia with uncertainty bands. Then run the same mask in L1/C1 and compare knee contact force against N1g and B24.

## PREREG and gate
Transform: median held landmark/residual ≤10 mm and 95th percentile ≤20 mm; at least 3 independent anatomical correspondences per segment, full rank. Force gain: ≥10% lower person-median RMSE than both N1g and B24 on ≥2 people; otherwise negative/UNKNOWN. Countertest: time shift ≥0,2 s, person swap, permuted correspondences. Require all three to lose the registration gate.

## Shared rules
- Write only under `results/CX-FLUOROLINK/`, plus large intermediate results in `external_media`. Read `tasks/NIGHT_PREAMBLE.md`, and relevant A rows in `notes/RESULTS_INDEX.md`. `~/projects/bodytwin` and other sessions' workspaces are read-only.
- Before the first computation: `PREREG.md` and `PREREG.sha256` with hypothesis, data/selection, numerical gate, strongest baseline, countertests and error definition. Document `Builds on` with graph node and source files plus `Not redone`. Preserve negative results and `UNKNOWN`.
- Locally only reading/short probes with at most two threads. Numerics >60 s or >1 GB via `tasks/cloud_run.sh` on OVH, within the BodyTwin quota of 12 vCPU and complete before 22:45. Internal restricted model data/the collaborator/LHDL data must never go to Modal. No email, push or publication.
- Deliver `results/CX-FLUOROLINK/RESULTS.md` with the first line `# CX-FLUOROLINK`, `results.json`, runnable code, provenance/hashes and meaningful checks. Report counts of both valid and excluded units.
