# CX-ATTACHAREA — higher resolution where it matters most: muscle attachments as AREAS instead of points → shrink the largest uncertainty in the chain

Why: attachment uncertainty dominates.
- The what-if shape effect (A369/A907: +15° CCD ≈ +5 % hip force) lies INSIDE N2b's attachment band (A363: right hip −8.5/+15.9 % at σ 9.5 mm).
- Earlier: the attachment error at 8–19 mm dominates (H2/AT1).
- N7c/N7b have attachment regions and certificates.

## Build on (requirement: "Builds on" with file paths; no new geometry base)
- `results/N7c/geomgr` (regions, attachment transfer, certificate v2) and `results/N7b` (b/c/e/h; regions survive editing).
- The N2b batched solver and bands (`results/N2b`, `results/BT-N50`).
- The emulator (BT-XF4); CX-SOLVER2; `bodytwin_core/`.
- The graph (grep attachment/insertion/footprint in `~/projects/bodytwin/data/MECHANISM_ANCHOR_GRAPH.json`) and `~/projects/bodytwin/scripts/msk` (grep attachment, insertion, footprint).
- Literature with DOI on attachment areas (footprints) for gluteus medius, adductors and hamstrings, if accessible.

## Tasks (PREREG.md + sha256 first)
1. Represent the attachments of the hip-crossing muscle groups as areas (a footprint on the bone surface, as a distribution) instead of points. The line of action = the effective force point integrated over the area, or several sub-elements.
2. Propagate through the chain on the collaborator's box lift: does the right-hip band shrink compared with the point model plus the σ 9.5 mm point error? How much of the band is area geometry, and how much is the attachment's location on the bone?
3. Does the what-if effect (CCD/AV) now become larger than the band, i.e. usable for decisions?
4. Counter-tests: a random area of the same size (placebo); the point model at the area's centroid.
5. Use N2b batched (GPU/CPU) for the ensemble. Internal the collaborator data runs locally or on OVH (`tasks/cloud_run.sh`, BodyTwin ≤ 12 vCPU, finish by 23:59).

Write in `results/CX-ATTACHAREA/`. `RESULTS.md` starting with `# CX-ATTACHAREA`, plus results.json and pytest.
