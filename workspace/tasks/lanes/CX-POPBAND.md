# CX-POPBAND — Calibrated population product for load bands

Builds on `results/L3` (1000 synthetic individuals and 23 geometries, but left hip band coverage 0,816), `results/N2b` (batched p=2/min-max), `results/CX-STRENGTH-POP`, and `results/CX-KNEENET` (insufficient coverage of the LOPO bands). Do not rebuild L3's pure scaling benchmark.

## Task
1. Build a population API/batch with explicit shared and individual latents: geometry, body weight, strength, attachments, solver/model uncertainty. Return per-person point forces, 90/95% bands, validity rate and source provenance for hip, knee, L5/S1 and GH.
2. Calibrate bands on a held-out person, not just a synthetic draw. Test the collaborator's box lifting and Grand Challenge walking separately; refrain from transferring bands between regimes without a test.
3. Profile 100/1000 people with identical numerical tolerance and report per-individual cost; run the heavy batch on OVH.

## PREREG and gate
KKT ≤1e-10 per accepted step; 90% band coverage 85–95% per person and 95% bands ≥90% on every evaluated joint; median width at most 1,25× unconditional band. At least 100 steg/s CPU batch; report if the ceiling falls. Countertests: permuted person geometry, removed correlation, B24/N1g and a conformal null band of the same width.

## Shared rules
- Write only under `results/CX-POPBAND/`, and large intermediate results in `external_media`. Read `tasks/NIGHT_PREAMBLE.md`, and relevant A rows in `notes/RESULTS_INDEX.md`. `~/projects/bodytwin` and other sessions' workspaces are read-only.
- Before the first computation: `PREREG.md` and `PREREG.sha256` with hypothesis, data/selection, numerical gate, strongest baseline, countertest and error definition. Document `Builds on` with graph node and source files and `Not redone`. Preserve negative results and `UNKNOWN`.
- Locally only reading/short tests with at most two threads. Numerics >60 s or >1 GB via `tasks/cloud_run.sh` on OVH, within the BodyTwin quota of 12 vCPU and finished before 22:45. Internal restricted model data/the collaborator/LHDL data must never go to Modal. No email, push or publishing.
- Deliver `results/CX-POPBAND/RESULTS.md` with first line `# CX-POPBAND`, `results.json`, runnable code, provenance/hashes and meaningful checks. Report the counts of both valid and dropped units.
