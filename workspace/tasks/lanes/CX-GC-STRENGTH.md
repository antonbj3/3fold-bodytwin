# CX-GC-STRENGTH — T3: does knee force in gait become more correct if the model's Fmax is calibrated against the person's measured strength? (Grand Challenge)

Facit: Grand Challenge competition zips `external_media*Competition-latest.zip` (gc1–gc6). L1 found isometric hip/knee strength trials (gc4–gc6), MaxEMG (gc1–4, gc6) and Biodex strength data (not extracted). Extract only what is needed to `external_media`. Map gc number → person (JW/DM/SC/PS) via the zips' own metadata and results/L1.
Chain:
1. Measured maximum knee extension/flexion moment (isometric/Biodex, the correct joint angle) per person.
2. The model's predicted maximum moment at the same angle from TLEM Fmax × moment arm (results/L1/code/tlem_arch.json, L1's geometry) under M0 (generic) and M1 (mass).
3. Calibrate the muscle-group Fmax (quadriceps, hamstrings, gastrocnemius) so that the predicted max moment = measured.
4. Rerun L1's V0/S1b chain for gait with calibrated Fmax. Does knee force RMSE fall against implant force, and against N1g?
Criterion: calibrated Fmax lowers V0's median RMSE by ≥ 10 % in the persons that have strength data, and the counter-test (another person's strength) does not give the same gain. Also report whether calibrated Fmax gets closer to N1g. Report honestly if the strength data do not cover a person.

## Common
- PREREG.md + PREREG.sha256 BEFORE the first run: facit, criteria with numbers, counter-tests (permuted inputs), and what counts as a fall. Symmetric skepticism.
- The lane runner sandbox has network access (public data may be downloaded; save it on external_media and record the URL + sha256). Locally: nice, 2 threads, ≤ 60 s per test; heavier work via Modal (`tasks/modal_run.py`) for public data, or OVH (`tasks/cloud_run.sh`, BodyTwin ≤ 8 vCPU) for internal data.
- Write only in `results/CX-GC-STRENGTH/`. `~/projects/bodytwin` is read-only. `RESULTS.md` starting with `# CX-GC-STRENGTH`, plus results.json and code with pytest.
