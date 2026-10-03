# CX-SPINEARM — Individual L5/S1 moment arm and back load in the collaborator's lift

Builds on `results/N2b` (full box lift, 141 step, GH/L5/S1-band), `results/BT-B63` and A188 (the reference model-VBR comparison regime problem), `results/CX-PELVIS` as well as existing `validate_spine_force.py`/`model_disc_pressure.py` in the source project. Do not attempt to reinterpret OrthoLoad activity as the same lift without matching.

## Task
1. Inventory available individual measurements for L5/S1 joint center, pelvic tilt, erector/multifidus attachments and extensor arm; build an explicit measurement and coordinate chain for the collaborator and 3 VSD shapes if identifiable.
2. Put moment arm/pelvis in N2b/N40b system and calculate the same 141 step with unchanged load and p=2. Distribute contributions from moment arm, attachment, strength and pelvic inclination.
3. Compare with load matched VBR-reference if available; otherwise mark external force validation UNKNOWN and report only internal power plus minimum required VBR-data.

## PREREG and gate
141/141 KKT ≤1e-10, zero perturbation ≤0,1 N; explained L5/S1-variation ≥20% and stable sign in ≥3 forms. External criterion: ≥20% lower median absolute percent error than weight×OrthoLoad median of ≥3 activity-matched subjects. countertest: mirrored and person-swapped torque arm as well as B24 and N1g where GRF is present. The uncertainty band must accompany each force value.

## Gemensamma regler
- Skriv endast under `results/CX-SPINEARM/`, samt stora mellanresultat i `/media/anton/sdc1-tmp/bodytwin/CX-SPINEARM/`Read it. `tasks/NIGHT_PREAMBLE.md`, and relevant A-rows in `notes/RESULTS_INDEX.md`. `~/projects/bodytwin` and other sessions' workspaces are read-only.
- Before the first computation: `PREREG.md` and `PREREG.sha256` with hypothesis, data/selection, numerical gate, strongest baseline, countertests and error definition. Document `Builds on` with graph node and source files plus `Not redone`. Preserve negative results and `UNKNOWN`.
- Deliver `results/CX-SPINEARM/RESULTS.md` with first line `# CX-SPINEARM`, `results.json`, executable code, provenance/hashes and meaningful checks. Report both number of valid and missing units.
