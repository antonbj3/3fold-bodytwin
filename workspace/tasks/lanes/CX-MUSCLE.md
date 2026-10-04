# CX-MUSCLE — the individual's muscle parameters from their own body (shape, volume, composition) instead of body mass


## Build on (read FIRST — much already exists in the original repo, read-only)
- `~/projects/bodytwin/docs/`:
  - MECHANISM_FMAX_PCSA_VALIDATION (knee prime movers over-strong 12/12, correction 1.515→1.255× OrthoLoad, hip unaffected)
  - FMAX_FIX_CROSS_SUBJECT, FMAX_CORRECTION_WAVEFORM
  - FIBER_LENGTH_CORRECTION, FIBER_LENGTH_CROSS_SUBJECT
  - TENDON_SLACK_SENSITIVITY, MUSCLE_AUDIT
  - BODY_COMPOSITION_FROM_VIDEO, GRIP_STRENGTH, MUSCLE_HYPERTROPHY, JOINT_FORCE_SCORECARD, CROSS_SUBJECT
- `~/projects/bodytwin/scripts/msk/subject_specific_scaling.py`, `body_composition_from_video.py`, and grep in scripts/msk for pcsa|fmax|fiber|slack|handsfield|ward.
- Graph `~/projects/bodytwin/data/MECHANISM_ANCHOR_GRAPH.json`: MSK-MUSCLE-ARCHITECTURE, MODEL-MUSCLE-ARCHITECTURE-PCSA, HOLE-MUSCLE-SPECIFIC-TENSION-CROSSMODEL-DISAGREEMENT-CONTRADICTION, MSK-TENDON-SPRING-SAFETY-FACTOR, AUTO-BODY-TWIN-MUSCLE-TENDON-MISMATCH-CELL-1 and neighbours. Read claim/notes, including corrections.
- The workspace:
  - `results/BT-B5`, `BT-B26`, `BT-B79`, `BT-N33`, `BT-N57` (A299), `BT-N55` (A300 EMG);
  - `results/CX-PARAM/PARAM_MAP.md` (P4/P5);
  - `results/L1` (coordinator lane, calibrated recruitment — do NOT redo it; read its PREREG and coordinate: you own the muscle PARAMETERS, L1 owns the recruitment CRITERION).
- Write "Builds on" and "Not redone" in PREREG.

## Question
Can the individual's muscle parameters be derived from the individual's own body so that joint force and muscle activity become more correct than with the mass-scaling baseline, measured against in vivo data? The parameters are Fmax/PCSA, optimal fibre length, and tendon slack length. The sources are muscle volume from segment shape/composition (as in N33) and bone geometry from geomgr, with population regressions (Handsfield, Ward) as priors.

## Tasks
1. PREREG.md + PREREG.sha256 first. Include:
   - facit: Grand Challenge knee force (JW/DM/SC/PS), OrthoLoad hip, OpenCap EMG, and any strength data you find (grep for dynamometer/MVC/isokinetic/grip in data maps: `tasks/index/DATASETS.json`, `~/projects/bodytwin/docs/DATASET_MAP.md`);
   - null models B24, N1/N1g (preamble §2) and the LengthMass baseline and LengthMassFat.
2. A parameter chain M0–M4, same everything else, LOSO/LOPO:
   - M0: generic;
   - M1: the LengthMass baseline;
   - M2: LengthMassFat;
   - M3: volume regression (Handsfield) from height/mass;
   - M4: individual volume from segment shape + layer model (N33 method) → PCSA with Ward fibre length, tendon slack from the geometry.
   Report the joint-force error and EMG correlation per model.
3. Uncertainty: the parameter band per muscle and how much of the joint-force error it explains (variance share). Identify which parameters are identifiable from the data we have, and which are not.
4. Symmetric skepticism: if M4 does not beat M1 and N1g, say so plainly and name the measurement that would be required.

## Resources and rules
- Local: nice, 2 threads, ≤ 60 s per test.
- Heavy runs go on Modal (`modal run tasks/modal_run.py ... --cpu 8 --mem-gb 16`). If the sandbox has no network, write shard scripts + `RUN_ON_CLOUD.md`.
- Write only under `results/CX-MUSCLE/`. `~/projects/bodytwin` is read-only. restricted model data and the collaborator's data stay internal. No emails or pushes.
- `RESULTS.md` starting with `# CX-MUSCLE`: table model → error per facit → vs N1g/B24 → verdict. Plus `results.json` and code with pytest.

## Handover from L1 (read first): `results/L1/FOR_CX_MUSCLE.md`
- TLEM architecture for 166 elements: `results/L1/code/tlem_arch.json`. Slack lengths do not fit straight-line paths (soleus ~2x L0, 27–56/166 outside range).
- STRENGTH FACIT: Grand Challenge has isometric hip/knee strength trials (gc4–gc6), MaxEMG (gc1–4, gc6) and Biodex strength data in the competition zips on external_media (not extracted). Extracted EMG: external_mount Use measured strength as a direct facit for M0–M4 (predicted vs measured max moment per person).
- OpenCap muscle parameters for 6 persons: `results/L1/n55x/*.npz`.
