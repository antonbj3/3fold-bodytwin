# CX-MUSCLE-CT — T2: does the individual's body shape predict their muscle size better than the mass laws? (VSD CT, 30 people)

Reference: muscle tissue in the thigh (and calf if the CT covers it) per person, from VSD CT. Segment muscle and fat by HU (muscle ~10–100 HU, fat −190…−30 HU). Measure the muscle volume and the CSA at 25/50/75 % of femur length, with the bone masked out. Look for existing thigh work: H5b/BT-N33/BT-B79 in results/ and `tasks/index/DATASETS.json` ("VSD postmortem-CT for H5b"). Reuse their paths, masks and hull_gated. Height/weight per person from VSD metadata.
Models (LOSO):
- M-a the LengthMass baseline (muscle size ∝ mass);
- M-b LengthMassFat (BMI fat proxy);
- M-c Handsfield regression height×mass;
- M-d skin shape + layer model (N33 method: segment volume from the skin surface minus bone minus a subcutaneous fat layer);
- M-e M-d + fat thickness from a skin-fold proxy.
Report the error in muscle volume/CSA in %, per person and in median. Counter-test: permuted persons. Criterion: M-d or M-e beats M-a and M-b by ≥ 20 % median absolute error. Also report whether postmortem CT biases the HU thresholds (state it; do not correct it silently).

## Common
- PREREG.md + PREREG.sha256 BEFORE the first run: facit, criteria with numbers, counter-tests (permuted inputs), and what counts as a fall. Symmetric skepticism.
- The lane runner sandbox has network access (public data may be downloaded; save it on external_media and record the URL + sha256). Locally: nice, 2 threads, ≤ 60 s per test; heavier work via Modal (`tasks/modal_run.py`) for public data, or OVH (`tasks/cloud_run.sh`, BodyTwin ≤ 8 vCPU) for internal data.
- Write only in `results/CX-MUSCLE-CT/`. `~/projects/bodytwin` is read-only. `RESULTS.md` starting with `# CX-MUSCLE-CT`, plus results.json and code with pytest.
