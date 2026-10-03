# CX-ANKLETERM — test the coordinator's frozen prediction points 2–3: the parameter-free knee law WITH an ankle/gastrocnemius term, on all Grand Challenge gait data with implant force

## Frozen prediction
`results/CX-BTHYPO/COORDINATOR_PREDICTION.md` (sha256 in `.sha256`, commit d137f81 — do NOT change it).
- Point 2: the knee-moment-only law F≈|GRF|+|M_flex|/r_q underestimates the SECOND peak by 0.5–1.0 BW (gastrocnemius).
- Point 3: with the ankle term F≈|GRF|+|M_knee|/r_q+β|M_ankle|/r_achilles, it lands within ~10 % of N1g without fitting.

## Data (local; F-8's jw_lungef1 knee_forces/eTibia is FORBIDDEN)
- The BT-BH-001…022 packets' inputs (now real files: markers, GRF, implant force, id_prep, n12b_scored per trial) for JW, DM, SC and PS across the gait variants.
- The net knee and ankle moments: id_prep (L1 prep b) or your own ID from markers+GRF (CX-LUNGEID's pipeline).
- r_q(θ): literature with DOI (patellar tendon arm vs angle) + a band. r_achilles: literature (~4–5 cm) + a band. β: the gastrocnemius share of the plantarflexor moment from PCSA ratios in the literature (NOT fitted to the implant force).

## Tasks (PREREG.md + sha256 first)
1. Per person × gait variant:
   - N1p_knee (without the ankle term) and N1p_ankle (with it), against the measured implant force: total RMSE (BW), the first and second peak errors, and correlation;
   - N1g (LOPO) and B24 on the same frames.
2. Score the frozen prediction point by point: point 1 (first peak within 15 %), point 2 (the second peak underestimated by 0.5–1.0 BW without the ankle), point 3 (with the ankle within 10 % of N1g), and point 5 if other activities exist. Report the numbers against the frozen intervals, with no retrospective adjustment.
3. Residual decomposition (the branch rule, `results/CX-BTHYPO/FRAMING_NOTE.md`): does the residual follow knee angle, EMG, activity or person, and which node expands.
4. Sensitivity: r_q/r_achilles/β over the literature bands. Counter-test: a time-shifted/permuted M_ankle.

Write in `results/CX-ANKLETERM/`. `RESULTS.md` starting with `# CX-ANKLETERM`, plus results.json and pytest.
