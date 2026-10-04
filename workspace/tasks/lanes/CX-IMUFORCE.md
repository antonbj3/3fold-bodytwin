# CX-IMUFORCE — Physical sparse sensors to joint force with uncertainty

Builds on `results/BT-B160` (three virtual sensors, 312,28 N RMSE and ground truth as a model proxy), `results/BT-IM1`, `results/CX-KNEENET` (N1g superior; correlated edges), and OpenCap/Grand Challenge sensor metadata. This is a test of observability and sensor transfer, not a rerun of the virtual 6-person analysis.

## Task
1. Inventory physically measured IMU channels, calibration, synchronisation, GRF and independent joint-force ground truth per person. Document separately which body sites and activities actually exist.
2. Build a sensor-noise/placement model and a nested LOPO joint-force model with intervals, sensor budget 1/2/3. Only measured channels may count as a physical sensor gain.
3. If invasive ground truth and real IMU never coexist in the inputs: perform a data-identifiability and domain-shift study, with clear UNKNOWN for physical joint-force performance and a concrete acquisition contract.

## PREREG and gate
Primary gate on compatible ground truth: person-median RMSE at least 20% lower than both N1g and B24 on ≥3/4 people; calibrated 90% interval coverage 85–95% per person, median width ≤2 BW. Control: time shift, incorrectly mounted IMU, body weight/GRF-only, sensor-label permutation. Also report the previous 300 N gate exactly.

## Gemensamma regler
- Skriv endast under `results/CX-IMUFORCE/`, samt stora mellanresultat i `external_media`Read it. `tasks/NIGHT_PREAMBLE.md`, and relevant A-rows in `notes/RESULTS_INDEX.md`. `~/projects/bodytwin` and other sessions' workspaces are read-only.
- Before the first computation: `PREREG.md` and `PREREG.sha256` with hypothesis, data/selection, numerical gate, strongest baseline, countertests and error definition. Document `Builds on` with graph node and source files plus `Not redone`. Preserve negative results and `UNKNOWN`.
- Deliver `results/CX-IMUFORCE/RESULTS.md` with the first line `# CX-IMUFORCE`, `results.json`, executable code, provenance/hashes and meaningful checks. Report both the number of valid and excluded units.
