# CX-IMUFORCE — Physical sparse sensors to joint force with uncertainty

Builds on `results/BT-B160` (three virtual sensors, 312,28 N RMSE and ground truth as a model proxy), `results/BT-IM1`, `results/CX-KNEENET` (N1g superior; correlated edges), and OpenCap/Grand Challenge sensor metadata. This is a test of observability and sensor transfer, not a rerun of the virtual 6-person analysis.

## Task
1. Inventory physically measured IMU channels, calibration, synchronisation, GRF and independent joint-force ground truth per person. Document separately which body sites and activities actually exist.
2. Build a sensor-noise/placement model and a nested LOPO joint-force model with intervals, sensor budget 1/2/3. Only measured channels may count as a physical sensor gain.
3. If invasive ground truth and real IMU never coexist in the inputs: perform a data-identifiability and domain-shift study, with clear UNKNOWN for physical joint-force performance and a concrete acquisition contract.

## PREREG and gate
Primary gate on compatible ground truth: person-median RMSE at least 20% lower than both N1g and B24 on ≥3/4 people; calibrated 90% interval coverage 85–95% per person, median width ≤2 BW. Control: time shift, incorrectly mounted IMU, body weight/GRF-only, sensor-label permutation. Also report the previous 300 N gate exactly.

## Shared rules
- Write only under `results/CX-IMUFORCE/`, and large intermediate results in `external_media`. Read `tasks/NIGHT_PREAMBLE.md`, and relevant A rows in `notes/RESULTS_INDEX.md`. `~/projects/bodytwin` and other sessions' workspaces are read-only.
- Before the first calculation: `PREREG.md` and `PREREG.sha256` with hypothesis, data/selection, numerical gate, strongest baseline, counter-test and error definition. Document `Builds on` with graph node and source files and `Not redone`. Preserve negative results and `UNKNOWN`.
- Locally only reading/short tests with at most two threads. Numerics >60 s or >1 GB through `tasks/cloud_run.sh` on OVH, within the BodyTwin quota 12 vCPU and complete before 22:45. Internal restricted model data/the collaborator/LHDL data must never go to Modal. No email, push or publishing.
- Deliver `results/CX-IMUFORCE/RESULTS.md` with the first line `# CX-IMUFORCE`, `results.json`, executable code, provenance/hashes and meaningful checks. Report both the number of valid and excluded units.
