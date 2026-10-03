# Reproduktion R2

Threads2 is set in the scripts. Choose **new** output/file names for reruns; existing destination files are rejected. PREREG/design must be verified before a new run. Large field data is on games-240, not shared_data.

```bash
python results/LANE_SURGICAL_SYNTHESIS/diagnose_r2.py --out results/LANE_SURGICAL_SYNTHESIS/r2/attribution_REPLAY
python results/LANE_SURGICAL_SYNTHESIS/tests_r2.py --out results/LANE_SURGICAL_SYNTHESIS/r2/tests_REPLAY
```

Ensemble: `plan_measurements_r2.py batch --start 0 --stop 32 --out <new lane directory> --store <new games-240 lane directory>` runs original sample IDs to exclusive new destinations. Reuse the frozen R2 design and prior. Run batch intervals0:32,32:96,96:160,160:224,224:256 and then summary and refinement with the same new --out/--store. Do not mix priors or results between rounds.

`report_r2.py` requires a complete ensemble, refinement and OK log and writes exclusive report files. OriginalR1 code is verified against CODE_MANIFEST_R1.json. Result binding uses the lane's local graph, kindreview; canonicalimport is performed by the coordinator.

Nonlinear information test: nonlinear_information_r2.py prior → fit → build --start 0 --stop 32 → build --start 32 --stop 64 → evaluate. V2 uses independent training/measurement noise. Its fixed exclusive destinations require a new lane-local replay copy with new OUT/STORE for reruns; existing data must not be overwritten.
