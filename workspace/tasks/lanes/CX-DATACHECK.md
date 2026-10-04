# CX-DATACHECK — settle two flagged findings and audit the The swarm measurements that will be used further (model hierarchy: The swarm < Sol)

Source: `notes/RESULTS_INDEX.md` A914 and `results/BT-PG-*`.

## 1. PG-027: the Biodex moment vs the implant compression, r ≈ 0 (`results/BT-PG-027`)
A knee extension moment should give tibiofemoral compression.
- Check the time sync between the Biodex file and eTibia/eKnee in the Grand Challenge archives: the sample rate, the start time, and the README's statement about synchronisation. Cross-correlate over ±2 s.
- Check the units and sign of each channel.
- Check that the same trial/session is being compared (file names; gc1 vs gc4).
- Verdict: a sync/data error (corrected r) or a real non-correlation (with an explanation).
- NOTE: the eTibia data for the JW4 lunge is Field's sealed F-8 facit. You may ONLY use the isokinetic trials, not jw_lungef1.

## 2. PG-037: mid-thigh CSA / volume-normalised mean CSA = 2.14 (`results/BT-PG-037`, source CX-MUSCLE-CT2)
Check whether the volume refers to 25–75 % of the femur (mid-thigh) while the length is the whole femur. That would give an artefactual ratio of ~2. Recompute with consistent definitions. If it is a definition error, list which PG-031…040 results are affected and correct their conclusions in `results/CX-DATACHECK/CORRECTIONS.md` (the coordinator writes inline corrections in the register).

## 3. Independent audit of the The swarm measurements that will be used further
Recompute the key number with your OWN code from raw data, per packet. Verdict: HOLDS / DEVIATES (number) / A MATTER OF DEFINITION.
- PG-011–020: fluoro rollback, hysteresis, ML/AP, coupling; the rollback's lack of benefit for patellar compression.
- PG-021–030: the isokinetic curves, the passive correction, cycle variation.
- PG-031–040: CT muscle geometry.

Write in `results/CX-DATACHECK/`. `RESULTS.md` starting with `# CX-DATACHECK`, plus results.json. lane_runner has full permissions in the workspace; the zips are read-only.
