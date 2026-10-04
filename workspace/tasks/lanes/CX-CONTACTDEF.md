# CX-CONTACTDEF — check the definitions in P-CONTACT-COMP and P-GRF-PEAK: why is k = implant peak / GRF peak ≈ 0.50 when ~2–3 is expected?

Source: `results/BT-DA-001-{a,b,c}` (pooled k = 0.5046 over 61 trials; GC3/SC, GC4/JW, GC5/PS, GC6/DM), `results/DATAMATRIX_TABLES/P-CONTACT-COMP.csv`, `P-GRF-PEAK.csv`, `tasks/datamatrix_compute.py`, and the raw files under `results/BT-DM-*/inputs/raw.csv` (eTibia: Time,Fx,Fy,Fz,Tx,Ty,Tz,…).
Expected: tibiofemoral contact force is ~2–3 × body weight in gait, and the GRF peak is ~1.1–1.2 BW, so k ≈ 2–3. k ≈ 0.5 indicates a definition error.

Check:
1. Which eTibia channel/combination is contact_peak_N? Is it the total axial force (Fz) in N, with the correct sign and unit, possibly per load cell or scaled? Compare with the Grand Challenge README/Data Description on eTibia/eKnee channels (units, the axial axis, the sign convention).
2. Is peak_resultant_N the GRF from ONE plate (the implanted leg) or the sum of several? In the calibration trials (1legstand, squats), which plate is the implanted leg?
3. Are both taken from the same trial and time window?
4. Correct the definition in `tasks/datamatrix_compute.py` if it is wrong. Recompute P-CONTACT-COMP (and P-GRF-* if affected), and recompute DA-001's pooled k as a check (expected 2–3 in gait).
5. List which BT-DA-* results are affected, and how their conclusions change, in `results/CX-CONTACTDEF/AFFECTED.md` (the coordinator writes the inline corrections).
6. EXCEPTION: jw_lungef1 (F-8) must not be read.

Write in `results/CX-CONTACTDEF/`. `RESULTS.md` starting with `# CX-CONTACTDEF`.
