# CX-KANAT — can N1g's k be derived from ANATOMY per person, with no fitting? (the elegant version of "simple wins")

Background: N1g = k·‖GRF‖ beats everything (0.396 BW). k varies 2.1–2.9 between persons (L1 LOPO, CX-CONTACTDEF, A1390). The parameter-free law (A1391/A1398) tried to be time-resolved and lost. A SINGLE anatomical ratio per person has not been tested.
Read: `results/L1/code/n12_io.py` and `n12_model.py` (per-person geometry, static trials, segment lengths), `results/CX-LAW2`, `results/CX-ANKLETERM`, `results/CX-KNEECENTER/geometry.json`, `results/CX-RQANGLE` (literature arms).

## Tasks (PREREG.md + sha256 BEFORE comparing against the implant force)
1. **Freeze 2–3 candidate formulas for k_person from anatomy/kinematics only** (no implant):
   - e.g. k = 1 + E[|r_GRF,knee|]/r_q + 0.29·E[|r_GRF,ankle|]/0.045, where r_GRF,knee is GRF's median lever arm about the knee in the person's stance phase (from markers+GRF), and r_q is the literature patellar/quadriceps arm scaled by the person's tibial plateau width/segment length;
   - OR a pure anthropometric formula: foot length/shank length/tibial plateau width from the static trial.
2. **Score:** k_formula vs k_measured per person (4 persons; plus the sessions JW1/JW4, DM2/DM6 as separate points = 6 points). Then N1g with k_formula against N1g with L1-LOPO k in RMSE per person on the L1 mask.
   Frozen criterion: |k_formula − k_meas| ≤ 10 % for ≥ 4/6 sessions AND RMSE within 5 % of LOPO-N1g for ≥ 3/4 persons.
3. With 4–6 points this is thin: report it honestly as a pilot, with a leave-one-out on the formula choice.

Deliver RESULTS.md starting with `# CX-KANAT`, results.json, the script, and pytest. Light computation: run locally with 2 threads. Internal data stays local. lane_runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
