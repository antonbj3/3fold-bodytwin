# CX-BOOKKEEP4 — book all finished but unbooked results (≈ 100 BT/FX + CX lanes) and summarise what is worth knowing

1. For each results/BT-*, results/FX-*, results/CX-* that has a RESULTS.md but whose id is missing from `notes/RESULTS_INDEX.md` (skip meta lanes: CX-BOOKKEEP*, CX-SWARMGEN*, CX-PLANNER, CX-PACKETFACTORY, CX-DSWAVE, CX-CLOUDLANES — one row in total for these), write ONE row in exactly the same format as the last rows. Numbering starts at A917 (A907–A916 are reserved). Append only (>>). Only numbers from the job's own files; mark "not reviewed"; mark a packaging error/UNKNOWN when the job says so. Mark The swarm-reproduces-The swarm (BT-Q/BT-S) as "reproduction" with HOLDS/DEVIATES/UNKNOWN.
2. CX lanes not yet booked (check with grep in the register) each get a full row.
3. `results/CX-BOOKKEEP4/RESULTS.md` starting with `# CX-BOOKKEEP4`, ≤ 30 lines, for the coordinator:
   (a) the 5–8 most important new findings with A-numbers;
   (b) reproductions that DEVIATED (the original result may be wrong — list id + deviation);
   (c) recurring packaging errors;
   (d) three recommended next steps with the highest value.
