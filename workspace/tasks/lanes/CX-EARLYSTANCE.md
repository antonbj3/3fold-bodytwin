# CX-EARLYSTANCE — why does the muscle model overestimate knee force in early stance (+1.02 BW JW, +0.67 BW DM) regardless of recruitment?

L1 (A312, `results/L1/RESULTS.md`, `scores.json`, `code/`) showed that no recruitment criterion removes the early-stance excess. The main candidate is the inverse knee extension moment. Build on L1's code and data (N12b chain, c′ geometry, N39 scoring, COMAK frames), plus BT-N12/N12b/N39/N57 and the graph's KNEE-CELL/COMAK nodes (`~/projects/bodytwin/data/MECHANISM_ANCHOR_GRAPH.json`, grep `knee_jw_transplant`, COMAK in `~/projects/bodytwin/scripts/msk`). Do not redo L1.

## Tasks (PREREG.md + sha256 first)
1. Decompose early stance (0–20 % of stance) for JW1 and DM6. Compare our knee extension moment against COMAK's inverse dynamics and against the moment from GRF × lever arm (Grand Challenge's own ID if available). Quantify the moment error in N·m and %.
2. Attribute the knee force difference to separate sources, each computed alone:
   (a) moment error from ID (segment inertia, CoP/GRF filtering, knee centre);
   (b) co-contraction of hamstrings and gastrocnemius;
   (c) patellar tendon lever arm/knee geometry;
   (d) the contact model.
   Swap one part at a time against COMAK's value. Which one closes ≥ 50 % of the excess?
3. Null models N1g and B24 on the same windows. The criterion for "explained" is that a single swapped part reduces the median early-stance error by ≥ 50 % in both persons. Counter-test: swap with a permuted/time-shifted COMAK value.
4. `results/CX-EARLYSTANCE/RESULTS.md` starting with `# CX-EARLYSTANCE`, plus `results.json` and code with pytest. Local ≤ 60 s per test; heavier runs via Modal (`tasks/modal_run.py`) or write RUN_ON_CLOUD.md. Write only in your directory. `~/projects/bodytwin` is read-only.
