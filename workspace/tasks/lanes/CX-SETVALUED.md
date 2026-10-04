# CX-SETVALUED — the constructive consequence of the gauge finding: joint loads as a SET (min/max over everything the data allow), not a point chosen by a cost function

## Background (A1384, A1385)
- Kinematics+GRF determine only the net joint moments (BODYGRAPH: rank 17, gauge 75; the collaborator's system: nullity 1,126, no contact coordinate determined).
- the reference model/OpenSim/our models present ONE contact force, chosen by the recruitment criterion.
- We have a fast exact min/max over the feasible set: the batched solver (A363: batched min/max = HiGHS ≤ 3.5e-7, 25× N14 bisection) and CX-SOLVER2.
- Field's U405 (running) asks which measurements CLOSE the gap. This lane asks how BIG the gap is, and what physiological constraints do.

## Tasks (PREREG.md + sha256 first)
1. **the collaborator's box lift** (bodytwin_core, 141 steps): per step, min/max of hip and knee contact force (and the medial/lateral split if the model has it) over {muscle forces ≥ 0, equilibrium}, then with physiological constraints added one at a time:
   (a) activation ≤ 1 (Fmax);
   (b) the force–length/velocity capacity at the actual fibre length;
   (c) co-activation limits from literature or EMG where available;
   (d) min/max muscle stress σ ≤ σ_max.
   Report the width of the set per constraint, and where the reference model's own value (the collaborator's h5) and our min-effort lie inside it.
2. **Grand Challenge gait with implant force** (L1/N12b chain, 4 persons; F-8's jw_lungef1 excluded):
   - Does the measured implant force lie inside the set? How wide is the set in BW?
   - Which constraints make the set narrow enough (< 0.3 BW) to still contain the measured value?
   - Compare with N1g: does N1g lie inside the set, and where?
3. **Product**: `bodytwin_core.solver.solve_bounds(...)` gives [min, max] per step and per joint, with a test. the collaborator's data stays local/OVH.
4. Every outcome expands a node (framing A1386 / FRAMING_NOTE): a wide set means the next measurement comes from U405; a narrow set that holds means a certified interval twin.

Write in `results/CX-SETVALUED/`. `RESULTS.md` starting with `# CX-SETVALUED`, plus results.json and pytest. lane_runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
