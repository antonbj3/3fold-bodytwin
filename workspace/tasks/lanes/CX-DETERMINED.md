# CX-DETERMINED — "the determined twin": joint load split into what the data DETERMINE, what they do NOT determine (as a set), and what has to be measured to close it

Today's findings are gathered into ONE capability in bodytwin_core. the reference model gives one number, chosen by the recruitment criterion (A1384, A1385, A1396). We give three things.

## Build on (read first; build no new models)
- The parameter-free law: A1391 (CX-ANKLETERM) and A1398 (CX-LAW2), with F_law=‖GRF‖+|M_knee|/r_q+0.29|M_ankle|/0.045. It reproduces N1g's level without fitting.
- The set: A1392 (CX-SETVALUED), `bodytwin_core.solver.solve_bounds`.
- The gauge and determinacy: A1384 (BODYGRAPH), A1396 (U414 on the collaborator), and A1394/A1400 (U413/BAYESOED: the minimal tendon triple).
- The null models: `bodytwin_core.nullmodels` (N1g, B24).
- The geometry priors: A1401.

## Build `bodytwin_core.determined` with a `knee_load(trial)` / `joint_load(system)` API that returns per time step:
1. **Determined:** the net joint moments + external load (from ID), and the law's total contact F_law with its uncertainty from the r_q/r_achilles/β literature bands. State that this is determined by the data up to the lever-arm uncertainty.
2. **Not determined:** [min, max] of the total and (where the model has it) the medial/lateral contact force over the feasible set, under physiological constraints (Fmax, one-sided contact), via solve_bounds. Also mark where the reference model's/min-effort's value would lie, and N1g as the reference.
3. **What to measure:** the minimal sensor set that closes the output (from BAYESOED/U413), with the expected posterior width per cost level, and a plain-text statement of whether the current data are enough.
4. A test on the Grand Challenge (JW, DM, SC, PS; gait + JW4's new activities; NOT jw_lungef1):
   - does the measured implant force lie inside the set, and where?
   - does F_law lie within its uncertainty band of the measured force?
   - coverage per PERSON (the person is the unit).
5. `python -m bodytwin_core.demo --determined` gives a table per trial and one figure (PNG): the determined band, the undetermined set, the measured value, and N1g.

PREREG.md + sha256 first (criteria: coverage of the measured force in the set ≥ 90 % of frames per person; F_law inside its band in ≥ 50 % of frames — report the numbers and do not tune afterwards). Every outcome is a node that expands (FRAMING_NOTE). Write in `results/CX-DETERMINED/`. `RESULTS.md` starting with `# CX-DETERMINED`, plus pytest. lane runner has full permissions in the workspace; internal the collaborator data stays local.
