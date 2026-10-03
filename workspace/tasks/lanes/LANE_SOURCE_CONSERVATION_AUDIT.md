# LANE_SOURCE_CONSERVATION_AUDIT

Why: the BodyTwin swarm (The swarm/swarm_worker/swarm_worker) builds hundreds of jobs on the source families' models in `tasks/free48/sources/<KEY>/` (43 families: Q005…Q168, MITOSTRESS, IMMUNITY, SOLBENCH, BIORESP, SURG_*). The graph lane found on 1/10 that Q115's `age_transition` silently lost the last age bin's cells — the inventory did not close (witness: 1000 cells, dt 1 h, hazard 0,12/h → 22,78 new, 120 shed, 880 disappeared). Such an error in a source model poisons every job that inherits it. Results directory `results/LANE_SOURCE_CONSERVATION_AUDIT/`.

## Uppgift

For each source family's executable model (`*_model.py` and other .py in sources/<KEY>/):
1. Identify the conserved quantities (mass, cells, fluid volume, charge, energy, amount of substance) and their ports (inflows/outflows, sources, sinks, boundaries, cutoffs).
2. Write an **inventory witness**: run the model's own functions on a constructed stress state (like G2's 1000-cells-in-the-last-bin) and check that Δinventory = the sum of reported ports to machine precision. Test limiting cases: last bin, empty source, zero/maximum parameters, long times, changed time steps.
3. Check units in every port sum (dimensional analysis) and that time-step dependence disappears in the limit dt→0.
4. Classify each finding: LEAK (unreported flow), UNIT ERROR, TIME-STEP ARTIFACT, CUTOFF THAT IS DECLARED AND NEGLIGIBLE (like Q115 at baseline), or OK. State size at the model's own default parameters and in which parameter region it becomes large.

## Leverans

- `AUDIT_TABLE.json`: family → conserved quantity → witness → result (balance residual, class, size at baseline and worst realistic case) → file:line.
- One executable test script per family (`witness_<KEY>.py`) that goes green when the model is conservative.
- For real leaks: a minimal **proposal** for correction (diff) in the results directory — do not change source files yourself; the coordinator applies it.
- RESULTS.md: which families are safe to build on, which have errors and which job types inherit them.

The witnesses are code that runs, not reading. Strongest control: G2's Q115 witness must be reproduced as a positive control (your method should find the known leak in the old version, `tasks/free48/Q115_model.py.bak_20261001_cutoff`). No internal data.
