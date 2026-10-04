# CX-WHATIF — parametric what-if: anatomical parameters as design variables → joint force, with derivatives and error bounds

Goal: the collaborator plans osteotomies with a femur population (2025) and wants to see how we parameterize shapes. In the reference model, every parameter value requires a new model and a new run. Show that our chain answers "what happens to the hip force if CCD changes −5° or anteversion changes +10°?"
- for any individual;
- with the derivative for each parameter;
- with an error bound;
- orders of magnitude faster than rebuilding per value.

## Build on (read first, do not redo)
- The geometry manager `results/N7c/geomgr` (ops.py: shape editing of CCD, AV and length while other measures are held, A183/A194).
- The exact fast solver `results/CX-SOLVER2/cx_solver.py` (the collaborator's box lift, 141/141 KKT ≤ 1e-10, A295).
- The patella fix (A278).
- The pelvis in the population chain (`results/CX-PELVIS`, A282).
- The emulator per active set with one-sided derivatives (`results/BT-XF4`, A303: 0.146 %).
- The Sobol/HJC sensitivity (A213).
- CCD/AV definitions (A158, A267).
- Grep `~/projects/bodytwin/data/MECHANISM_ANCHOR_GRAPH.json` and `~/projects/bodytwin/scripts/msk/` for osteotomy/what-if/sensitivity before designing. Name what you reuse in PREREG under "Builds on".

## Tasks
1. PREREG.md + PREREG.sha256 before the first run: hypotheses, criteria with numbers, counter-tests.
2. Parameter → force chain for at least 3 individuals (TLEM/z001 plus 2 VSD) and parameters CCD, AV, femur length, and HJC offset AP/ML. The output is the hip reaction peak and the time profile through the box lift, plus the gmed/adductor moment arms.
3. Derivatives: analytic or implicit, via KKT sensitivity for a fixed active set and one-sided derivatives at active-set changes (the XF4 approach). Compare against central finite differences with full rebuild. Criterion: relative error ≤ 1 % where the active set is constant; report active-set changes explicitly.
4. What-if curves: CCD −15…+15°, AV −15…+25°, in 1° steps. Compare three methods for accuracy and time:
   (a) full rebuild + solve per value;
   (b) first-order extrapolation with the derivative;
   (c) the emulator per active set.
   Report the valid interval per method with an error bound. Report the speed ratio.
5. Counter-test and null model:
   - Permuted parameters: the derivative must lose its predictive power.
   - The B24 null model, body weight × OrthoLoad median, does not change with CCD. Report where the shape effect is larger than the model's own uncertainty band (A194/N7c band) and where it is not.
   - Symmetric skepticism: if the effect is below the uncertainty, say so plainly.
6. Clinical anchor: compare the sign and magnitude of dF/dCCD and dF/dAV against published varus/valgus osteotomy directions, e.g. a reduced CCD angle increases the abductor lever and reduces hip force. Use only sources you can cite with a DOI; otherwise mark UNKNOWN.

## Resources
- Local: nice, 2 threads, < 60 s per test.
- Anything heavier (full rebuild over 3 × 60 parameter values) goes on Modal: `modal run tasks/modal_run.py --job-id CXW-<n> --src <dir> --cmd "<cmd>" --cpu 8 --mem-gb 16`, max 3 h. Split into shards.
- If the sandbox has no network, do everything that can run locally at ≤ 60 s per test, write shard scripts plus `RUN_ON_CLOUD.md` with exact commands, and state it in RESULTS.
- Write only under `results/CX-WHATIF/` (large intermediate results in `external_media`). `~/projects/bodytwin` is read-only.
- restricted model data and the collaborator's data stay internal. No emails or pushes.

## Deliverables
- `RESULTS.md` starting with `# CX-WHATIF`: table parameter → dF/dp per individual → FD agreement → valid interval → time (a)/(b)/(c) → effect vs uncertainty.
- `results.json`, the code, and a pytest with at least 3 tests: derivative vs FD, zero parameter change gives an identical force, permutation test.
