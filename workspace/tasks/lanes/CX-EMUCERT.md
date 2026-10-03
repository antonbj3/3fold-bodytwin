# CX-EMUCERT — faster where it is enough: the emulator per active set with a certificate inside the solver, and an exact solution as the fallback

Why:
- BT-XF4's emulator (A303, audited HOLDS: max 0.146 %) and warm start (B176: 1.5×, 141/141 KKT) exist, but only as separate results.
- The solver in `bodytwin_core` (A910) always runs the full QP.
- Principle (Field D3 / plan B3): cheap everywhere; a certificate decides where an exact solve is needed.

## Build on (requirement: "Builds on" with file paths)
- `results/BT-XF4` (emulator per active set, one-sided derivatives), B176 (warm start).
- `results/CX-SOLVER2/cx_solver.py`; `results/CX-SOLVER-ROBUST` (adaptive N2b, KKT/dual-gap checks).
- `results/CX-WHATIF2` (implicit KKT derivative).
- `bodytwin_core/solver`.
- Field's certificate principle: `romi_collab/build/U381` (certify_refine: cheap → band → INSIDE/OUTSIDE/UNCERTAIN → refine; the UNCERTAIN state gets a fixed rule), read-only.

## Tasks (PREREG.md + sha256 first)
1. The solver mode "certified":
   - per time step, predict with the emulator (active-set classifier + one-sided derivatives);
   - check KKT residual and dual gap on the prediction;
   - if the certificate holds (bound ≤ tolerance): accept; otherwise, an exact solve with warm start.
2. Measure on the collaborator's 141 steps and on N2b's population, 1,000 individuals × 141 steps (or a representative subset):
   - the share of accepted steps;
   - the maximum error against an exact solution;
   - time per step and in total.
   Criterion: max relative force error ≤ 0.5 % in 100 % of steps (the exact fallback guarantees it), and a speed-up ≥ 3× on the population.
3. Counter-tests:
   - deliberately perturbed emulator coefficients must be caught by the certificate (falls back, never accepts);
   - an active-set change mid-step.
4. Integrate as an option in `bodytwin_core/solver` (mode="certified"), with a test. The default stays exact.

Internal the collaborator data runs locally or on OVH. Write in `results/CX-EMUCERT/`. `RESULTS.md` starting with `# CX-EMUCERT`, plus pytest.
