# CX-SOLVER-ROBUST — make our fast solvers robust to scaling and ill-conditioning, and fix A144's suboptimal exoskeleton QP

Leads:
- Cloud W4-A363-ALGO: an independent batched active-set NNLS passed 1000/1000, but with column scaling 1e-6…1e6 only 31/100. A norm dominated by large columns can hide errors in small ones.
- W3-BATCH-CERT: batched FISTA reached 96.5 % (FAIL at cond 1e8), and a single projected-gradient stop gave false certificates.
- W4-A295-EXACT: the preregistered dual certificate proved ≤ 1e-8 gap in only 4/100, even though the true gap was small in 100/100 (a better dual point gave 100/100).
- A1269: A144's exoskeleton cap-QP is suboptimal in 95–99/141 steps; its own reduction is 17.1/30.3/36.5 % vs the source's 13.4/15.1/15.1 %.

## Tasks (PREREG.md + sha256 first)
1. Our solvers:
   - `results/CX-SOLVER2/cx_solver.py` (A295);
   - N2b's batched solver (`results/N2b`, A363);
   - the emulator (BT-XF4, A303).
   Stress-test them with the cloud's scaling/conditioning generators (copy them from `results/CLOUD-W4-A363-ALGO`, `CLOUD-W3-BATCH-CERT`, `CLOUD-W4-A295-EXACT`). Report the pass rate per scale/condition class.
2. Fixes:
   - column/row equilibration (Ruiz, or diagonal scaling);
   - a scale-invariant stop (relative per-column KKT + dual gap);
   - a better dual point for the certificate (the cloud's post-hoc variant, now preregistered).
   Criteria:
   - ≥ 99 % on the stress generator;
   - no regression on the collaborator's 141 steps (KKT ≤ 1e-10, forces identical to ≤ 1e-8);
   - time increase ≤ 20 %.
3. A144: repair the cap-QP with an optimality check (KKT/dual gap per step), then rerun. Write a correction `[CORRECTED … by CX-SOLVER-ROBUST]` inline in A144's row in `notes/RESULTS_INDEX.md` ONLY once the new frozen audit is done.
4. Counter-test: permuted scaling must not change the solution (invariance).

Resources: lane runner has full permissions in the workspace. the collaborator's data is internal (local/OVH only, via `tasks/cloud_run.sh`, BodyTwin ≤ 12 vCPU, finish by 23:59). Write in `results/CX-SOLVER-ROBUST/`. `RESULTS.md` starting with `# CX-SOLVER-ROBUST`, plus results.json and pytest.
