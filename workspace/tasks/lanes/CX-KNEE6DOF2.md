# CX-KNEE6DOF2 — step 2: a physical starting pose + a muscle-path operator for the knee's secondary DOF, then the first FDK solution

Background: A1788 CX-KNEE6DOF step 1 (read results/CX-KNEE6DOF: knee6dof.py, PREREG, results.json).
- Already built: ligaments, mesh contact + a 5×5 Jacobian, and the ID target screw.
- Two things are missing:
  - (1) The generic implant pose gives 17–27 mm penetration and ~430 BW, so it is physically unreasonable as a starting point.
  - (2) A pose-dependent 5-row muscle-path operator (the muscle forces' contribution to the secondary DOF, including the patella equilibrium).

## Tasks (PREREG.md + sha256 first; nothing fitted against implant force)
1. **A starting pose with no implant fitting:** per frame, search the secondary DOF for a pose with minimal contact penetration, near zero contact force under a small preload (e.g. 0.1 BW), and ligaments at their reference strain. That is the reference model's standard "initial condition". Report the distribution of penetration/force at the starting pose.
2. **Muscle-path operator:** move the tibia a small step in each of the 5 secondary DOF and compute the change in length of every L1 path (virtual work: moment arm = −dL/dq). Build it from the same path geometry as L1 (n12_model), with the quadriceps via the patella solver. Check it against A's knee row for flexion (it must reproduce it).
3. **FDK solution:** Newton/LM over the 5 DOF + recruitment (min-max + tie-break) per frame, warm start, on PS5 normal gait + right turn (the whole cycle). Report:
   - convergence per phase;
   - total RMSE against the implant force, stance/swing, in Marra's window, against Marra et al.'s FDK 0.26/0.29 and N1g (0.49 in swing);
   - time per cycle.
   Frozen criterion: convergence ≥ 90 % AND the swing RMSE < 0.49 BW.
4. Deliver `knee6dof_solve(trial)` + pytest (a synthetic case with a known solution).
If one sub-step does not hold: deliver the rest, and say exactly where it stops.

GPU allowed (a single process, NEVER `nvidia-smi -q`). 2 CPU threads under bigmem.lock, at most 300 MB of intermediate files on external_media Deliver RESULTS.md starting with `# CX-KNEE6DOF2` and results.json. Internal data stays local. lane_runner has full permissions in the workspace; `~/projects/bodytwin` and romi_collab are read-only.

## Restart 14:00
The previous run was killed (SIGKILL 13:46, system crash in another session). Continue from existing files in results/CX-KNEE6DOF2/ (start_pose.py, muscle_paths.py, validate_paths.py, solver.py, finalize.py, PREREG unchanged); do not rerun what is already done. Keep memory < 6 GB.
