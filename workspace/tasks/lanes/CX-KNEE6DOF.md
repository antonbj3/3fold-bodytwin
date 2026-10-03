# CX-KNEE6DOF — step 1 of the 6-DOF knee (the gap to Marra et al.'s FDK, A1786): quasi-static equilibrium in the knee's 5 secondary DOF, PS5 whole cycle

Background:
- A1786 CX-FDK2: L1 is a 1-DOF knee. To model swing/low load the way Marra et al.'s FDK does, we need equations for the 5 secondary DOF (AP/proximal/lateral translation, internal rotation, varus):
  - ligament forces,
  - contact force + Jacobian from the implant meshes,
  - muscle forces projected onto those DOF.
- There are building blocks:
  - the generic implant pose + MC draws (results/CX-FDK2, fdk.py, seed 1786);
  - the whole-cycle operator (results/CX-FULLCYCLE, /media/anton/sdc1-tmp/bodytwin/CX-FULLCYCLE/);
  - the PS5 STL (research/ANYBODY_PERFORMANCE_20260924);
  - Field's COMAK code (romi_collab/build/U434/code/comak.py and U458/U483 MovingContact; READ-ONLY, copy the functions with a source note);
  - the contact law 9.3 N/mm³ (Marra et al.) and BEM-calibrated alternatives.

## Build (PREREG.md + sha256 first; freeze all parameters BEFORE any implant scoring)
1. **Ligaments** as nonlinear springs (Blankevoort 1991 / Smith 2016 COMAK parameters: ACL-free for a PS implant, PCL-resected for PS; MCL, LCL, posterior capsule). Attachments scaled from the literature to PS's marker geometry, with the uncertainty stated.
2. **Contact:** femoral-component mesh against the insert mesh at the generic pose, with force + Jacobian w.r.t. the 5 DOF (finite difference or analytic).
3. **Muscles:** the L1 elements' force vectors projected onto the 5 DOF (from the path geometry; the quadriceps via the patella, as in L1).
4. **Per frame:** solve jointly for the 5 DOF + recruitment (min-max activation + quadratic tie-break; the reference model's default) so that knee equilibrium in all 6 DOF holds. Newton/LM with a warm start. Measure the time.
5. **Report:**
   - convergence per phase;
   - the total contact against the implant force (RMSE in Marra's window, stance/swing) against Marra et al.'s FDK 0.26/0.29 BW and N1g;
   - the spread over the MC poses.
   Frozen criterion (step 1): convergence ≥ 90 % AND the swing RMSE < N1g's swing RMSE (0.49 BW).

GPU allowed (a single process, NEVER `nvidia-smi -q`). 2 CPU threads under bigmem.lock, at most 300 MB of intermediate files on /media/anton/sdc1-tmp (check df). Deliver RESULTS.md starting with `# CX-KNEE6DOF`, results.json, `knee6dof_step(trial, frame, pose)` + pytest. Internal data stays local. lane runner has full permissions in the workspace; `~/projects/bodytwin` and romi_collab are read-only. If something is missing for a sub-step: report exactly what, and deliver the parts that work.
