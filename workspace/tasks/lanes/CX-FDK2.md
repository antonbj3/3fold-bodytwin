# CX-FDK2 — a knee-local FDK loop on PS5 over the WHOLE cycle with a GENERIC implant pose (without CT registration), scored against Marra et al.'s FDK

Background:
- A1784 CX-FULLCYCLE: the IK/ID/L1 operator exists for the whole cycle on PS5 (297 frames, LP feasible everywhere; /media/anton/sdc1-tmp/bodytwin/CX-FULLCYCLE/). N1g in swing gives 0.49 BW.
- A1785 CX-IMPLREG: a verified implant↔marker registration is impossible from the archive.
- A1781: Marra et al.'s FDK 0.26/0.29 BW (whole cycle) against our best 0.340/0.344.
- A1782 CX-FDK1: code for the FDK structure, the recruitment criterion (min-max + quadratic tie-break) and the contact law (9.3 N/mm³ as in Marra et al.) are already frozen there. Reuse them.

## The generic pose (freeze it BEFORE any scoring against implant force)
Place the femoral component and the tibial insert from marker/functional geometry ONLY:
- the joint line at the knee centre (the epicondyle midpoint / SCoRE);
- the insert's normal along the tibial axis;
- the femoral component's flexion axis along the transepicondylar axis (or SARA);
- the size from the PS5 STL.
Propagate the uncertainty: ±5 mm translation, ±3° rotation, as Monte Carlo over the pose (at least 20 draws, fixed seed). Report the result distribution, not a single number.

## Tasks (PREREG.md + sha256 first)
1. FDK per frame: 5 secondary DOF + recruitment + contact until the residual is below the tolerance. Warm start. The whole cycle on normal gait + right turn.
2. Frozen criteria:
   - convergence ≥ 95 %;
   - the MEDIAN total RMSE over the pose draws ≤ 0.29 BW in normal gait (the value reported by Marra et al.; the target is 0.26);
   - time per cycle.
   Also report stance vs swing, medial/lateral, and the spread over the pose draws (how much of the error the pose uncertainty explains).
3. Deliver `fdk_knee(trial, pose)` + pytest.
If it does not converge: which phase/DOF, and what that points to.

GPU allowed (a single process, NEVER `nvidia-smi -q`). 2 CPU threads under bigmem.lock, at most 300 MB of intermediate files on /media/anton/sdc1-tmp (check df first). Deliver RESULTS.md starting with `# CX-FDK2` and results.json. Internal data stays local. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
