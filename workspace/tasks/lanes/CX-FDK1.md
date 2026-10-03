# CX-FDK1 — step 1: comparison with the published reference: a KNEE-LOCAL FDK loop on PS5, scored in Marra's window against Marra et al.'s FDK (0.26/0.29 BW)

Background:
- A1781 CX-ANYBODYWIN: over the whole gait cycle Marra et al.'s FDK has lower reported error (0.26/0.29 vs our best 0.340/0.344 BW). The difference sits in swing/low-GRF phases, where our static/LP chain has no model of secondary knee kinematics and contact.
- Building blocks:
  - L1's operator + the quadriceps-arm geometry (results/CX-QUADARM/geometry.py);
  - the contact law calibrated against BEM (Field U434: 0.488×Smith, 9.3 µs/call, romi_collab/build/U434/code/comak.py, read-only);
  - the direct GPU contact (research/ANYBODY_PERFORMANCE_20260924, 0.41–0.58 ms);
  - PS5's implant geometry (the same meshes as in the performance study);
  - the recruitment QP (bodytwin_core/contact_band_batch.py, analyze.qp).

## Scope (keep it SMALL; this is step 1)
Per frame, iterate: secondary knee DOF (tibiofemoral translation/rotation, 5 DOF; flexion prescribed from IK) → contact force from the geometry + contact law → equilibrium of the knee's secondary DOF together with the muscle recruitment (min stress-2 or min-max, the reference model's default) → a new pose, until the residual is below the tolerance.
Warm start from the previous frame. PS5 normal gait + right turn; the whole cycle in Marra's window (results/CX-ANYBODYWIN defines the window).

## Frozen (PREREG.md + sha256 first)
- Recruitment criterion and contact law chosen BEFORE any scoring against implant force (the min-max polynomial used by Marra et al. as the default).
- Criteria:
  - (a) convergence in ≥ 95 % of frames;
  - (b) total RMSE ≤ 0.29 BW on normal gait (the value reported by Marra et al., with 0.26 as the target);
  - (c) timing per gait cycle against Marra's 4.5/13.6 min, with the hardware stated.
Report the medial/lateral split and the secondary kinematics too (against fluoro if it exists for PS5).

Deliver RESULTS.md starting with `# CX-FDK1`, results.json, `fdk_knee(trial)` in bodytwin_core-compatible form, and pytest. The GPU is allowed (a single process; NEVER run `nvidia-smi -q`). 2 CPU threads under bigmem.lock, at most 300 MB of intermediate files on /media/anton/sdc1-tmp (the disk is almost full). Internal data stays local. lane runner has full permissions in the workspace; `~/projects/bodytwin` and romi_collab are read-only. If it does not converge: report where (phase, DOF), and which node that points to.
