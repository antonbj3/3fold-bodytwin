# CX-MINCONTACT — does the body choose the MINIMAL knee contact? Clean test of the frozen hypothesis H-MIN

Read first:
- `results/CX-MINCONTACT/FROZEN_HYPOTHESIS.md` (P1–P3 are frozen; `FROZEN_HYPOTHESIS.sha256` is created by the coordinator; do NOT change it).
- `results/CX-INVERSEOC/RESULTS.md`, together with `locations.json`, `analyze.py`, `score.py`, `non_gait.py`, `non_gait_results.json` (these hold the mask, the N1g definition and the LP bounds per frame).
- `results/CX-SETVALUED` (LP machinery).
- `results/CX-DETERMINED` (A1406).
- `results/CX-ANKLETERM`/`CX-LAW2` (the parameter-free law).
- `results/L1` (muscle chain).

## Tasks
1. **P1 on INVERSEOC's exact common mask (11,474 frames) with L1's N1g definition.** Compare these predictors: lo+c (c by LOPO), lo+s̄·width (s̄ by LOPO), the parameter-free law, N1g, and the metabolic proxy. Report per person, per phase bin, and for early stance.
2. **P2, JW4 non-gait.** Take the LP bounds from non_gait.py/non_gait_results.json if they exist; compute them if they are missing. Compare lo+c against N1g per activity. Do not use jw_lungef1.
3. **P3.** Locate the frames where meas < lo, by person and by phase/knee angle. For each such frame, find which constraint is binding in the LP min solution: which muscle is at F0, which is at 0, and whether the contact is at 0 on one side. That constraint is the node the deviation points to.
4. **Mechanism (the most important one):** at the LP minimum, which muscles carry the knee moment? The hypothesis is that the minimum-contact solution uses muscles with the largest moment arm and bi-articular action (gastrocnemius/hamstrings/rectus). Compare the relationship between F_min and the parameter-free law (quadriceps only, fixed r_q) frame by frame. Does lo ≈ F_law? If so, the law is explained as a minimal-contact law. Also consider c as the co-contraction level: is it constant, or does it follow phase?
5. **A counter-test:** time-shift or permute the measured force relative to lo. The advantage must disappear.

PREREG.md + sha256 for tasks 1–5 BEFORE computing; the frozen P1–P3 are the criteria. Write in `results/CX-MINCONTACT/`. Deliver RESULTS.md starting with `# CX-MINCONTACT`, results.json, a script that reproduces the numbers, and pytest. Internal the collaborator data stays local; heavy LP computations go through the shared queue (bigmem.lock) or `tasks/cloud_run.sh` OVH (≤ 12 vCPU, finish by 07:30). Every outcome is a node that expands (no verdict words). lane_runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.

## Addendum (A1415, Field U441)
The hip shows the same pattern: the OrthoLoad cohort's measured force is BELOW gait2392's Fmin in 10/10 deciles (s −0,10). Knee + hip indicate that generic models' minimum is overestimated. Beyond P3: compute per person the minimum scaling of knee extensor/flexor moment arms (one factor per muscle group) that brings measured force ≥ lo in ≥ 95 % of frames, and compare with the literature's moment-arm spread (is the scaling physiologically plausible?). Report whether lo+c with that scaling (LOPO) changes P1.

## Addendum 2 (Field U447 runs the same operator analysis)
Field U447 (romi_collab/build/U447/) computes the minimum-norm shift Δf from LP-min to measured force and which operator component (moment arm/muscle limit/C0) must change least for the frames below Fmin. SKIP the moment-arm scaling in Addendum 1 (duplication). P3 (frozen) remains; put the effort into P1, P2 and the mechanism (task 4).
