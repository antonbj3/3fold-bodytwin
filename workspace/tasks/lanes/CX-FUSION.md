# CX-FUSION — sensor fusion: hard set + H-MIN prior + measured EMG → posterior for knee contact, per frame

Anton: "sensor fusion was the first word that struck me". The pieces already exist separately; this lane fuses them into ONE estimator and tests it against the implant force.

## Pieces (read first; do not rebuild what exists)
- **Hard constraint:** the L1 feasible set per frame (lo/hi), from `results/CX-INVERSEOC/locations.json` + analyze.py and the LP in `results/CX-SETVALUED`.
- **Prior:** `results/CX-MINCONTACT/FROZEN_HYPOTHESIS.md` (H-MIN: the measured force lies near the lower edge; A1410: median s 0.005–0.085). While CX-MINCONTACT is running, use only INVERSEOC's data; do not wait for it.
- **Measurement:**
  - the Grand Challenge's MEASURED EMG (the EMG columns in the gait files; the GC sessions in `results/L1/code/n12_io.py`). What has been done with EMG before: `results/CX-BEATN1G*`, `CX-IMUFORCE`, `CX-BAYESOED` (A1400, the minimal tendon triple), plus Field U405/U413/U419 (A1389/A1394/A1397: EMG ratios + stiffness make the contact structurally determined in gait2392).
  - the law/moments from ID (`CX-ANKLETERM`/`CX-LAW2`).
- **Field hip (A1409):** EMG shrinks the width by 25–56 %, and the effort cap by 77–94 %.

## Tasks (PREREG.md + sha256 BEFORE the comparison against the implant force)
1. Build the estimator: posterior over muscle forces f inside the LP set.
   - Prior: effort-based, concentrated towards min contact; the precision parameter is fitted LOPO over persons.
   - Likelihood: measured normalised EMG → activation → force, with an unknown gain per muscle and person (a hierarchical/marginalised gain, NOT fitted against the implant force).
   Output: posterior mean + 90 % interval for the knee contact per frame. Sampling method of your choice (hit-and-run in the polytope, or a Laplace/QP approximation), with the timing reported.
2. Test against the measured implant force (DM/SC/PS/JW gait; JW4's non-gait; not jw_lungef1):
   - RMSE/peak/early stance against N1g, the law, and lo+c;
   - coverage of the 90 % interval, per person.
3. Value per sensor: the posterior width and error with and without EMG; EMG per muscle group (leave-one-channel-out). Which channel reduces the most? Compare against A1400's minimal triple.
4. Counter-test: permuted EMG (between trials) must remove the EMG gain.

Deliver RESULTS.md starting with `# CX-FUSION`, results.json, code in bodytwin_core-compatible form (a `fuse_knee(trial)` function), and pytest. Internal data stays local; heavy computations go through the shared queue or `tasks/cloud_run.sh` OVH (≤ 12 vCPU, finish by 07:30). Every outcome is a node that expands. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.

## Addition (A1412, Field U440)
The body's gauge is actuator redundancy (956 out of 1 133 null directions = muscles with support on ≤2 segments), not graph cycles; the compliance path is ill-conditioned at 1 mm pose noise. This strengthens the case that fusion should rely on the recruitment prior + EMG, not on joint stiffness. Feel free to report the posterior width broken down by muscle group (which redundant groups EMG actually closes).
