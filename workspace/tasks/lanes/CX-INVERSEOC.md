# CX-INVERSEOC — what does the BODY choose inside the undetermined space? Inverse optimal control against in vivo implant force

## Question
- The contact force is not determined by kinematics+GRF (A1384, A1385, A1396). Every model, the reference model included, picks a point with a cost function.
- We have both things together: the exact geometry of the feasible set (the CX-SETVALUED/solve_bounds machinery, L1's muscle chain) AND the measured implant force per frame (Grand Challenge, 4 persons, 108 gait trials + JW4's stairs/step-up/chair rise/squat; NOT jw_lungef1).
- So turn the question around: WHERE in the feasible set does the measured force lie, frame by frame? What recruitment objective puts the solution there?

## Build on (read first)
- `results/L1` (the muscle chain, A, b, F0, the contact projection; code/, prep/*.npz).
- `results/CX-SETVALUED` (LP min/max over the feasible set).
- `results/CX-BEATN1G`/`CX-BEATN1G2` (H1: the feasible interval; C1: certified clipping — JW failed because of the DM-transplant geometry).
- `results/CX-EARLYSTANCE` (the early-stance overestimate).
- `results/CX-ANKLETERM` / `CX-LAW2` (the parameter-free law).
- `results/CX-GC-STRENGTH*` / `CX-SLACK` (the operating range).
- The graph session's input: `results/CX-BAYESOED/GRAPH_INPUT.md` (person as the unit, prior share).
- The eTibia units (A1390): lbf → ×4.4482 N/lbf.

## Tasks (PREREG.md + sha256 BEFORE analysis: fix the candidate objectives and the criteria)
1. **Location in the set:** per frame, compute [F_min, F_max] of the contact force (L1's system, Fmax constraints). Normalise the measured force to its position s = (F_meas − F_min)/(F_max − F_min).
   - Is s stable within a person or activity, or does it follow phase/angle/speed?
   - Is the measured force outside the set, which would be a model error: geometry, e.g. JW's transplant?
   Report the distribution of s per person, phase and activity.
2. **Candidate objectives** (inverse OC, fixed before looking at the fit):
   - (a) minimum muscle stress Σ(F/PCSA)^p for p = 2, 3 (the reference model/min-effort);
   - (b) minimum fatigue (min–max activation);
   - (c) minimum metabolic cost (Umberger/Bhargava-like);
   - (d) **a stability constraint:** minimum knee stiffness from co-contraction, k_joint ≥ k_req (with k_req as ONE parameter per activity group, fitted LOPO), plus min effort.
   For each: the predicted contact force vs measured (RMSE, peak, early stance), LOPO over persons, compared with N1g and the parameter-free law.
3. **Identifiability:** can the objectives be distinguished with the data (does the difference in predicted contact exceed the measurement noise)? Report which pairs are separable.
4. Counter-tests:
   - a permuted/time-shifted measured force must destroy the stability of s;
   - (d) with k_req = 0 must reduce to min effort.
5. Every outcome is a node that expands (FRAMING_NOTE).
   - If (d) wins: that is a physiological law for co-contraction, and a muscle model that beats N1g for the RIGHT reason.
   - If none wins: report what s shows about the body's strategy.

Resources: local (L1's data is internal) or OVH (`tasks/cloud_run.sh`, BodyTwin ≤ 12 vCPU, finish by 07:30). Write in `results/CX-INVERSEOC/`. `RESULTS.md` starting with `# CX-INVERSEOC`, plus results.json and pytest. lane_runner has full permissions in the workspace.

6. **Addendum 22:10:** read results/CX-INVERSEOC/EPS_OPTIMAL_NOTE.md (ε-suboptimal set: how close to optimal is the body, ε* per frame).
