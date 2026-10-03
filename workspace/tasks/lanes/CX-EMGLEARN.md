# CX-EMGLEARN — learn the EMG→force gain ACROSS persons (LOPO) and test on the held-out person: does real EMG then give precision?

Background: A1427 CX-FUSION forbade itself from calibrating the EMG gain against the implant force, and EMG added nothing (0.544 vs N1g 0.396; permuted EMG as good). Learning the gain on the 3 OTHER persons' implant force is legitimate supervised learning, and it has not been tested.

Read: `results/CX-FUSION/` (fuse.py, analyze.py: the EMG mapping to L1 muscle groups, λ-prior); `results/CX-INVERSEOC` (LP bounds); `results/L1/prep/*.npz` (emg 15 channels, A, b, F0, cj, C0, meas_tot); `results/CX-EPSBAND` (ε-band, eps_band.py).

## Tasks (PREREG.md + sha256 BEFORE the held-out person is scored)
1. **Model:** muscle force f inside the LP set (equilibrium, box, med/lat ≥ 0). Per frame the estimator = argmin E(f) + w·Σ_g (a_g(f) − G_g·emg_g)², where a_g(f) is the group activation (force/F0 summed per EMG group). G_g is the gain per channel group, and w is a weight.
   - Learn G_g (shared across persons, or with ONE scale factor per person set from that person's MaxEMG normalisation, not from the implant) and w on 3 persons by minimising the implant RMSE.
   - Freeze them, and score the 4th (LOPO ×4). Also an implant-free variant, as a baseline.
2. **Report:** RMSE/peak/early stance per held-out person against N1g (L1 LOPO), lo+c, and FUSION.
   - Permuted EMG (between trials) with the same learned G must lose clearly.
   - Which channel groups carry the gain (leave-one-group-out)?
3. **Frozen criterion:** beats N1g for ≥ 3/4 held-out persons, AND permuted EMG is ≥ 10 % worse.
4. Report separately the frames where meas < lo (the model-error cap). Also give the result on frames inside the set.

Deliver RESULTS.md starting with `# CX-EMGLEARN`, results.json, the script, and pytest. Run the LP/QP under bigmem.lock with 2 threads or on OVH via `tasks/cloud_run.sh` (finish by 07:30). Internal data stays local. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only. Every outcome is a node that expands.

## VARIANT B (the coordinator, 23:40; add it to PREREG before scoring): blind, implant-FREE gain calibration against ID moments (EMG-driven, CEINMS-like; free-energy framing)
Idea: the body's natural stride-to-stride pose variation changes A (the moment arms). A1439 says pose changes are exactly what make the hidden observable. If the gain G per channel group is constant across strides for a person, G can be calibrated WITHOUT the implant: minimise the moment residual ‖A f(EMG; G) − b‖ over all of that person's frames/strides. Here f(EMG;G) = the EMG-driven forces for the recorded groups + a minimal correction for the unmeasured muscles, within the box/set.
- Report the identifiability of G: the condition number/profile likelihood across strides. Does it improve with more pose variation? Compare normal gait against a mix of gait modes.
- Then predict the contact force (C0 + cj·f) and score it against the implant, as in A: per person against N1g, A (LOPO implant-learned), and permuted EMG.
- Frozen criterion: B beats N1g for ≥ 2/4 persons, OR B lands within 10 % of A's RMSE (which would show the implant is not needed for the calibration).
Reference: CEINMS/EMG-assisted (Lloyd & Besier 2003; Sartori 2012; Gerus 2013 on the Grand Challenge). Cite them; novelty is not a criterion.
