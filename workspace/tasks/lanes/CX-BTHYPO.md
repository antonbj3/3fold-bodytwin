# CX-BTHYPO — test the BREAKTHROUGH HYPOTHESIS broadly with Space The swarm: a parameter-free physical law for knee contact force

## The hypothesis (the synthesis of today's results)
- N1g = k(activity)·|GRF| beats all muscle models (A312, A324, A358, A361–A368, A913).
- Field's U400 (A1382): axial force determines only F_med+F_lat; the medial/lateral split is determined by the frontal-plane (adduction) moment.
- U401 (running): which forces are determined by motion+load and which only by the cost function.
- F-8: contact from measured kinematics is ill-posed.

**H:**
- Total knee contact ≈ the external load + the muscle force that is REQUIRED to balance the net knee moment: F_tot ≈ |GRF_axial| + |M_ext_flex|/r_ext(θ), with r_ext the extensor/flexor moment arm by knee angle.
- The medial/lateral split ≈ ½F_tot ± M_add/d_cond, with d_cond the condyle distance.
- k(activity) then follows from geometry + the moment ratio; it is NOT fitted.
- Muscle models fail because they add co-contraction that is not determined by motion/load.
Falsifier: the parameter-free law is not ≥ as good as N1g (median RMSE, LOPO) on implant force in ≥ 3/4 persons, or the medial/lateral split does not follow M_add.

## Data (local/OVH; do NOT open Field's sealed F-8 facit jw_lungef1 knee_forces/eTibia)
- Grand Challenge: the implant force (total + medial/lateral where the sensors provide it) for all persons and activities.
- Net knee moments from ID: `results/L1/prep/*.npz` (b = net moments, 113 gait trials), and compute more via CX-LUNGEID's ID pipeline for non-gait activities (stairs, step-up, chair rise, squat, etc.).
- `results/DATAMATRIX_TABLES/` (GRF peak/impulse, P-CONTACT-COMP).
- r_ext: literature curves with DOI (patellar tendon arm vs angle) + JW's geometry (A366); d_cond from the implant geometry.

## Build ~30 The swarm packets `results/BT-BH-001…` (queue proposal only; the coordinator reviews)
Each packet = one person × one activity group (or one held-out person), with raw data + the ID moments in inputs. Each computes:
1. N1p (the parameter-free law), N1g (LOPO k) and B24 against the measured implant force: RMSE, peak error, correlation;
2. the medial/lateral split from M_add against the measured medial/lateral force (where it exists);
3. a sensitivity to r_ext (literature band) and d_cond;
4. a counter-test: M_ext time-shifted/permuted → the gain must disappear.
Also 3–5 packets that pool across persons (LOPO) and one that tests activities N1g has never seen (stairs/step-up/chair rise).

Write `results/CX-BTHYPO/QUEUE_PROPOSAL.txt` + `SAMPLE.md` (2 full BRIEFs). Do NOT queue. Also write `results/CX-BTHYPO/HYPOTHESIS.md` (the hypothesis, the falsifier, and which registry rows it builds on). lane runner has full permissions in the workspace.
