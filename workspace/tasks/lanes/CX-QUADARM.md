# CX-QUADARM — is L1's QUADRICEPS arm the model error? A frozen confirmation of A1449's single-swap finding

Background: A1449 CX-L1ARMS2, pre-registered single swaps. Replacing ONLY the quadriceps/patella chain's sagittal knee arm in A[3] with Rajagopal 2016 → frames under lo 22.55 % → 9.71 %. The full knee swap only gives 18.3 % (the other swaps counteract). L1 uses TLEM 2.0 + its own patella solver (rho 0.4–1.3, results/L1/code/n12_model.py). The finding came out of an 8-way decomposition, so it must be confirmed with frozen tests before it is used.

## Tasks (PREREG.md + sha256 FIRST; freeze the criteria below)
1. **Mechanism.** Plot L1's effective quadriceps arm against Rajagopal's and against the literature over knee angle, per person. Literature: Herzog–Read, Im et al. EQma (in results/CX-RQANGLE), Buford. Which one is off, and by how much? Is the error in the patella solver (rho, the patellar tendon arm) or in the TLEM paths?
2. **Frozen confirmation 1, other activities (not used in the decomposition):** JW4 non-gait (11 trials; operators in /media/anton/sdc1-tmp/bodytwin/CX-INVERSEOC/non_gait/), and DM/SC/PS non-gait if CX-NONGAITOPS has delivered them. Criterion: the share under lo falls ≥ 40 % relatively with the quadriceps-arm swap only.
3. **Frozen confirmation 2, a literature curve instead of Rajagopal:** the Im et al. EQma curve (female/male) as the quadriceps arm, on the L1 gait mask. Criterion: the share under lo ≤ 13 %.
4. **Consequences** (report; no criterion):
   - ε=10 % coverage/width (eps_band.py);
   - lo+c and INVERSEOC's objectives (stress-2, metabolic proxy) against N1g with the new arm;
   - JW separately.
   Does the new arm move ANY model-based predictor closer to N1g?
5. Deliver `knee_geometry(model='l1_quadarm_rajagopal')` (+ an 'im_eqma' variant) + pytest.

Deliver RESULTS.md starting with `# CX-QUADARM` and results.json. Reuse CX-L1ARMS2/extract_consistent.py and geometry.py. Run the LP under bigmem.lock with 2 threads. Internal data stays local; no jw_lungef1. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only. Every outcome is a node that expands.
