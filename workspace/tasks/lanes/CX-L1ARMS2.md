# CX-L1ARMS2 — locate WHICH muscle groups, projections and phases carry L1ARMS's improvement, then a CONSISTENT geometry (arm + projection + frontal arm from the same path)

Background: A1444 CX-L1ARMS swapped L1's (TLEM 2.0) 39 knee-crossing elements for Rajagopal 2016 (the knee row of A + cj). The frames below LP-min went from 22.55 % to 18.22 %. The frontal arm Mkx stayed from L1, so it is a hybrid. the proof lane's night priority (research/COLLABORATOR_PRESENTATION_20260924/NIGHT_PRIORITIES.md §1): decompose, then test a consistent geometric model. Do not swap only a sagittal curve and interpret it as measured anatomy.

## Tasks (PREREG.md + sha256 BEFORE scoring)
1. **Decomposition with no implant fitting:** swap ONE muscle group at a time (quadriceps/patella chain, hamstrings, gastroc, gracilis/sartorius/TFL) and ONE quantity at a time (the moment arm in A vs the projection cj). Report per group/quantity the change in the share of frames below lo, per person and phase bin (5 bins). Which carries the 4.3 pp?
2. **Consistent geometry:** from the SAME Rajagopal path per muscle, take the sagittal arm, the frontal arm (the Mkx equivalent: the moment about the AP axis at the knee, for the medial/lateral split), and the projection cj, all in the same frame. Build the knee rows as a consistent whole; everything else stays as in L1. Also report the hip/ankle rows for bi-articular muscles, if they change consistency.
3. **Score:** the share below lo, meas−lo per person/phase, and the ε-band coverage at 10 % (CX-EPSBAND's eps_band.py) on the INVERSEOC mask. Compare L1, the L1ARMS hybrid, and the consistent geometry.
   Frozen criterion: the share below lo ≤ 11 % (half of 22.55) OR coverage at ε=10 % ≥ 75 %.
4. If CX-EMGLEARN's variant B (moment-calibrated EMG, implant-free) is in its RESULTS: run it on the consistent geometry as well, so that EMG and geometry are tested separately AND together.
5. Deliver `knee_geometry(model='rajagopal_consistent')` as a callable change + pytest.

Deliver RESULTS.md starting with `# CX-L1ARMS2`, results.json, and the scripts. Reuse CX-L1ARMS's OpenSim 4.6 extraction and inventory.json. Run the LP under bigmem.lock with 2 threads. Internal data stays local. lane_runner has full permissions in the workspace; `~/projects/bodytwin` is read-only. Every outcome is a node that expands.
