# CX-EPSBAND — a certified knee-contact band: ε-optimal set ∩ LP set ∩ C0 correction, with hamstrings as the free axis (BC-E24)

## Background (read these; do not rebuild them)
- A1420 `results/CX-INVERSEOC/epsilon.py` + `epsilon_results.json`: inside the set, the measured force lies at median 0.9–3.8 % above minimal stress-2 effort. 68–91 % of frames are within ε=10 %.
- A1418 `results/CX-MINCONTACT`: at the LP minimum, rectus femoris carries 54–66 % of the moment; meas−lo is phase-dependent.
- Field U447 (`the public staging tree/3fold-motion-engine/_private/romi_collab/build/U447/RESULTS.md`, read-only):
  - Δf from the LP min to the measured force is median 105 N, 88 % of it IN the contact gauge, carried by HAMSTRINGS (energy 0.45; gastroc 0.03, quad 0.08).
  - Frames below Fmin: a local C0 lowering (median −312 N) repairs 92/106; a contact-arm change of 3.9 % also works; widening F0 repairs only 13/106.
- L1 operator: `results/L1/prep/*.npz` (A, b, F0, cj, C0, medial/lateral), `results/CX-SETVALUED/run_gc.py`.
- CX-FUSION is running in parallel (posterior with EMG). Do not duplicate it; this lane is the constraint/certificate variant with no EMG.

## Tasks (PREREG.md + sha256 BEFORE comparing against the implant force)
1. **C0 correction without facit per frame.** Find out what C0 physically consists of in L1 (external load, gravity, or the contact normal/tibial-axis projection). Then build a correction from INDEPENDENT data:
   - the contact-normal direction from the tibial axis in the static trial / the implant geometry, or a person-level tilt;
   - OR ONE parameter per person fitted LOPO-free on other persons' data.
   Report how much of the −312 N follows from the geometry alone.
2. **The band:** per frame, [min, max] of the contact over {f : LP equilibrium, 0≤f≤F0, med/lat≥0, stress-2(f) ≤ (1+ε)·min} with the corrected C0. Use ε fixed at 5 %, 10 % and 20 % (frozen). Hamstrings is the free direction: report the band's endpoints in terms of hamstring force.
3. **Test on DM/SC/PS/JW gait (the INVERSEOC mask) + JW4 non-gait** (not jw_lungef1):
   - coverage of the measured force per person (target ≥ 90 % at ε=10 %);
   - median band width in BW (target ≤ 0.8 BW);
   - the band midpoint's RMSE against N1g (L1 LOPO 2.1153).
   Also run without the C0 correction, to isolate its effect.
4. **Counter-test:** the time-shifted measured force must give clearly lower coverage.

Deliver RESULTS.md starting with `# CX-EPSBAND`, results.json, code (`eps_band(trial, eps)` in bodytwin_core-compatible form), and pytest. Internal data stays local; LP/QP runs go under bigmem.lock or on OVH via `tasks/cloud_run.sh` (≤ 12 vCPU, finish by 07:30). Every outcome is a node that expands (no verdict words). lane runner has full permissions in the workspace; `~/projects/bodytwin` and romi_collab are read-only.

## Addition (Field): U447's minimum-norm-C0 per frame (COMPARISON only, not model/training)
the public staging tree/3fold-motion-engine/_private/romi_collab/build/U447/raw/repair_run1.json (nyckel rows), sampled_run1.json (muskelprojektioner, gaugebas), summary.json. Compare your independent geometric C0-correction to these; never use them to adapt.

## Addition (Field FB_HUNT_EPS_BAND_WIDTH, gait2392, unreviewed)
Band width in stance ~√ε: 199 N (1 %), ~600 N (10 %), 881 N (25 %); swing wider. Report your own width–ε curve (log-log slope) on L1 and compare; report stance and swing separately.

## Addition (Field U460)
Field certifierar alla 565 rutor (LP-band + ε 5/10/20 %both: med/lat-varianterna) → the public staging tree/3fold-motion-engine/_private/romi_collab/build/U460/raw/coverage.json. If it exists when scoring: compare your coverage/width per person and ε mot den (comparison only), redovisa avvikelser.
