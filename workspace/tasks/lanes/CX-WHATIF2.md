# CX-WHATIF2 — finish the parametric what-if: collect the OVH shards, what-if curves, derivatives, and the field lane's inertia derivatives

Read first:
- `results/CX-WHATIF/RESULTS.md` (A323: implicit KKT derivative vs FD 1.3e-9; the curves were missing), `RUN_ON_CLOUD.md`, `collect_cloud.py`, `FIELD_DERIVS_AVAILABLE.md`.
- `tasks/lanes/CX-WHATIF.md` (the original brief and criteria).

## State
- All 12 shards (identity, vsd_z001, vsd_z009 × CCD −15…15, AV −15…25) have run on OVH with EXIT 0. Output is in `/media/anton/sdc1-tmp/bodytwin/CX-WHATIF/cloud_bundle/modal_out/CXW-*/`.
- The field lane's femur mass-property derivatives are ACCEPTED after independent review:
  - `the public staging tree/3fold-motion-engine/_private/romi_collab/build/FX_FEMUR_MASSPROP_DERIV_20260924/RESULTS.md`;
  - z001 in `.../FD_REVIEW_FEMUR_Z001_20260924/RESULTS.md`.
  - Caveat: at ±10°, compute M(u) directly.

## Tasks
1. Run `collect_cloud.py` (adapted to the directory above). Require 141 KKT-approved steps per point and report any gaps.
2. What-if curves: hip reaction peak vs CCD and AV per individual. Compare the three methods:
   (a) full rebuild;
   (b) first-order extrapolation from the implicit derivative;
   (c) the emulator (BT-XF4).
   Report the valid interval (where (b)/(c) stay within 1 %/3 %) and the speed ratio (time per point).
3. Compare the direction and size against Kainz 2020 (DOI 10.1371/journal.pone.0235966, gait, not lifting — state the difference).
4. The field lane's inertia derivatives: the bone's contribution to the thigh segment's inertia → ID → hip force. Quantify the size of the effect (bone vs whole segment, cf. A181/N33). If it is < 1 % of the force effect, say so plainly.
5. The effect against uncertainty: compare the shape effect with N2b's 95 % band (A363: right hip −8.5/+15.9 % at σ 9.5 mm attachment error) and with the N7c band. Where is the what-if effect larger than the model's own uncertainty?
6. Update `results/CX-WHATIF2/RESULTS.md` (starting with `# CX-WHATIF2`), results.json, figures (PNG, what-if curves with bands), and pytest.

