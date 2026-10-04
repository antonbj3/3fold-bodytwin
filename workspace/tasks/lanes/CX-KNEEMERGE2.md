# CX-KNEEMERGE2 — the COMAK knee with Field's CLOUD-F-CART-KNEE (two layers + elliptical contact): the blocker from A1380 is lifted

Read: `results/CX-KNEEMERGE/RESULTS.md` + `MOTION_PAPERS_NOTE.md` (A1380: U384 was not valid in COMAK, because it has two elastic implant layers), and `../3fold-motion-engine/_private/romi_collab/build/cloud_results/CLOUD-F-CART-KNEE/` (RESULTS.md, code/, adapter.regime_check; A1393).
- Validation: 64/64 within 5 % vs FE for two layers and elliptical contact.
- Regime: ρ 1–4, t_i/R_e 0.03–0.375, ν 0.1–0.49, E1/E2 0.1–10, t1/t2 0.5–2, δ/(t1+t2) 0.02–0.30, β_eff ≥ −0.27.
- DM (5+5 mm, E 46.3, ν 0.46) lies within the regime for t/ν/E/t1:t2.
## Tasks (PREREG.md + sha256 first; build on CX-KNEEMERGE's code, do not build a new knee)
1. Per frame of JW/DM gait (the KNEE-CELL setup): fit R_x/R_y per contact patch from the implant meshes and pose (medial/lateral), and compute δ. Run adapter.regime_check → the share of frames inside the regime.
2. Automatic mode: CART-KNEE inside the regime, the original Smith2018 EF outside. Compare against the original: total and medial/lateral contact force, time per step.
3. Against measured implant force (eTibia in lbf → ×4.4482 N/lbf; the medial/lateral split from the four load cells PM/AM/AL/PL according to CX-LUNGESYNC/build.py). Use gait only; F-8's jw_lungef1 is excluded.
4. Report whether the medial/lateral split becomes more correct against eTibia, and how much of the split comes from the contact law vs the external KAM (A1382/A1389 context).
Write in `results/CX-KNEEMERGE2/`. `RESULTS.md` starting with `# CX-KNEEMERGE2`, plus pytest.
