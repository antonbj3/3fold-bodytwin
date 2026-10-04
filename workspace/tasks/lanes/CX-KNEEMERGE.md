# CX-KNEEMERGE — ONE knee model: the existing high-resolution knee (KNEE-CELL/COMAK) with Field's fast cartilage contact as a switchable partial computation + force-controlled formulation

Anton (24/9): "so it's partial calculations inside the high-resolution one that can make it faster?" Yes. The existing knee stays the model. The new contact becomes a switchable sub-layer, used where it is valid.

## Build on (read first — do NOT build a new knee)
- The existing knee: `~/projects/bodytwin/scripts/msk/` (add_knee_ligaments.py, knee_ligament_config.json, knee_jw_transplant/, static_opt_knee.py, comak_*, knee_d1_*, residual_decomposition_knee_lts.py, nonuniform_pcsa_comak_correction.py), the graph node KNEE-CELL, and today's `results/L1`, `results/N12b`, `results/CX-EARLYSTANCE`.
- Write "Builds on" with file paths and what is REUSED. That is a requirement: a new knee from scratch counts as failed.
- The new module: `results/CX-CONNECTIVE/knee.py` plus the U384 adapter. Only the CONTACT LAYER and the force-controlled formulation are taken from it.
- Field's regime limits: `results/CX-CONNECTIVE/U384_REGIME_FROM_FIELD.md` (one layer, sphere, frictionless; the knee's two layers/elliptical contact/sliding is outside) and `GAUGE_RULE_FROM_FIELD.md`. Field's CLOUD-F-CART-KNEE (the knee case) may come in during the run: grep romi_collab/build/cloud_results/CLOUD-F-CART-KNEE.

## Tasks (PREREG.md + sha256 first)
1. Find where the existing knee's contact is computed (elastic foundation on meshes, or equivalent). Make it switchable:
   (a) the original;
   (b) U384 per compartment where the regime indicator (a/ℓ, t/R, ν, δ/R) is within the validated range;
   (c) automatic: (b) where valid, otherwise (a).
2. Same inputs, JW/DM Grand Challenge gait (the KNEE-CELL setup): compare (a), (b), (c) on contact force, medial/lateral split and time per step. Criterion: (c) within 5 % of (a) on medial/lateral and total force, with a speed-up reported honestly; state the fraction of steps where (b) is valid.
3. The force-controlled formulation, where the existing knee is displacement-controlled: implement it as an option and test it on the same data.
4. Counter-tests: (b) outside its regime must be flagged/rejected, never used silently.

Resources: lane_runner has full permissions in the workspace. `~/projects/bodytwin` is read-only (copy with source references). Internal data stays local or on OVH. Write in `results/CX-KNEEMERGE/`. `RESULTS.md` starting with `# CX-KNEEMERGE`, plus pytest.
