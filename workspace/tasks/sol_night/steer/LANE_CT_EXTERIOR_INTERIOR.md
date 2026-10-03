# Steering LANE_CT_EXTERIOR_INTERIOR — round 2 (the coordinator, 1/10 04:00)

Review of round 1: important partial result — in 300 manually VAT-annotated people (AATTCT-IDS), pixel area gives R² 0,64 but **dimensionless shape → VAT fraction only R² 0,055** (verified MANUAL_ANCHOR_STATISTICS_R1.json). The node's ceiling ~0,65 thus appears to be body size, not shape. CT pilot (n=29): shape minus waist R² −0,04. And your scale stress test (shuffled labels R² −0,005 → 0,81 through a shared coordinate scale) is a trap that must remain in every future comparison.

## Round 2: do the size-versus-shape decomposition on physical CT

TotalSegmentator has physical voxel size in the NIfTI header (1,5 mm). HU threshold −190…−30 inside the abdominal wall at L3 is an established VAT method validated in publications (cite validation studies against manual segmentation and state their agreement as label uncertainty) — use it as a label with declared uncertainty rather than waiting for manual masks.

1. Scale from 29 to so many CT med hel L3 in the picture that the budget is sufficient (screena fler av de 1 228; stream from zip, delete extracts after features; hold /mnt/games-240-anvthe ending under ~40 GB totalt).
2. Preregister and report three numbers with 95 % CI: R²(size → VAT_L3 cm²), R²(size + dimensionless shape → VAT_L3), and R²(dimensionless shape → VAT fraction of the cross-section). The third is the key one: compare with AATTCT 0,055.
3. Controls as in round 1: waist/body area, shuffled labels (should give ≈ 0 with physical scale), cases without contrast, pathology mix.
4. Write in NODE_PROPOSALS a concrete proposal for a new formulation of MSK-VISCERAL-FAT-IDENTIFIABILITY (size part vs shape part) — in the results folder, not in the graph.
