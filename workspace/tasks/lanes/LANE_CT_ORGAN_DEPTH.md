# LANE_CT_ORGAN_DEPTH

The second half of Anton's CT seed (X bookmarks batch2 B1 / batch3 TOP5 #3): `ORGAN-DEPTH-VISCERAL` lacks an organ-imaging cohort (spleen volume varies 112–215 cm³ between cohorts). Results folder `results/LANE_CT_ORGAN_DEPTH/`.

## Underlag

TotalSegmentator v2 is already downloaded and verified on `/mnt/games-240/research/bodytwin_ct/` (read `results/LANE_CT_EXTERIOR_INTERIOR/` — pipelinen, screening av 1 228 CT, 249 kompletta L3- profiles, the shell trap and the physical storleks/form-faktoriseringen)Organ masks are in the dataset.

## Uppgift

1. Organ volumes (liver, spleen, kidneys, pancreas, heart where in the image) and organ depths (distance from skin to the organ's surface/centroid along AP and laterally) in physical units for all CTs where the organ is fully in the image.
2. Variance explained by body size s (same size mode as EXTERIOR_INTERIOR) and by dimensionless shape, per organ, with 95 % CI and shuffled-label control (should give ≈ 0).
3. Compare the distributions with published reference cohorts (e.g. the spleen volumes 112–215 cm³) and report the cohort's pathology mix.
4. NODE_PROPOSALS: proposed numbers for ORGAN-DEPTH-VISCERAL (the results folder, not the graph).

Game Data-240 (aldrig / eller /mnt/shared_data). Stream from zip, delete extracts after features. Public data, no internal.
