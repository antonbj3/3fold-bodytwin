# LANE_SURGICAL_RESPONSE

Anton’s direction 30/9: high resolution where the scalpel opens — materials, cells, tissues, bonds. This lane owns **the body’s response after the incision over time**: hemostasis → ischemic edge → inflammation and cell recruitment → matrix deposition → scar mechanics. It consumes ports from `LANE_SURGICAL_INCISION` (gap, injury zone, severed vessels) and `LANE_SURGICAL_BINDINGS` (collagen, cell injury thresholds). Results directory `results/LANE_SURGICAL_RESPONSE/`.

## Desired capability

Given the incision geometry and injury zone: predict bleeding until hemostasis, the oxygen/perfusion field at the wound edge (which cells become ischemic), inflammatory and migrating cell populations over time, collagen deposition and scar stiffness/strength at day 7/14/28/90 — with high spatial resolution only at the wound edge and a coarse body tissue field around it. The same kinetic parameters should predict several observables (hemostasis time, wound strength over time, edge necrosis).

## Existing code (read-only, build on it)

`~/projects/bodytwin/scripts/msk/{coagulation_hemostasis.py, platelet_hemostasis.py, wound_healing_cascade.py, acute_phase_inflammation.py}` and `data/msk_smoketest/wound_healing_cascade/`; `results/BT-HX-Q036/` (ischemia/reperfusion: minimal state vector q, o, A, C, R), `results/BT-HX-Q033/`; `tasks/free48/sources/SURG_HEMOSTASIS/`, `SURG_HEALING/`. Graph nodes: WOUND-HEALING-CASCADE, MODEL-COAGULATION-THROMBIN-GENERATION, MODEL-COAGULATION-FIBRINOLYSIS, MODEL-PLATELET-ACTIVATION-AGGREGATION, HEMATO-HEMOSTATIC-BALANCE, MODEL-NLRP3-INFLAMMASOME. Graph goals: `BT-CTX-SURG-HEALING` and `BT-CTX-SURG-HEMOSTASIS` (kind define/review).

## References, data boundary, control

Published measurements: bleeding time, wound break strength over time (klassiska tensile-strength-kurvor), oxygen stress at wound edge, cell kinetics (neutrofiler/makrofager/fibroblaster). Pre-register hit targets. No internal data (Nothing from `BT-DAT-Q034`, `BT-DAT-Q050`, the collaborator, Grand Challenge, restricted model data, external musculoskeletal solver). Starkaste kontroll: befintliga fenomenologiska kurvor/modeller calibrated per observable with the same data, full cost. The disks are almost full: large intermediate data on `external_mount<LANE>/` bara om < 1 GB, otherwise not at all.
