# LANE_SURGICAL_BINDINGS

Anton's direction 30/9: high resolution where the scalpel opens — materials, cells, tissues, bonds. This lane owns **the bond and cell scale** and delivers ports to `LANE_SURGICAL_INCISION` (which owns cutting mechanics, gap and bleeding). Results folder `results/LANE_SURGICAL_BINDINGS/`.

## Desired capability

From molecular bonds to tissue fracture properties and damage zone, without free fitting per scale: collagen triple helix and cross-links (enzymatic/AGE), fibril sliding and recruitment, fibre-network orientation → fracture toughness G_c, failure strain and anisotropy for dermis, fascia and vessel wall; cell-membrane rupture (strain, pressure, temperature) → cell-damage zone at the blade; thermal denaturation in electrosurgery (Arrhenius damage) → coagulation zone. The same parameters must predict several observables (toughness, damage zone, denaturation zone) — that is the cross-prediction a control calibrated per observable does not give.

## Existing code (read-only, build on it)

`~/projects/bodytwin/scripts/msk/{collagen_triple_helix_thermal_stability.py, tendon_collagen_hierarchical_mechanics.py, myofascial_transmission.py, skin_pulp_mechanics.py, bone_fracture_toughness_lefm.py}` and corresponding `data/msk_results/*.json`; `tasks/free48/sources/SURG_COLLAGEN/`, `SURG_INCISION/`; `results/BT-HX-Q049/` (thin layers), `results/BT-HX-Q036/` (ischaemic state vector). Graph nodes: MODEL-COLLAGEN-TRIPLE-HELIX-THERMAL-STABILITY, MODEL-TENDON-COLLAGEN-MECHANICS, SKIN-DERMAL-STIFFNESS-COLLAGEN-COMPOSITION, MSK-FASCIA-NETWORK. Graph target: `BT-CTX-SURG-COLLAGEN` (kind define/review).

## Portar till INCISION-lanen

Write `results/LANE_SURGICAL_BINDINGS/PORTS.json`: per layer G_c (J/m²), failure strain, anisotropy ratio, cell-damage threshold, denaturation kinetics (A, Ea), with uncertainty and source. Read `results/LANE_SURGICAL_INCISION/` for what it needs.

## References, data boundary, control

Published measurements (brottseghet hud/dermis, kollagen-denaturering DSC/Arrhenius, cell survival vs elongation/temperature). Pre-register hit targets. No internal data (Nothing from `BT-DAT-Q034`, `BT-DAT-Q050`, the collaborator, Grand Challenge, restricted model data, external musculoskeletal solver). Strongest control: phenomenological tissue model calibrated separately per observable with the same data, full cost.
