# LANE_FIBRIL_NETWORK_FRACTURE

The level below the entire surgical chain: a discrete network of fibrils and crosslinks that can break. Results folder `results/LANE_FIBRIL_NETWORK_FRACTURE/`. Form C in two steps, with B as method in steps 2.

## Why this lane now

Three lanes have stopped on the same kind of obstacle, and all three are pointing down towards the bindings:

1. `results/LANE_SKIN_TOUGHNESS_GAP/`: skin tear toughness 20–30 kJ/m² is matched in order of magnitude only by millimeter-scale fiber bridging (19–25 kJ/m²), but the mode splitting is unexplained and the fibril inventory ceiling is about 7 kJ/m².
2. `results/LANE_SURGICAL_INCISION/` R4: an efficient "seghet" does not move between tools (error 22–55 %) but stays within the same tool (≤ 2,8 %). So it is a tool-specific lump, not a material number.
3. `results/LANE_NEEDLE_FRICTION_ANCHOR/` R1–R3 (read them): the rate ratio 2,95 (1→80 mm/s) excluded viscous film (gives exactly 80 at all film thicknesses), poroelastic pore pressure (1,03–1,38) and a frozen transferred bond renewal law (1,455, under the gate 1,475). The direction difference UNI/NO 2,67 was dropped by elastic memory (1,001). Remaining candidate according to the lane's own diagnosis: **fiber recruitment and axial deformation in the dermis** — thus the material, not the interface.

Furthermore, Field's PROOF_LANE_AVLANKNING measured that the wavy surface at the layer boundary only gives 1,53× area, so the area of ​​the interlacing cannot explain the gap of ~10⁴ between the interface toughness and the skin tear toughness. Field owns the continuum level and consumes our Γ per layer (arranged 1/10 13:35).

## Capability

A bond-dissolved network that, from the fibrils' and crosslinks' own properties, gives: (a) tension hardening of the network, (b) tear energy per layer and its anisotropy, (c) the velocity dependence of force at large deformations via fiber recruitment. The same object must supply all three, otherwise it is a curve fitting per question.

## Do this, in order

1. **Step 1, mesh must pass a published collagen mesh before skin.** Build a discrete mesh (Mikado or Voronoi) with stretching and bending, without cross-link breakage. Reference: Licup et al. 2015 PNAS doi 10.1073/pnas.1504258112, strain-controlled solidification G(σ). Check the source and cite the figure. If the mesh does not handle the solidification deposit within factor 2: report it and fix the mesh before step 2.
2. **Step 2, breakage and anisotropy as held-out targets.** Add crosslink breakage with parameters from molecular and fibril level only (Gautieri 2011 doi 10.1021/nl103943u; Depalle 2015 JMBBM 52:1–13 — verify locator). Measure the J–R curve at a notch. Targets held: Annaidh 2012 JMBBM doi 10.1016/j.jmbbm.2011.08.016 (human skin, parallel to perpendicular Langer's lines) and Pissarenko's tear values ​​20–30 kJ/m². The length of the bridge zone should be millimeter scale **without being trimmed**.
3. **Step 3, rate dependence out of recruitment.** Fibers in dermis are curved (radius of curvature 6,56 µm, opening angle 32°, fibril diameter 82 ± 14 nm, fiber diameter 2,23 ± 0,96 µm — frozen in `results/LANE_SURGICAL_BINDINGS/SOURCES_R2.json`). The recruitment of straight fibers takes time towards the surrounding matrix. Predict the rate ratio 2,95 over 1→80 mm/s from the recruitment kinetics, with gate 1,475–5,9. Retrieve held-out points from `results/LANE_NEEDLE_FRICTION_ANCHOR/r1/CONTROL_COMPONENT_V1.json` (Barnett Fig11, Field digitization, semi-external). **Don't fart against 2,95.**

## Strongest control and falsifier

- **Control:** a continuum hyperelastic model with the same information (eg Ogden or Gasser–Ogden–Holzap error calibrated on the same tensile data), plus INCISION R4's per-tool toughness. The gain should be a guarantee the continuum model cannot provide: Γ per layer and mode splitting from bonds, not from a calibrated parameter.
- **Forger:** (a) step 1 misses solidification deposit by more than factor 2; (b) J is less than 7 kJ/m² or no anisotropy occurs; (c) brewing occurs only after per-sample weights; (d) step 3 ends up outside 1,475–5,9. In each case: name the term carrying the error and change a carrying condition in the same round. Don't end on a negative outcome.
- **Forbidden:** to calibrate bond parameters against tear toughness or against 2,95; to report a win against our own lineup as a result.

## Deliverable

`PORT.json` to Field with Γ per layer (interval), anisotropy ratio and stress hardening G(σ), as well as which lower level quantities each number is built from. Narrow follow-ups to the swarm in FOLLOWUPS.json with external_referent complete: locator, compared_quantity, our_value, command or excerpt.

No internal data. Anchors are held-out data, never model input. Everything PENDING_INDEPENDENT_REVIEW.
