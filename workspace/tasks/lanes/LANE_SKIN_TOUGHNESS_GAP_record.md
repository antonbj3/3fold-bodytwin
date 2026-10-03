# LANE_SKIN_TOUGHNESS_GAP

Uppkommen ur nattens kirurgiska jakt (1/10). Resultatmapp `results/LANE_SKIN_TOUGHNESS_GAP/`.

## Luckan

- Tear toughness in skin: Pissarenko 2020, pig skin, mode III 20,6 ± 2,15 and mode I 30,4 ± 4,9 kJ/m², with per-sample rows including orientation (long/trans) and thickness (results/LANE_SURGICAL_INCISION/SOURCE_TARGETS_R2.json).
- Needle/penetration toughness: effective fracture–geometry product ≈ 1,8 kJ/m² (results/LANE_SURGICAL_INCISION/r3/geometry/needle_geometry_results.json).
- Sharp blade, bond level: Γ_cut ≈ 150–380 J/m² (synthetic scenario, results/LANE_SURGICAL_BINDINGS/PORTS_R4_FINAL.json).
- Microscopic budget: published interface energies (DEJ 3–6, stratum corneum 3–15, aortic media ~100–180 J/m²) and fibril/fibre inventory give at most ≈ 7 kJ/m² (results/LANE_SURGICAL_BINDINGS/INTERFACE_BUDGET_RESULTS_R4.json). Volume hysteresis requires an unreasonable zone 44–66 mm (BINDINGS r3).

**No known mechanism explains 20–30 kJ/m².** The capability: a mechanistic explanation that predicts tear toughness, penetration toughness and cutting toughness from the same tissue parameters — and thus how a scalpel incision, a needle and a tear injury differ.

## Candidate mechanisms (test, judge using independent data)

1. **Large-scale fibre bridging**: intact fibres span the crack wake over mm distances before breaking; G ≈ ∫σ_bridge(δ)dδ — order of magnitude MPa × mm = kJ/m². Predicts strong orientation dependence (long vs trans) and thickness dependence.
2. **Poroelastic dissipation**: skin is ~60–70 % water; fluid redistribution at the tip dissipates energy (as in hydrogels). Predicts rate dependence with characteristic time ~ L²/(k·E); dermal permeability is published.
3. **Viscoelastic bulk dissipation** in a macroscopic zone: predicts different rate and temperature dependence than (2).

## Do this

Break down into the smallest constituents and preregister for each mechanism which observables it predicts without fitting (orientation, thickness, rate, mode I/III ratio, penetration versus tearing). Acquire independent data: the Pissarenko rows, skin tensile curves and failure strain, fibre length/orientation (Yang 2015, Pissarenko Table 1), dermal permeability, tear toughness at different rates if available. Judge the mechanisms against data they are not fitted to; several may contribute — quantify shares with uncertainty. Deliver PORTS.json for INCISION/RESPONSE.

## Data boundary and control

Only published literature and our code. No internal data. Strongest control: phenomenological cohesive toughness calibrated per test (tearing, penetration). Gain = a mechanism (or combination) that predicts both and the orientation/thickness dependence without being fitted to them.
