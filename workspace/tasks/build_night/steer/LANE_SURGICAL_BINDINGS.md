# Steering LANE_SURGICAL_BINDINGS — round 4 (coordinator, 1/10 01:45)

Review of round 3: clearly falsified — with measured partial-cycle dissipation w_d ≈ 0,46 MJ/m³, Pissarenko's tear toughness requires a process zone h ≈ 44–66 mm (verified CHECKPOINT_R3.json), against skin thickness ~1,5 mm. Volume hysteresis therefore does not explain tear toughness; the remaining ~20 kJ/m² must lie in an unmeasured fracture step (delamination/interface work, fiber fracture after reorientation). Sharp blade: Γ_cut ≈ 150–380 J/m² (bond/contact level) holds as a port.

## Last round in this form — make it data-driven

1. Search for published **interface/delamination energies** for skin and dermis: dermal–epidermal junction peel (peel energy per width), interlaminar/delamination toughness in dermis or similar collagen laminates (tendon, fascia, amnion, pericardium), fiber extraction/pullout work for collagen fibers. Freeze the numbers that exist.
2. Test: G_tear ≈ Γ_cut + (number of delaminating interfaces per fracture surface) × G_interface. The number of interfaces must come from the microstructure (lamellar density/fiber layer structure in dermis), not be fitted. Falsifier: if it requires more interfaces than the dermis has, the delamination picture is wrong.
3. Finish the lane: write final PORTS.json (Γ_cut per layer with uncertainty, best G_tear explanation and its status, transition radius, thermal change) and an exact **measurement specification** for what is missing (which sample, which measurement, which resolution) so that the next step can be an actual measurement or data collection.

No new synthetic network/LP machine in this round.
