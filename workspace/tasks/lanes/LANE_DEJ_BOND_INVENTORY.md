# LANE_DEJ_BOND_INVENTORY

Derive a bound for the dermo-epidermal junction's fracture toughness from the molecular inventory, where no measurement exists. Result folder `results/LANE_DEJ_BOND_INVENTORY/`. Form A.

## Why this lane, and why now

An acquisition job today asked whether any published measurement gives the DEJ's separation work as a function of mode mixity. The answer was **NO**, with `negative_result: true`. The consequence is stronger than "uncertain": the ratio k = Γ_II/Γ_I at the human DEJ is **not just undetermined but entirely unbounded by measurement**, and a parallel track showed that k decides the SIGN of our prediction of the age dependence of blister times. Another lane also showed that geometry cannot give a lower bound, because corrugation, interweaving and friction all ADD dissipation in shear and therefore can only raise Γ_II.

So the measurement route is closed for now and the geometry route is structurally closed. What remains is to calculate the bound **from below**, from the bonds that actually hold the interface. It is precisely the kind of task Sol shall have: a derivation where no measurement exists, not cataloging.

## What is already decided and shall not be redone

- Our own energy balance over a blister cavity gives an **UPPER** bound for Γ_i of 17,8–80 J/m², median 28. It is an upper bound because it ignores stored elastic strain energy, and it does not resolve mode mixity.
- **Peel force per width is not fracture energy** but a work balance. The low published DEJ values 1,47 and 2,84 N/m also concern CULTURED skin and do not identify the native interface.
- Elastic straightening and rotation of fibers is **recoverable** and must never be accounted as dissipation. It was my own error in earlier direction and it was rejected.
- The wavy surface gives only approximately 1,53× extra area, so area cannot explain the gap between interface toughness and the skin's tear toughness 20–30 kJ/m².
- A fibril-network lane has frozen microstructure numbers in `results/LANE_SURGICAL_BINDINGS/SOURCES_R2.json`: fibril diameter 82 ± 14 nm, fiber diameter 2,23 ± 0,96 µm, curvature radius 6,56 µm, opening angle 32°. Use them, do not measure them again.

## Do this

1. **Inventory the bonds that hold the DEJ**, with a locator per item: the collagen IV network, laminin-332, nidogen, perlecan, integrin α6β4 and collagen VII's anchoring fibrils. For each: areal density (count per µm²) and the energy required to separate it, with the type of dispersion named. Accept adjacent systems and **declare cross-species transfer explicitly** — better a declared transfer than an omitted number.
2. **Calculate Γ_i from below** as the sum of bond energy times areal density, and compare with our own upper bound 17,8–80 J/m². Two outcomes are both informative: if the sum calculated from below is far BELOW the upper bound, the interface carries more than its bonds, and then there is a dissipation mechanism beyond bond rupture that shall be named. If it is close the inventory is sufficient.
3. **The bound on mode mixity, which is the actual goal.** Separation work in opening breaks bonds normally; in shear the same bonds may break after sliding and rebinding. Derive what the inventory allows for the ratio k = Γ_II/Γ_I, and state whether it gives an **upper** bound, which the geometry route demonstrably cannot. An upper bound on k is what decides the sign of our age prediction, and no other route to it is known today.
4. **Forbid yourself to reach 20–30 kJ/m².** That number is the skin's tear toughness at millimeter scale and belongs to another level. If your inventory comes close something has gone wrong with the area convention.

## Strongest control and falsifier

- **Control:** our own energy balance over the blister cavity, thus 17,8–80 J/m² as an upper bound without mode resolution. The gain shall be a bound on the RATIO, not a better point value for Γ_i.
- **Falsifier:** if the bond inventory cannot give any bound on k without additional unmeasured quantities, say exactly which quantity is missing and which measurement would freeze it. It is a fully adequate outcome and better than a bound that rests on an assumption.
- **Forbidden:** accounting elastic straightening as dissipation; comparing a bond-level number directly against 20–30 kJ/m²; treating cultured skin as native DEJ; converting N/m to J/m² without a work balance.

## Delivery

`PORT.json` with Γ_i from below as an interval, the bound on k with its direction (upper or lower) and what it rests on, and a list of missing quantities with the measurement that would freeze each one. Narrow follow-ups to the swarm in FOLLOWUPS.json with external_referent complete.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
