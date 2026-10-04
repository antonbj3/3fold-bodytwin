# Steering LANE_DEJ_BOND_INVENTORY — round 1 (coordinator, 2/10 19:35)

New lane. It takes over the slot from LANE_NEEDLE_FRICTION_ANCHOR, which closed after nine rounds with an orderable measurement (Cread, n=7) and two remaining mechanisms.

- **Obstacle:** k = Γ_II/Γ_I at human DEJ is unbound. The measurement path is closed — an acquisition job today found that no published measurement gives separation work against mode mixity, `negative_result: true`. The geometry path is structurally closed, because corrugation, interweaving and friction all add shear dissipation and can therefore only raise Γ_II.
- **Changed operation:** calculate the bound from below, from the bond inventory, instead of waiting for a measurement or looking for a geometric bound.
- **Strongest control:** our own energy balance over the blister cavity, 17,8–80 J/m² as an upper bound without mode resolution. The gain should be a bound on the RATIO.
- **Falsifier:** if the inventory gives no bound on k without additional unmeasured quantities — name the quantity and measurement that freezes it.

## Why the bound's direction is the whole point
A lower bound on k cannot falsify us; an **upper** one can. If any path gives k_max below 0,148, the mechanism dies for every rete-ridge wavelength we calculated (our switches lie at 0,148 / 0,227 / 0,288 / 0,334 / 0,368 for 100–300 µm). Geometry demonstrably cannot give that direction. The bond inventory can, if shear breaks the same bonds after sliding and rebinding — so aim there, and state clearly which direction your bound has.

## Two error classes that cost us today
1. **The wrong regime gives correct arithmetic on the wrong question.** One of today's biggest claims fell because a requirement was calculated in a saturated regime and compared against a value in a linear one. Before comparing two numbers: say which level and area convention each is in.
2. **Predictions are not measurements.** Another lane got an assignment from me based on my reading a field with *conditional* in its name as a measured value. Mark each number you use as MEASURED or DERIVED, and never mix them in the same ratio without saying so.

## If the file or path is missing where you run
The frozen microstructure numbers are in `results/LANE_SURGICAL_BINDINGS/SOURCES_R2.json`. If `local_path` does not exist on your host: write `missing_prerequisite` in the first paragraph of RESULTS.md and use the numbers as they stand in the lane brief, marked as quoted from the brief. Three jobs today drew the wrong conclusion by silently reconstructing missing inputs.

Everything PENDING_INDEPENDENT_REVIEW. No internal data.
