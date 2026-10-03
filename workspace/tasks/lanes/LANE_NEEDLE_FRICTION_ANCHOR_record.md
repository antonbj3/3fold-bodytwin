# LANE_NEEDLE_FRICTION_ANCHOR

The surgical chain's sharpest open question (1/10 13:15). Results directory `results/LANE_NEEDLE_FRICTION_ANCHOR/`. Form C (ground truth that does not come from us), with B as the method.

## Status

- `results/LANE_SURGICAL_INCISION/` R3–R4: held-out needle geometry gives 11 % central error. Effective "needle toughness" does not transfer between tools (error 22–55 %) but holds within the same tool across velocities (≤ 2,8 %). It is thus a tool-specific lump, not a material constant. MEASUREMENT_SPEC_R4.md describes the measurement that would resolve the lump.
- The field lane (SOL_FALT_NALKRAFT_20261001, port in `results/LANE_SURGICAL_SYNTHESIS/EXTERNAL_PORT_PROVIDERS.json`): force = fracture energy × incision perimeter + friction + deformation, with shared tool parameters. Maximum error 31 %, and 26 % with velocity-dependent friction, against our per-tool control 25 %. The supporting assumption is that contact and friction are the same on reinsertion. Barnett 2016 assumes it, and none of us has tested it.

## Capability

Split the needle's insertion force into fracture, friction and deformation with **friction measured independently**: retraction force, reinsertion into a prepunctured hole or same-hole passage. Fracture energy then becomes a residual term that can be compared with an independent fracture measurement, and tool dependence should either disappear or be localized to a named geometry term.

## Do this

1. **Data hunt first**: published force–depth curves or tables with insertion, retraction and/or prepunctured holes in skin (preferably human or pig; otherwise named tissue, marked as cross-species). Known starting points to check, not assume: Okamura, Simone & O'Leary 2004 (IEEE TBME, prepunctured liver), Mahvash & Dupont 2010, Barnett 2016 and needle insertion reviews (Abolhassani 2007, van Gerwen 2012). Freeze numbers with table or figure, species, site, needle geometry and velocity. Check errata.
2. **Smallest constituents**: friction force = contact pressure × contact area × friction coefficient (or per length). Contact pressure comes from the hole's elastic recovery around the needle (cavity in hyperelastic material). Derive friction per length from tissue stiffness and hole radius/needle radius, instead of fitting it. Compare with measured retraction or prepuncture force.
3. **Preregister** before comparison: (a) friction per length from stiffness and geometry against measured; (b) fracture energy = insertion force − friction − deformation, against independent fracture toughness for the same tissue (`results/LANE_SKIN_TOUGHNESS_GAP/SOURCES_R1_FINAL.json`, mode I/III); (c) whether tool dependence in INCISION R4 (error 22–55 %) shrinks when friction is removed.
4. **Strongest control**: INCISION's per-tool toughness and Field's decomposition with fitted friction. Gain = an independent friction anchor that makes the fracture term shared across tools. Reject the claim if the friction anchor does not exist for skin: say it exactly and run the nearest tissue as a cross-species case.

No internal data. Anchors are held-out data, never model input. Everything PENDING_INDEPENDENT_REVIEW.
