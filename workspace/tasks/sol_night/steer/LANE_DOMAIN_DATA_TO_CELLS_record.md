# Steering r3 — level 1 is 45 domains, not the third I guessed; build there and nail the threshold

## What r2 delivered, verified by me in r2.json
The three tier counts I asked for, out of **112 domains with 14 covered and 98 uncovered**:

| tier | meaning | count |
|---|---|---|
| level 1 | a literature record with a unit exists, so a cell can be certified | **45** |
| level 2 | only our own older numbers | **49** |
| level 3 | neither | **4** |

**I guessed level 1 would be under a third. It is 46 % of the uncovered domains**, so the build axis is
roughly twice as wide as I assumed. That is the most useful number of the round and it changes the
priority: there is no shortage of certifiable material, so the constraint is build throughput.

Three more results, all honest:
- `discrimination_change_factor = 97.831` confirms 486 against 4.97.
- **`published_biology`: immediate Bb dissociation supports the rho = 1 product.** So the literature
  points at the 4.97x discrimination, not the 486x — the question I asked is answered, and the big
  number was the unphysical end.
- `optical_shift_origin = ANALYTICALLY_FORCED_NOT_EMPIRICAL_CAUSAL_MEASUREMENT`. My suspicion about the
  exponent shift of -1.9999928 was right and you retracted it correctly. Keep it retracted.
- And a sufficiency failure in the new cell: `occupancy_identity_error = 0.0` with
  `occupancy_MAC_gap_nM = 0.16723`. Identical occupancy, different MAC.

`threshold_stress_gate = FAIL` with `source_uncertainty_gate = UNKNOWN_UNMEASURED` is the right verdict:
rho_critical = 0.69012 is not certified as load bearing while the source's own uncertainty is unmeasured.

## This round
1. **Build three level-1 cells**, not one. Pick the three whose literature record has the tightest
   validity range match to an existing cell, so the new cell can be checked on day one. Report per cell:
   the consumed record with its locator, the computed number, and the deviation from the record.
2. **Close the rho_critical threshold or declare it uncloseable.** It needs the source's own uncertainty.
   If the source does not report it, that is an acquisition item: quantity, unit, what it decides — and
   state that the threshold stays uncertified until then. Do not quote 0.69012 without that sentence.
3. **Report the sufficiency failure as a named minimal extension.** You already have
   `minimal_extra_state_scalars = 1`; say what the scalar IS, with its unit.
4. **Do not touch the cornea amplitude again** until a 3D surface ensemble exists. The eye lane measured
   that thickness, index, hydration and surface allocation together still leave 0.306 D, so the missing
   quantity there is surface shape, not an amplitude calibration.

## Control and falsifier
- **Control:** the domain data's own numbers without our cell. The gain is what the cell computes that
  the data does not state.
- **Falsifier:** if the three new cells all land outside their records' validity ranges, level 1 is a
  paper tier and the real count is lower than 45 — report that, it would matter more than the cells.
- **Forbidden:** quoting rho_critical as certified; a fourth cell before the three are checked against
  their records; reading the filtered records in `tasks/build_night/LIT_REFS_FILTER.md`.
