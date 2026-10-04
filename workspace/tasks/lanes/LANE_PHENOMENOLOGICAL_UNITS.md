# LANE_PHENOMENOLOGICAL_UNITS — which cells CANNOT receive an external measurement?

Result directory `results/LANE_PHENOMENOLOGICAL_UNITS/`.

## The finding that makes the lane necessary
The immune lane retrieved an external anchor — 33,18 million resident macrophages per mL from a published
source — and could not use it. The reason is structural and not a missing number:
**`source_M_unit = PHENOMENOLOGICAL M-unit; no physical count/volume normalizer`**, which gave
`density_to_M0_calibration_gate = FAIL_MISSING_FUNCTIONAL_NORMALIZATION` and an identified set
for functional M of the entire **[0, +∞)**. A model whose state variable lacks physical normalization
cannot receive a measurement of the same thing, however good the measurement is.

It is almost certainly not unique to the immune cell, and the question is measurable: **how many of the 43 cells
have at least one state variable or parameter without physical normalization?** Every such one is a cell
where an external reference cannot get in — thus a silent blocker for seed 1, which is about making more
questions askable.

## Do this
1. Go through `tasks/free48/sources/*/` and for every cell list its state variables and
   parameters with declared unit. Classify each: **PHYSICAL** (unit a measurement can have),
   **NORMALIZED** (dimensionless against a declared reference that is itself physical), or
   **PHENOMENOLOGICAL** (no route to a measurable unit).
2. For every PHENOMENOLOGICAL one: say which normalization would make it receptive, with unit. It is
   a specification, not a measurement — and it is cheap to write and costly to lack.
3. **Report the distribution before making any statement about the system.** The fraction of phenomenological variables
   per cell, and the number of cells with at least one. A single case is not a system statement.
4. Connect to tonight's measurement: 13,3 % of 5 984 swarm reports make an explicit claim about
   structural unidentifiability against 0,8 % complaining about parameter uncertainty. A phenomenological unit is
   an instance of the same thing — a direction the observables cannot reach. Say whether the two measures point to
   the same cells.

## Control and falsifier
- **Control:** the assumption that a cell that runs and passes its dimension check can also receive a
  measurement. The immune cell shows the two are not the same thing.
- **Falsifier:** if fewer than 10 % of the variables are phenomenological, the immune case is an exception and
  not a pattern — report that and close the lane.
- **Forbidden:** proposing a normalization constant that makes an existing anchor fit.

## Deliverable
`PORT.json`: table cell → variable → class → missing normalization with unit, plus the distribution.
