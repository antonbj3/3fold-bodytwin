# LANE_DEPTH_TEMPERATURE — resolve the depth the surgeon cannot measure

Result directory `results/LANE_DEPTH_TEMPERATURE/`.

## Why this lane exists, measured
`LANE_LASER_SURGERY` r18 ran the sufficiency test on the one thermal quantity that is actually measured
during an operation, the surface temperature, and it failed:

| quantity | value |
|---|---|
| surface temperature history identity error | **0.0 K** — identical, not approximately |
| temperature gap at target depth | **0.3806 K** |
| numerical control, finest finite volume | 4.33e-05 K |

Two cases whose surface history is identical to machine precision differ by 0.38 K where the damage
happens, and the numerics are three orders of magnitude cleaner than the gap, so it is physical. The
laser chain separately sits **26.5 % above** a published ablation threshold of 1.15 J/cm², and its own
diagnosis was a mismatch between the source convention and the observation definition.

So the thermal state at depth is a stress point by the measured definition: its value decides a verdict
we already hold, and the quantity we can observe does not determine it.

## Do this
1. **Build the depth-resolved thermal field as the primary output**, not surface temperature with a
   depth correction. Report the temperature at the depth where the damage criterion is evaluated, with
   the grid spacing that was needed and the convergence order.
2. **Find the coarsest sufficient depth resolution.** The resolution atlas measured that of 108
   quantities only 3 require finer resolution, so do not assume this one does — prove it. Report the
   coarsest spacing at which the damage verdict does not change, and the factor against what the laser
   lane runs today.
3. **Then close the surface-to-depth question honestly.** Given the surface history alone, what is the
   tightest interval on the depth temperature? If that interval is wider than the damage criterion's
   margin, the surgical readout cannot be surface temperature and must be something else — say what,
   with its unit.
4. **Use the thermal-confound records rather than arguing from the model.** The literature collection has
   25 records whose property names are exactly this question, 13 of them PMID-asserted against the XML.
   Read `tasks/build_night/LIT_REFS_FILTER.md` first; two files are excluded entirely and the rest are read
   with a row filter, and you must report how many records the filter removed.

## Control and falsifier
- **Control:** the surface-temperature readout as the laser lane uses it today, on the same cases. Equally
  informed by construction — same delivered dose, same tissue, only the readout differs.
- **Falsifier:** if the depth temperature stays inside the damage criterion's margin across every
  admissible surface history, then surface temperature IS sufficient in practice and the 0.38 K gap is
  below the decision threshold. That would be a reassuring result and must be reported plainly, not
  softened — it would also mean this lane can close after one round.
- **Forbidden:** reporting a depth temperature without the grid spacing it was computed on; treating the
  published 1.15 J/cm² threshold as validated for our tissue without its measurement definition;
  consuming a filtered record.

## Delivery
`PORT.json` with the depth temperature per case, the convergence table and the coarsest sufficient
spacing, the tightest depth interval implied by surface history alone, the comparison against the damage
margin, and the count of records the filter removed.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
