# Shared collagen reference contract R3

PENDING_INDEPENDENT_REVIEW. `shared_reference_r3` is a synthetic forward model.

| Port/state | Unit and reference | Operation/source | Status |
|---|---|---|---|
| C and trace(Q_U+Q_I+Q_M) | intact collagen mass per fixed registered reference volume, dimensionless | Same initial0, B and lambda in every cell; RESPONSE FV and R3 tensor operator | SYNTHETIC |
| B | reference mass/day | k_deposit F h(1−C); R3≈5/day is inherited early strength calibration | SYNTHETIC |
| lambda | 1/day | .008 macrophage; equal proportional removal of U/I/M and D/H/A | SYNTHETIC |
| U→I→M | reference mass/day | k1 h U and k2 h I; .1/.05 per day, internal fluxes cancel | SYNTHETIC |
| D/H/A | mol/mol same reference collagen | New collagen gets s=.102 precursor-site equivalents; D→2H with eta or A; k from R4 cell analog | SYNTHETIC transfer |
| D+2H+A | mol site equivalents/mol reference collagen | s C: chemical site equivalents are different units from collagen mass | DERIVED under declared closures |
| Strength proxy | %intact | 150 bridge eᵀQ_Me; isotropic75% reference | SYNTHETIC |
| Chemical strength diagnostic | %intact | 75 C clip((H/C−.003)/.040); scales from R3, no external HP curve | SYNTHETIC, no native law |
| Native traction/rupture and joint posterior | Pa / unknown | Requires matched connected chemistry×age×angle/anchors | UNKNOWN/null |

`birth_scale` and `turnover_U_I` are explicitly unused historical ports in R3. No extra .03/day loss is added. C is not projected separately from marks; an inadmissible state is rejected with a requirement for smaller dt. R3 updates the same O2/FV state as before; the U/I/M partition does not govern oxygen beyond its shared sum. The isolated fixedC10 reservoir is preserved only as an older diagnostic and is not a second collagen amount in the feedback chain.

The R3 checkpoint stores reference version, full FVstate, Q, chemical marks, born/removed/initial, hazard/earlystate, grid, clock and laws/config. The R1 checkpoint belongs to a different reference and is rejected by the new solver operator. No history conversion is assumed; the entire prefix has been rebuilt and the cost reported.
