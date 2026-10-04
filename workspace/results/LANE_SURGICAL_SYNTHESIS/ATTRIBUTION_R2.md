# Misattribution for the surgical chain, R2

PENDING_INDEPENDENT_REVIEW. The same nine conserved targets are used; RMSE is counted on the seven lines with `independent_day=true`, with replicas remaining. ΔRMSE = 43,703696 − new RMSE; positive number is improvement. No parameter has been customized in R2.

The isotropic R1 observer is exactly `S=75 A m`, where `A=trace(U+I+M)` and `m=trace(M)/A`. R3 uses `S=75 C_R3 x_HP`. R3's HP curve is external input and C_R3 a model state that inherits force calibration from RESPONSE R1. The D/H/A reservoir in synthesis R1 has no power consumer.

| Ersatt port | RMSE pp | ΔRMSE pp |
|---|---:|---:|
| None: original chain | 43.703696 | +0.000000 |
| Quantity A only → R3 C | 13.000632 | +30.703064 |
| Enbart mognad m → R3 HP-normalisering | 39.456781 | +4.246914 |
| Both A and m → R3 | 8.306367 | +35.397329 |
| Quantity A → chain's own FV C | 29.560901 | +14.142795 |
| m → isolated D/H/A-HP with R3 scale | 47.268908 | -3.565212 |
| Bridge → R3 referens1 | 43.703696 | +0.000000 |
| Riktning → isotrop R3 referens | 43.703696 | +0.000000 |
| Tidig state/hazard → R3 | 45.697008 | -1.993312 |
| All six upstream R3 ports | 60.007651 | -16.303955 |
| gap → R3 | 52.892706 | -9.189010 |
| biological_width → R3 | 43.004371 | +0.699325 |
| perfusion_width → R3 | 45.734780 | -2.031084 |
| face_conductance → R3 | 45.366986 | -1.663290 |
| k_deposit → R3 | 59.827729 | -16.124033 |
| k_mature_legacy → R3 | 43.703696 | +0.000000 |

The quantity gate carries 30.926739 pp of the total 35.397329 pp two-factor Shapley enhancement; ripe berry 4.470590 pp. The interaction in RMSE is -0.447351 pp. Shapley here is descriptive misattribution, not a biological causal effect. Absolute port changes cannot be equated with independently acquired information.

The exact obstacle is in the coupling of the fixtures. FV C is formed with `k_deposit F gate (1-C)` and lost with `.008 macrophage C`. U/I/M forms simultaneously with `birth_scale F gate (1-C_FV)` and further loses `turnover_UI (U+I)`. Internal maturation fluxes cancel each other out, but the total U/I/M pool is never forced to FV C. When FV C fills, the U/I/M's birthgate closes even if its own inventory is empty.

Dag90 is the chain's FV C=0,640170, U/I/M sum=0,313012 and R3 C=0,998320. Just importing R3's k_deposit≈5/day to the chain drives C to≈1 and the U/I/M sum till0,009881, which degrades RMSE till59,827729 pp. All upstream R3 ports simultaneously reproduce R3's FV C at the saved times, but the remaining U/I/M observer provides endast0,620360% strength dag90. This distinguishes a misconnected reference set from an early O2 error.

R3 calibration ancestry: `surgical_chain/vendor/response_r1.py:main` uses `least_squares` vs older styrka3% dag7 och20% dag21 and get k_deposit≈5/day. R3 did not do **new** strengthfit; its quantity port is nevertheless not independent force-free information. Synthesis's k_deposit=0,1 is another synthetic default choice. The statement "same chemistry input" thus does not describe this comparison.

The port change is a run diagnostic reconstruction, not a changed biological model. Native traction/fracture law remains null, both proxy branches retain previous FAIL and strongest equally informed port/FV-/cohort control is TIE. Data: [ATTRIBUTION.json](r2/attribution_v1/ATTRIBUTION.json), [raw curves](r2/attribution_v1/ATTRIBUTION_CURVES.npz), [inventory diagnosis](r2/attribution_v1/INVENTORY_DIAGNOSIS.json).
