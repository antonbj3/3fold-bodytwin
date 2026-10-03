BT-HX-Q168

## What was built

A four-state, first-principles synthetic tracer model separates precursor/microbial TMA production, epithelial transport, portal FMO3 conversion, and systemic TMAO loss. With concentrations in `µM`, the core fluxes are

- `J_m = Vmax_microbe*g/(K_microbe+g)` (`µM h⁻¹`)
- `J_b = P_epithelial*A_epithelial*(l-p)*l/(K_barrier+l)` (`µmol h⁻¹`, Fick-type effective transport with saturation)
- `J_f = Vmax_FMO3*p/(K_FMO3+p)` (`µM h⁻¹`)

The balances, urine observation `Q_urine_dot = V_s*k_renal*o`, parameter table, units, and provenance are in `model.py`; all parameter values are synthetic assumptions except the published challenge dose and reference anchor.

## Prediction versus published reference

Koeth et al. (2013), *Nature Medicine*, DOI `10.1038/nm.3145`, reported median fasting plasma TMAO **4.6 µM** in `n=2,595`; the page consulted was https://pmc.ncbi.nlm.nih.gov/articles/PMC3650111/ (Figure 4f and accompanying text). The same paper specifies the 250 mg d3-(methyl)-L-carnitine challenge and serial plasma/24 h urine collection. The source was found and is **verified**; it is not marked `OVERIFIERAD`.

| Quantity | Model | Reference | Difference/status |
|---|---:|---:|---|
| steady-state plasma TMAO | 4.52349 µM | 4.6 µM | absolute error 0.07651 µM; ratio 0.98337; factor-2 interval 2.3–9.2 µM |
| tracer TMAO AUC, 0–24 h | 76.1538 µM·h | no matching published AUC found | synthetic output, not a measured comparison |
| tracer TMAO peak | 4.32932 µM | no matching published peak found | synthetic output |
| 24 h urinary tracer loss | 19.0384 µmol | no matching published amount found | synthetic output |

The numerical reference criterion is **UPFYLLT**. The AUC/peak values are not presented as validation because the reference estimand is total fasting TMAO, not the synthetic labelled-tracer AUC.

## Counterfactual mechanisms

| Model condition | Tracer AUC (`µM·h`) | Peak (`µM`) | 24 h urine (`µmol`) |
|---|---:|---:|---:|
| full model | 76.1538 | 4.32932 | 19.0384 |
| microbial capacity = 0 | 0 | 0 | 0 |
| epithelial permeability = 0 | 0 | 0 | 0 |
| host FMO3 capacity = 0 | 0 | 0 | 0 |
| taxonomy-only operator | 0 | 0 | 0 |

These are structural nulls in the synthetic model, not biological measurements. They show which link is necessary for this particular tracer output; they do not establish an intervention effect.

## Sensitivity (`±50%`)

Primary outcome: labelled plasma-TMAO AUC. The largest full-matrix effects were:

| Rank | Parameter | −50% change | +50% change |
|---:|---|---:|---:|
| 1 | systemic volume `V_s` | +100.0% | −33.3% |
| 2 | microbial capacity `Vmax_microbe` | −49.2% | +43.5% |
| 3 | microbial saturation `K_microbe` | +33.8% | −19.0% |

Among direct mechanism knobs, the order was `Vmax_microbe` (49.2% maximum absolute AUC change), `Vmax_FMO3` (25.8%), and `k_renal` (24.4%). `P_epithelial` had 6.0% maximum change in this parameterisation. `P_epithelial` and `A_epithelial` had identical effects, so the present observations cannot separate permeability from exchange area without geometry.

## Information gain and cost

The model uses time points `0, 0.5, 1, 2, 4, 8, 12, 24 h`. The full observation vector is plasma TMAO, lumen TMA, portal TMA, and cumulative urine loss. The local Fisher comparison for `Vmax_microbe`, `P_epithelial`, and `Vmax_FMO3` gave:

- plasma-only log-determinant: `11.1564`
- full multi-site log-determinant: `25.6236`
- information gain: `14.4672 nats`; determinant ratio `1.91875×10^6`
- assumed relative cost: plasma-only `1.0`; full `9.5`; gain per full-cost unit `1.52286`

This is a model-based identifiability calculation with assumed measurement noise, not an empirical information gain. The numerical information-separation criterion is **UPFYLLT**; biological identification remains **UNKNOWN**.

## Verification

- `python3 model.py` completed and wrote `results.json`.
- `python3 test_model.py` passed 6 tests.
- maximum mass-balance residual: `1.75014×10⁻16`
- analytical four-state linear-chain error: `4.59188×10⁻13`
- dimensional checks: all declared flux and pool units passed.
- preregistration digest: `PREREG.sha256` contains the digest of the frozen `PREREG.md`.

## Next resolution step

Run a preregistered, paired longitudinal human tracer study using the same participants and held-out participants for taxonomy, labelled precursor, TMA/TMAO, and barrier measurements. Collect serial plasma and 24 h urine, paired stool/lumen samples, and a validated portal/organ-accessible surrogate; measure epithelial area or recover it from segment-specific permeability geometry. Prespecify the perturbation (microbial suppression, barrier perturbation, or host FMO3 perturbation), isotope recovery, sampling time, assay limits, renal function, diet, medication, and correlated error. The current model identifies the required measurements and the remaining confounding; it does not identify a real microbial, barrier, or host effect.

All numerical outputs above are recorded in `results.json`.
