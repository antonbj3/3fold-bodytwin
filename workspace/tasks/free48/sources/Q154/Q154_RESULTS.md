BT-HX-Q154

# Resultat

## Primary outcome

The model is an executable, time-resolved mechanistic model for D-glucose in a brain capillary endothelial cell. It has separate pools for luminal glucose, free endothelial glucose, receptor-bound glucose, vesicular glucose, abluminal glucose, endothelial lactate and brain lactate. Passive diffusion, luminal transporter influx, luminal efflux, receptor binding, endocytosis, vesicle release, abluminal transport and endothelial metabolism are calculated as different directed fluxes. The run window is 3600 s and the adaptive Radau solver can be tested with finer time steps than the explicit 1-sampling.

The prediction `J_u` is the mean of the total luminal flux to the endothelium during the last 600 s. The result is **0,4432469692 µmol g⁻¹ min⁻¹**. The Hertz reference is **0,46 µmol g⁻¹ min⁻¹** (Hertz et al., 1981, *Journal of Clinical Investigation*, DOI `10.1172/JCI110073`; verified abstract page: https://pmc.ncbi.nlm.nih.gov/articles/PMC370607/). Relative error is **0,0364196321**. The frozen criterion `|fel| ≤ 0,30` is therefore **met**.

This is a numerical order-of-magnitude result, not an independent validation of parameter values. The passive null model gives **0,0020030086 µmol g⁻¹ min⁻¹** and does not meet the criterion; it thus cannot alone reproduce the reference's order of magnitude under the assumed geometry and permeability values.

## Modellens utdata

| Quantity | Value | Unit |
|---|---:|---|
| Gross intact glucose to the brain side | 0,2317585024 | µmol g⁻¹ min⁻¹ |
| Net intact glucose to the brain side | 0,1854021467 | µmol g⁻¹ min⁻¹ |
| Lactate to the brain side, gross | 0,0482665479 | µmol g⁻¹ min⁻¹ |
| Lactate to the brain side, net | 0,0482665479 | µmol g⁻¹ min⁻¹ |
| Luminal efflux fraction | 0,0738930453 | 1 |
| Endothelial retention at the end | 2,8731099128e-9 | mol |
| Free endothelial glucose | 2,8701206291e-9 | mol |
| Receptor-bound glucose | 9,9642791307e-13 | mol |
| Vesicular glucose | 1,9928558262e-12 | mol |
| Luminal glucose at the end | 4,1842265329 | mM |
| Maximum mass-balance residual | 2,4815418377e-24 | mol s⁻¹ |
| Maximum charge residual | 0,0 | mol eq s⁻¹ |

The luminal efflux fraction is defined as `J_luminal_efflux/(J_luminal_in+J_luminal_efflux)`, where reverse passive flux, transporter efflux and receptor unbinding count as efflux. Transmembrane abluminal transport is not included in this definition.

## The paths

Mean for the last 600 s, in µmol g⁻¹ min⁻¹:

- Luminal passive diffusion in: **0,0015768674**.
- Luminal transporter influx: **0,4416671125**.
- Receptor-mediated endocytosis: **2,9892837408e-6**.
- Luminal transporter efflux: **0,0353602048**.
- Luminal receptor unbinding: **5,9785674816e-6**.
- Abluminal passive glucose flux: **0,0018887277**.
- Abluminal transporter flux: **0,2298673833**.
- Vesicle release: **2,3914269929e-6**.
- Endothelial metabolism: **0,2224840323**.
- Vesicle substrate lost through the described degradation: **5,9785674822e-7**.

The receptor path is thus an explicit, separately accounted micro-path and not a synonym for permeability. It is small at baseline but can change independently through `r_total_pmol_per_g`, `k_endocytosis_s`, `k_release_s` and `k_vesicle_degradation_s`.

Lactate is treated as `lactate^-` with an equal, oppositely directed `H+` cotransport in the charge ledger. The proton pool and pH-buffer dynamics are not solved, but the electrical net charge per membrane is explicitly zero. All pools and both mass balances are non-negative in the checked run.

## Sensitivity

For primary `J_u`, ±50 % gave the following largest relative changes:

| Parameter | −50 % | +50 % |
|---|---:|---:|
| `vmax_luminal_in_umol_min` | 0,2305056191 (−47,996 %) | 0,6481319140 (+46,224 %) |
| `c_lumen_source_mM` | 0,3343578342 (−24,566 %) | 0,4917392894 (+10,940 %) |
| `km_luminal_in_mM` | 0,5066946925 (+14,314 %) | 0,3955884764 (−10,752 %) |

For secondary gross brain-side glucose, the largest parameters are `vmax_luminal_in_umol_min`, `vmax_metabolism_umol_min` and `c_lumen_source_mM`. On halving and a 1,5 factor for metabolism, brain-side glucose becomes **0,3124094085** and **0,1794288737 µmol g⁻¹ min⁻¹** respectively. This shows that transporter capacity governs primary influx while metabolism and abluminal allocation govern what reaches the brain side. Sensitivity is model sensitivity, not a measurement interval.

## What is and is not identified

Numerical acceptance is **YES**, but causal identification of passive, transporter, efflux, transcytosis and metabolism is **UNKNOWN**. `results.json` contains modelled directed fluxes but no independent substrate-specific inhibitor, tracing or metabolite measurements. In particular, glucose efflux capacity and receptor identity are assumptions, not published measured data in this run. No new measured data has been fabricated.

## Next resolution step

1. **Data:** measure independent influx and efflux with substrate-specific tracers, simultaneous arterial/venous glucose and lactate, and perturb one path at a time with a verified inhibitor or receptor ligand.
2. **Geometry:** replace the network-equivalent 1-g parcel with 3D capillary geometry, measured endothelial area and thickness, luminal mixing volume and blood flow; keep two membranes and the vesicle pools as separate objects.
3. **Measurement:** connect endothelial GLUT/efflux expression, vesicles and metabolite profile to time-resolved flux densities. An independent directed flux and metabolite profile is necessary to distinguish active transport from passive permeability.

## Reproduction and files

- Equations, parameter source/implicit assumption, unit contract and dimensional check: `model.py`.
- Analytical and numerical tests: `test_model.py`.
- Preregistration and frozen hash: `PREREG.md`, `PREREG.sha256`.
- All reported numbers and sensitivity cases: `results.json`.
- Commands run: `python3 model.py` and `python3 test_model.py`; the test suite passed, the dimensional check passed and analytical equilibrium gave zero flux.
