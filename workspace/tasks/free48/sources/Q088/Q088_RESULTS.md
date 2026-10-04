BT-HX-Q088

# First executable mechanistic model

## What was built

`model.py` is a reduced cylindrical model of a skeletal muscle fiber/territory with three coupled states: a radially resolved oxygen field, ATP/PCr energy and a normalized ionic disturbance `I`. The oxygen field uses implicit radial diffusion with Robin conditions (`k_transfer=1.0e-5 m/s`), Michaelis-Menten oxygen consumption and volume weighting. ATP and PCr fluxes are driven by oxygen factor, substrate factor, energy demand and fast CK transfers. The ion pump is ATP- and oxygen-dependent. The preregistered functional outcome is `F=0.50*PCr/PCr_rest+0.30*ATP/ATP_rest+0.20*(1-I)`. It is an equivalent surface-supplied cylinder, not a literal central-capillary geometry.

All assumptions, units and sources are in `model.py` (`PARAMETER_TABLE`) and frozen rules are in `PREREG.md` + `PREREG.sha256`. No measurement data have been created or used. The model run is reproducible with:

```text
python3 model.py --output results.json
python3 -m unittest -v test_model.py
```

## Reference anchors

- **Source, Layec et al. 2013**, DOI `10.1152/japplphysiol.00257.2013`, Table 3: PCr-tau `33 +/- 21 s` (free flow) and `27 +/- 10 s` (reactive hyperemia); `V_ATP` `28.7 +/- 13.3` and `41.2 +/- 13.6 mM/min` respectively. Table 4: reoxygenation mean response time `70 +/- 15` and `24 +/- 15 s` respectively. Methods: `ATP_rest=8.2 mM`.
- **Source, Heskamp et al. 2021**, DOI `10.1113/JP280771`, Table 2: distal/proximal `k_PCr=0.44 +/- 0.26` and `1.50 +/- 0.57 min^-1` respectively; `V_PCr=5.2 +/- 3.1` and `23.3 +/- 8.9 mM/min` respectively; `k_O2Hb=5.4 +/- 3.8` and `7.8 +/- 4.4 min^-1` respectively.
- **Source, Kushmerick et al. 1992**, DOI `10.1073/pnas.89.16.7521`, abstract/first page: `ATP=8 mM`, `PCr=32 mM`, total creatine `39 mM` for fast fibers.
- **Source, Piiper & Scheid 1986**, DOI `10.1016/0034-5687(86)90118-0`, abstract: the Krogh cylinder model is more suitable for skeletal muscle; capillary/fiber ratio approximately `2` and capillary/fiber radius ratio approximately `0.1`.
- **Source, Segal & Faulkner 1985**, DOI `10.1152/ajpcell.1985.248.3.C265`, abstract: calculated critical oxygen diffusion radius `1.19 mm` at `20 grader C` and `0.51 mm` at `40 grader C`. This is context for a whole-muscle piece, not a claimed cell value.

The reference values are scale and direction anchors; they are not a trained model and do not demonstrate validation on the target person.

## Frozen result

| Condition | pO2 (kPa) | work (mM/s) | R (mikrometer) | T90_F (s) | V_PCr0 (mM/min) | PCr at end of exercise (mM) | O2, volume mean (mol/m3) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Control | 13.3 | 0.60 | 20 | 50.41 | 27.69 | 10.42 | 0.3063 |
| Model reference run | 13.3 | 0.60 | 25 | 53.52 | 27.37 | 9.18 | 0.2847 |
| Joint stress | 2.0 | 0.90 | 40 | 250.56 | 9.56 | ~0 | 0.0271 |
| Null model, same stress | 2.0 | 0.90 | 40 | 72.11 | 32.70 | ~0 | — |

`T90_F` is 200.15 s longer in the stress corner than in the control; `V_PCr0` is 34.5 % of the control value. The preregistered prediction (at least 90 s longer and at least 50 % lower rate) passes. The reference run's `PCr-tau=28.83 s` is close to published 27–33 s, and `V_PCr0=27.37 mM/min` is within the published overview range `5.2–41.2 mM/min`; this is a scale comparison, not calibration.

The combination is explicit: at `2 kPa`, `0.90 mM/s` and `R=40 mikrometer`, `V_PCr0=9.56` and `T90_PCr=233.28 s`, compared with `12.17 mM/min` and `184.02 s` at the same oxygen deficit and work but `R=20 mikrometer`. Under normoxia, `0.60 mM/s`, an increase from `R=20` to `R=40` gives `V_PCr0=27.69 -> 26.22` and `T90_F=50.41 -> 66.07 s`. At `R=40`, increasing work from `0.60` to `0.90 mM/s` under normoxia increases `T90_F=66.07 -> 87.76 s`; rate saturates at `26.22 mM/min`.

## Sensitivity, ±50 percent

| Parameter | lower (−50 %) | nominal | higher (+50 %) |
|---|---:|---:|---:|
| pO2 (kPa) | 6.65: V=22.22, T90=110.35 s | 13.3: V=27.37, T90=53.52 s | 19.95: V=29.26, T90=36.05 s |
| work (mM/s) | 0.30: V=0, T90=46.57 s | 0.60: V=27.37, T90=53.52 s | 0.90: V=27.37, T90=84.16 s |
| R (mikrometer) | 12.5: V=28.13, T90=46.14 s | 25: V=27.37, T90=53.52 s | 37.5: V=26.43, T90=63.75 s |

At 0.30 mM/s no PCr depletion occurs under normoxia; the zero rate therefore means no initial resynthesis flux, not an error. It expresses that this state does not test the same recovery load.

## Controls, limitations and next steps

`results.json` contains the entire time course, radial oxygen profiles, grid, sensitivities and controls. The unit check passes: `0.399 mol/m3` at the boundary, `0.119 mol O2 m^-3 s^-1` as oxygen consumption, `0.8 s^-1` as boundary-transfer and a dimensionless function score. Four tests pass, including the analytical zero-consumption limiting case. The anoxia, sham and null controls pass; all frozen mechanical criteria are `true`.

What is not solved: `T90_ion` did not reach 90 % within 600 s in the control or stress run and is therefore reported as `>600 s`; the composite functional outcome can still pass 0.90 because ATP/PCr dominates. Radius, Robin coefficient, energy demand and substrate factor are assumptions, not measured person parameters. No pH, Ca2+, ROS, substrate or capillary network model is included. The grid's low-work state has no PCr depletion and shall not be interpreted as the same recovery test.

The next resolution step is matched longitudinal measurements of oxygen/perfusion, PCr or ATP, an ion-balance proxy and function under crossed oxygen deficit and work levels, with actually measured radius or cross-sectional area. The parameters shall be hierarchically estimated per person/fiber, geometry and capillary density shall be measured, and the preregistered functional outcome shall be validated without subjective post-adjustment. Until then the combinations above are model predictions, not measured recovery differences.
