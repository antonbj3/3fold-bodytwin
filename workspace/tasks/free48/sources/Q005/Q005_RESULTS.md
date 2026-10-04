BT-HX-Q005

# Result

## Core results

The model is a first-principle chain `f_u -> Cu -> OCT2 uptake -> MATE-efflux -> tubular secretion -> renal clearanc`. The main outcome is creatinine clearance (`CL_Cr`, mL/min); the independent secondary renal readout is metformin clearance (`CL_Met`, L/h). No internal data or new measurements are used.

The nominal model prediction is `CL_Cr = 140.616 mL/min` and `CL_Met = 29.585 L/h` (`results.json:nominal`). The filtration-only null model gives `94.99 mL/min` and `5.6994 L/h` (`results.json:null_model`) respectively. The model's GFR reference is `94.99 mL/min/1.73 m^2`; the source baseline for CrCl is `152.3 mL/min` and the metformin renal clearance `30.0 L/h` (`results.json:reference_anchors`). Reference and unit check passes (`results.json:reference_check`, `results.json:unit_check`).

## Which molecular parameter?

The leading molecular candidate is **fraction unbound, `f_u`**: score `1.9300 mL/min per kostnadsenhet`, RMSE reduction `9.397 %` and primary variance proportion `14.545 %` (`results.json:uncertainty.candidates`, entry `f_u`; `results.json:uncertainty.leading_molecular_parameter`). However, it does not meet the pre-frozen 10% threshold. Therefore, it is correct to report **no confirmed winner**, not to upgrade `f_u` to a validated answer.

`GFR` meets the criterion (`15.851 %` reduction, score `1.0852`) but is a physiological, not molecular, parameter-driven value. `Jmax_OCT2` reaches `4.389 %`; `K_OCT2` gives `-2.158 %` in the frozen RMSE rank. The negative K ranking is due to the high K/J correlation and the fact that a free K is not identical to a free transporter expression. `P_MATE` does not affect the steady-state clearance in the linear efflux model, but halving it increases the Cr cell concentration from `316.85` to `633.70 microM` (`results.json:sensitivity.records`).

## Sensitivity and correlation

The ±50 % run shows the direction in the mechanism: `K_OCT2` down `50 %` raises the Cr clearance to `163.833 mL/min`, while up lowers it to `129.110`; `Jmax_OCT2` down/up gives `117.803/163.429 mL/min`; GFR down/up gives `93.121/188.111 mL/min` (`results.json:sensitivity.records`). `f_u` is physically bound to at most 1, so the plus-50 percentage drop is capped at nominal 1.0.

The two outputs are not independent in the latent sense: the calculated correlation between `CL_Cr` and `CL_Met` is `0.585363` (`results.json:uncertainty.latent_output_correlation`). K and J have log correlation `0.70` (`results.json:base_parameters`), which is intended to capture that these cannot be optimized independently of each other.

## Sources and reproduction

- Severance et al., DOI `10.1124/jpet.117.242552`, Table 1: metformin-OCT2 `Jmax = 1046 +/- 171 pmol cm^-2 min^-1`, `Ktapp = 518 +/- 45 microM`.
- Topletz-Erickson et al., DOI `10.1002/jcph.1750`, Table 2/Table 3/Figure 3A: metformin renal clearance, GFR and CrCl.
- Mathialagan et al., DOI `10.1124/jpet.123.001890`, abstract: OCT2/MATE and approximately 30% active secretion of creatinine clearance.

`PREREG.md` and `PREREG.sha256` are created before the first run. `python3 -m unittest -v test_model.py` passed five tests, including three analytical limit cases. No K02/K08/K10/K11 nodes or internal BodyTwin files were available in the job's accessible workspace.

## Next resolution step

1. Measure or validate `f_u` and OCT2/MATE kinetics in the relevant human kidney cell/plasma context; especially if you want to lower the residual error that causes `f_u` to end up just below the threshold.
2. Replace the assumptions `C_Cr = 100 microM`, total area `2.0e5 cm^2`, and the substrate factors with validated data.
3. If the organ outcome is to be damage rather than transport clearance, replace the steady-state clearance with a time-dependent cell and damage chain; then `P_MATE` can become controlling.
4. Repeat pre-registration with real failure model and independent organ outcomes; this model cannot specify a universal molecular parameter.
