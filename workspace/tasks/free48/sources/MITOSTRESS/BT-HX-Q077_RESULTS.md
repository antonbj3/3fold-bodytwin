BT-HX-Q077

## Result

A runnable finite-volume crista model was built in `model.py`, with `test_model.py` covering the equal-boundary no-sink analytical limit, unit conversion, matched area/volume/site counts, topology direction, and diffusion-time scaling. The model solves steady IMS ADP diffusion with Robin CJ exchange, a low-substrate mass-action ANT flux, matrix ADP/ATP balance, and ATP-synthase occupancy. All output is in `results.json`.

At the frozen defaults (`L=0.9 um`, matched `A_IM` and `V_IMS`), one-CJ versus two-CJ ATP output is `1.4223183584034947e-19` versus `1.4501088856720255e-19 mol/s` (`42.827006799189874` versus `43.66379913421842 molecules ATP/ms/um2`). Thus `R_flux=1.0195388937394612`; the preregistered `R_flux>=1.10` criterion fails. Maximum fractional ADP depletion is `0.9997656984506117` versus `0.9997079398178961`. Maximum diffusion time is `0.0810000000000001 s` versus `0.019274241522903045 s`; the two-CJ topology passes the preregistered diffusion-time check. The branched topology gives `1.44332757392832e-19 mol/s` and lowers the maximum distance further.

The equal-boundary no-sink error is `4.163336342344337e-17 mol/m3`; the test suite passes. The frozen criterion therefore has one failed condition, recorded as `flux_ratio_ge_1_10=false`; no post-hoc parameter retuning was used.

## Sensitivity

For `D`, ANT `kcat`, and total CJ area, respectively, the two-CJ/one-CJ flux ratios over `-50% / baseline / +50%` are `1.0105273528872998 / 1.0195388937394612 / 1.0318276117750695`, `1.0283598744046196 / 1.0195388937394612 / 1.015539285626353`, and `1.008095321695957 / 1.0195388937394612 / 1.037033374597828`.

## Source, limits, next step

The external anchor is verified: Adams et al., DOI `10.3390/cells14040257`, Figure 4, reports `0.60` versus `0.90 J/JMAX` and a `17.0%` connectivity-associated increase; these are literature anchors, not fitted data. The next resolution step is to replace assumed low-substrate ANT and effective diffusivity with measured or independently calibrated kinetic/transport parameters, then add finite saturation and 3-D CJ geometry before validation against ADP/ATP transients.
