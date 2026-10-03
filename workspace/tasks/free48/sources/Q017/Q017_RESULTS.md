BT-HX-Q017

## Resultat

**Builds on:** no aborted implementation was found in the directory. `PREREG.md` and `PREREG.sha256` were frozen before the run. The model is based on the nodes `INT-CARDIORENAL-PRESSURE-AXIS`, `ORG-KIDNEY-NEPHRON`, `MODEL-RAAS-BLOOD-PRESSURE`, `MODEL-ARTERIAL-WINDKESSEL`, `SYS-BAROREFLEX-AUTONOMIC` and `scripts/msk/renal_filtration.py`, `nephron_transport_cell.py`, `raas.py`, `baroreflex.py`, `venous_return.py` and `fluid_compartments.py`. The segmented nephron model, full RAAS-kaskaden and pulse wave model were not rebuilt.

**Primary Test:** protocol is 3 L/3 h, which matches Kumar et al. 2004, DOI `10.1186/cc2844`, Table 1 (n=24; cardiac index +14,7 ± 2,4%). The model's predicted relative CO change is **+0,8367347**. The reference range is 0,097–0,197; absolute error is 0,6897347. **Frozen criterion: NOT met.**

At the end of the infusion, the full model reports: CO 9,0000 L/min (numerical limit 9,0), MAP 45,9184 mmHg, venous pressure 18,8130 mmHg, renal blood flow 0,4356 L/min, GFR 0,000 mL/min, urine 0,000 mL/min and ECF 20,8210 L. Sympathetic, Ang II and natriuretic activity is −1,5, −0,5 and 2,0. The relative final error of the mass balance is 6,11×10⁻¹⁵; analytical pressure–flow boundary case and unit control passes. The state minimum value for MAP was 34,5054 mmHg, below the declared 40-mmHg calculation limit; the boundary is not a hard state projection and is therefore marked as a model constraint.

**Other checks:** 5/5 tests pass; `py_compile` passes; `PREREG.sha256` is verified. The renal pressure sweep 80–180 mmHg gives maximum GFR-avvikelse 0,2701, above the frozen limit 0,20. The zero filtration control gives GFR=urine=0. Neural/hormonal ablation gives at the end of the infusion MAP 183,6735 mmHg and GFR 181,4207 mL/min; volume–pressure ablation gives CO 4,3600 L/min and MAP 22,2449 mmHg. Afterload and renal resistance steps alone produce small CO changes (−3,12% and −1,17%, respectively). Sensitivity: ±50% preload coefficient changes relative CO by up to 32,19%; ±50% myogenic coefficient changes renal flow by up to about 12%; ±50% RAAS-time changes the reported final numbers by about 1,15%. The primary CO-kurvan is therefore saturated and not identifiable for all three parameters in this case.

**Source Status:** Davies & Shock 1950 (DOI `10.1172/JCI102286`, Table II), Kumar 2004, Bikia 2024 (DOI `10.1038/s41598-024-56137-8`, Table 2), Maas 2012 and Bhave & Neilson 2011 are verified published sources. Other values ​​are assumptions or reduced coefficients. No internal or measured transient data has been used; this is not a clinical claim.

**Next resolution step:** measure synchronously MAP/PCWP, CO, renal flow, GFR, urine and volume during a defined perturbation. Then add and preregister a measured or independently sourced pressure–venous return coupling, pump limit, and renal venous congestion; do not calibrate the current run by the result.
