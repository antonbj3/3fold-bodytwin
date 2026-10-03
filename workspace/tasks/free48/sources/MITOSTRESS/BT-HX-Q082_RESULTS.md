BT-HX-Q082

## Built and sourced
A first executable, inventory-based moving-control-volume model is provided in `model.py:326`. It uses `V_L=A_r x_m`, `V_R=A_r(L-x_m)`, pairs each internal flow and shifts the substance in the swept volume with the concentration of the donor. diffusion and positive ion transport are given by the Nernst equation; deformation gives `C_m=c_A A_m` and the surface charge gives `V_m=(q_b A_m+F N_i^M)/C_m`. Parameter table and unit check are in `model.py:597`; source and assumptions in `PREREG.md`.

Pérez-Mitta & MacKinnon (2023), DOI `10.1073/pnas.2221541120`, Fig. 3B, report specific capacitance `0.3–0.5 microF/cm^2`; the model uses the center value `0.004 F/m^2`. Other parameters are explicit assumptions, not measurement data.

## Results
All frozen criteria in `results.json` are true. The main run yields mass residual `4.3368086899420177e-16`, charge residual `1.1293772630057336e-16`, minimum concentration `0.2 mol/m^3`, minimum inventory `0.0 mol`, membrane transfer `5.055916812565041e-06 mol` neutral substance, and `4.0075468050517423e-07 mol` positive ion. The end surface is `0.011205825438169701 V`; surface bound numbers are `1.5648902929053098e-08 mol` and `4.999e-15 mol` respectively.

The analytical limiting case (`model.py:523`) has relative concentration deviation `3.3750779948604757e-15`. Placebo without swept volume flow (`model.py:556`) falls with deviation `0.4285714285714288`, thus demonstrating that the criterion is not trivial. The marginal flow run has volume residual `0.0 m^3/s` and preserves both balances.

## Sensitivity
For `k_s`, membrane transfer from `2.628774938967258e-06` to `7.297419981983366e-06 mol` increases by factor `0.5` to `1.5`. For `k_i`, ion transfer increases from `2.0263126397462614e-07` to `5.944727013560453e-07 mol`. For `N_max_i`, the final voltage changes from `0.01060291271908485` to `0.011808738157254549 V`; no negative inventories arise.

## Next step
No measurement data is invented. The next resolution step is a spatial Poisson–Nernst–Planck model with counterions, the moving grid itself and independently measured permeabilities/surface tension; only then can dilution and drift be separated from physical redistribution effects.
