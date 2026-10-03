BT-HX-Q009

## Resultat

A runnable first-principles model is in `model.py`: total plasma curve → free/bound plasma → passive tissue exchange → reversible tissue binding and free clearance. It builds on `MET-PHARMACOKINETICS-ADME` and the verified primary source Gill et al., DOI `10.1093/jac/dkac055`, Table 2; source data, derivations and hypotheses are separated in `results.json`.

With an identical total plasma curve, `f_u,p=0.40` versus `0.20` gives a free tissue AUC ratio of **2.0009** (predicted **2.00**). The case `f_u,p=0.05` gives **0.25** of the reference's free tissue AUC. Equality in total plasma is therefore not enough: plasma friction is a direct multiplicative input, while transport and clearance filter the temporal shape and level.

The reference study reports tissue penetration **0.54** and **0.66**, respectively; the model's unbound nominal value is **0.920**. This is not calibration or validation and empirical validity is `UNKNOWN`. Sensitivity to `±50 %` in `f_u,p`, `k_transport` and `k_elimination` is in `results.json`; the checks give ratio **1.0** for the same binding and **0.0** at zero transport. `test_model.py` passes, including the analytical limiting case.

## Next step

Measure matched total and free plasma, extracellular free tissue and intracellular target tissue; estimate finite binding capacity, time-dependent protein binding, `PS/V_t` and clearance; validate on held-out individuals. The source's abstract/Table-2 discrepancy for free plasma AUC is documented and is not used as calibration.
