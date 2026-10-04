# PREREG CLOUD-DXAINCREMENT

Train on 1999–2000; evaluate once on 2001–02, age ≥50. Use identical complete-case people for every model. Primary outcome is knee-extension peak force in Newton, as recorded in MSX. Comparator: ordinary least squares with intercept, body mass, height, sex and age. Primary augmentation: add bilateral DXA leg lean mass as one feature. Secondary: add leg lean × height instead. Never choose a model by looking at holdout RMSE. Use grouped training cross-validation only for optional hyperparameters.

Primary criterion: DXA augmentation yields at least 5% lower holdout RMSE than the strong comparator, with paired bootstrap 95% CI lower bound above zero, and the sign agrees in all five DXA imputations. If units/field meaning cannot be confirmed, report criterion UNKNOWN. Subgroup estimates are exploratory; no threshold. Record departures, leakage checks, imputation handling and unverified source facts. This preregistration is frozen.
