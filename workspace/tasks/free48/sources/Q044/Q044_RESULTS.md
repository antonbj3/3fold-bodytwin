BT-HX-Q044

Builds on `inputs/QUESTION.md`, `inputs/NIGHT_PREAMBLE.md` and the empty interrupted session. `model.py` solves the quasistatic Poisson equation in five concentric tissue layers with a harmonic Green solution, potential at electrodes, average/mastoid reference and a fixed depth grid for the inverse error. `PREREG.md` was frozen before the run; the digest is in `PREREG.sha256`.

The run is synthetic: fixed tangential dipole source, 20 mm below the CSF/skull boundary, 1 nA and 6 mm separation. The baseline topography norm is 1,1449×10⁻⁹ V and the reference solver recovers the source at 20 mm (`results.json:baseline`). The null model gives RDM 0; the homogeneous spherical limit gives relative error 1,25×10⁻⁷; three tests pass (`results.json:null_model`, `verification`; `test_model.py`).

±50 % skull conductivity (0,010→0,005/0,015 S/m) gives RDM 0,0939/0,0604 and inverted depth 29/16 mm (9/4 mm error). ±50 % skull thickness gives RDM 0,0254/0,1103 and depth 18/12 mm (2/8 mm error). ±50 % CSF thickness gives RDM 0,0050/0,0176 and depth 20/19 mm. All values are in `results.json:sensitivity`; the computation is run with `python3 run_sensitivity.py`.

Reference value: Vorwerk, Wolters & Baumgarten (2024), DOI `10.3389/fnhum.2024.1335212`, Table 1 (standard conductivities) and Figure 2/7 (skull Sobol and localisation sensitivity). No measured data or individual anatomy is claimed. Average and mastoid reference gave RDM 0,9387, so the reference choice is a separate uncertainty.

The next step requires an individual tetrahedral MRI mesh, measured or calibrated tissue/electrode impedances and a preregistered comparison against EEG or a phantom.
