BT-HX-Q084

## Results

No model or results file existed from the interrupted session; the work therefore builds on `inputs/QUESTION.md`, `inputs/NIGHT_PREAMBLE.md` and the frozen `PREREG.md`. The primary source is Destrian et al., PNAS 123(4), e2519599123 (2026), DOI `10.1073/pnas.2519599123`: Fig. 3C states `16–20 µm² s⁻¹`, Fig. 3D about `20 %` low-porosity effect and Fig. 4D–G up to `70 %` nano-obstacle effect.

`model.py` implements reflecting–absorbing finite-volume diffusion, available volume, tortuosity, Kozeny–Carman-like permeability, hydrodynamic hindrance `H` and reversible binding. The parameter table, equations and unit check are in the model; `results.json` contains all the numbers. The analytical empty-space case gave `2.5 s` against `2.5 s` and relative error `0`.

At the frozen nominal profile, the mechanistic passage time is `38.4363 s`, adjusted diffusivity `36.3435 s`, thus `+5.758 %` (the explanatory model is slower). Mean `D` is `9.0668 µm² s⁻¹` (`0.4533 Dα0`). At `10 s`, total bound probability is `0.4491`, spatial CV `0.8284` and profile L1 against the adjusted model `0.08286`.

Frozen criteria: spatial signature and analytical sanity pass; `≥10 %` passage separation fails. The conclusion is therefore not that crowding necessarily changes passage time by as much, but that this particular heterogeneous model lies below the threshold. ±50 % sensitivity gave 5.60–6.65 % passage contrast and 22.10–58.90 s mechanistic time for the three control parameters.

The next step is a measured local obstacle field and independent calibration of `K_nano`, binding kinetics and geometry; sharper bottlenecks and validation against independent FRAP/FCS are then tested.
