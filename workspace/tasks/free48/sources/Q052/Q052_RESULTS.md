BT-HX-Q052

## Results

I built on `Q052` in `inputs/QUESTION.md` and `BRIEF.md`/`inputs/NIGHT_PREAMBLE.md`; no external BodyTwin anchors could be read in the session. `model.py` couples a process state `s` to the same pore geometry `(φ,d)` and from there to effective stiffness, Darcy permeability, diffusion/advection and a hypothesis-based response window.

The nominal case (`s=0.5`, `φ=0.80`, `d=800 µm`) gives:

- `E_eff = 3.30 GPa` — in Chao et al. (2021), DOI `10.3389/fbioe.2021.779854`, figure 8/9: 2.6–4.0 GPa.
- `K = 1.87e-8 m²` — Chao et al. figure 13A: verified reference value.
- `D_eff = 1.322e-9 m²/s`, `c_mean = 0.99886`, `R = 0.48876`; `R` is a dimensionless hypothesis, not a measurement.
- The process sweep's correlation between `s` and mean concentration is `0.9503`; randomized placebo gives median absolute correlation `0.2843`.

±50 % sensitivity: the tortuosity parameter changes `K` by about −8.5/+9.8 %; the mechanics correction changes `E` from 1.65 to 4.95 GPa; consumption affects `c_mean` little. All five tests, analytical limits, unit checks and frozen criteria passed. Full provenance and all numbers are in `results.json`; the hash is `PREREG.sha256`.

## What next

Calibrate `s→(φ,d)` and the response law against the same sample's micro-CT, compression, flow and biological measurements; use an independent K holdout and replace 1D homogenization with 3D flow. The current tissue response should therefore not be interpreted as validated cell or tissue growth.
