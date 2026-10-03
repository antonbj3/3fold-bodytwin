BT-HX-Q026

## Builds on

- `inputs/QUESTION.md`, Q026: spatial signal/receptor measurements with time anchors.
- Harish RK et al., *Development* 150, dev201559 (2023), DOI `10.1242/dev.201559`, primary source. Fig. 3C reports `D_fast = 55 um^2 s^-1`, `D_slow = 4 um^2 s^-1`, and 93%/7% fast/slow components. Fig. 4C reports source concentration about `8 nM`, a `0.8 C0` target domain about `40–70 um`, and a `0.4 C0` domain about `70–140 um`.
- The existing job directory contained `agent.log`, `ALLOW_WEB`, `BRIEF.md`, and `inputs/`; no earlier model output was present. No internal or invented measurement data were used.

## First runnable mechanism

`model.py` solves the steady form of

`dL/dt = D d2L/dx2 - k_clear L + s0 exp(-x^2/(2 sigma^2))`

with no-flux outer boundaries. The ligand field is in nM, distance in um, time in s. Receptor readout uses the fast-equilibrium law `theta = L/(K_D + L)`. The raw source amplitude is arbitrary; the linear solution is rescaled to the source reference `C0 = 8 nM` for occupancy and threshold ratios, without changing spatial crossings. The slow/HSPG-bound pool is deliberately excluded from this first model.

Frozen parameters and units are in `PREREG.md` and `model.py`: `D=55 um^2 s^-1` (verified), `k_clear=0.0015 s^-1` (assumption), `sigma=25 um` (assumption), `K_D=2 nM` (assumption), `s0=1 nM s^-1` (normalization), domain half-width `500 um` (assumption). `unit_check()` gives `D/dx^2 = 0.55 s^-1` at `dx=10 um`, matching the clearance-rate unit.

## Execution and spatial result

`python3 -m unittest -v test_model.py` ran 4 tests, all passed, including the analytic `D=0` limit `L=S/k_clear`. `python3 model.py` generated valid `results.json`.

At `dx=10 um` (`results.json:coarse_metrics`):

- `x_0.8 = 64.48 um`, inside the source window `40–70 um`.
- `x_0.4 = 203.81 um`, outside the source window `70–140 um`.
- The `0.8 C0` to `0.4 C0` span is `139.34 um` (`13.93` grid cells).
- `L(100 um) = 5.48 nM`; receptor half-occupancy distance is `459.52 um` for the assumed `K_D=2 nM`.
- Halving the grid to `5 um` changes `x_0.8` by `0.116%` and `x_0.4` by `0.051%` (`results.json:grid_convergence`).

Operational resolution for this modeled Fgf8a-like gradient is therefore a spatial grid no coarser than `10 um`, with a threshold span of about `140 um`; pooling over the whole span would erase the two reported response levels. The comparison is mechanistic, not a tissue-wide measurement.

## Sensitivity, plus/minus 50%

Values below are copied from `results.json:sensitivity`.

| controlling parameter | x_0.8 at -50% / +50% (um) | x_0.4 at -50% / +50% (um) | derived change |
|---|---:|---:|---|
| `D` | 51.33 / 75.54 | 146.26 / 260.61 | threshold span 94.93–185.06 um |
| `k_clear` | 85.95 / 56.15 | 328.50 / 166.48 | threshold span 110.33–242.54 um |
| `K_D` | 64.48 / 64.48 | 203.81 / 203.81 | half-occupancy distance is `null` at `K_D=1 nM` within the 500-um domain and 351.80 um at `K_D=3 nM`; baseline is 459.52 um |

`null` means no crossing was found before the finite boundary, not a measured zero.

## Comparison and next resolution step

The primary `x_0.8` comparison is numerically stable, while `x_0.4 > 140 um` is a pre-registered falsifier for the reduced pure free-diffusion description. The next model needs the measured slow/HSPG-bound pool, source geometry, and time-resolved ligand plus receptor measurements; then compare the two-pool solution on `10 um` and `5 um` grids. Until those data exist, `K_D`, clearance, and source width remain assumptions, not BodyTwin or human observations.
