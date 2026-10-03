# BT-HX-Q026 — preregistration

## Scope and hypothesis

This is a bounded, first-principles model of one paracrine field: secreted Fgf8a-like ligand released at a tissue margin and read by receptor-bearing cells. It is not a claim that one ligand represents all paracrine signals.

Hypothesis: a free-ligand reaction–diffusion field is spatially resolvable when the cell-state transition induced by ligand concentration is narrower than the field's reporting grid. The primary predicted quantity is the first distance from the source at which the extracellular concentration reaches 0.8 of the source concentration, `x_0.8` (micrometres). The secondary diagnostic is `x_0.4`, the corresponding 0.4-level distance.

## Frozen model before execution

Use a one-dimensional, steady, extracellular field `L(x,t)` in micrometres, with concentration in nM and time in seconds:

`dL/dt = D d2L/dx2 - k_clear L + S(x)`

with a localized, normalized source `S(x) = s0 exp(-x^2/(2 sigma^2))`, `s0` in nM s^-1 and `sigma` in micrometres. Use a finite domain with homogeneous no-flux outer boundaries. The source amplitude is only a normalization and is not interpreted as a measured secretion rate. Receptor access is represented by a fast-equilibrium occupancy law, `theta = L/(K_D + L)`, where `theta` is dimensionless and `K_D` is in nM. This reduced receptor law is a readout, not a fitted receptor-trafficking model.

The following are frozen implementation choices: `D = 55 um^2 s^-1`, `k_clear = 0.0015 s^-1`, `sigma = 25 um`, `K_D = 2 nM`, domain half-width `500 um`, source amplitude `s0 = 1 nM s^-1`, and candidate grid spacing `dx = 10 um`. `D` is the fast-component value reported in the primary source. `k_clear`, `sigma`, `K_D`, domain size, and `s0` are explicit assumptions, not measurements. The source paper's 7% slow component is recorded but excluded from this first model.

## Reference value and source check

Primary source: Harish RK, Gupta M, Zöller D, Hartmann H, Gheisari A, Machate A, Hans S, Brand M. “Real-time monitoring of an endogenous Fgf8a gradient attests to its role as a morphogen during zebrafish gastrulation.” *Development* 150, dev201559 (2023). DOI: `10.1242/dev.201559`.

Verified quantitative anchors from the primary source:

- Fig. 3C: fast diffusion coefficient `55 um^2 s^-1`; slow coefficient `4 um^2 s^-1`; approximately 93% fast and 7% slow in the fitted measurements.
- Fig. 4C: source-domain extracellular concentration approximately `8 nM`; the 0.8-level is associated with a target-expression domain approximately `40–70 um` from the margin, and the 0.4-level with approximately `70–140 um`. These are reported spatial windows from the experiment, not a claimed exact exponential fit.
- The source explicitly describes Fgf8a as propagating through extracellular space and forming a graded distribution. This makes it a primary experimental anchor for a spatial paracrine-gradient model.

If the source value cannot be located in the publication, it must be marked `OVERIFIERAD`; it was located in Fig. 3C and Fig. 4C, so the values above are not marked `OVERIFIERAD`.

## Frozen decision criteria

1. `x_0.8` is computed by linear interpolation of the first crossing of `L/L(0) = 0.8`.
2. Primary external check: `40 um <= x_0.8 <= 70 um`. This is the pre-registered comparison to the primary-source spatial window.
3. Grid check: with `dx = 10 um`, halving `dx` changes `x_0.8` and `x_0.4` by no more than 5%, and the interval between the two thresholds contains at least four grid cells.
4. A pure free-diffusion model is considered mechanistically incomplete if `x_0.4 > 140 um`, because that misses the primary-source 0.4-level window. This is a falsifier/diagnostic, not a parameter-tuning instruction.
5. Report receptor half-occupancy distance (`L = K_D`) as a derived readout, while keeping `K_D` visibly an assumption.

## Counterfactual and failure rules

The null/placebo is the same PDE with `S(x) = 0`; it must remain at the initial condition. A non-positive concentration, a non-finite value, a missing threshold crossing, or a grid-convergence failure is a model failure. No measured BodyTwin or human data are used. The slow/HSPG-bound Fgf8a pool is deliberately not included in this first run; it is the next model refinement, not an unreported success.

## Builds on

- `inputs/QUESTION.md`, Q026: spatial signal/receptor measurements with time anchors are the required input class.
- Harish et al. 2023, DOI `10.1242/dev.201559`, Fig. 3C–4C: quantitative Fgf8a diffusion, concentration, and spatial target-expression anchors.
- Existing files in this job directory: only `agent.log`, `ALLOW_WEB`, `BRIEF.md`, and `inputs/`; no prior model output was present.

## Not redone

