# BT-HX-Q033 — preregistration

Date: 2026-09-25
Status: frozen before the first run of `model.py`

## Scope and hypothesis

This is a first, dimensionally checked mechanistic model of a local elliptical incision in a loaded soft-tissue layer. The model must simultaneously give:

1. the force required to appose the incision edges,
2. the opening/deformation and local stress,
3. the exposed wound area,
4. outward fluid and evaporative boundary fluxes.

The hypothesis is that the exponential opening gives a faster-than-linear force increase with incision width, and that the exposed area gives two parallel transport paths: pressure-driven fluid transport according to Darcy and evaporative barrier loss. Increased natural prestress must increase the mechanics but not change the geometric area; the barrier-damage factor must affect transport but not force.

The model is a reduced-order first-principles model. It is not clinical validation and makes no claim that a particular incision gives a particular measurement.

## Predicted quantity

The default case is a full-thickness elliptical opening with geometry 46 × 13 mm, thus half-length `a = 23 mm` and half-width `b = 6,5 mm`.

Predictions before the run:

- `F_close > 0` and is of the order of a few newtons.
- `F_close`, `epsilon_max` and `A_exposed` increase when `a` or `b` increases.
- `epsilon_max = 2b/L_ref` and `A_exposed = pi*a*b`.
- `J_evap_cut/J_evap_intact = barrier_damage_multiplier`.
- `Q_total_cut` is positive for a positive pressure difference and zero when `b = 0`.

## Published reference values

### Mechanics

Lackmann F, Rohwedder T, Maron A, Stegen L, Brunnberg M, Brunnberg L, Burger M, Böttcher P. “Quantification of skin wound tension using a newly designed wound tensiometer.” *Tierarztliche Praxis Ausgabe K: Kleintiere/Heimtiere*. 2023;51(6):386–393. DOI: `10.1055/a-2150-0587`.

Primary source: ex vivo dog cadavers, 19 dogs. Table 1 states for experienced surgeons in experiment 1 an upper recommended tension of **5,4 N** and a lowest tension where simple apposition was advised against of **6,0 N**. For experiment 2, pooled fresh-cadaver data, the same table states **4,2 N** and **6,0 N** respectively. Figure 7 shows the relationship between measured tension and the surgeon's decision. The method description states elliptical full-thickness skin wounds and that incision size was increased stepwise.

The reference value is verified in Table 1 and Figure 7. It is **not** `OVERIFIERAD`: the source and table/figure reference are checked. The value is not a human-specific or incision-exact calibration.

### Transport analogue

Barthe M, Clerbaux LA, Thénot JP, Braud VM, Osman-Ponchet H. “Systematic characterization of the barrier function of diverse ex vivo models of damaged human skin.” *Frontiers in Medicine*. 2024;11:1481645. DOI: `10.3389/fmed.2024.1481645`.

Primary source: ex vivo human skin. Figure 2 shows that mechanical dermabrasion increased transepidermal water loss **5,9 times** relative to the control; the measurement is in **g m^-2 h^-1** and Figure 2 states `p = 0,004`, `n = 9`. Figure 4 also shows a **133-fold** increase in Lucifer Yellow in the receptor fluid after dermabrasion, as a separate traceability expression.

The reference value is verified in Figure 2 and Figure 4. It is **not** `OVERIFIERAD`. It is a barrier-damage analogue from human skin, not a measurement of surgical wound fluid; it is therefore used only as a justified assumption for an effective evaporative flux factor.

## Frozen criteria

The criteria are fixed before the first run:

1. `unit_check()` passes all declared SI unit checks.
2. The default case's total `F_close` lies between **0,5 N and 6,0 N**. This is a broad sanity interval spanning the published ex vivo forces; no individual square is claimed to be calibrated.
3. For positive `a`, `b`, `E` and `sigma0`, `F_close` and `A_exposed` are strictly increasing in `b`; `A_exposed` is strictly increasing in `a`.
4. Placebo/negative control `b = 0` gives `F_incremental = 0`, `A_exposed = 0` and `Q_total_cut = 0` within numerical tolerance `1e-12` in SI units.
5. With default parameters, `J_evap_cut/J_evap_intact = 5,9`; this is the frozen transport analogue, not an invented measurement.
6. The sensitivity run for `modulus_pa`, `prestress_pa` and `barrier_damage_multiplier` at `0,5x` and `1,5x` must be finite, lack negative area and show the expected direction for at least force and evaporative flux fraction respectively.

## What counts as error

- SI unit errors or mixing Pa, N and m³/s.
- Negative area, negative incr. force or negative absolute flux.
- Non-zero placebo effect at `b = 0`.
- Monotonicity error in geometry.
- Calling the source analogues measurements on the current incision.
- Reporting a parameter search, fitting or new training after the preregistration is frozen.

## Builds on

- `BRIEF.md`, `inputs/QUESTION.md` and `inputs/NIGHT_PREAMBLE.md`.
- The interrupted session's `agent.log` was read; no previous `PREREG.md`, `model.py`, `test_model.py` or results files remained in the working directory.
- Public primary sources: DOI `10.1055/a-2150-0587` and DOI `10.3389/fmed.2024.1481645`.

## Not redone

- No surgical decision-support model, no patient risk and no prognosis.
- No claimed measurement of a particular incision size, pressure, fluid or tissue.
- No reaction/biology, blood circulation, electrical transport or drug diffusion beyond the explicit Darcy and evaporative boundary model.
- No 3D-FEM, mesh refinement, time-dependent healing or individual calibration.

## Freezing

This file and `PREREG.sha256` must exist before `model.py` is run. The hash applies to this exact text before any later edits.
