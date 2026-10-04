# BT-HX-Q107 — preregistration

## Hypothesis and reference

A two-compartment model (fundus as storage reservoir and antrum as mixing/contraction chamber) with a time-dependent pyloric opening and separate particle transport should be able to distinguish liquid from solid emptying. The prediction is half-emptying time (`T50`) and volume/retention curves for a standard meal. The geometry is thought to provide the most information when the antral pressure and the pyloric opening are large enough to alternately block and release outflow; a pure scalar `k_empt` should not be able to reproduce both a fast liquid phase and a slow, retained solid phase.

Reference: total gastric volume in mL, residual solid mass in g, pyloric flow in mL/s, antrum and fundus pressure in Pa, and `T50` in minutes. `Q_in` is zero after initial meal loading; mass and volume conservation is checked at each step.

## Referensankare

- Collins, P. J. (1983), *Gut*, normal-subject scintigraphic gastric-emptying reference: approximately `T50 = 20 min` for liquid meal and `T50 = 90 min` for solid meal. **OVERIFIERAD, out of memory**. This is a broad literature anchor, not a measurement from this session and not a source dataset used for calibration.
- The reference values ​​are only used as external direction. No partial volumes, pressures, flows, food properties or individual data are found.

## Frozen acceptance criteria

1. Numerical mass/volume conservation: maximum relative deviation in final state and cumulative balance `<= 1e-8` .
2. The floating `T50` of the model must be within factor 2 of 20 min, thus `[10, 40] min`.
3. The fixed `T50` of the model must be within factor 2 of 90 min, i.e. `[45, 180] min`.
4. A claim that the geometry improves *the prediction* on held-out data is only accepted if the same measurement curves, meal, volume and time axis are used and the geometry model reduces RMSE against at least one reference `T50` by at least 10 percent compared to scalar `k_empt`. Since no such measurements are found in this catalog, this criterion is `UNKNOWN`, not approved.
5. An approved solid food claim must also show residual solids over time; a better average alone is not enough.

## Controls and counter-tests

- Zero model: one room model with `V(t)=V0 exp(-k_empt t)` and separate fixed `M(t)=M0 exp(-k_s t)`.
- Diagnostic control: fixed analytical limit case without activity, with constant pressure and closed pylorus, should reproduce mass conservation.
- Diagnostic check: the orifice law should halve the flow when the pressure difference is quadrupled.
- No calibration against the target references is done in this run.
- Sensitivity ranking: each parameter is varied alone with `+50%` and `-50%`; reported as absolute change in both `T50` values.

## The frozen assumptions of the model

- SI-units internal; mL, mmHg and min are only converted on output.
- Fundus and antrum are compressible wall volumes with separate resting volumes, passive pressures and active contractions.
- Electrical activity is approximated by a sinusoidal slow-wave with 3/min and a positive contraction envelope function.
- The pylorus is a time-dependent orifice-like opening; liquid uses pressure flow, while solid particles require mobilization and size selection.
- No contact, fracture, microstructure, or 3D geometry parameters are identified; they are replaced by explicit assumptions and marked as such in `model.py`.

## Building on

- `inputs/QUESTION.md` (core mechanism, parameters and control tests).
- `inputs/NIGHT_PREAMBLE.md` (rules of conduct, zero model and pre-registration requirements).
- No previous implementation or results file was in `/opt/agents/jobs/BT-HX-Q107` at startup.

## Not redone

No internet access, no access to or modification of BodyTwin source reps, no assumptions about the collaborator/the reference model, no fabricated measurements, no LOSO sweeps and no empirical validation. The model is a controlled first run mechanistic model, not a substitute for 3D contact/tissue solution.
