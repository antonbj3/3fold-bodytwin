# PREREG — BT-HX-Q013

## Frozen before the first run

Status: `FROZEN_BEFORE_MODEL_EXECUTION`  
Model version: `BT-HX-Q013-0.1`  
Date: 2026-09-25

This is a mechanism model, not an empirical effect model. No BodyTwin, patient or other internal measured data are used. No parameters are adjusted after the model has been run.

## Question and hypothesis

How can a combination of a pharmacokinetic interaction (changed exposure for the victim drug) and changed cell sensitivity (changed coupling from local exposure to effect) be distinguished?

Hypothesis: if the exposure increase is entirely pharmacokinetic, the effect of the measured local exposure follows an unchanged cell model. If the combination gives an effect outside that model at the same local exposure, there is a cell-sensitivity component.

## Predicted quantity and primary reference

Primary model quantity:

`R_local = AUC(C_local, combination) / AUC(C_local, drug alone)`.

Secondary quantities:

- `R_plasma = AUC(C_plasma, combination) / AUC(C_plasma, drug alone)`.
- `E_PK(C_local)` = computed effect with the combination's local exposure but unchanged cell sensitivity.
- `delta_E_PD = E_kombination - E_PK` with the same local exposure.
- `gamma` = the combination's sensitivity factor; `gamma = 1` means unchanged cell sensitivity.

The first numerical state is an F-mediated PK interaction: `R_plasma = 1,27` and therefore, under the model's linear and stationary conditions, `R_local = 1,27`. The value is a reference calibration, not a measured local exposure.

## Reference value

Primary source: Tachibana M, Horwitz S, Jacobsen E, et al. “Effect of Valemetostat on the Pharmacokinetics of Midazolam and Digoxin: A Phase 1 Drug–Drug Interaction Study in Patients With Non-Hodgkin Lymphoma.” *Clinical and Translational Science*. 2025;18(11):e70330. DOI: `10.1111/cts.70330`. PMID: `41150706`.

Table 3, comparison of PK parameters:

- drug: digoxin 0,25 mg as reported, alone compared with digoxin 0,25 mg + valemetostat 200 mg once daily;
- `AUC_last` GMR: **1,27**;
- 90-percent CI: **1,06–1,52**;
- `n = 16` for AUC_last;
- `Cmax` GMR: 1,30 (90-percent CI 1,07–1,57).

The source measures systemic plasma exposure, not independent cell exposure. Therefore `R_local = 1,27` is a model prediction under the assumption that the F change transfers linearly to the cell compartment. `AUC_inf` is not used as calibration because the source gives only `n = 3` for that comparison.

A separate older primary source, Tsunoda SM et al., *Clin Pharmacol Ther* 1999;66:461–471, DOI `10.1016/S0009-9236(99)70009-3`, PMID `10579473`, reports in its public abstract a 5-fold intravenous and 16-fold oral AUC increase for midazolam with ketoconazole. It is used as a mechanistic quality check, not as a calibration value for the executable model.

## Frozen criterion

1. **Calibration:** with `bioavailability_multiplier = 1,27`, `clearance_multiplier = 1` and unchanged absorption rate, the numerical `R_plasma` must be `1,27 ± 0,01`; it must lie within the source's 90-percent CI `[1,06, 1,52]`.
2. **PK-only:** if the combination's local exposure curve is used and `gamma = 1`, the measured/predicted residual `|delta_E_PD| / Emax <= 0,02` must count as PK-explained in this first pass.
3. **Cell sensitivity:** at the same local exposure, a deviation `|delta_E_PD| / Emax > 0,05` must count as a sensitivity component. Between `0,02` and `0,05`, the result is classified as `UNRESOLVED`.
4. **Exposure:** plasma may be used as a proxy only to compute `R_plasma`; local cell concentration or an independent local measure is required for cell-sensitivity classification.
5. **Countertest:** gamma must be tested at `1,0` and `1,5`; the PK parameter must be tested at `0,5` and `1,5` times the nominal value. No parameter may be chosen after the residual criterion has been evaluated.

## The model's physical structure

### Exposure

The victim drug is given as extraction into a central compartment and, for the oral route, a first-order absorption compartment:

`dA_g/dt = -ka A_g`  
`dA_p/dt = F ka A_g - CL A_p/V_p`  
`C_p = A_p/V_p`.

The interaction is factorized into two observable/composable PK branches:

`F_eff = F * F_mult`  
`CL_eff = CL * CL_mult`.

For a linear system with the same `ka`, `V_p` and complete sieving law, `AUC_inf` is proportional to `F/CL`, so

`R_plasma = F_mult / CL_mult`.

The reference case is frozen at `F_mult = 1,27`, `CL_mult = 1`. This is an identification assumption, not a claim that the source can distinguish these two branches.

### Local exposure

The cell compartment uses mass balance across a local volume segment:

`dC_local/dt = k_in C_p - k_out C_local`.

`k_in` and `k_out` have unit `h^-1`; the concentrations have the same normalized concentration unit. `C_local` is therefore not equated with plasma exposure.

### Cell sensitivity and effect

The receptor/target-site saturation part is a Hill occupation:

`O = C_local^h / (EC50^h + C_local^h)`.

The effect is

`E = E0 + gamma Emax O`.

`gamma` is a signal amplification factor after activation. `gamma = 1` is unchanged sensitivity; `gamma > 1` is increased sensitivity. `EC50` and `Emax` are kept separate from the PK parameters so that sensitivity cannot be confused with exposure.

## Assumptions that are not measured data

- `F`, `ka`, `CL`, `V_p`, `k_in`, `k_out`, `EC50`, `Emax`, `h` and `E0` are normalized demonstration values.
- The only published calibration parameter is `R_plasma = 1,27` from the primary source above.
- No measured cell exposure or cell sensitivity exists in the job package. `gamma = 1,5` is a sensitivity scenario, not an estimate from study participants.
- The source study's patient population and administration form are not used to claim cell response.

## Builds on

- `inputs/QUESTION.md`: Q013's contract — separate exposure and effect in individual and combined courses, and local/independent exposure measurement.
- `inputs/NIGHT_PREAMBLE.md`: requirements for frozen registration, sources, provenance, `UNKNOWN` and the next resolution step.
- `K03`, `K08`, `K11`: specified inputs in the brief; no node files or external BodyTwin sources were available in this limited job directory, and no code has been copied from them.

## Not redone

- No external solver runtime, no Batch/LOSO sweep, no cloud run and no internal data.
- No patient model, no measured data and no claimed local exposure are created.
- No separation between absorption and clearance changes is calibrated from a single AUC ratio; that mechanism is left `UNKNOWN`.
- No cell type, receptor or disease is chosen without measured evidence; the model is therefore drug- and cell-type agnostic.

## What counts as an error

- calling a plasma AUC a local cell-concentration measure;
- using the source's `AUC_inf` with `n = 3` as the main calibration;
- changing `gamma`, `EC50`, `F_mult` or `CL_mult` after the criteria have been frozen;
- reporting a synthetic scenario as an observed effect;
- reporting `PK-only` or `cell sensitivity` when local exposure is missing.
