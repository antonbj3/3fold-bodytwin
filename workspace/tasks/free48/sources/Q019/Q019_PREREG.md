# BT-HX-Q019 — Preregistration

This file is frozen before the first numerical run of the bounded Q019 model.

## Question and estimand

Q019 asks which variation between nephron segments changes the organ-level observable response over filtration, transport, and excretion time. The model represents one superficial nephron and scales its terminal output to one kidney. The independently observable organ outputs are terminal urine flow, urinary NaCl-equivalent excretion, urinary urea excretion, and terminal urine osmolality. Segment states are latent model variables and are not treated as measurements.

For a controlling parameter `p`, the frozen local sensitivity is

`S_p = (Y(+50% p) - Y(-50% p)) / (2 Y(0))`

where `Y` is one independently measured organ output. The ranking is by `abs(S_p)`, with ties below `0.02` treated as unresolved.

## Prediction frozen before execution

With the default antidiuretic state (`avp = 1`) and all other parameters held fixed:

1. A ±50% change in collecting-system water permeability (`CD_Lp`) will be the largest controller of terminal urine flow and urine osmolality.
2. A ±50% change in thick-ascending-limb active NaCl reabsorption (`TAL_vmax`) will be the largest controller of urinary NaCl-equivalent excretion.
3. A ±50% change in proximal-tubule water permeability (`PCT_Lp`) will be a weaker controller of both outputs because proximal transport is upstream and partly buffered by the axial flow.

These are model predictions, not measurements. The expected winner is re-evaluated from the computed sensitivities rather than hard-coded into the result.

## Reference value looked up

Primary research source: Layton AT, Layton HE. *A computational model of epithelial solute and water transport along a human nephron*. PLOS Computational Biology. 2019;15(2):e1006108. DOI: `10.1371/journal.pcbi.1006108`. PMCID: `PMC6405173`.

The source was checked in its full text. The following values are transcribed, not invented:

- Table 1, “Nephron segment lengths and luminal diameters”: proximal tubule `1.7 cm`, descending limb `0.32 cm`, medullary thick ascending limb `0.5 cm`, cortical thick ascending limb `0.5 cm`, distal convoluted tubule `0.2 cm`, connecting tubule `0.4 cm`, cortical collecting duct `0.4 cm`, outer-medullary collecting duct `0.5 cm`, inner-medullary collecting duct `1.2 cm`; corresponding listed diameters are `37, 26, 26, 26, 20, 24, 45, 45, 50 micrometre`.
- Table 2, “Interstitial concentrations”: cortex osmolality `285.9 mosm/(kg H2O)`, outer-medullary/inner-medullary boundary `688.5 mosm/(kg H2O)`, papillary tip `737.9 mosm/(kg H2O)`.
- Table 3, “Baseline”: urine flow `0.62 ml/min` per kidney, Na+ `73 mM`, K+ `72 mM`, Cl− `69 mM`, urea `280 mM`, pH `6.2`.
- The article’s model-parameter text states single-nephron GFR `100 nl/min` and approximately one million nephrons per human kidney.

The reference values are anchors for plausibility and reporting, not fitted observations from the present run. The present model is a deliberately reduced first-principles model and is not expected to reproduce every electrolyte in the source model.

## Mechanistic model frozen for this run

Each segment is a well-mixed axial compartment. Time is in minutes, volume is in mL per nephron, solute amount is in micromol per nephron, and concentration is in mM. Glomerular filtration is a net-pressure filtration law,

`J_f = K_f * [(P_gc - P_bs) - (pi_gc - pi_bs)]`.

The segment states obey mass conservation,

`dV_i/dt = J_in,i - J_out,i - J_water,i`,

`dM_k,i/dt = J_in,i C_k,in - M_k,i/tau_i - J_k,i`,

with `J_out,i = V_i/tau_i`. Water uses the osmotic part of the Kedem–Katchalsky law,

`J_water,i = L_p,i A_i (C_interstitial,i - C_lumen,i)`,

where `1 mM` is converted to `19.37 mmHg` at 37 degrees C. NaCl-equivalent and urea transport use saturable active reabsorption,

`J_k,i = Vmax_k,i C_k,i / (K_m,k,i + C_k,i)`,

with an explicit assumption for any paracellular term. Concentrations are bounded by the compartment states and axial outflow is removed from the same state, so the terminal output is an independent measurement of the simulated organ response.

The source article supplies geometry, filtration scale, and interstitial anchors. The reduced transport coefficients, residence times, and lumped NaCl representation are explicit model assumptions; they are not claimed to be newly measured human values.

## Controls and countertests

- Null perturbation: run the model twice with the same nominal parameters and require identical terminal observables to numerical precision.
- No-transport limit: set every segment water and solute permeability/active flux to zero; the terminal flow must equal the filtration input and total water must be conserved.
- Strongest-baseline check: compare the default model with a placebo that changes only an output label and with the zero-transport analytical limit; a segment perturbation is not accepted if its effect is smaller than the numerical floor or reverses the sign without a physical explanation.
- Sensitivity is one-factor-at-a-time at exactly ±50%, with the same initial state, duration, solver tolerances, and AVP state.

## Frozen acceptance criteria

The run is accepted only if all of the following hold:

1. The unit check reports `PASS` for the filtration, water, concentration, and solute-flux conversions.
2. The maximum relative water-mass residual over the recorded trajectory is below `1e-6` and every terminal flow/concentration is finite and non-negative.
3. The zero-transport limit returns terminal aggregate flow equal to `100 ml/min` per kidney within `1e-8` relative error.
4. The CD water-permeability perturbation has the largest absolute normalized sensitivity for urine flow, and the TAL active-transport perturbation has the largest absolute normalized sensitivity for NaCl-equivalent excretion, unless the run reports a tie or a failed criterion rather than silently changing the prediction.
5. The null perturbation is below the solver’s numerical resolution and no result is described as a measured clinical value.

## What counts as an error

A negative compartment volume or amount, a unit mismatch, a mass-balance failure, a non-finite output, a post hoc change to the acceptance rule, or a reported value without a corresponding entry in `results.json` is an error. A sensitivity that fails to meet its preregistered ranking is a failed prediction, not grounds for retuning.

## Builds on

- Node-id: none; Q019 had no earlier implementation in this directory. The interrupted session had only inspected `BRIEF.md` and `inputs/`.
- Files read: `BRIEF.md`, `inputs/NIGHT_PREAMBLE.md`, `inputs/QUESTION.md`, and `ALLOW_WEB`.
- Published anchor used: Layton & Layton 2019, DOI `10.1371/journal.pcbi.1006108`, full text and Tables 1–3.

## Not redone

- No external solver runtime, external repository edits, subagents, or web-derived data files are imported.
- The model does not attempt to identify a person, diagnose disease, or infer a clinical effect from the synthetic perturbations.

## Execution

Python with NumPy and SciPy is used for a short, single-threaded deterministic run. The result is written to `results.json`; the interpretation and limitations are written to `RESULTS.md`.
