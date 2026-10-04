# BT-HX-Q168 — preregistration

## Freeze and scope

This preregistration is frozen before the first model run on 2026-09-25. The bounded deliverable is a runnable, synthetic, first-principles tracer model for the causal chain

`oral L-carnitine tracer → microbial TMA production → epithelial passage → portal TMA → hepatic FMO3/TMAO → systemic clearance`.

No internal BodyTwin data, human samples, or new intervention results are used. The model is a measurement-design object, not a biological-effect estimate.

## Causal estimand and observables

The perturbation `A` is a 250 mg oral dose of d3-(methyl)-L-carnitine, with serial sampling at the time points specified in the model. The primary synthetic outcome is the systemic labelled-TMAO exposure

`Y = ∫_0^24 h C_plasma(TMAO-d3)(t) dt`

with unit `µM·h`. The secondary outcome is the total steady-state plasma TMAO concentration `C_ss` with unit `µM`. The observation operator is

`y(t) = [C_lumen(TMA-d3)(t), C_portal(TMA-d3)(t), C_plasma(TMAO-d3)(t), Q_urine(t)]ᵀ`,

where `Q_urine(t) = V_s ∫_0^t k_renal C_plasma(s) ds` is the cumulative systemic-loss amount in `µmol`. A plasma-only operator is the first-order placebo for identifiability: it cannot assign a terminal curve uniquely to production, passage, or metabolism.

The pre-specified mechanism perturbations are: microbial conversion capacity, epithelial effective permeability, host FMO3 capacity, and systemic clearance. A taxonomy-only operator is a negative control because it contains no tracer-flux information.

## Published reference

The numerical scale anchor is the median fasting plasma TMAO concentration of **4.6 µM** reported in the GeneBank cohort (`n = 2,595`). The same paper reports the oral d3-(methyl)-L-carnitine challenge dose of **250 mg**, serial plasma sampling, and 24 h urine collection.

Source: Koeth RA, Wang Z, Levison BS, et al. (2013), *Intestinal microbiota metabolism of L-carnitine, a nutrient in red meat, promotes atherosclerosis*, **Nature Medicine** 19:576–585, DOI **10.1038/nm.3145**. Page consulted: https://pmc.ncbi.nlm.nih.gov/articles/PMC3650111/ (abstract, Figure 1 caption, Figure 4f, Methods). The value is **verified**; it is not labelled `UNVERIFIED`.

The source value is an observational total-TMAO median, whereas `Y` is a labelled-tracer AUC. The comparison is therefore a scale/anchor check, not independent validation of the tracer estimand. No exact human tracer AUC matching the synthetic parameterisation is asserted.

## Frozen acceptance criteria

1. **Reference scale:** `C_ss` must lie between 2.3 and 9.2 µM (within a factor of 2 of 4.6 µM). This is a numerical anchor check only.
2. **Conservation:** all compartment balances must close with relative residual `<1e-8` over the 24 h integration.
3. **Analytical test:** the linear first-order chain must agree with its closed-form solution to `<1e-8` maximum absolute error.
4. **Mechanistic separation:** the full multi-site Fisher-information determinant must exceed the plasma-only determinant by at least `1.5` times, with no parameter reduced to zero in the finite-difference design.
5. **Epistemic status:** biological identification remains `UNKNOWN` without independent, paired interventional observations. A numerical criterion pass does not upgrade that status.

## Counterfactual and negative controls

- `Vmax_microbe = 0`: labelled TMAO exposure must approach zero while precursor remains represented.
- `P_epithelial = 0`: labelled TMAO exposure must approach zero while lumen TMA can remain.
- `Vmax_FMO3 = 0`: labelled TMAO exposure must approach zero while portal TMA can remain.
- taxonomy-only operator: no sensitivity to any tracer-flux parameter; it is not an acceptable substitute for the paired operator.
- urine is a systemic-loss observation and cannot by itself identify the microbial production rate.

## Frozen parameter policy

All values not directly reported as the challenge dose or the 4.6 µM reference are synthetic assumptions, not measurements. They are listed with unit and provenance in `model.py`. The epithelial flux is a Fick-type flux with explicit area and effective permeability, followed by a saturating factor. This keeps geometry, transport, kinetics, and conservation visible rather than replacing them with a correlation.

## What counts as failure

- A terminal plasma curve is interpreted as a production or metabolism result without a compartment-specific contrast.
- Units, dose conversion, or mass balance are inconsistent.
- A reference value is silently treated as a tracer AUC match.
- A synthetic number is presented as a measured human value.
- A taxonomy-only improvement is claimed as a mechanistic improvement.

## Builds on

- `inputs/QUESTION.md`: Q168 measurement-operator and identifiability target.
- `inputs/NIGHT_PREAMBLE.md`: bounded-run, precision, and no-invented-data requirements.
- Koeth et al. 2013, DOI 10.1038/nm.3145: published human challenge design and 4.6 µM plasma-TMAO anchor.

## Not redone

- external solver runtime, restricted model data, OrthoLoad, GRF, joint-force models, and other unrelated prior nodes are not used.
- No public dataset is treated as an intervention result.
- No real antibiotic, permeability, FMO3, or dietary intervention is recommended or performed; those are future controlled comparisons.
