# BT-HX-Q021 — preregistration

## Scope and hypothesis

This is a bounded, first-principles regional lung model for gas exchange under a load-like change from rest to exercise. The model uses three or more lumped regions, a pressure/resistance perfusion description, a compliance/resistance ventilation description, hypoxic pulmonary vasoconstriction (HPV), an O2 content curve, and an alveolar gas balance. It is a mechanism test, not a clinical or population estimate.

Hypothesis: a single homogeneous compartment is adequate only when regional ventilation and perfusion are sufficiently uniform. When load changes their relative distributions, regional gas exchange produces a materially different mixed arterial O2 result from a homogeneous compartment with the same total ventilation and perfusion.

## Frozen prediction and reference values

Primary endpoint: the model-predicted alveolar-arterial O2 difference, `A-aDO2 = PAO2_global - PaO2_mixed`, in Torr, at rest and at `exercise_fraction = 1.0` (heavy-load analogue). The primary comparison is the heavy-load regional model with the heavy-load homogeneous null model.

Reference anchors (public primary sources; values were checked and are not OVERIFIERAD):

1. Hall ET et al. *The effect of supine exercise on the distribution of regional pulmonary blood flow measured using proton MRI.* Journal of Applied Physiology 2014;116:451-461. DOI `10.1152/japplphysiol.00659.2013`. Table 4 reports perfusion in `ml min^-1 ml^-1`: nondependent `2.9 ± 1.7` at rest and `4.0 ± 1.5` during exercise; middle `4.5 ± 1.7` and `5.7 ± 2.1`; dependent `4.2 ± 1.7` and `5.1 ± 1.6`. The derived nondependent/dependent ratios are `0.690` at rest and `0.784` during exercise.
2. Tedjasaputra V et al. *The heterogeneity of regional specific ventilation is unchanged following heavy exercise in athletes.* Journal of Applied Physiology 2013;115:126-135. DOI `10.1152/japplphysiol.00778.2012`. Table 3 reports `A-aDO2 = 6.3 ± 3.7 Torr` at rest and `23.3 ± 5.3 Torr` averaged during exercise, together with arterial blood gases. These are external observations, not calibration data.
3. Harf A, Pratt T, Hughes JM. *Regional distribution of VA/Q in man at rest and with exercise measured with krypton-81m.* Journal of Applied Physiology 1978;44:115-123. DOI `10.1152/jappl.1978.44.1.115`. The indexed abstract reports a `40–150%` apical blood-flow increase within 30 s of 50 W cycling in six subjects. This is contextual support, not a fitted value.

## Frozen run rules

- Use the default nine-region parameter set in `model.py`, one thread, and no internal observations.
- Run rest (`exercise_fraction=0`) and heavy-load (`exercise_fraction=1`) regional cases plus matched homogeneous null cases.
- Primary quantitative screen: heavy-load regional `A-aDO2` must be at least `10 Torr` above the regional rest value and must lie in the broad `10–40 Torr` observation-screening interval around the published heavy-exercise value. This is a preregistered plausibility screen, not a claim of validation.
- Mechanistic need criterion: the absolute heavy-load regional-minus-homogeneous difference in `A-aDO2` must be at least `5 Torr`. The model is declared regionally informative only if this threshold is met; a failure is reported as a failure.
- Distribution check: the heavy-load nondependent/dependent perfusion ratio must move toward the published exercise ratio `0.784` relative to rest `0.690`; the frozen tolerance is `0.20` absolute ratio units.
- The zero-heterogeneity/no-HPV case must be numerically identical to the homogeneous null case to `1e-12` for its gas outputs.
- Conservation and validity failures, non-finite values, negative flows, or a failed analytical test count as errors. No parameter may be changed after seeing the run to make a criterion pass.

## Countermodels and controls

- `homogeneous`: every region receives `1/N` of alveolar ventilation and perfusion, while total flows, metabolic inputs, and gas equations are identical to the regional case.
- `no_heterogeneity`: regional gradients set to zero and HPV set to zero; this is the numerical placebo.
- `no_ventilation_region`: a zero-alveolar-ventilation unit is tested as a shunt limit; its blood leaving the unit must equal mixed venous blood while its O2 exchange is zero.
- `uniform_flow_limits`: when perfusion is zero, the local O2 exchange must equal inspired alveolar O2; when ventilation is zero, the local O2 exchange must approach the venous content.

## Builds on

- `BRIEF.md`, `inputs/NIGHT_PREAMBLE.md`, and `inputs/QUESTION.md` in this job directory.
- The two open, public primary studies above for external anchors; no internal BodyTwin data were available or used.
- The equations are an explicit reduced model: hydrostatic pressure difference divided by regional vascular resistance for perfusion, pressure-driven compliance/resistance flow for ventilation, an O2-content/Hill relation, and an alveolar O2 mass balance.

## Not redone

