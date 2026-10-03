# BT-HX-Q023 — preregistration

## Scope and data boundary

This is a first executable, reduced-order model of the human hypothalamic–pituitary–adrenal (HPA) feedback axis after an upstream hormone intervention. It represents a CRH-like drive, pituitary ACTH, biologically available free cortisol, and a slower total-cortisol compartment. It does not claim that the supplied package contains an observed individual trace: the inputs contain no hormone measurements. All generated trajectories are model predictions, not measurements.

The intervention is a rectangular normalized drive pulse of amplitude 1 for 10 min at t = 0. The target is an ACTH-to-cortisol phase lag, a transient cortisol response, and delayed negative feedback visible in the subsequent ACTH trace.

## Builds on

- `inputs/QUESTION.md`, Q023: represent phase, amplification, and transient in one specified hormone axis with simultaneous hormones and intervention time.
- `inputs/NIGHT_PREAMBLE.md`: preregistration-before-run, source/assumption separation, and no invented measurement data.
- `Q023` and the referenced `K03/K07/K09` labels are not present as files in this bounded directory, so no additional BodyTwin claim or data is silently imported.
- `agent.log` confirms that the earlier session stopped after reading the brief; no prior model or result is available to reuse.

## Not redone


## Hypothesis

A delayed HPA feedback loop is identifiable only if the downstream signal is generated after the upstream signal and if cortisol inhibits ACTH after a nonzero transport/feedback delay. With source-anchored ACTH-to-cortisol lag near 3 min and free-cortisol half-life near 2.2 min, a feedback delay in the several-minute range should produce a measurable phase ordering and a transient; setting all delays to zero should not be accepted merely because it fits an instantaneous curve.

## Source-anchored reference values

1. **Dorin RI, Qiao Z, Qualls CR, Urban FK III. 2012. “Estimation of maximal cortisol secretion rate in healthy humans.” Journal of Clinical Endocrinology & Metabolism 97:1285–1293. DOI: 10.1210/jc.2011-2227.** Primary randomized cosyntropin study, n = 21 healthy adults, frequent sampling after 250 microgram ACTH1–24. The Results accompanying **Figure 2** report a mean intervention-to-cortisol response lag of **3.0 +/- 0.92 min**. **Figure 5** and the same Results section report free-cortisol half-life **2.2 +/- 1.3 min** (pooled placebo/DEX conditions; abstract reports 2.2 +/- 1.1 min). CSRmax is **0.44 +/- 0.13 nmol L^-1 s^-1**. These are reference anchors, not values fitted to an internal trace.
2. **Keenan DM, Roelfsema F, Veldhuis JD. 2004. “Endogenous ACTH concentration-dependent drive of pulsatile cortisol secretion in the human.” American Journal of Physiology—Endocrinology and Metabolism 287:E652–E661. DOI: 10.1152/ajpendo.00167.2004.** Primary study, 32 healthy adults sampled every 10 min for 24 h. Its abstract reports total-cortisol half-life **49 +/- 2.4 min**, ACTH half-life **20 +/- 1.3 min**, and free-cortisol half-lives **1.8 +/- 0.2 min** and **4.1 +/- 0.3 min**. The two free components motivate representing free cortisol separately from a slower total pool.
3. **Angst MS et al. 1998. “Pharmacokinetics, cortisol release, and hemodynamics after intravenous and subcutaneous injection of human corticotropin-releasing factor in humans.” Clinical Pharmacology & Therapeutics 64:499–510. DOI: 10.1016/S0009-9236(98)90133-3.** Primary human intervention study; the abstract reports hCRF elimination half-life **45 +/- 7 min** after a 10-min IV infusion and **37 +/- 10 min** after a 180-min IV infusion. This is used only to keep the upstream delay biologically motivated, not as a direct estimate of portal transit.

## Frozen prediction and acceptance rules

- **Prediction P1 (phase):** with the frozen default parameters, the free-cortisol response onset must occur after the ACTH response by approximately **3 min** (allowed model tolerance: 1–5 min because the delay is represented by a distributed first-order line).
- **Prediction P2 (transient):** free cortisol must peak before the slower total-cortisol proxy, and the ACTH trace must show a delayed suppression/recovery when feedback is enabled. A monotone constant solution is not an acceptable feedback demonstration.
- **Prediction P3 (gain):** increasing loop gain by 50% must increase feedback suppression or rebound contrast relative to the default; if a numerical trajectory is saturated, the direction is reported as saturated rather than hidden.
- **Prediction P4 (sensitivity):** vary `delay_ac`, `delay_fa`, and `feedback_gain` independently by +/-50%, and report peak times, phase lag, overshoot/rebound, and stability for every case.
- **Falsifiers:** wrong ordering (cortisol before ACTH), failure to return toward baseline after the pulse, numerical instability, or a unit-check failure. A parameter set that cannot exhibit the observed ordering is not rescued by tuning.
- **Interpretation boundary:** “supports the mechanism” means only that the reduced model reproduces the preregistered phase/transient criteria. It does not identify a person, diagnose a disorder, or establish a causal delay from the absent internal data.

## Frozen default parameterization

All concentrations are dimensionless relative to a reference state. Time is in minutes; rates are min^-1. The implementation defaults are `tau_r = 5 min`, `tau_a = 20 min`, `tau_f = 2.2 min`, `tau_total = 49 min`, `delay_ra = 2 min`, `delay_ac = 3 min`, `delay_fa = 8 min`, `feedback_gain = 0.9`, `feedback_scale = 0.5`, Hill exponent 2, and a 10-min normalized upstream pulse. Source values and assumptions are marked separately in the model parameter table. The delay representation is a cascade of first-order transport/filtration compartments, so reported delay is an effective distributed delay rather than a hidden pure delay.

## Analysis plan

Run the deterministic model once at the frozen defaults, run the analytical no-feedback/no-input limit test, then run the 3 x 3 x 3 sensitivity grid (three controlling parameters at 50%, 100%, and 150% of default). Save only scalar predictions and provenance to `results.json`; do not label any scalar as a measured value. The next resolution step is to fit the same state equations to preregistered, timestamped simultaneous CRH/ACTH/cortisol measurements with assay and compartment metadata.
