# BT-HX-Q100 preregistration

## Scope and frozen hypothesis

This is a bounded, no-internet first-principles comparison of a two-mass vocal-fold oscillator with and without a compliant layered tissue, one-sided contact, and a mucus-film surrogate. The primary predicted quantity is the fundamental frequency `f0` of the autonomous, symmetric phonation limit cycle at a fixed subglottal pressure `Psub = 800 Pa`. Secondary outputs are transglottal pressure, characteristic mucus-wave speed, contact fraction, and the displacement/area oscillation amplitude.

The preregistered mechanistic hypothesis is:

> At the same pressure, geometry, masses, air state, and measurement window, adding the thin surface layer, explicit one-sided contact, and mucus-film coupling will change the predicted `f0`, contact fraction, and mucosal-wave proxy relative to the calibrated two-mass reference. This is a mechanistic improvement only if it lowers the error to the frozen reference benchmark; extra detail alone is not counted as improvement.

The comparison is a prediction exercise, not a clinical or physiological validation. No measured pressure, geometry, video, or mucosal-wave data are used.

## Reference value

The order-of-magnitude benchmark is approximately `120 Hz`, from the classic two-mass study by N. Ishizaka and J. L. Flanagan, “Two-mass model of the vocal folds,” *Journal of the Acoustical Society of America* (1972), representative phonation condition. The numerical value and the precise pressure/geometry associated with it are **UNVERIFIED, FROM MEMORY**; no online verification was possible. It is therefore used only as a broad mechanical benchmark, not as a matched individual measurement.

## Frozen approval criterion

The primary model receives `PASS` for the broad benchmark only if, after the same fixed transient and measurement window for both variants:

1. `60 Hz <= f0_layered <= 240 Hz` (within a factor of two of the memory benchmark), and
2. the final-window fundamental is accompanied by a finite, non-degenerate oscillation with peak-to-peak displacement at least `5 um`.

The benchmark subcriterion is frozen before execution. A value outside the interval is `FAIL`, and a non-finite or non-oscillatory solution is `FAIL`/`UNKNOWN` rather than being retuned. The separate improvement flag is `PASS` only when `abs(f0_layered - 120 Hz) < abs(f0_baseline - 120 Hz)`; otherwise it is `NO`. Contact and mucus-wave outputs are reported as `UNKNOWN` against empirical references because no matched data are available.

## Fixed protocol

- `Psub = 800 Pa`; zero external acoustic pressure in this first run.
- The same two lateral generalized masses, rest area, area derivatives, air density, flow-resistance parameters, coupling stiffness, and measurement window are used in both variants.
- The layered variant changes only the tissue law, the contact indentation law, and the mucus-film force. No parameter is fitted to the frequency benchmark.
- Integration begins from a fixed small outward displacement, runs for `0.24 s`, and uses the final `0.08 s` for frequency and contact statistics.
- Frequency is the dominant positive Fourier peak of the mean-wall displacement after mean removal. Contact percentage is the fraction of final-window samples with raw glottal area at or below the contact threshold.
- Mucus-wave speed is a declared reduced proxy, `c_m = sqrt(T_sheet/(rho_m h_m))`, not a measured wave-tracking speed.
- Sensitivity is one-at-a-time `+50%` and `-50%` perturbation of exactly three preselected parameters: surface stiffness, contact stiffness, and subglottal pressure. The signs and rankings are reported; no parameter is selected after seeing a favorable result.

## Null models and falsification

- **Null/baseline:** homogeneous two-mass tissue law with the same pressure, geometry, masses, and air coupling; it is the comparator requested by the question.
- **Placebo:** run the same layered equations with the surface-layer correction, contact, and meniscus terms set to zero. This must reproduce the baseline within numerical tolerance.
- **Falsifier:** a layered result is not called an improvement if the placebo is not reproduced, if the frequency is not measurable, or if the result depends on an undocumented retuning.
- A useful mechanistic distinction requires a non-negligible change in at least one of `f0`, contact fraction, or mucus-wave proxy, while the placebo remains at baseline.

## Data and source separation

- `inputs/QUESTION.md` defines the requested mechanism and output quantities; it is a task specification, not a measurement source.
- `inputs/NIGHT_PREAMBLE.md` supplies the write/resource constraints and the requirement to use explicit provenance.
- The reference benchmark is a memory citation, explicitly unverified.
- All numerical predictions below are derived outputs of this local model, not measured data.
- No internet, external dataset, or the reference model/runtime component is used because the task is bounded to this directory.

## Builds on

- `BT-HX-Q100` and the question crosslinks `Q049`, `Q050`, `Q055`, `Q101`, and source entry `S02` are retained as scope references only.
- No prior code, numerical result, geometry, or reusable vocal-fold model exists in this task directory at inspection time; the only pre-existing files are `BRIEF.md`, an empty `ALLOW_WEB`, and `inputs/NIGHT_PREAMBLE.md` plus `inputs/QUESTION.md`.

## Not redone

- No empirical calibration to a subject or recording.
- No 3D CFD, acoustic radiation, turbulence closure, heat/mass transfer inside tissue, or 3D contact geometry.
- No claim that the mucus proxy is a measured mucosal-wave velocity.
- No claim that the memory benchmark identifies this model's exact operating point.
