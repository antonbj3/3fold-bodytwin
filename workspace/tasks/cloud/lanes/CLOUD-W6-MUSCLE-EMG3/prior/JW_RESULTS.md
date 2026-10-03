BT-PG-189
# BT-PG-189

## Status and scope

This is a single-subject, observational measurement package for JW. It is not a population estimate and does not compare JW with another person. The source metadata is *Grand Challenge Competition to Predict In Vivo Knee Loads* (DOI `10.1002/jor.22023`). The canonical machine-readable outputs are `results.json` and the single combined numerical table `measurements.csv`; both were produced by one run of `python3 extract.py` after the method was fixed in `PREREG.md`.

The local frozen-source integrity check completed successfully for all 6 source files (`python3 inputs/load_test.py`). No independent reference target or held-out test values were supplied, so `reference_test` is `NOT_PERFORMED`; this package makes no reference or model-validation claim.

## Combined package table

The table below is the requested one-table summary. Curve rows point to the complete 10-degree-bin curves below and in `measurements.csv`.

| output | value | unit | method | non-circular control |
|---|---:|---|---|---|
| 60°/s flexion curve peak | 71.046 at 59.400° | Nm | largest absolute qualified signed torque; full 10° curve is tabulated below | 3 full qualified flexion segments; bin IQR and n retained |
| 60°/s extension curve peak | 36.551 at 63.600° (signed −36.551) | Nm | largest absolute qualified signed torque; full 10° curve is tabulated below | 3 full qualified extension segments; bin IQR and n retained |
| 90°/s flexion curve peak | 52.861 at 61.650° | Nm | largest absolute qualified signed torque; full 10° curve is tabulated below | 2 full qualified flexion segments; bin IQR and n retained |
| 90°/s extension curve peak | 31.855 at 84.704° (signed −31.855) | Nm | largest absolute qualified signed torque; full 10° curve is tabulated below | 2 full qualified extension segments; bin IQR and n retained |
| 60/90 moment ratio, flexion | median 1.367; IQR 1.325–1.369 | dimensionless | median of 9 common 10° bins of `abs(T60)/abs(T90)` | peak ratio 1.344; paired bins and n retained |
| 60/90 moment ratio, extension | median 1.282; IQR 1.145–1.604 | dimensionless | median of 7 common 10° bins of `abs(T60)/abs(T90)` | peak ratio 1.147; paired bins and n retained |
| 60°/s implant-load delay, flexion | 9.40 (phase 1.49°) | ms | total load `PM+AM+AL+PL` versus torque; frequency-weighted phase, equivalent lag | 3 full segments; quadrant phases retained in `results.json` |
| 60°/s implant-load delay, extension | 412.88 (phase −111.22°) | ms | total load `PM+AM+AL+PL` versus torque; frequency-weighted phase, equivalent lag | 3 full segments; direction/quadrant controls retained |
| 90°/s implant-load delay, flexion | 66.32 (phase −25.38°) | ms | total load `PM+AM+AL+PL` versus torque; frequency-weighted phase, equivalent lag | 2 full segments; quadrant phases retained in `results.json` |
| 90°/s implant-load delay, extension | 179.88 (phase −74.57°) | ms | total load `PM+AM+AL+PL` versus torque; frequency-weighted phase, equivalent lag | 2 full segments; direction/quadrant controls retained |
| 60°/s EMG-to-torque phase, flexion | 148.48° (equivalent −412.79 ms) | deg (ms equivalent) | source-filtered vaslat EMG versus torque; 0.2–2 Hz weighted phase | 3 full segments; raw EMG correlation r=0.004 |
| 60°/s EMG-to-torque phase, extension | −179.38° (equivalent −641.78 ms) | deg (ms equivalent) | source-filtered vaslat EMG versus torque; 0.2–2 Hz weighted phase | 3 full segments; raw EMG correlation r=0.004 |
| 90°/s EMG-to-torque phase, flexion | 169.58° (equivalent −442.25 ms) | deg (ms equivalent) | source-filtered vaslat EMG versus torque; 0.2–2 Hz weighted phase | 2 full segments; raw EMG correlation r=−0.056 |
| 90°/s EMG-to-torque phase, extension | −169.10° (equivalent 447.28 ms) | deg (ms equivalent) | source-filtered vaslat EMG versus torque; 0.2–2 Hz weighted phase | 2 full segments; raw EMG correlation r=−0.056 |

## Synchronization, overlap, and units

| stream | 60°/s recording | 90°/s recording |
|---|---|---|
| Biodex torque/angle/velocity | 1,327 samples, 0.000–11.050 s, median Δt 0.008333 s | 773 samples, 0.000–6.433333 s, median Δt 0.008333 s |
| Implant force quadrants | 1,327 samples, 0.000–11.050 s, median Δt 0.008333 s | 773 samples, 0.000–6.433333 s, median Δt 0.008333 s |
| Dedicated EMG | 11,058 samples, 0.000–11.057 s, median Δt 0.001000 s | 6,442 samples, 0.000–6.441 s, median Δt 0.001000 s |
| Common synchronized interval | 0.000–11.050 s; 11.050 s; 1,327 Biodex samples | 0.000–6.433333 s; 6.433333 s; 773 Biodex samples |
| Biodex–force timestamp residual | median 0 s; maximum absolute 0 s | median 0 s; maximum absolute 0 s |

Signals were aligned by timestamp; linear interpolation was used only when evaluating the force or dedicated EMG stream at Biodex timestamps. Torque, angle, velocity, and force received no additional digital filter or smoothing. The EMG phase uses the source column `Filtered vaslat EMG`; it is linearly detrended within each full segment and projected over 0.2–2.0 Hz without a new waveform filter. The dedicated raw `vaslat` channel is a scale/filter control, not the primary phase signal.

Time is in s, angle in degrees, velocity in degrees/s, and torque in Nm. The four implant channels are compressive quadrants; their sum is treated as N, although the supplied force header does not state a unit. EMG remains in source/arbitrary units because no calibration is supplied.

## Repetitions and velocity mask

The inclusive mask was `abs(velocity) >= 48°/s` at 60°/s and `abs(velocity) >= 72°/s` at 90°/s. A qualified segment is a maximal contiguous run; a full traversal spans at least 80°.

- 60°/s: 1,091 qualified samples in 7 segments; 3 full flexion and 3 full extension traversals plus one 10.5° partial extension segment at the start.
- 90°/s: 563 qualified samples in 5 segments; 2 full flexion and 2 full extension traversals plus one 58.0° partial extension segment at the start.

Partial edge segments remain in the sample-weighted primary curves, as specified in `PREREG.md`; the full/partial status and contributing segment counts are retained per bin.

## Torque–angle curves

Flexion and extension are shown separately because flexion torque is positive and extension torque is negative. Each cell is the signed median torque in Nm followed by the sample count in brackets. The 80% velocity-masked criterion is applied before binning.

| angle bin (°) | 60°/s flexion | 60°/s extension | 90°/s flexion | 90°/s extension |
|---|---:|---:|---:|---:|
| 20–30 | 42.58 [58] | −2.48 [59] | 36.12 [27] | −1.33 [22] |
| 30–40 | 56.92 [60] | −10.13 [60] | 43.82 [27] | −4.51 [28] |
| 40–50 | 64.81 [59] | −20.71 [60] | 47.41 [26] | −9.68 [27] |
| 50–60 | 68.17 [61] | −29.97 [60] | 49.80 [27] | −16.97 [36] |
| 60–70 | 67.20 [60] | −35.17 [60] | 50.73 [27] | −24.37 [41] |
| 70–80 | 65.62 [60] | −34.59 [61] | 48.00 [26] | −26.98 [39] |
| 80–90 | 61.16 [58] | −33.29 [60] | 44.64 [28] | −27.99 [42] |
| 90–100 | 56.74 [60] | −24.94 [64] | 40.40 [26] | −23.50 [39] |
| 100–110 | 39.39 [57] | −11.02 [74] | 29.72 [24] | −10.02 [40] |
| 110–120 | — | — | — | −10.34 [7] |

The largest absolute qualified sample torque and the corresponding angle are:

| speed | direction | peak absolute torque (Nm) | signed torque (Nm) | peak angle (°) | largest absolute binned median (Nm; bin centre °) |
|---:|---|---:|---:|---:|---:|
| 60 | flexion | 71.046 | 71.046 | 59.400 | 68.168; 55 |
| 60 | extension | 36.551 | −36.551 | 63.600 | −35.166; 65 |
| 90 | flexion | 52.861 | 52.861 | 61.650 | 50.727; 65 |
| 90 | extension | 31.855 | −31.855 | 84.704 | −27.995; 85 |

The 60°/s flexion curve has 533 qualified samples and the extension curve 558; the 90°/s curves have 242 flexion and 321 extension samples. The 60/90 ratio uses only common angle bins with at least five samples in each curve and a 90°/s median magnitude of at least 5 Nm.

## Motprov: curves without the velocity mask

The following control uses every finite sample in each directional branch, including low-speed reversals and transition samples. It is not substituted for the criterion estimate. The complete IQR, bin n, and segment counts are in `measurements.csv`.

| angle bin (°) | 60°/s flexion | 60°/s extension | 90°/s flexion | 90°/s extension |
|---|---:|---:|---:|---:|
| 10–20 | 31.03 [65] | −1.62 [31] | 21.05 [60] | −0.09 [34] |
| 20–30 | 41.78 [66] | −2.37 [64] | 36.12 [27] | −0.36 [31] |
| 30–40 | 56.92 [60] | −10.13 [60] | 43.82 [27] | −4.51 [28] |
| 40–50 | 64.81 [59] | −20.71 [60] | 47.41 [26] | −9.68 [27] |
| 50–60 | 68.14 [62] | −29.97 [60] | 49.80 [27] | −16.97 [36] |
| 60–70 | 67.20 [60] | −35.17 [60] | 50.73 [27] | −24.37 [41] |
| 70–80 | 65.62 [60] | −34.59 [61] | 47.96 [27] | −26.98 [39] |
| 80–90 | 61.16 [58] | −33.29 [60] | 44.64 [28] | −27.99 [42] |
| 90–100 | 56.74 [60] | −24.94 [64] | 40.40 [26] | −23.50 [39] |
| 100–110 | 37.20 [66] | −9.17 [92] | 27.40 [30] | −10.02 [40] |
| 110–120 | 13.94 [27] | −0.68 [71] | 13.68 [50] | −6.04 [58] |

The unmasked curves add the edge bins and alter transition-bin medians, while the interior high-speed bins remain close to the masked estimates. This is the requested analytic counter-test, not an independent validation.

## Implant load and EMG phase controls

For JW, the four transducer channels are the compressive tibial quadrants. The primary scalar implant load is `PM + AM + AL + PL`; `GON` is retained as a sensitivity channel but is not added to that total. For each full qualified segment, phase is the frequency-weighted circular phase of the comparison signal relative to torque over 0.2–2.0 Hz after linear detrend. The equivalent lag is `-phase/(2*pi*f_weighted)`; it is a representation of the phase offset, not a causal transmission delay.

| speed | direction | implant phase (°) | equivalent lag (ms) | full segments | phase concentration |
|---:|---|---:|---:|---:|---:|
| 60 | flexion | 1.49 | 9.40 | 3 | 0.993 |
| 60 | extension | −111.22 | 412.88 | 3 | 0.900 |
| 90 | flexion | −25.38 | 66.32 | 2 | 0.999 |
| 90 | extension | −74.57 | 179.88 | 2 | 0.960 |

The direction-specific estimates and all quadrant controls are in `results.json`. The large extension phase spread and the differing quadrant phases mean these delay values are exploratory; they should not be read as a validated physiological lag.

For the source-filtered EMG channel, the phase estimates are near inversion rather than a small activation lag:

| speed | direction | EMG phase relative to torque (°) | equivalent lag (ms) | full segments | phase concentration |
|---:|---|---:|---:|---:|---:|
| 60 | flexion | 148.48 | −412.79 | 3 | 0.951 |
| 60 | extension | −179.38 | −641.78 | 3 | 0.999 |
| 90 | flexion | 169.58 | −442.25 | 2 | 1.000 |
| 90 | extension | −169.10 | 447.28 | 2 | 0.999 |

The dedicated raw `vaslat` channel has Pearson correlation 0.004 with torque at 60°/s and −0.056 at 90°/s over the common interval, so it is not used to replace the source-filtered channel. The near-180° filtered phase is reported as a phase inversion/offset only.

## Reproducibility and validation

- `PREREG.md` fixes the mask, bins, direction handling, implant-load definition, phase method, and control rules before the final run.
- `extract.py` is the executable extraction and uses only the Python standard library.
- `measurements.csv` is the combined numerical table; `results.json` contains the full curves, per-segment estimates, controls, synchronization metadata, and status fields.
- `python3 inputs/load_test.py` verified source integrity; `python3 -m py_compile extract.py` passed syntax validation.
- No held-out values, reference targets, model fitting, or population comparison were used. The result is one new observation, not model validation.
