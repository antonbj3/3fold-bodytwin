BT-PG-179

## Results
Status: OBSERVATION for DM; no population estimate and no independent reference validation. PASS is not claimed.

### Coherent measurement package
| quantity | value | unit | method | non-circular control |
|---|---:|---|---|---|
| 60°/s torque–angle curve | peak 50.189 at 64.94° (10°-bin 60–70°) | Nm; angle ° | 10° bin median of `knee torque (Nm)` for `velocity >= 0.8 × 60°/s`; all valid samples pooled | 5 velocity groups; raw peak 62.608 Nm; curve without mask below |
| 90°/s torque–angle curve | peak 25.958 at 44.84° (10°-bin 40–50°) | Nm; angle ° | 10° bin median of `knee torque (Nm)` for `velocity >= 0.8 × 90°/s`; all valid samples pooled | 3 velocity groups; raw peak 27.710 Nm |
| torque ratio 60 versus 90°/s | median 1.559; peak ratio 1.933 | dimensionless | Paired 10°-binmedian: M60/M90 over 9 shared bins | Ratio of medians of groupwise peak torque 1.922; per-bin ratios in `results.json`; no other individual is compared |
| torque → implant load | 60°/s: -133.3; 90°/s: -295.8 (median per velocity; positive = load after torque) | ms | Exploratory centered cross-correlation, five-channel norm PM/AM/AL/PL/GON, ±500 ms, complete groups | Proxy because of local column schema; status 60°/s: unstable; 90°/s: unstable; correlation median 60°/s: 0.468; 90°/s: 0.440; group spread and component check in `results.json` |
| EMG → torque | 60°/s: 116.7; 90°/s: 58.3 (median per velocity; positive = EMG after torque) | ms | Supplied filtered `Filtered vaslat EMG` against torque, same 120 Hz time axes, ±500 ms | Correlation median 60°/s: 0.850; 90°/s: 0.840; raw 1 kHz-interpolated EMG control 60°/s: 66.7; 90°/s: -187.5 |

### 60°/s primary curve
Number of velocity groups: 5; complete groups: 5; samples: 759. Intervals: 0.750–1.908; 4.108–5.400; 7.483–8.717; 10.642–11.942; 13.942–15.242.
| 10°-bin (°) | median angle (°) | median torque (Nm) | n | repetition groups | median velocity (°/s) |
|---:|---:|---:|---:|---:|---:|
| 20–30 | 28.02 | 27.431 | 40 | 5 | 59.40 |
| 30–40 | 35.05 | 33.157 | 100 | 5 | 59.56 |
| 40–50 | 45.00 | 40.479 | 101 | 5 | 59.60 |
| 50–60 | 54.97 | 48.373 | 100 | 5 | 59.77 |
| 60–70 | 64.94 | 50.189 | 100 | 5 | 59.86 |
| 70–80 | 75.07 | 41.089 | 101 | 5 | 59.88 |
| 80–90 | 84.91 | 26.725 | 100 | 5 | 59.85 |
| 90–100 | 94.23 | 19.000 | 91 | 5 | 59.57 |
| 100–110 | 102.07 | 15.082 | 26 | 3 | 58.45 |

### 90°/s primary curve
Number of velocity groups: 3; complete groups: 2; samples: 285. Intervals: 0.000–0.475 (partial); 2.183–3.125; 4.767–5.700.
| 10°-bin (°) | median angle (°) | median torque (Nm) | n | repetition groups | median velocity (°/s) |
|---:|---:|---:|---:|---:|---:|
| 20–30 | 25.78 | 20.574 | 34 | 3 | 89.77 |
| 30–40 | 35.02 | 24.849 | 39 | 3 | 89.70 |
| 40–50 | 44.84 | 25.958 | 41 | 3 | 89.79 |
| 50–60 | 55.13 | 25.259 | 40 | 3 | 89.72 |
| 60–70 | 64.17 | 24.830 | 34 | 3 | 89.45 |
| 70–80 | 75.13 | 17.938 | 27 | 2 | 86.73 |
| 80–90 | 84.82 | 11.728 | 28 | 2 | 86.50 |
| 90–100 | 95.00 | 14.722 | 28 | 2 | 86.50 |
| 100–110 | 102.47 | 13.508 | 14 | 2 | 82.30 |

### Control without velocity mask
Same positive registered branch, but `velocity > 0` without the 80% threshold. This is a control and does not replace the primary curve.
- 60°/s: peak 50.189 Nm at 64.94°; 1073 samples.
- 90°/s: peak 25.958 Nm at 44.84°; 386 samples.

### Data, filters and limitations
- Biodex–knee overlap: 60°/s 0.000000–17.008333 s (2042 samples); 90°/s 0.000000–7.116667 s (855 samples). The time axes are identical; EMG is 1 kHz and used only as a control after interpolation to the Biodex time axis.
- Torque, angle and velocity are used as supplied; no extra smoothing is applied. `Filtered vaslat EMG` is the existing filter chain in the Biodex file.
- The knee-force header is shorter than the rows. Only PM, AM, AL, PL and GON are used; three extra numerical fields after them (index 6–8) are ignored because their name mapping is uncertain. The force unit is assumed to be N but is not explicitly stated in the local file.
- The negative velocity branch and raw EMG result are in `results.json` as sensitivity analysis. No population comparison or training/test fitting is performed.

### Reproduction
Run `python3 inputs/load_test.py` and then `python3 analyze.py` from the package root. The latter writes `results.json`, `curve_table.csv` and this report from the frozen inputs.
