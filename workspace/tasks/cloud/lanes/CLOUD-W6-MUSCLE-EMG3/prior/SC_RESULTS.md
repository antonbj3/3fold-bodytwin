BT-PG-199

## Resultat

All four requested quantities were generated in the same run of `extract.py`, with the same hash-validated sources and separate results for 60 and 90°/s. This is a descriptive observation from SC, not a population estimate.

### Sammanfattande tabell

| Measured quantity | 60°/s | 90°/s | Unit and method | Non-circular control |
|---|---:|---:|---|---|
| Quadriceps–hamstrings co-activation | 0.187307 | 0.196585 | Unitless Q/H; absolute amount, 4th-order 6 Hz zero-phase low-pass, pooled RMS in common 10° bins 10–110° | All 10 bins and complete cycles are reported; 0,2 s countertest are reported below |
| Medial–lateral vastus ratio | 0.048468 | 0.037358 | Unitless `vasmed`/`vaslat` RMS-kvot; same envelope and bin | Bin and cycle values ​​are reported; no left/right repeat exists |
| Activation drift early→late | 1.144592 | 1.178359 | Dimensionless `RMS_sen/RMS_tidig`; last/first complete cycle over five EMG channels | Separate Q, H and channel ratios and all cycle RMS values are reported |
| EMG–medial contact time course | UNKNOWN | UNKNOWN | No medial contact series or defined onset marker is present in the frozen files | `GON/FP` and force columns were not used as surrogates; no circular imputed time |

### Primary RMS-integration

| Hastighet | Antal EMG-samples i 10–110° | Quadriceps RMS | Hamstrings RMS | Q/H | Medial RMS | Lateral RMS | Medial/lateral |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 60°/s | 9 701 | 0.073524 | 0.392532 | 0.187307 | 0.005347 | 0.110327 | 0.048468 |
| 90°/s | 10 426 | 0.069662 | 0.354359 | 0.196585 | 0.003997 | 0.107002 | 0.037358 |

RMS values are given in the raw unit scale of the source; the source does not contain a calibration unit. The ratios are therefore dimensionless within this observation.

## Metod

- Primary sources: `inputs/emg/SC_isokin{60,90}_emg.csv` and `inputs/biodex/SC_isokin{60,90}_biodex.csv`. Biodex files prefiltered EMG-columns were not used.
- Timestamp interval: actual overlap only, 0–10.4167 s (10 417 EMG-samples) at 60°/s and 0–11.6750 s (11 676) at 90°/s samples. No extrapolation.
- Synchronization: Biodex `knee flexion (deg)` linearly interpolated to valid EMG timestamps.
- Filter: `abs(original EMG)` followed by fourth-order Butterworth low-pass at 6 Hz, zero-phase `sosfiltfilt`, 1 kHz.
- Groups: quadriceps = `vasmed`, `vaslat`, `rf`; hamstrings = `semimem`, `bifem`. For each group, the mean channel effect was used and then RMS over samples.
- Primary angle integration: `[10,20), [20,30), ..., [100,110)`. The bin ratio and the pooled total are reported. `0–10°` and `110–120°` are available as sensitivity data in `results.json`.
- Cycles: turnaround detection on the Biodex angle with five-point Savitzky–Golay smoothing, 40° prominence and minimum 0,5 s distance. A complete cycle is low turn → high turn → low turn with at least 80° span. This gives 2 complete cycles at 60°/s and 3 at 90°/s.

## Numeriska 10°-bin

| Hastighet | Bin | Samples | Q/H | Medial/lateral |
|---:|---:|---:|---:|---:|
| 60°/s | 10–20° | 1 060 | 0.355837 | 0.036650 |
| 60°/s | 20–30° | 1 012 | 0.395037 | 0.040674 |
| 60°/s | 30–40° | 1 002 | 0.177754 | 0.061200 |
| 60°/s | 40–50° | 1 003 | 0.172130 | 0.061689 |
| 60°/s | 50–60° | 1 005 | 0.213438 | 0.038788 |
| 60°/s | 60–70° | 912 | 0.152500 | 0.047617 |
| 60°/s | 70–80° | 844 | 0.094282 | 0.063824 |
| 60°/s | 80–90° | 832 | 0.124838 | 0.043425 |
| 60°/s | 90–100° | 839 | 0.095176 | 0.071161 |
| 60°/s | 100–110° | 1 192 | 0.036826 | 0.183727 |
| 90°/s | 10–20° | 1 645 | 0.418644 | 0.021796 |
| 90°/s | 20–30° | 1 026 | 0.387729 | 0.024150 |
| 90°/s | 30–40° | 966 | 0.302955 | 0.034127 |
| 90°/s | 40–50° | 934 | 0.219862 | 0.040437 |
| 90°/s | 50–60° | 920 | 0.159347 | 0.041533 |
| 90°/s | 60–70° | 911 | 0.144632 | 0.048698 |
| 90°/s | 70–80° | 909 | 0.112846 | 0.057011 |
| 90°/s | 80–90° | 950 | 0.101867 | 0.060108 |
| 90°/s | 90–100° | 1 023 | 0.061414 | 0.113585 |
| 90°/s | 100–110° | 1 142 | 0.059809 | 0.136102 |

## Cycles and activation operation

| Hastighet | Cykel | Tid (s) | Samples | Q/H | Medial/lateral | RMS fem kanaler |
|---:|---:|---:|---:|---:|---:|---:|
| 60°/s | 1 | 1.0750–4.6917 | 3 617 | 0.226974 | 0.047642 | 0.201834 |
| 60°/s | 2 | 4.6917–8.3667 | 3 675 | 0.264761 | 0.031428 | 0.231018 |
| 90°/s | 1 | 2.0500–4.8250 | 2 775 | 0.212820 | 0.037625 | 0.205493 |
| 90°/s | 2 | 4.8250–7.4500 | 2 625 | 0.188885 | 0.041056 | 0.228400 |
| 90°/s | 3 | 7.4500–10.3583 | 2 909 | 0.167182 | 0.037401 | 0.242145 |

The drift value in the main table is the last/first complete cycle RMS over all five channels. Control values:

| Hastighet | Q sen/tidig | H sen/tidig | `bifem` | `rf` | `semimem` | `vaslat` | `vasmed` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 60°/s | 1.318201 | 1.130066 | 1.130192 | 1.258574 | 0.975614 | 1.338180 | 0.882754 |
| 90°/s | 0.937151 | 1.192979 | 1.193139 | 0.824142 | 1.020275 | 0.975903 | 0.970085 |

## 0,2 s-motprov

The hamstrings envelope was replaced with `e_h(t−0.2 s)`; quadriceps, angle bins and the shared time window were otherwise unchanged. The shifted run is compared with an allocated Q/H run on exactly the same samples.

| Rate | Samples in countertest window | Allocated Q/H | Hamstrings +0,2 s Q/H | Absolute change | Relative change |
|---:|---:|---:|---:|---:|---:|
| 60°/s | 9 501 | 0.184744 | 0.186356 | 0.001613 | +0.87 % |
| 90°/s | 10 226 | 0.202663 | 0.215392 | 0.012729 | +6.28 % |

## Channel noise and data quality

`P05 envelope` is the 5th percentile of the rectified 6 Hz envelope and is used as a conservative activation noise floor proxy. `>6 Hz` is the RMS fraction of `raw − LP6(raw)`; it is a band/quality measure, not a calibrated electrode noise estimate. All channels had 0 % missing or non-finite samples.

| Rate | Channel | Raw RMS | Raw SD | P05 envelope | Share >6 Hz |
|---:|---|---:|---:|---:|---:|
| 60°/s | `semimem` | 0.018927 | 0.018928 | 0.003371 | 1.000000 |
| 60°/s | `bifem` | 0.679595 | 0.679628 | 0.062270 | 1.000000 |
| 60°/s | `vasmed` | 0.008083 | 0.008083 | 0.000797 | 0.999992 |
| 60°/s | `vaslat` | 0.144901 | 0.144908 | 0.003891 | 1.000004 |
| 60°/s | `rf` | 0.080031 | 0.080035 | 0.004098 | 0.999884 |
| 90°/s | `semimem` | 0.019251 | 0.019252 | 0.003473 | 1.000004 |
| 90°/s | `bifem` | 0.624870 | 0.624896 | 0.069552 | 1.000278 |
| 90°/s | `vasmed` | 0.005311 | 0.005311 | 0.001008 | 1.000231 |
| 90°/s | `vaslat` | 0.127647 | 0.127651 | 0.003827 | 1.000002 |
| 90°/s | `rf` | 0.066911 | 0.066914 | 0.004775 | 1.000000 |

The large difference in the raw scale of the channels means that the ratios should not be interpreted as calibrated physiological force measures.

## Reference and file status

- `python3 inputs/load_test.py`: PASS for six frozen files' hash, scheme and numeric basic requirements.
- OIndependent reference sample: UNKNOWN. There are no independent held-out labels or a reference target in the package; therefore nothing is reported PASS.
- Executable extraction code: `extract.py`; machine results: `results.json`; pre-registration: `PREREG.md`.
