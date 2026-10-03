BT-DAT-Q058

# BT-DAT-Q058 — public measured datasets that constrain BT-HX-Q058

Status: **completed**. Continuation of an interrupted session; all downloading, preregistration and
calculation had already taken place. This passage is the verification of the existing work plus the
final report. No model parameters have been changed (PREREG N3).

## 1. What was built on

| Underlag | Roll |
|---|---|
| `inputs/Q058_QUESTION.md` | Q058 — switching between detail levels without losing biological memory |
| `inputs/Q058_model.py` | the model being constrained, `PARAMETER_TABLE` rows 51–157 |
| `inputs/Q058_PREREG.md` | the model task's own frozen criteria and `UNKNOWN` rule |
| `inputs/NIGHT_PREAMBLE.md` | work order, writable only here |
| `PREREG.md` + `PREREG.sha256` | frozen **before** download; `sha256sum -c` → `PREREG.md: OK` |

Model verified by import, not visual reading: **15** parameter rows, of which **9**
`assumption` rows (`tau_a, tau_x, gamma_h, tau_f, k_f, tau_s, k_s, w_f, w_s`) and 3
with literature sources (`tau_h=22.0 s`, `c0=0.050 s`, `beta_ct=0.379`; the remaining 3 are
normalisation and two frozen perturbation values). **Zero geometry/anatomy parameters** — no
length, mass, moment arm or anatomical attachment exists. Every geometry/anatomy dataset is therefore
`NOT_APPLICABLE` for this model, and "coordinate frame" is reported as `NOT_APPLICABLE` instead
of fabricated. It is a property of the input, not a void.

## 2. Nyckeltal med fil

All numbers below are traceable to `results.json`.

**Sample** — `samples/zenodo_5189275_emg_fatigue.zip`, 10 603 743 B (≤ 50 000 000 B, C1 ✓),
md5 `a954b4ef…b864` matches the downloaded file, CC-BY-4.0, 15 healthy subjects.
Loads: `24085 × 8` per person, 200 Hz, 8-bit ADC counts (`−128..+127`), 8 channels,
~120,2 s per person. **Unit:** 8-bit ADC counts, the gain is not published →
microvolts `UNKNOWN`. **Coordinate frame:** `NOT_APPLICABLE` (EMG channels carry no anatomical frame;
channel–muscle mapping is missing in the record).

| Measured quantity | Result | Class |
|---|---|---|
| **M1** on-transient | No load onset exists in the sample — 6 kg is already applied at the first sample. Plateau/first-0,5 s = 0,938–1,127 in all 15, all below threshold 1,15. Floor according to C3 is 2/200 Hz = 0,010 s | `NOT_APPLICABLE` 15/15 → `tau_a` **`NOT_IDENTIFIED`** |
| **M2** fatigue | Median-channel ratio last 20 s / first 20 s = **1,399 ± 0,286** (median ± SD across 15 people), range 1,076–2,326. 90,0 % of 120 channel comparisons > 1 (p05 0,918 / p95 2,879). Aggregated log-slope **+0,00352 ± 0,00208 /s** | `BOUND_ONLY` |
| **M3** spread | 15/15 people have the majority of channels rising (5–8 of 8) | `CONSTRAINED` |

M2 is a **rise**, not a decline — that is the finding itself, not a pipeline error. R²
for log-linear fitting is 0,05–0,88, so a single-exponential `tau` is not acceptable
for most people and is therefore not reported (`aggregate_tau_s = null` for all 15).

**N1 sham window: NOT EXECUTED.** PREREG required a window with `drive=0`; the record contains a single
120 s segment with 6 kg held throughout, so an unloaded window does not exist. A replacement
control was performed instead: channel-versus-channel contrast within the same segment. Falling channels exist
(ratio down to 0,491 in sub1), which upper-bounds recording drift. The M2 metric is thus not
invalid, but its sign is not robust to drift.

**N2 unit check: PARTIAL.** The frozen µV band [0,5, 5000] could **not** be applied, because
the published unit is ADC counts and the gain is missing. What could be checked:
1,18 % of the samples are clipped at full scale (1135 at ≥ +127, 1125 at ≤ −128). The sample is therefore used
for shape and unit checks, **not** for absolute quantification.

**Catalogue** — `DATA_SOURCES.json`, 7 candidates, 12 constraints:
`CONSTRAINED` 2, `BOUND_ONLY` 5, `NOT_IDENTIFIED` 2, `NOT_APPLICABLE` 3.

| # | Dataset | Licence | Size | People | Constrains |
|---|---|---|---|---|---|
| 1 | [10.5281/zenodo.5189275](https://doi.org/10.5281/zenodo.5189275) sEMG during 6 kg isometry | CC-BY-4.0 | 10,6 MB | 15 | `tau_f,tau_s` (BOUND_ONLY); between-person spread (CONSTRAINED); `tau_a,tau_x` (NOT_IDENTIFIED) |
| 2 | [10.5281/zenodo.232032](https://doi.org/10.5281/zenodo.232032) EMG vs ultrasound onset | CC-BY-4.0 | 15,6 GB | 10 | `tau_a` (CONSTRAINED); `gamma_h` (BOUND_ONLY) |
| 3 | [10.5281/zenodo.1213604](https://doi.org/10.5281/zenodo.1213604) M-wave/MEP, hypnosis | CC-BY-4.0 | 878 kB | 13 | `gamma_h` (BOUND_ONLY); `c0,beta_ct` (NOT_APPLICABLE, C5) |
| 4 | [10.5061/dryad.jq8bg5c](https://doi.org/10.5061/dryad.jq8bg5c) MEF/MDF at 80 % MVC | CC0-1.0 | 1,05 MB | 13 | `tau_f,tau_s` (BOUND_ONLY); model outputs (NOT_IDENTIFIED) |
| 5 | [10.5281/zenodo.4641292](https://doi.org/10.5281/zenodo.4641292) MEP during perturbation | CC-BY-4.0 | 78,4 MB | not stated | `NOT_APPLICABLE` — not downloadable under 50 MB |
| 6 | [10.5281/zenodo.5534422](https://doi.org/10.5281/zenodo.5534422) MMG per quadriceps head | CC-BY-4.0 | 101 kB | 25 | `gamma_h,w_f,w_s` (BOUND_ONLY); the table → C6 |
| 7 | [10.5061/dryad.326qs26](https://doi.org/10.5061/dryad.326qs26) motor-unit pool *simulation* | CC0-1.0 | 203 MB | 0 | `NOT_APPLICABLE` — model outputs, not measured data |

## 3. What failed

- **No force-domain measurement was found.** C4 prohibits EMG from setting `gamma_h, k_f, k_s, w_f, w_s`,
  because log amplitude is activation and the model's `load = clip(x·Phi(…), 0, 1)` is force-dimensioned.
  No fatigue time constant could therefore receive a point value — only a side constraint.
- **`tau_a` and `tau_x` could not be identified**, in line with C3: the floor 0,010 s is next to
  `tau_a = 0,025 s` and the sample also contains no load onset at all.
- **N1 could not be run** as preregistered (see above). This is a deviation from PREREG, not
  an optional substitution.
- **No detail switches are documented publicly.** What Q058 actually asks — is there a
  documented switch between detail levels with a subsequent perturbation — is entirely missing in the
  public material found.
- **N2 could not be applied as frozen** (ADC counts, not µV).

## 4. Korrigering i denna session

`DATA_SOURCES.json` rank 4 stated "8 sheets" for the Dryad file. Verification against
`pandas.read_excel` gave **17 sheets**: 13 person sheets ("1".."13", 55×23 except "4", which is 54×23) plus
4 analysis sheets (63×23, 60×20, 62×23, 60×20). The field is corrected. `results.json` already had 17 and was
correct. The 13 person sheets independently confirm the cohort size 13.

## 5.Verification performed here

- `analyze_sample.py` was requeued from scratch → `results.json` **byte-identical**. No numbers are
  handwritten.
- All 6 records were checked live against the Zenodo API: licence and total size match for all.
- `md5` and `sha256` for all files in `samples/` match; no sample exceeds 50 MB; all five
  files `loads: true` with a declared reader.
- Model parameters verified by importing `DEFAULT_PARAMETERS`, not by reading.

## 6. Next steps

1. **Acquire 10.5281/zenodo.232032** (per-subject DICOM ~3,2 MB, well below 50 MB). The only dataset
   whose measured quantity *is* activation delay. Requires a DICOM reader with pixel spacing to
   give `tau_a` an actual bound; the file here was read only to file metadata.
2. **Search specifically for the force domain**: held isometric force in F/F_max with documented switching.
   Without this, `k_f, k_s, w_f, w_s` remain unidentifiable — this is the largest gap.
3. **Evoked-twitch time series** for `c0` and `tau_x` (C5 requires a registered twitch).
4. Record 5 (78,4 MB) can be downloaded if C1 is raised; it is structurally closest to the model's perturbation test.

## 7. Bounds

No internal BodyTwin/external musculoskeletal solver data, no LOSO sweep, no cloud run, no post hoc adjustment.
No verdict words. `UNKNOWN` means a measurement is missing and has never been filled with literature values.
EOF marking: no claim about the number of subjects is made without a source in the dataset's own
documentation or archive structure.
