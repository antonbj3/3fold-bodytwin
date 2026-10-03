BT-DAT-Q019

## What it built on

`inputs/` (Q019 = nephron segment → organ-level filtration/transport/excretion). The model in `inputs/Q019_model.py` has filtration scale, segment geometry, residence times, water permeabilities, `vmax`/`Km`, interstitium and AVP state as **assumptions**; only Layton & Layton (2019) Table 1–3 is an external source. `PREREG.md` (sha256 `afe42e80…0b499ec0`) was frozen before download.

## Nyckeltal med fil

- **Flow scale (K5).** Model: `13.2391 ml/min/kidney` (`inputs/Q019_results.json`). MIMIC-IV demo, 137 stays: median `1.24`, p99 `3.63`, p99.9 `5.36`, max `5.59 ml/min` → the model at the **100th percentile**. eICU demo, 1803 stays: median `1.17`, p99.9 `32.60`; with recording window `12.28` and p90 `2.97` → flag set in both cases (full/occupied window gives different answers, reported, threshold unchanged). Source: `samples/mimic4_demo_2.2_icu_outputevents.csv.gz`, `samples/eicu_crd_demo_2.0_intakeOutput.csv.gz`.
- **Na excretion (K6).** Urinary Na 40 mEq/L × flow, 55 paired stays: max `465.9 µmol/min`. The model `3400.5 µmol/min/kidney` = **14,6×** the highest measured value.
- **Osmolality (K6b).** Urine osmolality (itemid 51093, n=69): median `352`, max `1074 mOsm/kg`; the model `522.9 mM` lies **inside** the distribution (81st percentile). The only terminal output the open data do not contradict.
- **Geometry (K7).** TCIA C4KC-KiTS (CC BY 3.0, 300 patients) loaded: CT 62 slices, `0.9238 mm` in-plane, `3.0 mm` z, MONOCHROME2/HU, `ImagePositionPatient` given. The SEG object is **internally inconsistent** (sum over slices `547.5 mL` vs in-plane union `31.7 mL`; z-grid is not on CT-grid) → kidney volume **not usable**, reported as `criterion_failure_not_evaluable`.
- **Controls.** Null-run difference `0.0`; parser check `100466 = 100467−1` rows; unit placebo with 1000× upwards reverses the flag (downwards does not — the prereg variant is reported as non-discriminating).

## What failed

No open measured dataset was found for `lp_cm_s`, `vmax`, `Km` or `residence_min` per segment; the organ-level anchor is missing because of the SEG inconsistency. MIMIC-III demo returns 404 (registered as S4).

## Next step

Acquire kidney geometry from a source with verified per-slice segmentation (KiTS21-NIfTi), and fit `residence_min`/`vmax_na_umol_min` against the paired urinary Na/osmolality pairs in the MIMIC-IV demo instead of a single literature anchor.
