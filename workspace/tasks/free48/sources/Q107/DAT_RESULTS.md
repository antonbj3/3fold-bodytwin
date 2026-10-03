BT-DAT-Q107

## Outcome

Six public measured resources were found that constrain the Q107 gastric model, spanning all four
requested classes (anatomy/geometry, parameters, time series, cohorts), and **8 of 10** frozen
acceptance criteria in `PREREG.md` pass. Two fail, and they fail for the same substantive reason:
**the model's most sensitive parameter has no public measured data at all.** All numbers below are
in `results.json`; all source metadata is in `DATA_SOURCES.json`.

## Built on

`inputs/Q107_model.py` (39 parameters, **0** tagged `measured:` — every one is `assumption:` or
`control:`), `inputs/Q107_RESULTS.md` (its own ±50 % sensitivity ranking), and `inputs/Q107_PREREG.md`.
`PREREG.md` + `PREREG.sha256` were frozen before any fetch; sha256
`e5d077746579fe42fe9e48abcb6731299cc648a41f1de6075c0ead62aecdace4`.

## Sources (URL / licence / size / subjects / what it constrains)

| id | licence | size | subj | constrains |
|---|---|---|---|---|
| **MRE-114** [zenodo 13839321](https://zenodo.org/records/13839321) | CC BY 4.0 | 441.5 MB zip, 228 members | 114 | `fundus_rest_volume_m3`, `antrum_rest_volume_m3`, radii, `transfer_area_m2` |
| **TOTALSEG-V2** [zenodo 8367088](https://zenodo.org/records/8367088) | CC BY 4.0 | 23 587 MB zip | 1228 | same volumes/radii, fasting-state; `stomach` is a TS class |
| **VHP** [NLM Visible Human](https://www.nlm.nih.gov/research/visible/visible_human.html) | US public domain | per-section sets | **2** | wall thickness, wall modulus prior, radii |
| **GASTEMPT** [CRAN](https://cran.r-project.org/web/packages/gastempt/index.html) | GPL ≥3 | source tarball | n/s | `scalar_*_k_per_min`, meal definition |
| **UWMADISON-GI** [Kaggle](https://www.kaggle.com/competitions/uw-madison-gi-tract-image-segmentation/data) | Competition Data | per-case | 85 | radii, between-day size variability |
| **NIDDK-GPR3** [study 196](https://repository.niddk.nih.gov/study/196) | research use permitted | **0 rows** | 406 | solid-meal protocol, `scalar_solid_k_per_min` |

## Download proof (A5) — PASS

One sample, **4.887 MB** of a 441 MB archive, pulled with HTTP range requests only
(`extract_zip_range.py` reads the ZIP central directory from the tail; nothing was silently
truncated). Verified by `load_sample.py`, which decodes the NIfTI header directly because
`nibabel` is absent here.

- `samples/100_data.nii.gz` — shape `(360, 384, 29)`, `float32`, range `[0.0, 2205.0]`, 4 834 833 B,
  sha256 `b30928155af8a26bf72cf5489d71e7636416b3060059b824cbca5e3107d63e3e`
- `samples/100_label.nii.gz` — same shape, labels `0–10`, 51 800 B,
  sha256 `0afb3d922b66efd78f71617e957f0363e0ebc23a42440d092fd397749242084b`
- **Units:** `xyzt_units=2` → **mm** in both files; voxel `1.1719 × 1.1719 × 4.800 mm`
  → 6.5919 mm³/voxel
- **Frame:** `sform_code=1`, patient-scanner anatomical, affine
  `diag(−1.1719, −1.1719, +4.800)` with zero origin — the negative i/j signs are the LPS
  convention, not a defect
- **Measured (source):** stomach **399.58 mL**, duodenum **108.57 mL**, whole labelled GI 3464.49 mL
  (derivation: voxel count × voxel volume). Label bounding box `i 71–297, j 7–354, k 2–28`.

## What fell

- **A3 FAIL, A4 FAIL.** Of the top-4 sensitivity items in `Q107_RESULTS.md`, **0/3** have measured
  public data. `pylorus_radius_m` (largest effect, max |ΔT50/T50| = 2.88), `pylorus_span`, and
  `mobilization_rate_per_s` (0.284), plus `particle_diameter_m`/`particle_cutoff_m`, are all
  **UNKNOWN**. Coverage is 6 of 10 parameter groups, not 8.
- **NIDDK-GPR3 exposes 0 dataset files** ("Datasets (0)"). It is counted for its *cohort class and
  fully specified 4-h egg-white protocol*; it contributes **no data rows today**.
- **VHP n = 2.** Public-domain wall thickness can falsify the frozen spherical geometry; it cannot
  constrain a population parameter.
- **PhysioNet has no GI-motility, EGG or gastric slow-wave database** (searched gastric, antral,
  electrogastro, manometry, slow wave). Recorded as a deliberate negative: `slow_wave_rate_hz` and
  `antrum_active_stress_pa` have no public raw EGG corpus.
- **Two traps caught and rejected** (`DATA_SOURCES.json/not_counted_candidates`): ICES *Stomach
  Content* Database is fish and invertebrate stomachs, a lexical false positive; SIRF Zenodo 4940072
  is phantom data — simulated, and it would otherwise have fit under the 50 MB cap.
- **UW-Madison is one dataset, not two.** Lee et al. 2024 (107 patients, 467 serial scans) declares
  "Data Availability: Kaggle UW Madison" — it is a publication on the same data and must not be
  double-counted.

## Honest caveat on the one number that looks like a win

The 399.58 mL stomach volume is **not** a resting volume. These subjects drank 1600–2000 mL of 2.5 %
mannitol and received 10 mg raceanisodamine IM 10 min pre-scan: the stomach is maximally distended
*and* pharmacologically relaxed. It bounds the compliant/distended state only. TotalSegmentator CT
(unsedated, clinical fasting) is the better source for the rest volumes — but at 23.6 GB it was not
fetched. **Fundus vs antrum is not separable in any public label scheme found**, and the real
stomach is a J-shaped viscus while Q107 freezes it as two spheres.

## Next step

Extract per-subject stomach masks from `TOTALSEG-V2` by the same range method for a real fasting
volume distribution; then treat the missing pylorus as the real research problem — it needs
prospective MRI/manometry of pyloric diameter, and no public dataset will supply it. Do not let
Q107 claim geometric specificity before beating `GASTEMPT`-style bulk-retention curves.

**Provenance:** *source* = fetched file or cited publisher page. *derivation* = arithmetic in
`load_sample.py` / `build_results.py`. *hypothesis* = the coverage judgement above, which is
interpretation and is labelled as such. No measured value was invented.
