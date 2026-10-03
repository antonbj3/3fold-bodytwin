BT-DAT-Q115

## What was built on
`BRIEF.md`, `inputs/NIGHT_PREAMBLE.md`, `inputs/Q115_QUESTION.md` and the model in `inputs/Q115_model.py` (readable, reparameterized: **some parameter was not changed**). `PREREG.md` was frozen with `PREREG.sha256` before first metadata call and download; `sha256sum -c` gives OK, and `inputs_before.sha256` shows that `inputs/` is byte-identical after the run (CT5).

## Published measured datasets that constrain the model
Seven sources were verified in-session via Zenodo-API (license, size, file count, date from the record, not from the paper): `DATA_SOURCES.json`.

- **zenodo:22834925** (CC BY, 22,1 GB) per-cell attribute for 32 051 epithelial cells in 7 ROI (3 pcs 12 h, 4 pcs 72 h) → limits `n_cells`, `mature_cell_area_mm2`.
- **zenodo:20200485** (CC BY, 0,12 MB) TEER, 42 samples/126 values, 6 arms, 24 h → limits `g_mature`/`g_young`.
- **zenodo:3403670** (CC BY-NC, 514 MB) villus axis zoning, 3 mice → limits `maturation_tau_h`, `transport_birth_fraction`.
- **zenodo:6768583** (CC BY, 194 MB) 840 organoid images, 23 066 annotations, classes "1–2 crypts"/"≥3" → `crypt_count`.
- **zenodo:4311473** (CC BY, 405 MB) Ilastik organoid segmentation, 4 replicates × 5 time points → `crypt_depth_mm`.
- **zenodo:7544194** (CC BY, 6,68 GB) 3D organoid microscopy → villus/crypt geometry; **over 50 MB, not downloaded**.
- **zenodo:19541885** (CC BY, 7,9 MB) villus length/width, crypt depth, **guinea pig** → `villus_height_mm` (system mismatch flagged).

Rejected: `20185659` (license = null), `895554`/`826475` (only PDF), celltrackingchallenge.net (read: 3D) series are carcinoma/embryo, no intestine).

## Sample (≤ 50 MB, total 8 049 858 B)
`samples/22834925_sprinkled.zip` (5 503 250 B) — ZIP, 38 members, 7 × `cell_attributes.csv` with the columns `label, area, center_x, center_y`; 6 979 cells in ROI 12 h/1; median cell area 26 664; **units: `not_stated_in_file`**, coordinate frame 2D image coordinates without pixel size → absolute mm² **not done**. `samples/20200485_Processed Data Figure-1.xlsx` (12 491 B) — TEER in % of control: control 100,08 %, TNF-α+IFN-γ 84,26 % (**−15,81 %**), +BSW 108,80 %, +BWW 102,02 %; no absolute ohm·cm². Placebo CT1: `samples/895554_article.pdf` = `%PDF-1.4`, no numerical table.

## What fell
No single open source provides µm-calibrated crypt–tip length, migration rate in mm/h, absolute conductance in mS/cm² or cell age distribution in cells. The model reference (1 885 cells, 26 cells/h, 65 h) remains literature source. Upstream FAIL (mechanistic ≥10 %) and UNKNOWN (empirical superiority) remain unchanged.

## Next step
Retrieve `zenodo:7544194` (or a µm-calibrated 3D sample) for villus geometry; search open pulse-chase/EdU line tracking for the age distribution; get absolute TEER in ohm·cm² from the `20200485` authors or the Ussing campaign. Command: `python3 fetch_and_check.py && python3 prov_stats.py`.
