BT-DAT-Q156

## Results
Two PUBLIC measured datasets verified with URL, license, size, format, subjects and parameter. C1 claimed three; there were two. It is reported, not filled out.

**S1 — Harrison et al. (UCL), Dryad `10.5061/dryad.121hs31`** (`zenodo.org/records/4973308`), **CC0**, 7 411 757 B, 7z → 28 NIfTI-1 float64 (`magic n+1`, h348), rat: MDDW n=10, DTI n=6, ECG-gated n=5, dobutamine n=6. All 28 headers parsed, all values ​​finite. **Not time-resolved**: part 3/4 has exactly 2 volumes = 2 states, not a time series. Limits PVS- liquid signal; `eps_pvs`/`D_pvs` requires the b values ​​from the article.

**S2 — Mortensen, Mendeley `10.17632/78rvmj2wh8.1`**, **CC BY 4.0**, 5 227 422 599 B / 76 files. Six files ≤ 50 MB downloaded, **all sha256 matches Mendeley-API**. Limits pial PVS cross-sectional area [mm²] (`calculate_pvs_csa.m`: volume/length, pixdim in mm) and CSF→brain-tracer time [mmol/L = mol/m³ exact model unit] (`dce_mri_calc_concentration.py`, 3 baseline + 42 frames à 4 min 5 s, Gd-DTPA 0,5 kDa). The series are 201–216 MB per file → **was not downloaded** (50 MB cap).

## What failed
- C1: 2 of 3. C1b (time-resolved) **not met** — no time-resolved payload verified.
- `eps_isf=0,20`, `D_isf=1e-10` (3 kD), `L=250 µm` and `k_h=1,9e-15`: **olimited** by something downloaded.
- **600 The s criterion is untestable** with public data: S2 has 245 s/frame and no step response.
- OpenNeuro unavailable (DNS-wrong); no public deposition found for Fultz 2019.

## Unit/frame findings (most important)
1. **S1 declares no units at all** (`xyzt_units=0`, pixdim all 1,0) → voxel size and b-values are completely missing.
2. **S2 declares *meter* while pixdim is *millimeter*** (0,08/0,08/0,11 and 0,09115/0,09115/0,156) → 1000× in length, 1e6× in volume if header is taken on Tro. Silent error.
3. `t2.nii.gz` has invalid `qform_code=-31288`; read via `sform` (RAS, x reverse, Bruker/Paravision).

## Avskrivna artefakter (C5)
Zenodo 16996736 (Gan, code, CC-BY-4,0, 356 MB), 13739055 (Sun, 9 PNG-fluorescenspaneler, EEG-mat 60–225 MB), 3241364 (finite-element-** simulation**, confirmed i record text), 12518349 (22,9 GB, not inspected → UNKNOWN).

## Building on
`inputs/Q156_QUESTION.md`, `Q156_PREREG.md` (BT-HX-Q156), `Q156_model.py`, `Q156_RESULTS.md`. PREREG.md + sha256 `515b9d8f…` frozen before any inspection; `DATA_SOURCES.json` carries all five fields per source.

## Next step
Either (a) get size exceptions and retrieve **one** `dce.nii.gz` (201 MB) to fit the arrival curve and normalize CSA → real `eps_pvs`, or (b) rewrite the Q156 criterion from 600 s to hourly exchange/clearance that public data carries. A 600 s test requires a sub-minute animal two-photon/fiber optic tracer array with declared geometry — no localized public repository provides it.
