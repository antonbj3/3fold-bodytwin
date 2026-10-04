BT-DAT-Q036

# Results — public measured dataset constraining the Q036 model

## Core answer

Six public, measured datasets can constrain the frozen Q036 model
(`inputs/Q036_model.py`, unchanged), one of which is downloaded and verified:

| # | Dataset | License | Size | Format | Meters | Limits |
|---|---|---|---|---|---|---|
| DS01 | Zenodo 10.5281/zenodo.7668420 — *Activation of skeletal muscle … dual-filament* (fig 1/3/4/6) | CC-BY-4.0 | 1 824 082 B (1 file) | xlsx, 4 sheet | pCa sweep against T/T0 and **kPa force**, UV flash force/L time series | `C`, `F`, `k_function_stroke`, `tau_function` |
| DS02 | Dryad 10.5061/dryad.b8d0g = Zenodo 4967455 — *Anthr1 … mouse model of hindlimb ischemia* | CC0 | 104 301 B | xlsx, 333 r × 34 k | day −1…+28, ligated versus contralateral leg, + human block (Non-ischemic/Bypass/Amputation) | `D`, `q` (time axis + contrast) |
| DS03 | GEO **GSE47730** — miRNA profile, mouse thigh, 2 h ischemia → reperf 0 h/4 h/1 d/7 d | NCBI public domain | RAW.tar 3,4 MB + filelist 1,6 kB | Affy miRNA microwall (GPL17252) CEL archive | damage time series in target tissue | `D`, `Gly` |
| DS04 | GEO **GSE246962** — unilateral hindlimb I/R 4 h, RNA-seq, treatment arms | NCBI public domain | FPKM table 2,2 MB + count table 1,2 MB (gz) | gzipped text, FPKM + integer counts (GPL24247) | damage + intervention contrast (placebo type of model) | `D` |
| DS05 | Dryad 10.5061/dryad.sbcc2fr2q = Zenodo 4985221 — SR function versus cytosolic Ca²⁺ overload | CC0 | 77 705 223 B (1 file, **> 50 MB**, therefore cannot be the sample) | zip raw data | cytosolic + SR-Ca²⁺ in soleus/EDL/gastrocnemius | `C`, `k_ca_influx/clear/mito` |
| DS06 | Dryad 10.5061/dryad.pp56ck7 = Zenodo 4936109 — q-space fibre properties, TA versus soleus | CC0 | 154 635 B | xlsx | fibril-CSA/fibre-type distribution | geometry: `tau_flow`, `k_edema_clear` (compartment structural scale) |

`URL`, `licens`, `storlek`, `format`, `subjects` and `constrains` by post
is machine-readable in `DATA_SOURCES.json`; same structures plus
criteria output values are in `results.json`. No figure from any dataset has
entered into the model — the connection is symbolic (dataset → model symbol).

## Verified sample

`curl` against `https://zenodo.org/api/records/7668420/files/Data%20source%20file_dualfilament_regulation.xlsx/content`
→ **HTTP 200, 1 824 082 B** (= 3,6 % of the 50 MB boundary), saved as
`samples/dualfilament_7668420.xlsx`,
sha256 `b16ec2f0883daca16c33d9d8667ce5adf054803e9a197e2b43dfa3366b667f3b`.

- **Form:** magic bytes `50 4b 03 04` = ZIP/OOXML package, 17 zip entries,
  4 worksheets — `Figure 1` (A2:AI43, 37 nonempty rows), `Figure 3`
  (A1:BL4563, 4563), `Figure 4` (B1:V18210, 4840), `Figure 6` (B2:H20902, 20901).
  Read with `openpyxl` 3.1.5; first row of data in each sheet has parsed numeric
  cells (see `results.json → sample_verification.sheets`).
- **Units:** read from file's own column headers, not assumed — 42 unit bearing
  headings, including `pCa` (calcium activity, −log10), `Force (kPa)`, `F(kPa)`,
  `T/T0`, `L/L0`, `SE_Force (kPa)`, `SE_DL %L0`, `time from UV(ms)`.
- **Coordinate frame:** `null`, reported as `null` with reason according to PREREG criterion 3.
  Header search found **0** hits on x/y/z, voxel, slice, ROI, origin or
  FOV; the measurement is an isolated muscle fiber against a UV flash time origin, so
  time series no spatial coordinates. No frame has been found.

## What fell (and why)

- **Best human ATP+feature anchor is not data published.** Zenodo 8122234
  (burosumab; 31P-MRS-ATP synthesis *and* muscle function test, human) consists of
  **four PDFs**. The whole point about measurable ATP synthesis + function in the same human
  found only in the article. The same error occurs in 21930313, 5930674, 14474015,
  20101221 and in Dataverse `doi:10.7910/DVN/PSIOBE` (single file = PDF report).
  18447885 / 18448175 / 18447847 has `files = 0`.
- **19472107 "31P-MRS reproducibility"** consists of 8 MATLAB script, 33 015 B,
  empty description — code without measurements.
- **MetaboLights / Metabolomics Workbench / figshare / OSF / Mendeley: zero hits**
  on skeletal muscle ischemia metabolomics. Workbench responded HTTP 200 with empty array
  for all 4 questions. MetaboLight's hits were kidney/heart. The consequence: **`A`, `L`,
  `H`, `Gly` has no metabolomic time series in this search.**
- **The wrong compartment checks were sent through the same table and rejected.**
  MTBLS15452 is "no reflow injury after myocardial I/R" — correct mechanism
  (`k_no_reflow`), but heart and 3,61 GiB. Zenodo 5001822 is brain, 6,4 GB.
  Dataverse `VKZNXV` is true skeletal muscle-GC–MS (40 netCDF, 5,06 GB) but without
  ischemia/low-flow/reperfusion/geometry linkage in the retrieved metadata →
  rejected. The filter therefore does not only capture "muscle".
- **GSE144270 was rejected for confounder**: cardiotoxin is given at the same time as I/R, so
  the transcript signal cannot be attributed to ischemia alone.
- **DS02 cannot be used as a sample**: the deposited xlsx file has zero
  unit rows and zero headers, so physical quantity behind each figure block is
  `UNKNOWN` out of the landfill. Admitted for time axis/cohort/contrast, not as flow anchor.

## Hatches (saying this bluntly)

- `M` (mPTP), `k_mptp_open`, `k_mptp_close`: **no public measurement data found.**
  GEO's 25 mPTP hits were brain/neuro, heart, and mitochondrial dynamics. mPTP
  in vivo in human skeletal muscle is not measurable by open-set assay; the anchors
  E2–E4 in `inputs/Q036_PREREG.md` is item level.
- `C` **during** ischemia/reperfusion: only non-ischemic Ca²⁺ datasets (torp,
  activation).
- `R` during skeletal muscle ischemia: none; transcriptomic ROS signatures are counted
  intentionally not (PREREG error condition).
- **`q` as quantitative time series: the biggest gap.** None public
  laser-Doppler/laser-speckle/CEUS/MRI perfusion course for skeletal muscle I/R
  was found, and DS02 is missing devices. The primary state in the model thus has
  still no measurable series — just what `inputs/Q036_RESULTS.md` points out as
  next resolution step.
- `E` (edema), `F` in acute human cohort, as well as `K01/K03/K09/K12` semantics: `UNKNOWN`.

## Reproduction

- `PREREG.md` was written **before** first request; `PREREG.sha256` =
  `9a1ae3a61668ac434c21e48f4171ca5da30139671696cc319e0f30a3cc6a1c7b`, unchanged i
  `results.json → prereg` (`sha256_recorded_before_any_query` == `sha256_now`).
- `python3 scratch/build_outputs.py` genererar `results.json` +
  `DATA_SOURCES.json` out of the actually retrieved answers and out
  `samples/dualfilament_7668420.xlsx`; no sample values are handwritten.
- Raw API response is in `scratch/raw/*.json`; the search and detail scripts i
  `scratch/query_geo.py`, `query_repos.py`, `query2.py`, `detail1-3.py`.
- **Continued from the aborted session:** just the network range from
  `agent.log` was reused. Its `scratch/bio.json` (GitHub-Pages 404-HTML),
  `zen.json` (MDPI-PDF Hit List), `bs1.json` (114 088 BioStyles Hits), `mt.json`,
  `mw.json`, `zen2.json` are garbage and **not** used as evidence; all metadata here
  are retrieved from API responses in this session.
- `inputs/` is untouched; no model parameter changed; no numeric above 60 s, none
  thread, no cloud run.

## Next step

1. The probing question for Q036 is not more datasets but a **flow progression with
   units**: either a LAS-X speckle/laser Doppler or CEUS time series for
   open repository skeletal muscle I/R, or the contact authors of DS02/DG
   b8d0g to get the `filelist`/unit lines. Without it anchors are missing for `q` and
   thus for `k_no_reflow`.
2. A `n` field must be read before using anything: DS01–DS06 are all missing
   group sizes in the retrieved records (written as `UNKNOWN`, not guessed).
3. The human cohort gap (VQI/TOPVAS/UK Biobank) is not public and should not
   suggested as calibration path; the next viable step is therefore DS01-formen off
   function coupling (kPa force vs. pCa) plus DS03-formen of damage time series.
