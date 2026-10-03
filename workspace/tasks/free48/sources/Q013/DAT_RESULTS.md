BT-DAT-Q013

# Public measured datasets that constrain the BT-HX-Q013 PK/PD model

Frozen preregistration: `PREREG.md`, sha256 `40228d8fcd24f861a20885fad6681aac85898d91368d3cd5f8ffa3c69c5d1674` (also in `PREREG.sha256`, written before any download or analysis).
Machine-readable result: `results.json` (sha256 `161f0fd0d36e938e8c5cf63d7454b9223ce5b1872e67df49c17849751bc007a0`). Source table: `DATA_SOURCES.json`. Sample loader: `dat_probe.py` -> `sample_report.json`.

## What was built on

`inputs/Q013_model.py` (BT-HX-Q013-0.1) is a first-principles PK/PD model: oral one-compartment PK, a driven local compartment `dC_local/dt = k_in C_p - k_out C_local`, Hill occupancy and `E = E0 + gamma Emax O`. Eleven parameters are `assumption`; the only published anchor is `R_plasma = 1.27` from Tachibana 2025 (PMID 41150706). The model itself reports `independent_local_exposure: UNKNOWN`, `cell_sensitivity_measurement: UNKNOWN` and no split of `F_mult` vs `CL_mult`. This job supplies public data for those gaps. The model was **not** re-run, no parameter was estimated (PREREG M3), and `inputs/Q013_results.json` is unchanged.

## Sources (all URLs HTTP-verified in this job, 2026-09-25)

| id | source | license | size | format | n | constrains | status |
|---|---|---|---:|---|---:|---|---|
| S1 | ClinicalTrials.gov NCT02732275 (the study the model anchors on) | US Gov, **not** machine-verified here | 17 880 B | JSON (API v2) | 100 est. enrolled; DDI sub-study 24 (15/16 PK-evaluable) | cohort and design only — `hasResults = false` | registry, no results |
| S2 | Dryad `10.5061/dryad.z612jm683` — nevirapine DDI on artemether-lumefantrine, CYP2B6-stratified | CC0-1.0 (API field) | 375 026 B, 3 files | XLSX | 60 (30 HIV+ / 30 HIV−) | B1: `R_plasma = F_mult/CL_mult` as a **signed** ratio, `ka`, `CL`, `V` | accepted |
| S3 | Zenodo `10.5281/zenodo.18363164` — lumefantrine in plasma **and** breast milk | CC-BY-4.0 (record field) | 8 506 B | CSV, 11 cols, 195 rows | 6 subjects; 78 plasma + 81 milk rows | B2: `k_in`, `k_out`, `R_local` (paired non-plasma compartment) | accepted, **downloaded** |
| S4 | PubChem AID 743139, Tox21 aromatase qHTS (concise) | public domain (NCBI policy, quoted in job) | 3 272 962 B | JSON/CSV table | 10 486 tested substances | B3: `ec50`, `emax`, `hill` at absolute uM scale | accepted |
| S5 | Zenodo `10.5281/zenodo.21757686` — MedPelvis3D pelvic anatomy | CC-BY-4.0 (record field) | 31 351 310 903 B, 100 files | NIfTI + STL + CSV + NPY | 99 CT cases (59 M / 40 F) | **nothing** — negative control | accepted as control |

Not usable, with reason: UCSF-FDA TransPortal (DNS did not resolve — would have been the P-gp/CYP3A mechanistic source for B1), NLM Visible Human (HTTP 403), BodyParts3D (HTTP 404), TDC curated ADME tables (no per-dataset license field in the fetched file; and F4(c) rejects curated tables without raw per-subject records).

## Sample: it downloads and it loads

`samples/Ojara_et_al_Lumefantrine_plasma_and_breast_milk_data.csv`, 8 506 bytes (limit 50 000 000), sha256 `6fccfd6d1dd720b743561f90db64ec13b9232fb71eb6bc58a1a2156397b2a5d3`, content-length confirmed in the HTTP header and re-checked with `stat`. Parsed with `dat_probe.py` (standard-library `csv` only): **195 data rows, 11 columns, 6 subjects, exit 0**.

- **Form:** long format, one observation per row; `CMT` 1 = 36 dosing rows (`DOSE = 480` mg, no concentration), 2 = 78 plasma rows, 3 = 81 breast-milk rows — the compartment legend was read out of the file itself, not from the abstract, and it reproduces the record's stated 78/81.
- **Units:** `PLASMA CONC (ug/L)` and `BREAST MILK (ug/L)` in µg/L; time in h; `DOSE` in mg (480); `WEIGHT` kg; `HEIGHT` cm; maternal age in years, infant age in months.
- **Coordinate frame:** no spatial frame (this is a time series). Two coexisting time bases: `TIME` (h since first dose, 0–333.5 h, 107 distinct) and `Time After Dose` (h since the most recent dose, 0–276.6 h, 57 distinct) — i.e. a 3-day multi-dose schedule, not a single-dose profile.
- **Data quality, not cleaned:** 1 of 195 rows has a blank subject id (excluded, not imputed); 5 concentration cells are exactly 0 (3 plasma, 2 milk); 3 compartment rows carry `.` in both columns; 234 cells are `.`; missing token is `.`. No censoring rule applied.
- **Descriptive only:** plasma 76 values, median 1128.7 µg/L, max 7084.6; milk 80 values, median 168.4 µg/L, max 3912.4; per-subject milk-Cmax/plasma-Cmax from 0.171 to 0.567, median 0.384. No `ka`, `CL`, `k_in`, `k_out`, `EC50`, `Emax` or `gamma` was fitted from this.

## What the data does to the model

- **B1 `F`/`CL`/`ka`/`V` — PARTIAL.** S2 has exactly the design needed to split `F_mult` from `CL_mult`: per-subject victim-drug plasma concentration-time with and without a perpetrator drug, and the perpetrator here is an inducer, so the sign is the opposite of valemetostat/digoxin. Still `UNKNOWN` for the model's own pair.
- **B2 `k_in`/`k_out` — PARTIAL.** S3 is the structural prototype the model lacks: the same subject and the same time base measured in plasma and in a non-plasma compartment. Numbers do not transfer (other drug, other tissue, fitting forbidden).
- **B3 `ec50`/`emax`/`hill` — PARTIAL.** S4 shows cell-based potency exists publicly with absolute units at n = 10 486, but on aromatase/CYP19A1, and the model's `EC50`/`Emax` are dimensionless.
- **B4 `gamma` — UNKNOWN.** No public factorial design (drug alone / PK-only / cell-only / combination with the same cell readout) was found. This is the quantity the model exists to separate, and it is the one with no public data.

## What fell, and the one finding that matters

1. **H1 falsified as stated.** There is no public measured source for B4.
2. **The anchor is not data.** NCT02732275 has `hasResults = false`; `R_plasma = 1.27` exists only as a published table value. A pipeline that needs it machine-readable has nothing to read.
3. **A falsification risk inside the same study.** Its midazolam arm went the other way — `AUClast` GMR **0,87** (90% CI 0,75–1,03, n = 15) while digoxin rose to **1,27** (1,06–1,52, n = 16), from one perpetrator. The frozen `F_mult = 1.27, CL_mult = 1` branch is therefore exposed to a within-study discordance, and because no per-subject data are public it cannot currently be tested. This is provenance, not a dataset, and is not counted as one.
4. **Geometry is irrelevant here (H4 supported).** MedPelvis3D carries full metadata — 99 cases, mm, patient-level LPS, NIfTI affines — and constrains nothing in a model with no geometric parameter. Listed so the category is answered, not because it helps.

## Not claimed

No parameter estimated, fitted or changed; model not re-run; no value interpolated, imputed or invented; the missing subject id left missing.

## Next step

Run the identifiability procedure on S2 first: signed per-subject DDI plasma data with a genotype split can separate `F_mult` from `CL_mult` and gives the test the valemetostat/digoxin pair currently lacks. S3 is the template for the missing `C_local` measurement. Neither was started here.
