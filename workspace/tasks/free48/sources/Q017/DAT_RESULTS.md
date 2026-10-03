BT-DAT-Q017

## Status

Building on the aborted session in the same directory: `PREREG.md` + `PREREG.sha256` were
already frozen and sha256-checked before this run (`sha256sum -c PREREG.sha256` → OK),
and the predecessor had downloaded `samples/eicu_crd_demo_2.0.1/` (8 files) and `samples/nhanes_2017_2018/`
(5 files). `scripts/load_sample.py` had only had time to write `patient` (line
1 in `load_sample.out` ). I completed the parsing, did not load a new file, did
not change a frozen limit value, and did not run `inputs/Q017_model.py` .

A single fix was made in `scripts/load_sample.py`: the eICU table lookup
`infusionDrug` was case sensitive, while the publisher's own archive named the file
`infusiondrug.csv.gz`. The file name on disk was not changed; only the lookup was done
case insensitive. Without this, `infusiondrug.csv.gz` (38 256 lines) would never have been opened.

## Prov: form, enheter, koordinatram

`scripts/load_sample.py` → `results_load_sample.json`. 13 files, **42,407 MB of 50 MB**
(42 407 332 bytes, summed over files in `samples/`; `du -sb` gives 42 419 620 bytes including
three directories). All 13 opened without error.

| source | file | rows × columns | format |
|---|---|---|---|
| eICU demo 2.0.1 | `vitalPeriodic.csv.gz` | 1 634 960 × 19 | gzip CSV |
| | `vitalAperiodic.csv.gz` | 274 088 × 13 | gzip CSV |
| | `intakeOutput.csv.gz` | 100 466 × 12 | gzip CSV |
| | `lab.csv.gz` | 434 660 × 10 | gzip CSV |
| | `patient.csv.gz` | 2 520 × 29 | gzip CSV |
| | `infusiondrug.csv.gz` | 38 256 × 9 | gzip CSV |
| NHANES 2017–2018 | `DEMO_J.xpt` | 9 254 × 46 | SAS XPORT |
| | `BPX_J.xpt`, `BMX_J.xpt` | 8 704 × 21 | SAS XPORT |
| | `BIOPRO_J.xpt` | 6 401 × 41 | SAS XPORT |
| | `DXX_J.xpt` | 5 114 × 93 | SAS XPORT |

**Integrity:** all six eICU data files match sha256 in the publisher's own `samples/eicu_crd_demo_2.0.1/SHA256SUMS.txt`
. NHANES does not publish checksums; local sha256 for all 13 files is
in `samples/SHA256SUMS_LOCAL.txt` (verified `OK` for each line).

**Tidsreferens (eICU):** `observationoffset` / `intakeoutputentryoffset` / `labresultoffset` =
minutes after ICU-intag. The `vitalPeriodic` lines lie at the 5 minute offset and are sparse
placed; no interpolation is done. **No spatial coordinate frame exists** — the values are
bedside measurements from Philips IntelliVue, not reconstructed anatomy. **NHANES:** no timeline
and no coordinate frame; one row = one survey visit, DXA-regionerna (DXX\*) is
body compartment, not coordinates.

**Units:** `systemicmean` /`pamean`/`noninvasivemean`/`cvp`/`paop` in mmHg, `cellvaluenumeric`
in mL, `labresult` according to `labmeasurenamesystem` (creatinine mg/dL, albumin
g/dL), NHANES `BPXSY*` /`BPXDI*` in mmHg, `BMXWT` in kg, `LBXSCR` in mg/dL, `LBXSAL`
in g/dL. The eICU unit strings `svr` /`svri`/`pvr` could **not** be verified verbatim
in this run (see below), so no conversion to model mmHg·s/mL is claimed.

## Kandidater — `DATA_SOURCES.json`

8 entries, of which **3 are counted as candidates**. All 8 have (a) a landing page that answered HTTP 200
in this run. Seven have (b) a verbatim copied license string; NHANES-posterna is missing one
license string on the source side and are provided with the field instead
`NO_LICENSE_STRING_ON_THE_SOURCE_PAGE` plus the operationally verified statement that the files
loads without login, DUA, application or payment — it is the only record where (b) is not
literally coated, and it's flagged. The five entries marked `counts_as_candidate: false`
lacks public file index or requires DUA (4 st) or is `not_applicable` (1 st); each one
carries a `missing_fields` list.

| entry | status | what it limits |
|---|---|---|
| eICU-CRD demo 2.0.1 | `downloaded_parsed` | MAP-kanal + P_a0, P_v0 (CVP/PAOP), G0 via eGFR, f_reabs0 via urin, π_gc via albumin, 3 L/3 h-protokollet, V0 via vikt |
| NHANES 2017–2018 | `downloaded_parsed` | P_a0 (vuxen-MAP), G0 (eGFR), π_gc (albumin), V0 (kroppsmassa) |
| MIMIC-IV demo 2.2 | `listing_verified_only` | basically everything; **the demo contains no measurement tables** (only LICENSE/README/SHA256SUMS/demo_subject_id.csv), so nothing was loaded |
| MIMIC-IV 3.1 | `metadata_only` | synchronous MAP + crystalloid volume + urine + creatinine per care session — the most important unused limitation |
| HiRID 1.1.1 | `metadata_only` | the only time-resolved candidate that can solve τ = 5–30 s |
| eICU-CRD full 2.0 | `metadata_only` | same channels, "over 200,000 admissions" against the demo's 2 520 |
| CARDIA | `metadata_only` | the only candidate with repeated measurements in a population cohort (n = 5 115) |
| OpenCap (SimTK) | `not_applicable` | Q017 model is lumpy and missing geometry parameter |

Verified HTTP 200: `physionet.org/content/eicu-crd-demo/2.0.1/` , `wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm`
, `physionet.org/content/mimic-iv-demo/2.2/` , `physionet.org/content/mimiciv/3.1/`
, `physionet.org/content/hirid/` (→ 1.1.1), `physionet.org/content/eicu-crd/2.0/`
, `biolincc.nhlbi.nih.gov/studies/cardia/` , `simtk.org/frs/?group_id=91`
. Nine rejected/incorrect URLs are logged with HTTP code in `DATA_SOURCES.json:rejected_url_checks`
(among others `mimic-iv/3.1` → 404, error sly; `ukbiobank…/costs`
→ 403; `zenodo.org/records/5137251` → 403; `eicu-crd.mit.edu/static/*.pdf`
→ 200 but HTML shell, not PDF; CDC termination pages
→ 403 or HTTP 200 with "Page Not Found" body).

## Frysade kontroller K1–K9 — `results.json:controls`

| id | quantity | value | limit | outcome | share of limit |
|---|---|---|---|---|---|
| K1 | proportion measured MAP outside model span [34,5054; 183,673] mmHg | 0,00192 (raw 0,05917), n = 504 210 | ≤ 0,02 | PASS (raw: FAIL) | 0,096 |
| K2 | \|median CVP − P_v0 = 5\| | 6,0 mmHg (median CVP 11,0 filtrerad, 12,0 raw) | ≤ 10 | PASS | 0,60 |
| K3 | median eGFR CKD-EPI 2021, eICU-demon | 78,40 mL/min (n = 2 371) | [60; 150] | PASS | pos. 0,20 |
| K4 | median measured urine / model's 0,9824 mL/min | 3,393× (3,333 mL/min = 200 mL/h) | ≥ 10× | **FAIL** | 0,34 |
| K5 | \|π_gc_est − 32\|, eICU-albumin 2,80 g/dL | 3,07 mmHg | ≤ 8 | PASS | 0,38 |
| K6 | percentil av 3 000 mL i per-patient max 3 h-intake | 87,0 (median max-3h 1 070 mL, n = 1 387) | ≥ 1:a percentil | PASS | — |
| K7 | \|medel vuxen-MAP NHANES − 100\| | 9,80 mmHg (medel 90,199; SBP 124 / DBP 72) | ≤ 12 | PASS | 0,82 |
| K8 | \|median eGFR NHANES − 122,8\| / 122,8 | 0,2085 (97,20 mL/min, n = 5 154) | ≤ 0,25 | PASS | 0,83 |
| K9 | V0/73 kg | 0,24233 (median kroppsmassa 78,9 kg; 73 kg = 37,9:e percentilen) | [0,18; 0,30] | PASS | pos. 0,52 |

**8 PASS, 1 FAIL.** K1 primary is the declared filtered statistics; the raw amalgamated
stat (5,92 %) fails and is driven almost entirely by `vitalPeriodic.pamean`, whose
demo distribution is artifact-laden (median 25 mmHg; 29 954 of 32 632 rows fall out
the declared 40–260 mmHg window). The two NIBP-kanalerna are clean (0,70 % and 0,31 % raw).
K1's NIBP-only execution is reported as `sensitivity_nibp_channels_only` and is explicitly
marked as declared after the merged outcome was known; it does not replace the primary.

## Motprov M1–M4

- **M1 channel switching pig (distinguishable).** Within-patient correlation between mean MAP and
  cumulative volume delivery on a common 60 minute patient and hourly grid:
  r = 0,1150 (n = 19 603 patient hours, 433 care sessions). Placebo with 1 000 permutations of
  time order within each patient (seed 20260925, frozen): median 0,0001, 95 % band
  [−0,0198; 0,0189]. Gap 0,1149 > 0,02 → the relation is **not** a time axis artifact. The band
  is too narrow to rule out a completely arbitrary connection between MAP and
  fluid balance in ICU-data; it just rules out that this particular correlation is artifactual.
- **M2 constant-/synthetic-data detector.** No `DEGENERATE` columns. En `ZERO_INFLATED`-flag
  was set to `vitalAperiodic.cardiacoutput` (median 0 L/min, 80,5 % of those 14 235 inaccurate
  values are zero, 94,8 % of column missing, 256 unique values): the column carries no control. Therefore, `cardiac_output_l_min` (Q0 = 4,9 L/min)
  **not** limited by any downloaded sample.
- **M3 self-execution.** `inputs/Q017_model.py` was never imported or executed.
  The only numbers from `inputs/Q017_results.json` are the eleven model anchors quoted
  with their exact JSON path in `results.json:model_anchors` (eg `/full_model/summary/baseline/urine_ml_min`
  = 0,9824). No synthetic values ​​are used.
- **M4 null group.** The null model "no fluid in, no urine out" (U = 0) against the measurement:
  median 3,333 mL/min, 90th percentile 11,667 mL/min (95th 16,667), 0 % of windows exactly zero.
  The urinary tract is not a constant, so K4 is a comparison between two non-degenerate series.

## Hypoteser

- **H1 (partially met).** Eight model quantities have a measured value in the two downloaded
  the samples. Synchronous requirement is met *within* the eICU demo: **1 039 care sessions (41,2 %) have
  MAP, intake, urine and serum creatinine on one and the same time axis**, 96 has all six channels,
  1 340 has MAP + intake. But the requirement "minst fem nedladdade dataset" does not hold: 2 downloaded,
  3 candidates in total. 50 MB-budgeten and the session time was the binding constraint.
- **H2 (not supported).** PREREG expected at least two controls to fall outside the measurement
  data. Only K4 fails. K2, K7 and K8 pass but use 60 %, 82 % and 83 % of their frozen limits,
  and K1 drops if the artifact filter is removed. Symmetric reading: the model's frozen
  values ​​lie *within* but not comfortably *within* the measurement data.
- **H3 (met).** Two datasets were loaded and parsed without login, DUA
  or application; The eICU demo is "Access Policy: Anyone can access the
  files, as long as they conform to the terms of the specified license."

## Source / derivation / assumption

**Source (measured, traceable to a line in a downloaded file):** all K1–K9 stats and all
M1–M4 statistics. Each entry in `results.json` carries `source` with filename and column name.
All six downloaded eICU data files match publisher `SHA256SUMS.txt`, and NHANES
BIOPRO_J's code table specifies 6 401 entries, which matches the parsed file.

**HDerivation (calculated from source data, not measured):** CKD-EPI 2021-eGFR from serum creatinine.
citation: Inker LA et al., "New Creatinine- and Cystatin C-Based Equations for GFR without Race",
*N Engl J Med* 2021;385(19):1737–1749, PMID 34554658 — journal, year, volume, issue, pages and PMID
verified in this run via NCBI E-utilities `esummary`. **The coefficient string
(142 / 0,7 / 0,9 / −0,241 / −0,294 / 0,9938 / 1,012) could not be machine verified from
a retrieved full-text document **: NKF-ASN-uppsatsen that reproduces the equation
(*J Am Soc Nephrol* 2021;32(12):2994–3015, PMID 34556489, PMC8822996, retrieved to
`web_evidence/ckdepi2021_nkj.html`) has the equation in a figure. The coefficients are therefore
taken out of the published equation and marked as assumption in `results.json`, not as
source verified string. No grid term; the female coefficient 1,012 is reported separately as
variant because eICU is de-identified without biological sex. Other derivations: NHANES-MAP
as (SBP + 2·DBP)/3, mean of reading 1–3 (reading 2–3 also reported: 90,253 mmHg);
urine flow as summed urine volume in a 60 minute window divided by 60 (semantics verified
empirically: `nettotal == intaketotal − outputtotal` on each line); max-3-hourly-intake as
sliding 180 minute window; all age and body mass reversals.

**Assumption (neither source nor derivation):** the artifact windows 40–260 mmHg (MAP), 0–30 mmHg (CVP),
2–35 mmHg (PAOP), declared before K1–K9 were calculated; oncotic conversion 6,2 mmHg per g/dL albumin (reproducing
~25 mmHg at 4,0 g/dL) with Donnan gain (1+f)/(1−f) = 1,6667 at f = 0,25, and the alternative 2,5
mmHg per g/dL; NHANES survey weights not applied; 73 kg as the assumed body mass of the model.

## What went wrong / what is unknown

- **K4 FAIL.** The model's `reabsorption_fraction0 = 0,992` gives 0,9824 mL/min
  = 59 mL/h resting value. Measured median in the eICU demo is 200 mL/h: ratio
  3,39× does not reach the ≥ 10× limit. The model's resting point urine is a non-ICU
  value, just as PREREG predicted — but the distance is 3,4×, not ≥ 10×.
- **K5 is coefficient and cohort dependent.** At 6,2 mmHg per g/dL, eICU albumin
  gives π_gc = 28,93 mmHg (PASS) while NHANES adult albumin 4,10 g/dL gives 42,37
  mmHg (FAIL, Δ = 10,37). With the 2,5 mmHg per g/dL option, both fail. Status
  in `results.json` is `COHORT_AND_COEFFICIENT_DEPENDENT` . The model's π_gc
  = 32 mmHg corresponds to 3,10 g/dL albumin during the declared conversion.
- **V0 remains arithmetically but not measured controlled.** 17 690 mL for 73 kg gives 0,24233, in the middle
  the window [0,18; 0,30]. But neither the eICU demo nor NHANES 2017–2018 contains measured
  total body water or extracellular volume per participant; DXA yields fat/lean-/bone compartments
  but not water. The field `independent_measured_ecf_fraction` is `null` with status
  `NOT_CONSTRAINED_BY_DOWNLOADED_SAMPLES`. 73 kg is at the 37,9th percentile of measured
  adult body mass (median 78,9 kg); V0 converted to that median mass becomes 19 120 mL.
- **Hcardiac output and system resistance could not be constrained.** `cardiacoutput` is
  zero diluted (M2) and no unit string for `svr`/`svri` could be verified verbatim:
  `eicu-crd.mit.edu/eicutables/*` is rendered client-side and `/static/*.pdf` returns
  HTML-skalet (18 829 B) with HTTP 200. Therefore nothing in `results.json` that converts
  these columns to the model's mmHg·s/mL. **Oknown**, not assumed.
- **K1 raw statistics fail** (5,92 % of 504 210 MAP values outside the model span),
  driven of `pamean`-artefakter i demot.
- **Synchronous feedback times are completely absent.** τ_S = 5 s, τ_N = 30 s, τ_a = 20 s and τ_A = 300 s
  are the model's own assumptions. No downloaded sample has the required time resolution;
  HiRID (1,1 s–1 min) and eICU-`vitalPeriodic` (5 min) are the only candidates, and HiRID requires
  DUA. This is the biggest remaining gap.
- **NHANES license string could not be quoted verbatim.** The CDC pages do not carry
  a license block; its terms-of-use pages respond 403 or HTTP 200 with "Page Not Found"
  body. The only verified thing is that the files are loaded without login, DUA, application
  or payment. Checked `NO_LICENSE_STRING_ON_THE_SOURCE_PAGE` in `DATA_SOURCES.json`
  ; treatment as a public-domain work is an assumption, not a quoted string.
- **CKD-EPI coefficients are not source verified** (see "Derivation" above). If K3 or K8 are
  to be used as sharp limits, the coefficient string must be read from the original publication
  first. The numbers are not affected by that character-level uncertainty — the medians are
  far from the edges of the window — but the "source" grading does not apply to them.
- **Demo ≠ population. The** eICU demo is 2 520 care sessions from a convenience sample, not
  a random selection of ICU posts.

## Next step

1. Apply for MIMIC-IV 3.1 and HiRID 1.1.1 (DUA + CITI). MIMIC-IV provides the synchronous
   pressure–volume–kidney series in a single care session as required by the question; HiRID is the only one
   identified datasets that can measure or exclude τ = 5–30 s.
2. Maintain a written oncotic conversion (6,2 vs. 2,5 mmHg per g/dL) and report both
   K5 variants until it is pre-registered.
3. Search a dataset of measured total body water per participant (NHANES has no general
   published TBW-fil; candidates to verify: NHANES 2017–2018 BIA or isotopic dilution
   in a cohort with DUA) before V0 may be called measured limited.
4. Do not calibrate `inputs/Q017_model.py` against these numbers. PREREG "What counts as an error" point 4
   forbids it, and the model's primary criterion failed already in `inputs/Q017_RESULTS.md`
   (r_CO,pred = 0,8367347 vs 0,147 ± 0,050). These numbers are a catalog of where
   frozen values have measurement support, not a calibration.

## Reproduktion

```
sha256sum -c PREREG.sha256
python3 scripts/load_sample.py      # -> results_load_sample.json  (shape, dtypes, units, coordinate frame)
python3 scripts/audit_samples.py     # -> results.json             (K1-K9, M1-M4, hypotheses)
python3 scripts/fetch_page.py <url> # -> stdout: HTTP status + verbatim license strings
```

All files are written under `results/BT-DAT-Q017/`. No file outside this directory is written.
