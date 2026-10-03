BT-DAT-Q012

# Public measured data that constrains the Q012 model (diet -> microbial TMA -> FMO3 -> TMAO -> renal clearance)

**Deliverables:** `DATA_SOURCES.json` (29 sources), `results.json`, `samples/` (1 sample, 2 311 B),
`load_sample.py` + `samples/load_sample_output.txt` (the printed demonstration), `PREREG.md` +
`PREREG.sha256` (`6f1bef01badcf9dea1befd477ca61273c343bf1bd1aca0f7407425491155efc3`).
All numbers below are in `results.json` / `DATA_SOURCES.json` and tagged **measured** (stated by the
repository), **derived** (computed here, formula in the JSON) or **UNKNOWN**.

## 1. What was found

The Metabolomics Workbench REST index holds **98 human blood** and **16 human urine** studies with a
measured TMAO feature (`GET /rest/metstat/;;;Human;Blood;;C01104` -> HTTP 200, 54 077 B, 101 analysis
records; `…Urine…` -> HTTP 200, 8 506 B). All 98 blood studies state **CC BY 4.0** in the repository
record itself, and all are retrievable by unauthenticated GET — so the plasma block of the model is
genuinely coverable by public measured data.

The single downloaded sample (prereg rule **S1**: measured human plasma TMAO at subject level, chosen
over every other candidate) is:

| field | value | tag |
|---|---|---|
| source | Metabolomics Workbench **ST000488 / AN000754** — "Sleep apnea and cardiovascular metabolites carnitine, tma, tmao, betain, choline in plasma" | measured |
| URL | `https://www.metabolomicsworkbench.org/rest/study/analysis_id/AN000754/datatable/txt` | measured |
| HTTP / bytes / sha256 | **200** / **2 311 B** (limit 50 MB) / `4f878a73bdd264236f87066cf1f39b1e302388c32f648542b937c20f3ba8e554` | measured |
| licence | **CC BY 4.0** (`license_url` in the study summary) | measured |
| format | TSV data table, 7 columns x 36 rows + companion mwtab (13 504 B) | measured |
| matrix / units | **plasma**, **micromolar** (verbatim from `MS_METABOLITE_DATA:UNITS` and the `/analysis` record; the TSV header itself carries no unit string and none was invented) | measured |
| subjects | **36** (36 distinct subject ids `s_01..s_36` in the mwtab; summary `number_of_samples`=36; the `/factors` endpoint returns 35 rows — repository inconsistency, reported, not resolved) | measured |
| coordinate frame | `N/A (scalar concentration data; no spatial reference in source)` — prereg A4, field present with an explicit value | derived |

Parse (python `csv`, delimiter `\t`) — shape **36 samples x 5 analytes = 180 values**, 0 missing:

| analyte (µM) | n | min | median | max | distinct |
|---|---:|---:|---:|---:|---:|
| tma | 36 | 8.640 | 13.955 | 20.710 | 36 |
| tmao | 36 | 0.910 | 2.645 | 54.590 | 35 |
| carnitine | 36 | 22.480 | 38.420 | 56.230 | 36 |
| choline | 36 | 4.650 | 8.100 | 13.980 | 36 |
| betaine | 36 | 23.070 | 38.270 | 52.960 | 36 |

Checks: **C1** TMAO non-degenerate (35 distinct values, not all-zero) PASS. **C2** parser idempotent
(two parses identical) PASS, plus a cross-file check: the TSV equals the mwtab re-keyed by sample id to
within 0.005 µM = the TSV's 2-decimal rounding. **A1–A6** all pass (`results.json:checks`).

## 2. Which model parameters this constrains — and which it does not

`C_ss` is measured directly (subject-level plasma TMAO, µM). Derived bounds (`results.json`):
median plasma TMA/(TMA+TMAO) = **0.840**, TMAO share = **0.160**, median TMA/TMAO = **5.27**
(range 0.19–17.6). These are **cross-sectional concentration ratios, not fluxes**: they bound `f_fmo3`
only under an extra assumption (one shared steady state, TMAO has no other sink) — recorded as
`separability: false`. Spearman(carnitine, TMA) = **0.980** (scipy cross-checked) while
Spearman(TMA, TMAO) = **0.205**: precursor and microbial product track each other almost perfectly in
this cohort while the hepatic product does not track its own substrate — suggestive of a saturating or
rate-limited hepatic step, but correlation is not flux and nothing was fitted.

**No parameter moved from ASSUMPTION to identified.** `f_colonic_escape` and `y_microbial_conversion`
stay non-identifiable: an observational cohort with no isotope tracer, no diet record and no portal
sampling constrains only their **product**, so the **G3 FAIL** in `inputs/Q012_RESULTS.md` is *not*
resolved. `r_host_umol_h` is not covered (no antibiotic arm). `precursor_intake_g_day` is not covered
(plasma carnitine/choline are endogenous pools, not ingested dose).

## 3. Negative results (recorded, not hidden)

- **FMO3 kinetics** (`fmo3_km_umol_l`, `fmo3_vmax_umol_h`) — **UNKNOWN**: SABIO-DB unreachable
  (curl exit 000, 0 bytes); BRENDA `result_download.php` returns HTTP 200 with a **0-byte** body
  (account required). Prereg criterion (b) fails.
- **Diet block** (`precursor_intake_g_day`) — **UNKNOWN**: the USDA FDC API answers HTTP 200 with the
  shared `DEMO_KEY` (461 hits for "carnitine"), but the three food records fetched (beef 85 % lean, SR
  Legacy, 114 nutrients; greenlemon pomelo drink, 12; spinach mature, Foundation, 42) contain **no
  L-carnitine amount**, and `fdc-datasets.html` returned **HTTP 403**, so even the licence text is
  unobserved.
- **Renal anchor** (`gfr_ml_min_per_kg`) — **UNKNOWN**: all three NHANES URLs returned either a CDC
  "Page Not Found" body served with HTTP 200 (20 905 B) or HTTP 404 (1 245 B). MW **ST002820**
  (AASK plasma cohort, N = 1 001, CC BY 4.0, retrievable) is the only renal-context candidate
  obtained, but its MW factors carry only `Study:AASKG1` — no per-sample GFR in the MW record.
- **TMA-gene carriage** (`y_microbial_conversion`) — **UNKNOWN**: Qiita 10317 returns a
  browser-shim HTML page and needs an account; redbiom API unreachable (curl 000).

## 4. Candidate sources listed but not downloaded (all HTTP-observed, sizes not invented)

`ST002820` AASK renal context N=1001 · `ST002027` 6-week flaxseed diet intervention N=356 ·
`ST003662` exercise time-course N=352 · `ST004611` FMT targeted plasma N=198 · `ST000974` N=366 ·
`ST003177` N=2492 (largest) · plus all 16 human urine TMAO studies (`ST000020` 88, `ST003032` 416,
`ST001047` 140, `ST000291` 45, …). Units for these are **UNKNOWN** (data tables not fetched), so they
are candidates, not constraints. 23 of 29 entries carry a repository-stated licence; 1 was downloaded.

## 5. Next step that would actually move a parameter

Only an **isotope-tracer** dataset (d3-choline/d3-carnitine → d3-TMA recovery with a controlled oral
dose and 24 h urine) can separate `f_colonic_escape` from `y_microbial_conversion`; the 98 blood
studies retrieved here contain none. Within the public repositories probed, no such table was found —
that absence, not the TMAO concentrations, is the binding data gap.

*Scope note (PREREG §6):* the model was not re-run, re-tuned or fitted; geometry/anatomy and
coordinate frame are N/A-with-reason for a pharmacokinetic concentration matrix.
