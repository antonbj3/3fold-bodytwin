BT-DAT-Q005

# Public measured data that constrains the BT-HX-Q005 renal model

**Job type:** data identification. No model re-fit; `inputs/` read-only. Every
number below is in `results.json` (provenance tags `source` / `derived` /
`hypothesis`) or in `DATA_SOURCES.json`. `PREREG.md` + `PREREG.sha256` were
frozen before the first download; 14/14 unit tests pass.

## What was found

8 public candidate sources; verdicts: **1 `CONSTRAINS`, 3 `CONSTRAINS-WEAK`,
1 `NO-CONSTRAINT` (placebo), 3 `OVERIFIERAD`**, plus 2 `DATA-GAP` parameters
(`results.json:parameter_constraint_map`, `DATA_SOURCES.json:sources`).

**Sample downloaded** — eICU-CRD **demo** v2.0 (open to anyone, ODbL 1.0, DOI
`10.13026/gxmm-es70`): 3 of 30 tables, 7,693,200 B ≤ 50 MB, sha256 per file in
`results.json:sample`. 1,841 patients / 2,520 unit stays, 2,414 adult. Form:
`patient` 2,520×31, `lab` 434,660×10, `intakeOutput` 100,466×12 gzip-CSV. Units:
creatinine in mg/dL (3 spellings, verified all-`MG/DL`), urine volume mL, age in
years (`> 89` right-censored, excluded). Reference frame: `patientunitstayid`,
all offsets in **minutes from ICU admission**, no absolute dates, de-identified.

## Constraints found (cohort, n = 2,091 adult stays)

| Model target | Empirical | Verdict |
|---|---|---|
| `gfr_mL_min` 94.99, prior CV 0.15 | CKD-EPI 2021 median **71.04** mL/min/1.73 m² (IQR 43.35–96.84); only **71.1 %** inside the model's frozen 0.5–2×94.99 window; cohort CV 0.471 = **3.19×** the prior | CONSTRAINS |
| `CL_Cr` 140.616 (primary output) | Cockcroft-Gault median **74.46** mL/min; at the model's `CL_Cr` + measured urine flow the implied urine creatinine has median **110.3 mg/dL** | CONSTRAINS-WEAK |
| `Cr_error_cv` 0.05 | within-24 h creatinine CV median **0.110** (IQR 0.061–0.194, n=1,572) — **2.2×** the frozen residual | CONSTRAINS |

Measured inputs only: serum creatinine (14,631 rows), measured urine volume,
age, sex, weight, height. CKD-EPI, Cockcroft-Gault, Du Bois de-indexing and the
mass balance are `derived`, not measured.

## What fell / is still unconstrained

The one sample that fits the budget touches **none** of the molecular
parameters: no albumin, no metformin, no uptake, efflux or accumulation series.
So the model's leading candidate `f_u` and `K_OCT2`, `Jmax_OCT2`, `P_MATE`,
`area_cm2` are **untouched** — `kappa_creatinine` (0.40) and
`jmax_factor_creatinine` (0.067) stay `DATA-GAP` derived factors, not measured.
ChEMBL 37 holds 7 metformin×OCT2/MATE records, only 2 numeric and both
inhibition IC50 (OCT2 1.7 mM, MATE1 250 µM) — not uptake Kₘ, so `K_OCT2`=518 µM
still rests on one in-vitro CHO study. Metrabase (the most on-target public
resource) failed DNS from this sandbox → `OVERIFIERAD`, unused. The empirical
secretion share here is 2.1 % vs the model's 32.4 %, but that is a
formula-vs-formula difference in an ICU cohort and is **not** evidence against
~30 % secretion; a verified 24-h urine creatinine reference interval is
`UNKNOWN`, so no pass/fail is claimed. Every constrained parameter rests on
exactly one source (N3). K02/K08/K10/K11: `UNKNOWN`, absent from this directory.

## Next step

Spend the 50 MB budget on an *eICU/MIMIC full* credentialed pull, or on
NHANES albumin (f_u) + kidney CT volume (`area_cm2`) — the two parameters no
bedside endpoint can reach.

**Repro:** `sha256sum -c PREREG.sha256 && python3 analyze_sample.py && python3 -m unittest -v test_analyze.py`
