BT-DAT-Q009

# RESULT — public measured datasets that can constrain BT-HX-Q009

## What was done
Inventory of public measured datasets that can place bounds on `inputs/Q009_model.py`
(BT-HX-Q009: total plasma → f_u,p → passive tissue exchange `k_transport = PS/V_t` →
reversible tissue binding `k_on,t/k_off,t/K_d,t/P_t,total` → `k_elimination`). The model
was not run, calibrated or validated. `PREREG.md` + `PREREG.sha256` were
frozen before the first download (sha256 for PREREG.md and for all six `inputs/` files
verified against the freeze).

Outcome: **3 qualifying sources** (requirement C1 required 6) → **C1 FAIL**, C2–C6 PASS.

## Qualifying sources (complete fields in `DATA_SOURCES.json`)

| # | Source | Licence | Size / format | Units | Units (subjects) | Constrains |
|---|---|---|---|---|---|---|
| S1 | NCATS ADME@NCATS, Harvard Dataverse, doi:10.7910/DVN/21LKWG v105.0 (RELEASED 2025-08-27) | CC0 1.0 (SPDX), certificate from the version API | 676 files / 25,3 GB total; 10 tab files downloaded = 1,78 MB, largest 0,85 MB; TSV | PPBR %, t½ h, Vss L/kg, hepatocyte CLint µL/min/10⁶ cells | 2 828 rows / 1 797 drugs; human filter 1 614 rows; t½ 665; Vss 1 130 | `f_u,p` directly; `k_elimination` via t½; `V_t`/`P_t,total` only via Vss |
| S2 | Gill et al. 2022, J Antimicrob Chemother, doi:10.1093/jac/dkac055, PMC9047675 | CC BY-NC 4.0 (in PMC-XML) | 62 686 B, JATS-XML with 3 machine-readable tables | mg/L, mg·h/L, h, L/h, L; AUC₀₋₂₄ h | n=6 infected DFI + n=6 healthy volunteers; soft-tissue microdialysis + matched plasma | `f_u,p`; only matched plasma+tissue pair; `k_transport` only as AUC ratio |
| S3 | Dynamic distribution of antibiotics in orthopaedically relevant tissues, PMC11582342 | CC BY 4.0 | 384 303 B, JATS-XML, 2 tables | Cmax mg/L and µg/mL; penetration dimensionless ratio; sampling window in h | per row: horse n=6, pig n=8, healthy men n=6 et al. | `k_transport` as tissue/plasma CMax ratio; indirect `P_t,total` |

## Sample downloaded and read with code
`samples/` contains 14 files, all below 50 MB (max 852 658 B). All read by
`scripts/parse_samples.py` and `scripts/build_results.py`; no number is handwritten.

- Shape: NCATS = TSV with 1 header row, e.g. `ppbr_az.tab` → `Drug_ID, Drug, Y, Species`.
- Units: from the dataset’s own documentation, see `source_values` in `results.json`.
- Time axis/coordinate frame: **missing in S1** (endpoint values per drug, no time).
  S2 has an explicit time frame, AUC₀₋₂₄ h. The unit test N2 is thus met for S2 and
  S3 but not for S1 — S1 gives non-empty, dimensioned columns but no time dimension.
- Check: `sha256` saved for each file (`scripts/ncats_adme_fetched.json`,
  `scripts/pmc_samples.json`).

## Key numbers (source = file in `samples/`)
- `ppbr_az.tab`, filter `Homo sapiens`: 1 614 rows, median PPBR 95,43 % → `f_u,p` median
  0,0457 (Q1 0,0142, Q3 0,1481). **The model’s reference `f_u,p` =   0,20 lies at only the 4,34th
  percentile** of the measured human values, i.e. in the most unbound quartile. 1 544 of
  1 614 rows have `f_u,p` ≥ 0,20. Source: `samples/ppbr_az.tab`.
- `half_life_obach.tab`: 665 drugs, median 4,2 h (min 0,065, max 1200) — `k_elimination`.
- `vdss_lombardo.tab`: 1 130 drugs, median 0,95 L/kg (min 0,01, max 700) — `V_t`.
- Gill table 2 (from `PMC9047675.xml`): total plasma AUC₀₋₂₄ 6,27 (1,38) and 14,06
  (3,40) mg·h/L respectively; free plasma AUC₀₋₂₄ 1,30 (0,37) and 2,78 (0,55) mg·h/L respectively; plasma-free
  fraction 0,21 (0,03) and 0,20 (0,02) respectively; tissue AUC 0,82 (0,38) and 1,37
  (0,48) mg·h/L respectively. **Derivation** (ratios recomputed from the table): free/total plasma AUC
  0,207 and 0,198 respectively; tissue/free plasma AUC 0,631 and 0,493 respectively. The paper’s own
  "tissue penetration" (0,66 / 0,54) is not the same ratio — it is defined differently.
- The model’s `REFERENCE_ANCHOR` in `inputs/Q009_model.py` agrees with the table for all
  numbers used; the documented internal discrepancy (1,13 in abstract versus 1,30 in table) concerns
  the abstract value, which is absent from the downloaded table and therefore not used here.

## What failed
- **C1 FAIL**: 3 of 6 qualifying sources, and three parameter classes completely lack qualifying
  material: `P_p`/`K_d,p` measured directly; `k_transport = PS/V_t` measured as a number; and
  `k_on,t`, `k_off,t`, `K_d,t`, `P_t,total` measured on tissue. This confirms
  the preregistered hypothesis that the gaps are largest precisely there.
- 11 candidates registered as `excluded` with reasons, including PAMPA and Caco-2 (N1:
  in vitro permeability without tissue), BBB/HIA (only 0/1 class labels), solubility (non-PK),
  hepatocyte CLint (in vitro, does not constrain in vivo `k_elimination`), linezolid microdialysis
  PMC1489800 (only 8 253 B abstract, no open full text/table), platinum microdialysis
  PMC6936607 (no named licence in the file), PK-DB (registration), TDC/PyTDC (N4: the same
  measurements as S1 repackaged) and a paywalled BJCP article.
- **Unverified candidates, left open**: BindingDB (all API calls 404, landing page
  286 B shell script — this is the most important remaining path to measured `K_d` against albumin
  and AAG and `k_on`/`k_off`) and DruMAP (host did not respond; published contents table
  lists 2 319 human `f_u,p`, 4 408 Caco-2 values, 1141 tissue concentrations, but `Kp,uu`
  only predicted).

## Physical validity
BT-HX-Q009 remains `UNKNOWN`. No measured value from these sources has been fed into the model,
and no calibration or validation is claimed.

## Next steps
1. BindingDB via documented download path (SDF/TSV per target protein HSA P02768 and
   AAG P04217) — the only realistic path to measured `K_d,p` and `k_on`/`k_off`.
2. Check DruMAP’s access conditions and download its raw data; distinguish measured and predicted
   fields (`Kp,uu` is predicted according to the published table).
3. Search specifically for an open machine-readable time series with free plasma **and** free tissue
   (microdialysis) — it is the only parameter direction that can test the model’s
   free-exposure statement empirically, not just aggregate AUC ratios.
4. If no such dataset exists openly, it should be stated as a lasting data-gap
   state for Q009, not as a search failure.

## Files
`PREREG.md` (frozen, sha256 in `PREREG.sha256`), `DATA_SOURCES.json`, `results.json`,
`samples/` (14 files), `scripts/fetch_ncats_adme.py`, `scripts/fetch_pmc_samples.py`,
`scripts/parse_samples.py`, `scripts/build_results.py`.
