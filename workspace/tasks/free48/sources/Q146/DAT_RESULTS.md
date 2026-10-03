BT-DAT-Q146

# Q146 — public measured datasets that can constrain the nose model

**Builds on:** `inputs/Q146_QUESTION.md`, `inputs/Q146_model.py`, `inputs/Q146_PREREG.md`,
`inputs/Q146_RESULTS.md`, `inputs/NIGHT_PREAMBLE.md` §1–3. PREREG frozen before download:
`PREREG.md` + `PREREG.sha256` = `a3d4a94df9e4ceea98a359cdea105444e0e75d7836d5756bf24d6d1b11c6949f`.
Final answer ≤ 300 words. No own measured data; no verdict words.

## Vad som gjordes

Seven sources verified (URL + HTTP response or provider API), each mapped to a
specific symbol in `inputs/Q146_model.py`. One sample downloaded and loaded. Complete
source list with licence/size/format/subjects/constraint: `DATA_SOURCES.json`.

## The sample (meets the requirement: ≤ 50 MB, loaded)

`samples/P_COT.xpt` — **NHANES 2017–March 2020, serum cotinine and hydroxycotinine**.
HTTP 200, **522 560 bytes**, sha256 `2c989209…`. US public domain.
**Shape:** SAS transport V5, table **13027 × 5** (row = participant, column = variable;
`SEQN, LBXCOT, LBDCOTLC, LBXHCOT, LBDHCOLC`).
**Units:** ng/mL in serum (confirmed against codebook P_COT.htm), ID-HPLC-APCI-MS/MS.
**Coordinate frame:** no spatial frame — cross-sectional biochemical cohort, no voxel array
and no time axis. Index frame: participant (SEQN) × analyte; measurement: NHANES MEC morning sample,
fasting, venous serum.

**Measured numbers** (the only measured data in this run):
cotinine n = 11 395, median 0,034 ng/mL, mean 41,5, p5 0,011, p95 302, max 1620 ng/mL.
NM = hydroxycotinine/cotinine, median **0,337** (IQR 0,221–0,488; n = 3365 with
cotinine ≥ 0,3 ng/mL). log10 SD of cotinine **1,54 dex = 35×**.

## Quantitative constraints (PREREG C6)

| Quantity | Model (`inputs/Q146_results.json`) | Publicly measured | Deviation |
|---|---|---|---|
| Cmax spray | 1,18 ng/mL | 2–12 ng/mL (D2, n = 30) | 1,7× below the low end |
| t_max spray | 77 min | 4–15 min (D2) | **5,1× above the high end** |
| F_abs spray | 0,440 | 0,53 ± 0,16 (D2) and 0,58 (R1) | −17 % against D2 |
| f_met | 0,10 (assumption) | not measurable from D1 | **UNKNOWN** |

The two independent public measurements of F (0,53 and 0,58) agree; **the model is
the outlier, not the data**. The time axis is dead: 77 min against 4–15 min for the same route and same dose.
D1 also shows that a single CYP2A6 number (f_met = 0,10) cannot represent a
population with 35× spread in cotinine.

## What failed

- **No public machine-readable data exists for the blocks that determine the model.** `DEPOSITION`
  (mass fraction per region), `chi_swallow` (swallowed fraction per region) and `k_epi`
  (1,6e-4–4,9e-4 m/s required according to the model itself) have **no open measured source**.
  They are tables in papers, not files. The assumptions remain.
- **D6 (Ličen et al. 2026, the source behind REGIONS/DEPOSITION) could not be verified:**
  mdpi.com responds HTTP 403, mdpi-res.com 404. The numbers 150/130/15 cm² and "<1 %" in
  `inputs/Q146_PREREG.md` are thus **not** reconfirmed here.
- **D1 cannot identify E_H.** A steady-state cohort is no nicotine mass balance.
  28,7 % of the rows (3744) are at the detection limit 0,011 ng/mL for both analytes
  and give NM = 1,0 exactly — informationally empty. Non-quotiented values are needed.
- **No public candidate over 50 MB was downloaded** (D4 224 MB, D5 337–532 MB,
  D7 multi-GB). D4 (n = 130, region masks in mm, RAS) is the only open source that can
  measure S_cm2 per region — but the size discipline makes it unavailable here.

## Next step

1. Acquire D4 within the size discipline (streaming subset, or request an exception from 50 MB)
   and measure S_cm2 per region; compare with 10/15/50/80/10 cm².
2. Turn D1 into a real mass balance: cotinine/3'-HC over time ⇒ E_H and f_met.
3. Search specifically for Ussing-chamber Papp for nicotine on human nasal mucosa — the only
   remaining freely set number that determines C3.
4. Update `inputs/Q146_RESULTS.md`: the t_max error is now **two independent public
   pieces of evidence**, not an internal spot check.

## Filer

`PREREG.md`, `PREREG.sha256`, `DATA_SOURCES.json`, `results.json`, `sample_report.json`,
`samples/P_COT.xpt`, `read_sample.py`, `build_results.py`, `zq.py` (API search tool),
`verify/` (two FDA PDFs, 1,3 MB, verification evidence), `p_cot_doc.html` (codebook).
