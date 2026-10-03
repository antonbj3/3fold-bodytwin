BT-DAT-Q140

Status: PASS with three declared gaps. 4 public measured datasets admitted (one per required
category), 1 sample downloaded and parsed, 5 candidates rejected with the failing clause. The model
was not re-run, not refit, not edited.

**Built on.** `inputs/NIGHT_PREAMBLE.md` (node discipline, prereg-before-computation), the renal
question in `inputs/Q140_QUESTION.md`, the parameter table in `inputs/Q140_model.py`, and the
frozen numbers in `inputs/Q140_results.json`. `PREREG.md` (sha256
`2b356dbff951515489f51320ee6be0d734270b105498a9739ae51d78d0975a88`, verified) was written before the
first byte was fetched. The earlier interrupted session had produced only `agent.log` — no
prereg, no data files, no sample — so nothing was reused numerically.

**Not reused, and why.** `NIGHT_PREAMBLE.md` §1/§4/§7 point at `~/projects/bodytwin` and
`/mnt/shared_data`. Both are absent in this sandbox (`No such file or directory`), so the
`MECHANISM_ANCHOR_GRAPH`, `scripts/msk/`, `DATASETS.json` register, OrthoLoad, TLEM, OpenCap, Keast
and the rest could not be opened. Recorded as unreachable, not scored on merit.

**Search priority was inherited, not chosen** (`P2`): from the upstream ±50% sweep in
`inputs/Q140_results.json`, on a baseline of 518.642 mL/min and 489.741 mg.

## What was found

| ID | Category | Dataset | Constrains |
|----|----------|---------|-----------|
| DS1 | parameters | UCSF-FDA Transportal, metformin in-vitro table (HTTP 200, 40 431 B HTML) | `mate_km_mg_ml`, `oct2_permeability_cm_s` — T1 |
| DS2 | cohort | NHANES 2017–Mar 2020, Albumin & Creatinine–Urine `P_ALB_CR` (open, US Gov public domain) | `gfr_ml_min`, `segment_water_reabsorption` — T2/T3 |
| DS3 | geometry/anatomy | BodyParts3D 4.0 OBJ meshes (CC BY-SA 2.1 JP; 142 903 898 B / 64 888 505 B archives) | `cell_volume_ml`, `segment_residence_min` — T5 |
| DS4 | time series | `Metformin_PK_data_Kuan2021.xlsx`, Figshare 13524335 (CC BY 4.0, 42 796 B) | `gfr_ml_min`, `mate_vmax_mg_min` + all clearance outputs — T1/T2/T4 |

T1 and T2 are the only groups that move the reported endpoint: ±50% on secretory capacity spans
465.772 mL/min (89.81 % of baseline), ±50% on GFR spans 100.000 mL/min (19.28 %). DS1 and DS2/DS4
cover exactly those. DS1 supplies measured Km in µM with the cell system and reference for each —
e.g. OCT2 990 µM (HEK293-OCT2, Kimura 2005), MATE1 780 µM (HEK293-MATE1, Tanihara 2007),
MATE2-K 1050 µM (Masuda 2006) — against a model that sets `mate_km_mg_ml` and
`oct2_permeability_cm_s` as bare assumptions.

**DS4 is the important one.** It is human metformin PK *stratified by kidney function* — the rare
public design that separates filtration from secretion instead of reporting one pooled clearance,
and precisely the "timed urine concentration and volume against free plasma" measurement the
upstream run named as missing. Its 42 796 B workbook could **not** be fetched: this sandbox's egress
proxy 403s `ndownloader.figshare.com` even with a browser User-Agent. Metadata is verified via
`api.figshare.com` (HTTP 200, license CC BY 4.0); its sheets, units, time grid and subject count are
UNKNOWN here and no value from it was used anywhere.

## The sample — `samples/P_ALB_CR_NHANES_2017_2020.xpt`

835 600 B (budget 52 428 800 B), sha256 `c493f6b7…6116bf`.

- **Form.** SAS Transport (XPORT) v5, header `HEADER RECORD*******LIBRARY HEADER RECORD!!!!!!!`;
  `pandas.read_sas(format='xport')` → **13 027 × 8**, all `float64`, columns `SEQN, URXUMA, URXUMS,
  URDUMALC, URXUCR, URXCRS, URDUCRLC, URDACT`. Non-null: 13 027 SEQN, 12 510 urine albumin,
  12 509 urine creatinine.
- **Units.** Taken from the CDC codebook, not inferred: urine albumin µg/mL and mg/L, urine
  creatinine mg/dL and µmol/L, albumin-creatinine ratio mg/g. Observed: `URXUCR` 3.54–739.0 mg/dL,
  `URXUMS` 0.21–16 070.0 mg/L, `URDACT` 0.27–11 676.92 mg/g. `URDUMALC`/`URDUCRLC` are comment
  codes, not measurements — their 5.4e-79 extremes are the SAS missing-value sentinel and were left
  uninterpreted.
- **Coordinate frame.** **Not applicable** — a subject-by-variable laboratory table has no origin,
  axis or handedness. This is the explicit "not applicable (dataset is not spatial)" branch of
  criterion A4c, not an omission.
- Method: urinary albumin by solid-phase fluorescent immunoassay (codebook cites Chavers 1984).
  Eligibility: all examined participants aged 3+.

## What fell

- **No open measured (iohexol/inulin) human GFR dataset found.** So `gfr_ml_min = 100.0 mL/min`
  stays an assumption with only a measured-serum-creatinate/derived-CKD-EPI bracket (DS2). UK
  Biobank, BioMe, All of Us, MIMIC-IV, eICU-CRD, CRIC/ARIC are all registration- or credential-gated.
  The biggest accessible route, the Diabetologia 2023 dataset (191 adults, 7476 plasma + 316 urine
  levels), states its data are *not* public and routes requests through Vivli — a larger and better
  stratified series than DS4, logged as a known gated gap.
- **No nephron-segment-resolved human geometry.** BodyParts3D stops at whole organ, so the
  per-segment split of `cell_volume_ml` and `segment_residence_min` is still open.
- **No measured segmental lumen pH series.** `segment_ph` is unconstrained, and its sensitivity is
  UNKNOWN because the upstream run never swept it.
- **T3 and T4 left deliberately unconstrained** (criterion P3). Water reabsorption moves 24 h urine
  by 0.188 mg (0.0384 %) and zero-to-double reabsorption by 0.203 mg (0.041 %) — non-identifiable.
  Manufacturing a "bounding" dataset for those would have been a failure, not a result.
- **Placebo P1 passed:** BioModels Zake2021 PBPK metformin models (open, kidney compartment, `Qgfr`,
  OCT/MATE-flavoured reactions, zero measured content) and the Graham 2011 review value were both
  admitted-looking and both rejected on I1.

## Next step

From an unrestricted host, pull DS4 (42 796 B) and read its sheets, units, dose and per-kidney-
function-stratum subject counts. That one file is what would let filtration, secretion and
reabsorption be *separated* rather than only bounded, and would move T4 off UNKNOWN.

Files: `PREREG.md` + `PREREG.sha256`, `DATA_SOURCES.json`, `results.json`, `samples/`. Every number in
`results.json` is tagged `source` / `derived` / `assumed`; model numbers are copied from
`inputs/Q140_results.json` (sha256 of all 7 `inputs/` files unchanged).
