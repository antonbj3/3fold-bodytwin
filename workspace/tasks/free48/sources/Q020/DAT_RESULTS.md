BT-DAT-Q020
# Which public measured datasets constrain the BT-HX-Q020 model?

**Role:** data availability audit. No new model, no modification of `inputs/`,
no verdict on Q020's mechanism. All numbers are in `results.json` / `DATA_SOURCES.json`.
PREREG frozen before download: `PREREG.sha256` = `6575aee933aadc6ff09eb49af4433196ba42e30955070c4379370e45a5c541b1`.

## 1. Dom

**31 of 34 parameters in `PARAMS` have no public, individually measured data source.** The model's
microvascular half (`K_f0`, `L_p`, `A_cap`, `sigma_g0`, `sigma_s0`, `D_ratio`, glycocalyx
parameters) and the entire lymph half (`K_lymph`, `P_pump`, `beta`, `sigma_L`, `tau_l`,
`f_p_lymph`) are **pure literature values**. That is precisely where Q020's mechanistic content lies.
Open data anchors **one** parameter (`c_plasma`), binds **three** water sums weakly, and **one
preregistered test failed** (D5, `Pc_cap`).

## 2. What was downloaded and read (sample, ≤ 50 MB, no login)

6 files, 35.7 MB total, largest individual 12.2 MB. All HTTP-verified 200/206.

| Fil | Form | Enheter (ur filens egen header/codebook) | Koordinatram | N |
|---|---|---|---|---|
| `samples/BIOPRO_J.xpt` | SAS XPORT v5, row = participant, key `SEQN` | `LBDSALSI` g/L, `LBXSAL` g/dL, `LBDSTPSI` g/L, `LBXSOSSI` mmol/kg | **N/A** – no geometry | 6401 (5905 with albumin) |
| `samples/BIX.xpt` | SAS XPORT v5, BIS 5 kHz–1 MHz + Cole model | `BIDECF` L, `BIDICF` L, `BIDTBW` L, `BIDFFM` kg, `BIXS*`/`BIXC*` ohm | **N/A** – body level, electrode geometry fixed in protocol | 5311 (4083 with volume) |
| `samples/BMX.xpt` | SAS XPORT v5 | `BMXWT` kg, `BMXHT` cm, `BMXBMI` kg/m² | **N/A** | 9282 |
| `samples/BPX.xpt` | SAS XPORT v5, oscillometry 3–4 readings/arm | `BPXSY*`/`BPXDI*` mmHg | **N/A** – cuff | 9282 (6368 used) |
| `samples/DEMO.xpt` | SAS XPORT v5 | `RIDAGEYR` years, `WTINT2YR` kg (weight **not used**) | **N/A** | 9965 |
| `samples/SC4002E0-PSG.edf` | EDF+ , int16, 2048 B header + 2830 data records of 30 s = 23.58 h | `EEG Fpz-Cz`/`Pz-Oz` uV, `EOG horizontal` uV, `EMG submental` uV, `Resp oro-nasal` –, **`Temp rectal` DegC** | **Signal positions, no rigid 3-D frame:** EEG in 10–20 (Fpz-Cz, Pz-Oz), EOG lateral, EMG submental, nasal thermistor oronasal, **rectal** core temperature. The file carries no landmarks. | 1 (out of a total 197) |

The `EDF` file uses the **blocked** header variant (field by field), not the usual
256 B/signal interleave; the parser autodetects the layout — reading it interleaved makes
the numbers silently wrong. Patient row: `X F X Female_33yr`, start `25.04.89 14.50.00`.

## 3. What measured data says about the model's assumptions

| Prereg | Parameter | Model value | Measured (unweighted) | Outcome |
|---|---|---|---|---|
| D3 | `c_plasma` | 40.0 g/L | albumin `LBDSALSI`, n=5905, mean 40.79, sd 3.45, **[34.0, 47.0]** g/L | **PASS** — 40 g/L is in the distribution |
| D4 | `V_isc0+V_plasma0` (ECF) | 0.200 L/kg | `BIDECF/BMXWT`, n=1985 adults, mean 0.214, **[0.164, 0.266]** | PASS (but near median) |
| D4 | `V_cell0` (ICF) | 0.400 L/kg | `BIDICF/BMXWT`, mean 0.293, **[0.201, 0.402]** | PASS **at the p97.5 edge**; model +39 % against median |
| D4 | sum (TBW) | 0.600 L/kg (=42 L @70 kg) | `BIDTBW/BMXWT`, mean 0.507, **[0.373, 0.656]** | PASS; model +19 % against median, at ~p90 |
| D5 | `Pc_cap` | 20.0 mmHg | MAP, n=6368, **[60.7, 114.0]** mmHg | **FAIL** |
| D6 | `T` | 310 K (37 °C) | `Temp rectal`, 84 900 samples/23.6 h, mean 37.053, **[36.34, 37.94]** °C | **PASS** |
| D7 | `c_isc0` | 15 g/L (0.38·c_p) | — | ** negative reference:Not testable** |
| D9 | `K_f0`, `beta` | 0.05 | — | PASS (countertest returns `no_open_dataset`, no numbers) |

Filter for D4: `BIAEXSTS==1`, `BIDFIT<=2` (the file's own quality requirement), `RIDAGEYR>=18`
→ n=1985. NHANES survey weights are **not** used (PREREG §4.2); the model is defined for a
70 kg reference person, so unweighted percentiles are the right comparison.

**Three findings that matter:**

1. **The model overallocates water to the cells.** 42 L for a 70 kg adult is at ~p90 in
   a measured BIA distribution of 1985 people (median 35.4 L/kg·70). The ICF parameter is the
   largest single deviation (+39 %). Every `V_isc0`/`V_cell0`-dependent output is thus
   systematically shifted — not by chance.
2. **The albumin amount in plasma cannot be confused with total protein without data catching it.**
   `LBDSTPSI` median 72 g/L against albumin 41 g/L in the same file; 1.75× error in `c_plasma` ->
   1.75× error in `pi_p` and thus the entire Starling branch. The method shift between cycles is
   also documented (−4.42 %, r=0.968) — uncorrected comparison against older cycles counts
   as an error according to PREREG §4.3.
3. **D5 failed, and the test was structurally wrong.** `Pc_cap` is *capillary* pressure; NHANES measures
   *brachial artery* pressure, about four times higher. Comparing 20 mmHg against the MAP distribution
   was doomed already in PREREG. I ran it anyway and report the error rather than replace it:
   the test is invalid, not the model. The capillary/arterial ratio lacks an open data source.

## 4. What does NOT exist openly (negative reference, the core of the assignment)

| Block | Parametrar | Status |
|---|---|---|
| Microvascular | `K_f0`, `L_p`, `A_cap`, `sigma_g0`, `sigma_s0`, `D_ratio`, `f_glyc`, `n_glyc`, `k_leak`, `E_min` | only Michel & Curry 1999 / Squire. **No open data.** |
| Interstitium | `c_isc0`, `a_stiff`, `V_k_ratio`, `R_mob_max`, `P_half_mob`, `P_min`, `P_max` | only Guyton 1965/1966 + subcutaneous PIV. **No open data.** |
| Lymph | `K_lymph`, `P_pump`, `beta`, `sigma_L`, `tau_l`, `f_p_lymph` | in-vitro pump studies. **No open data.** |
| Liver/cell | `tau_hep`, `c_prop`, `K_cell`, `c_cell`, `tau_cell`, `P_cell`, `M_body`, `n_Hill` | calibration in the file, no external data. |

**The single worst gap: `c_isc0` = 15 g/L and the ratio 0.38·c_p.** Interstitial
(subcutaneous wick) albumin is measured in half a dozen small human/perfusion studies published
as articles, without individual deposition. Consequence: the entire **retention branch** (`J_p`, `sigma_s`,
`R_ret`) rests on an unverified ratio. Plasma albumin (D3) is the only albumin concentration
in this audit that is individually measured in a public file.

**Geometry/anatomy: no hits.** LIDC-IDRI (1018 thorax CT, DICOM/LPS/mm, CC BY 3.0) was
the first candidate and is **dismissed** — lung nodule CT constrains no parameter in a
fluid model. The model's only anatomical input is four scalars (`V_plasma0`, `V_isc0`,
`V_cell0`, `A_cap` = 3000 m²). The first three are now measured via BIA; `A_cap` is a
morphometric estimate without an open data source.

**Best available candidates, blocked in this sandbox:**
* **OAI** (4796 participants, annual knee MRI with reader-graded joint fluid, i.e. measured swelling
  in mL on the same scale as the model's `dV_isc_L`) — the best path to validate the model's **outputs**,
  but NDA is required.
* **MIMIC-IV / eICU** — the only realistic source of per-patient fluid balance time series against
  which `t90_h`, `tau_l`, `tau_hep` can be fit. Credentialed (DUA + CITI).

## 5. Coverage

`1` direct + `3` weak out of `34` parameters = **11.8 %** with usable open data.
Repos probed with real HTTP: CDC/NCHS, PhysioNet, TCIA, GEO, Zenodo, Figshare,
BioLINCC, OAI/NDA. One probe failed and was logged as failed:
`DEMO_A.xpt` → HTTP 404 (wrong cycle suffix; 1999–2000 is called `DEMO.xpt`).

## 6. What the next step is

1. **Correct `V_cell0`/`V_isc0` against NHANES-BIX** (ICF 0.29 vs the model's 0.40 L/kg) and
   report how `dV_isc` and `R_ret` shift. It is half a day and it affects the results.
2. **Look for ways to bring `c_isc0` into physical existence**: contact the Kawai/Lehr material's
   authors for raw data, or run your own wick measurement. Without it, the retention branch is
   unverified in principle.
3. **Apply for OAI** (and MIMIC-IV if theCredentialed flow works) — it is the only thing that can be
   compared against the model's *outputs* in mL over time.
4. Leave `Pc_cap` as a literature value and rewrite D5 as a test of *MAP ranges* if anything,
   not capillary pressure.

## 7. Filer

`PREREG.md`, `PREREG.sha256`, `analyze_sample.py` (pipeline + countertests), `DATA_SOURCES.json`,
`results.json`, `samples/` (6 files, sha256 in `results.json`).
`inputs/` is read but untouched.
