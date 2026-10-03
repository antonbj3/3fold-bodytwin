BT-DAT-Q160

# RESULT — public measured dataset constraining the Q160 model

PREREF: `PREREG.md`, sha256 `6f70a4999e3aa6ea98429d3d4454deddd0262ca8b20afa1102c11d8842c7bf44`
(frozen 2026-09-26 09:46:44 UTC, **before the first download**). Total: `results.json`.
Candidate register: `DATA_SOURCES.json`. Sample files: `samples/` (11 files, 5,3 MB).
Verification: `verify_samples.py` (V1–V4) + `build_outputs.py` (V5–V8 + extraction).

This mission **did not** run the Q160 model and **did not** calibrate anything. The only calculation
 is to extract and check downloaded files.

---

## 1. Main result: Q160's primary reference is not open data

`Q160_model.py:372-406` hardcodes `REFERENCE` — Table 1's 119 `t90` pairs and Table 2's
`D_e`/`D_aq` for 10 solutes — from Takenaka S, Pitts B, Trivedi HM, Stewart PS (2009),
DOI 10.1128/AEM.02279-08, PMCID PMC2655469.

Europe PMC responds to `PMC2655469`: **`isOpenAccess: N`, `inEPMC: Y`, `hasSuppl: N`**, and
`/PMC2655469/supplementaryFiles` gives `errCode 0` / *"Article with id PMC2655469 is not open
access one"* (296 byte XML). **There is no machine-readable deposit of the numbers that Q160 uses
for reference.** `Q160_RESULTS.md` writes that R1 is "Verified (full text retrieved 2026-09-25)";
it is true for the full text but does not mean that the values are openly verifiable, and none of C1,
C2, C3 can be reproduced from a deposit. Registered as `X1` in `DATA_SOURCES.json`.

Corollary: **C2 and C3 are built on a non-reproducible reference.** Not "wrong" — but they are not
falsifiable by outsiders, and it should be in the RESULTS register.

## 2. What was actually downloaded and verified (V1–V8)

11 files, alla `http=200`, alla sha256 registrerade, alla magiska byte-signaturer OK, alla
`Content-Length` matched, **all ≤ 50 MB** (`results.json:V5_form_summary`). V5 form parsing
on 5/5 datasets: **0 error after two corrections by my own reader** (see §4). V8 Anchorage 5/5 PASS.

| ID | Dataset | License | What is measured | Status |
|----|---------|--------|--------------|--------|
| **S1** | Zenodo 10367138 — OCT structure, data to *Polymers* 14:4410 | CC-BY-4.0 | thickness, porosity, pore volume, biovolume, contour, wet weight | **usable** |
| **S3** | Zenodo 14767251 — *Ralstonia* biofilm rheology, 83 CSV | MIT | flow curves, amplitude/frequency sweep | **usable-with-caveat** |
| **S3b** | same entry — moisture content | MIT | wet/dry mass, water content | **usable-with-caveat** |
| **S5** | Zenodo 7627305 — interferometric colony profiles | CC0 | radial height profiles, 9 trunks, 2729 rows | **usable** |
| **S6** | Zenodo 5644626 — PIV metadata of movie | CC-BY-4.0 | water temperature, px/inch calibration | **usable-weak** |

### Measured numbers limiting Q160 (§3 in PREREG)

| Target | Q160's Value | Measured (Source) | Ratio |
|-----|--------------|----------------|-----------|
| **T8** `L` | 5,86·10⁻³ m (false upper limit) | **33,35–335,0 µm** (median 98,9), S1, 18 rows, 7–49 d | Q160 are **17–175× too thick** |
| **T8** `L` | — | `max_height` median **145,5 µm**, max 1651 µm; `mid_height` median 106,9 µm; S5, 2729 rows, 9 stems | confirm µm scale |
| **T1** `eps` | 0,65 (**ASSUMPTION**) | OCT porosity **2,69–21,54 %**, median 11,50 % (S1) | see discordance §3 |
| **T1** `eps` | — | moisture content **90,25 %w/w** (WT, n=9) vs **64,64 %w/w** (`epsB` mutant, n=8) (S3b) | see §3 |
| **T2** `r_m` | 27,1 nm (**ASSUMPTION**) | non-contiguous pore volume** 1911,68–2990,07 µm³ → volume equivalent sphere diameter **15,4–17,9 µm** (S1) | **4 orders of magnitude** larger |
| **T6** `tau_y` | 5 / 200 Pa (**ASSUMPTION**) | the deposit has **no `tau_y` column**; G′/G″ junction not extracted | **UNKNOWN** |
| **T7** `J_det` | 0 cells/(m²·s) | S6 does not contain shear stress values (calibration only) | **UNKNOWN** |
| **T4** `D_eff` against size | 2/10 within factor 2 | **no open deposit found** | **UNKNOWN** |
| **T3** sorption, **T9** `Y_EPS`, **T10** `K_o` | assumptions | no open deposit found | **UNKNOWN** |

## 3. Discordance between measured porosity and Q160's `eps` — reported, not replaced

S1 measures **2,69–21,54 %** pore volume. Q160 assumes `eps` = 0,65, i.e. **φ_cell + φ_EPS = 0,35**.
As numbers, they are ~5–25× apart. **I do not replace Q160's value**, for three reasons that can
 verify from the deposits themselves:

1. **Definition difference.** S1's porosity is a **OCT segmented fraction** of a *cyanotubular*
   biofilm on hard surfaces. Q160's `eps` is pore volume fraction in the Darcy/EPS blocking sense. S1
   **does not** publish voxel size, stack axis order, or segmentation threshold
   (`results.json:...S1...V7_frame.V7c_volumetric`) — without them it is impossible to determine whether
   low porosity is real or an undersegmentation. **UNKNOWN.**
2. **System difference.** Cyanobacteria on glass/epoxy/CNT ≠ odontological model biofilm. To
   transferring 0,11 as Q160's `eps` would be an assumption I must not make.
3. **Mass vs Volume Fraction.** S3b gives 90,25 %w/w water — a **mass**fraction. Without dryness
   matrix density it does not go to porosity. **UNKNOWN** (PREREG §3 T1 the line is just this).

What **is** robust from S3b and that Q160 does not have: the pure EPS contrast. The `epsB` mutant
loses 25,6 percentage units of water (90,25 → 64,64 %w/w) with otherwise the same cultivation. It is one
**experimental rearranged EPS→water retention**, not an assumed `phi_e`.

## 4. Two real finds from the form and unit check (V2/V5/V6)

**a) Encoding typo in S3.** 83 members have suffix `.csv` but are **UTF-16LE with BOM** (`ff fe`)
and are **tab-separated key/value blocks**, not comma tables. A standard reading gives
71 rows × 1 column with NUL bytes — a fake "FAIL form". After correction: e.g.
`WT GMI 7 - amplitude sweep.csv` → 35 rader, header `Project: / Test: 20230223 WT GMI7 Amp Sweep oil
/ Result: Amplitude sweep 1`, and the directory is `GMI1000` ⇒ the subject connection is internally consistent
(V8 PASS). **This is my reader that was wrong, not the file** — I'll say it straight.

**b) Unit defect in S5's own README.** `width: width of the colony. Units of milliliters.`
**A width cannot have volume units.** Observed range 1,925–18,89. The actual unit is
**UNKNOWN** (reasonably etc., but it's a guess and I'm not claiming it). Consequential errors: `hL`, `hR`,
`border_l`, `border_r`, `displacement` are **pixel index**, and the deposit does not provide px→µm scale,
so they **cannot** be used as lengths. `mid_height`/`max_height`/`avg_height`/`offset_*` is on the other hand
explicitly declared in µm and their magnitudes (median 106,9 µm) are incompatible with nm and mm ⇒ V6 PASS.

## 5. "Coordinate Frame" (V7) — the trap that had given a silent error in C3

PREREG §5 V7b requires each dataset to specify whether its length is **spherical cluster radius R**, **slab
half depth d/2** or **full thickness L**. It's not cosmetic: `Q160_model.py:t90_cluster()` I guess
**sphere of radius R**. Wrong interpretation gives a factor 2 silent in C3. Results:

| Dataset | V7a enhet + datumpunkt | V7b geometrisk tolkning |
|---------|------------------------|-------------------------|
| S1 | µm; zero level = **cup solid surface**; no profile, just mean ± SD per cup | **full thickness** of a laterally spread layer. As spherical R → ~2× error |
| S3 | no spatial datum (bulk rheology on scraped colony) | N/A — and that's *exactly* why S3 doesn't transfer to Q160's `tau_y` (see below) |
| S3b | none; mass fraction | N/A |
| S5 | µm; **media background is subtracted ⇒ z = 0 = agar surface** ⇒ full colony heights | `mid_height` is a **height**, not a cluster radius. Using it as spherical R is a **category error** |
| S6 | px/inch ⇒ 11,9–13,1 µm/px | PIV Field **above** the cinema, not in the pores |

**Symbol collision on `tau_y` (most important to forward).** S3's G′/G″ intersection is a
**oscillatory weak stress at 1 Hz–10 rad/s on a scratched colony**. Q160's `tau_y = 5 + 200·phi_EPS`
(Pa) is a **static yield stress in a porous matrix under bulk flow shear**. There are **two different ones
measured quantities with the same symbol**. Even if I had extracted the intersection it would not have been one
substitution for Q160's parameter. This cuts right through Q160's H2 ("the same yield stress governs
compaction and loosening") and should be named as definition unclear in Q160.

## 6. What did NOT succeed — is reported as error, not silent

- **Two of the top candidates could not be fetched from the sandbox.**
  - `10.6084/m9.figshare.11917890` (Jana et al. 2020, CC-BY) — the three small files had passed I4
    (7,8 MB), but `ndownloader.figshare.com`, `figshare.com/ndownloader`,
    `api.figshare.com/v2/file/download` and `ndownloader.figstatic.com` all answer **403**, with
    and without `Referer`/`User-Agent`. The article publishes yield stresses **BS 248 / CD 822 / PF 57 /
    PA 100 Pa** — that is, up to 165× Q160's `tau_y0` = 5 Pa. **Without deposit I can
    however, don't verify it, and I don't count it as verified.** (`M1`)
  - `10.5258/SOTON/D2439` (Snowdon et al., Soft Matter, CC-BY) — **401** on the landing side, **403**
    on all four direct file endpoints. Would have given yield stress 15±4 / 13 / 11 Pa **and**
    creep at 0,25–20 Pa, which **reworks Q160's `tau_b` = 0,5 Pa within the tested window**.
    Not verified. (`M2`)
- **Never downloaded 5 files that had passed I4** (`M4`: `Growth curve.zip` 155 kB,
  `cblaster output` 3,6 MB; `M5`: `monod_diffusion.csv` 462 kB, `all_bootstrap.csv` 741 kB,
  `column_profiles.zip` 3,4 MB, `bootstrap_trajectories.csv` 39,0 MB, `model_predictions.csv` 73 kB).
  The time budget ran out on the T4/T5 search. **Listed as `partly-available`, not omitted.**
- **T4, T5, T3, T9, T10 remains UNKNOWN.** Any open deposit of
  *effective diffusivity versus molecular size*, *penetration time per cluster with specified geometry*,
  *size-dependent sorption in EPS*, *EPS per biomass* or *oedometer module* not found
  (Zenodo, figshare, Europe PMC, Dryad). That T4 is empty is not a weakness in the audit — **it is
  same slot as §1: the field publishes its measurements in articles, not in deposits.**
- **`M6`** (Zenodo 4642554, CC0, rotating-disc rheometry, detachment of *S. gordonii*) is 10–12 MB
  and had passed I4, but is **AVI video without machine readable numbers ⇒ I3 excludes it.** Best T7 candidate
  which is still not a dataset.
- No file above 50 MB was downloaded. No measured figure is constructed. No number in `results.json`
  comes from `Q160_RESULTS.md` and is presented as a measurement.

## 7. Recommended next step (for Q160, not for me)

1. **Rewrite `REFERENCE` as non-reproducible** and mark C2/C3 as weaker than they look.
2. **Declare two different yield stresses** (`tau_y^statisk` for porous matrix, `tau_y^rheo` for G′/G″) in
   instead of a `tau_y`. §5 shows that the merge is a definition error, not an approximation.
3. **Search further on T6/T7 with working output** (human browser or a non-blocking client):
   figshare 11917890 and SOTON/D2439 are worth the download, the two individual numbers there
   (`tau_y` 57–822 Pa, `tau_b` window 0,25–20 Pa) directly touch Q160's tired unloading joint.
4. **Get the 5 remaining ≤50 MB files in `M4`/`M5`** — especially `monod_diffusion.csv`, which is
   the only tracer diffusion data set found.
5. Repeat the porosity measurement on **odontological** biofilm with published voxel size and
   segmentation rule before `eps` is changed from 0,65.

---

*Based on: `inputs/Q160_PREREG.md` §5/§6/§9, `inputs/Q160_RESULTS.md` §2–§5,
`inputs/Q160_model.py` (line 372–406, 494–517), `inputs/NIGHT_PREAMBLE.md`.
None `import *` from foreign code; `~/projects/bodytwin` is not in the sandbox, so someone
internal data cross-referencing against Q160's own records was impossible — that's a limitation, not a choice.*
