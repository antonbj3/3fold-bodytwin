BT-DAT-Q031

# Public measured datasets that constrain the Q031 model

## Kort svar

Nine public, measured datasets mapped; a sample was downloaded and read. All nine are freely available
without a login, and seven of the nine are CC0-1.0. `PREREG.md` was frozen before the search,
sha256 `57b79258f7f4a67ee2b36d0dd23add9dc0aa544cbc4ac9f2200d7c191435d97f` (`PREREG.sha256`).

The only expression that received a **direct measurement basis in the model's own unit**
and is still far from the model's value is the membrane capacitance. The test measures
specific membrane capacitance for **1994 individual cells** in 5 cell lines:

- upmeasured `C_specific = 2,2059 ± 0,5171 µF/cm² = 0,022059 ± 0,005171 F m⁻²` (pooled, n = 1994)
- the model: `membrane_capacitance_F_m2 = 1,0e-2 F m⁻²` , own label "Assumption"
- relative deviation `|Δ|/|modell| = 1,206` (121 %)

Per cell line (mean ± sd, F m⁻²): 95D `0,018629 ± 0,003125` , 95C `0,020343 ± 0,003556`
, H460 `0,021009 ± 0,003803` , A549 `0,024536 ± 0,005730` , H446 `0,025192 ± 0,005384`
. No cell line is close to the model's 1 µF/cm²; the range 0,80–5,55 µF/cm²
covers 1 µF/cm² but as a narrow lower corner, not as a typical value.

**Important caveat:** these are suspended, trypsinized tumor cells, not the MDCK layer
to which Q031's reference anchor applies. The measurement limits a *generic* membrane
capacitance value. It does not validate the Q031 model and is not counted as such.

## Provet som laddades ner

`samples/Pogoda_PNAS2018_DataSource_20160830.xlsx`, **257 130 B** (≤ 50 MB, frozen limit).
md5 `e408dd94073a9b46b3f48caa3f1ac1f8` = registered checksum in the repository;
sha256 `05f958f7270990b77ddafa81e22f8c055bd35e943e238155699586f83c6b4ef4`.

- **Form:** XLSX, 5 worksheets (95C, 95D, A549, H446, H460), 9 columns,
  shapes 290×9, 415×9, 442×9, 410×9, 437×9 → 1994 saturated cells.
- **Units:** read verbatim from the column headings in the file — `Cspecific membrane (μF/cm²)`
  , `σconductivity (S/m)` , `Einstantaneous (kPa)` , `L*(μm)` . `A1 kHz` / `A100 kHz` only
  indicates frequency in header, no unit → reported as UNKNOWN, not taken from article.
- **Coordinate frame:** none. The measurement is cell-resolved but not voxel-resolved:
  no x/y/z, no voxel pitch, no alignment between rows. Prints explicitly, not blank.
- **Recalculation (the only calculation):** `1 µF/cm² = 1e-6 F / 1e-4 m² = 1e-2 F m⁻²`
  , i.e. `value_F_m2 = visde_uF_cm2 × 1e-2` .

**Checking the loading:** the two cell lines whose summary statistics are verbatim in
the Dryad abstract were reproduced from the raw lines. H460 `2,1009 ± 0,3803` vs published
`2,10 ± 0,38` µF/cm² (0,04 %), H446 `2,5192 ± 0,5384` vs `2,52 ± 0,54` (0,03 %); σ H460
`0,9126` vs. `0,91` , H446 `0,8286` vs. `0,83` S/m. The deviations are rounding
level. This shows that the file was read correctly — nothing about the model.

## The Nine Sources (complete fields in `DATA_SOURCES.json`)

Each record: URL/DOI, license, total size, format, file list with size
and checksum, number of samples, which model expression it restricts.

| id | source | license | size | format | samples | limits |
|---|---|---|---:|---|---|---|
| S1 | Dryad `10.5061/dryad.3019k` (Zenodo 5025108) | CC0-1.0 | 67 216 758 B| 2 × XLSX | **1994 cells**, 5 lines | `membrane_capacitance_F_m2` (direct); Na/K/Cl conduction indirectly via impedance |
| S2 | Dryad `10.5061/dryad.vq83bk3rg` (Zenodo 4312130) | CC0-1.0 | 343 095 B | 6 × XLSX | UNKNOWN (4 lines × 2 treatments + 2 A549) | `membrane_capacitance_F_m2` (direct) |
| S3 | Zenodo `10.5281/zenodo.5829478` | CC0 | 6 830 865 B | 2 × CSV + README | 1696 volume time series, 4 replicates, HeLa | `initial_volume_m3` (direct); volume dispersion |
| S4 | Zenodo `10.5281/zenodo.10064853` | CC0 | 897 255 018 B| 4 × XLSX + 5 × MRC | UNKNOWN, 3 amphibian species, enterocytes 228–10 593 µm³ | `area_to_volume_initial_m_inv` , `initial_volume_m3` (direct) |
| S5 | Dryad `10.5061/dryad.r7sqv9sb6` (Zenodo 4459501) | CC0-1.0 | 2 040 054 B| 3 × XLSX | UNKNOWN, 7 ocean areas, plankton | `area_to_volume_initial_m_inv` (direct) |
| S6 | Dryad `10.5061/dryad.vt966` | CC0-1.0 | 1 831 424 B | XLSX | UNKNOWN, unicellular green algae | `area_to_volume_initial_m_inv` (direct) |
| S7 | Zenodo `10.5281/zenodo.4941338` | CC0 | 7 817 747 B | 7 × Origin `.opj` | UNKNOWN, manifold cells + TRPC3-ko | volume `V(t)` and `[Ca²⁺]i(t)` during osmotic step; challenges "channels constantly open" |
| S8 | Dryad `10.5061/dryad.31dj3` (Zenodo 4970083) | CC0-1.0 | 2 792 943 838 B | 39 × ZIP (clamp files) | UNKNOWN, cerebellar nucleus neuron in vivo | `membrane_capacitance_F_m2` (direct), resting (`initial_voltage_mV = -80,4`) |
| S9 | Dryad `10.5061/dryad.7m0cfxpt2` | CC0-1.0 | 58 716 061 B | 5 files | UNKNOWN, *E. coli*, 50 000 generations | `area_to_volume_initial_m_inv` (indirect, lower limit) |

Identifiers, titles, licenses, sizes, dates and the entire file lists are retrieved
from the repositories' public APIs when building `DATA_SOURCES.json` (Dryad `/api/v2/datasets/{doi}`
and `/api/v2/versions/{id}/files` , Zenodo `/api/records/{id}` ).

**Coverage of model expressions:** 3 out of 15 expressions have a direct measurement
basis (`membrane_capacitance_F_m2`, `initial_volume_m3` , `area_to_volume_initial_m_inv`
), 6 have a direct or indirect basis. Counter defined in `results.json → parameter_coverage`
; simulated sources and literature values ​​are not counted.

## What cannot be decided — UNKNOWN

- **No public dataset found** for `pump_density_mol_m2_s` , `pump_na/k_half_saturation_mM`
  , `pump_baseline_activation` , `initial_na_mM` , `initial_k_mM` , `initial_cl_mM` , `initial_fixed_anion_mM`
  , `water_area_fraction` or `hydraulic_permeability_m_s_Pa` . The water permeability
  thus remains as a single quantity converted from the literature without a measurement
  basis, and the bearing surface factor 0,025 as a pure assumption.
- **No public dataset was found with membrane potential, ion abundances and cell volume measured
  simultaneously in the same cell under the same osmotic step.** That is the Q031 question itself.
  S7 is closest (volume + intracellular Ca²⁺ during osmosis), S8 closest on the electrical side
  (Vm + capacitance), but none of them have all three. This is a gap in the public data supply,
  not a weakness in the search that is reported without being filled with estimates.
- Number of samples for S2, S4–S9 is UNKNOWN because the frozen rule allows
  a single sample. Estimated numbers do not count as measurement.
- `cell-dimensions-database.com` (Loeffler) could not be reached on check and has been excluded
  rather than quoted from memory. The Figshare search returned HTTP 404. The Allen Cell Types
  API responded but could not be verified within the time budget and is therefore not listed.
- The units of `A1 kHz` / `A100 kHz` are not declared in the file and have not been taken from the article.

## What the next step is

1. Get S2 (343 kB) and compare against S1 to see if ~2,2 µF/cm² is technology proprietary; that
   determines whether 1 µF/cm² should be raised or whether MDCK-like epithelium is at the lower end of S1.
2. Get S4's two small XLSX files (17 kB + 19 kB) for a measured SA/V and replace
   `area_to_volume_initial_m_inv` with a range instead of a single cell layer value.
3. Search targeted by Na/K-ATPase turnover per area unit and for bearing geometry; the two expressions
   which Q031's own result points out as most assumption dependent.
4. Redo the download for S7 if `.opj` can be read; it is the only way to an osmotic
   time series of simultaneous ion signal found.

## Filer

`PREREG.md` + `PREREG.sha256` (frozen before the search) · `DATA_SOURCES.json` · `results.json`
· `results_sample_load.json` · `samples/Pogoda_PNAS2018_DataSource_20160830.xlsx`
· `scripts/{dryad,search,load_sample,build_data_sources,build_results}.py`

No simulation was run. `inputs/` is untouched. Each number above is in `results.json` with source.
