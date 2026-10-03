BT-DAT-Q127

## Results — public measured datasets that constrain the Q127 model

**Builds on** `PREREG.md` (SHA-256 `ff1dd3b76eb7ad3fd75490fe6839e279e697362766b076a020042ad6ef3d31e5`, frozen before first download), `inputs/Q127_QUESTION.md`, `inputs/Q127_PREREG.md` (frozen parameter table, 32 rows) and `inputs/Q127_results.json`. The model was **not** rerun. `inputs/` is untouched.

### Inventariet

9 datasets, 6 with full credit + 3 with half → **7.5 ≥ 6** (C1). Measured per parameter: `A0, M0, S0, a_base, a_load, fiber_length, k_cycle, r_nuc` = **8 ≥ 6**; of the target group `{M0, A0, k_loss, J_fuse/y_fuse, fiber_length}`, **3** are hit; of `{S0, a_base, a_load, k_cycle}`, **3** are hit (C2).

| id | dataset | license | size | format | n | constrains |
|---|---|---|---|---|---|---|
| DS01 | Vézina 2021, snow bunting pectoralis, per-fiber | CC0-1.0 | 11 538 B | XLSX 15×14 | 15 birds | M0, A0, r_nuc |
| DS02 | Okafor/Wu 2023, satellite cells per fiber over time | CC0-1.0 | 12 213 B | XLSX 38×26 | fiber as replicate, 25/fiber | S0, a_base, a_load, k_cycle |
| DS03 | Hodson-Tole 2020, human in-vivo US thickness | CC0-1.0 | 204 843 B + 568 B ReadMe | XLSX 7 sheets × 99×20 | 9 people | A0, k_turn, k_inc (scale check) |
| DS04 | Lavin 2019, human biopsy WGCNA | CC-BY-4.0 | 10 977 B | XLSX 45 modules | missing in the file | M0 (fiber-type mixture) |
| DS05 | Lassche 2020, human per-fiber force, FSHD | CC0-1.0 | 18 594 B | DOCX | 26 people, 547 fibers | A0, r_nuc |
| DS06 | Katzke 2022, fiber length from CT | CC0-1.0 | 4 031 257 851 B | ZIP | nonhuman | fiber_length |
| DS07 | Son/Ward/Lieber 2024, 896 human muscle lengths | **not specified** | 1,8 MB PDF | article + PDF | 34 cadavers, 896 muscles | fiber_length, A0 |
| DS08 | Human Protein Atlas, skeletal muscle | CC BY 4.0 | per-gene TSV/JSON | TSV/JSON | Karlsson 2021 | S0, D_star, M0 (indirectly) |
| DS09 | GTEx V8–V11, Muscle – Skeletal | GTEx license | e.g. 20,4 MB GCT | GCT/TSV/H5AD | 838 donors | M0, A0 |

### Downloaded sample — loads and verified (C3, C4)

5 filer i `samples/`, alla ≤ 50 MB (max **204 843 B**), SHA-256 i `samples/SHA256SUMS.txt`:

- **DS01** (11 538 B, sha `91324276…`): 1 sheet, 15 rows, 14 columns. Units: `Fiber diameter (um)` = µm (in the column header), `Number of nuclei per fiber` = count/fiber. **Measured: 77,2514 ± 13,16 nuclei/fiber (n = 14, 58,08–102,29); diameter 35,4347 ± 4,7738 µm (28,28–43,25).** Derivation: `π/4·d²` = 628,14–1469,10 µm². The model's `A0 = 4029 µm²` corresponds to d = **71,62 µm**, outside the dataset's range. The `Myonuclear domain` column has **no unit** and fails internal consistency (6,141–25,296 µm²/nucleus derived against 6469,7–10394,2 reported) → excluded from all numbers.
- **DS02** (12 213 B, sha `8e02e2f9…`): 1 sheet `Raw Counts`, 38 rows × 26 columns (13 WT + 13 KO), 0h/24h/48h/72h. **Measured (WT 0h): median 5,0 cells/fiber, 2,0–8,0, n = 25 fibers. WT 72h: median 29,5 (21–43, n = 18).** Derivation: +24,5 cells/fiber over 72 h = 0,340278 cells/(fiber·h). Units documented, coordinate frame = fiber section, no axial coordinates. **Not a fusion rate and no nuclear count.**
- **DS03** (204 843 B, sha `795246a8…` + ReadMe `43980d03…`): 7 sheets, each 99×20 float. ReadMe specifies **mm**, 1–100 % of the gait cycle, 20 measurement locations distal→proximal. **Measured: largest absolute thickness change 3,698658 mm within one gait cycle (peak-to-peak up to 4,49491 mm), 9 people.** Derivation: no day rate can be formed — the file gives change in mm, not absolute thickness.
- **DS04** (10 977 B, sha `1670d47d…`): 45 WGCNA modules, 79–42 492 unique transcripts, kME 0,094–0,692. **No fiber geometry in the file** → no numerical ceiling.

### What was tested and rejected (C-P)

- **R1 Cumming et al. 2024** (`10.1113/JP285675`) — rejected as a *data source*: open article but no machine-readable per-fiber data; values are means in text and figures (Wiley gave HTTP 403 here). Consequence according to falsifier F1: human `M0` and `A0` are downgraded to `literature_only`, and the entire reference `R_excess = 0,8888889` inherits that status.
- **R2 Xenium human muscle** (zenodo 20237760) — rejected: "files are restricted", 671 GB. Not publicly accessible.

### Parameter coverage (C5) — 32 frozen rows

`publicly_constrained` 8 (`M0, A0, r_nuc, S0, a_base, a_load, k_cycle, fiber_length`) · `literature_only` 1 (`k_loss`) · `unconstrained` 23.

### Nyckeltal

- The model's `R_excess = 1,3542475` against reference `0,8888889`, relative error **52,35 %**, `criterion_met = false` (source: `inputs/Q127_results.json`).
- The model's fusion flow **0,0085089 / 0,0022092 / 0,0081554** events/fiber/day (source: same file).
- `S0 = 0,1` vs measured resting pool 2,0–8,0 cells/fiber → **20×–80×** (source: DS02 above).
- `k_loss = 0,00035 day⁻¹` has **no** public measurement; the sensitivity range for `R_excess` is 1,295–1,420 over 3× in `k_loss`.

### Slutsats

1. The two human values the entire model rests on (`M0 = 2,4`, `A0 = 4029`) have **no public machine-readable source** — they are literature means from a single study. All quantitative acceptance criteria in `inputs/Q127_PREREG.md` are therefore `literature_only`.
2. The four parameters creating the retention result — `k_fuse`, `y_fuse`, `k_loss` and thus `J_fuse` — lack publicly measured values from any species. No public datasets report a fusion rate per fiber and day.
3. Per-fiber nuclear counting exists publicly, but only in nonhuman (DS01, snow bunting) and with a unitless column that cannot be used. **Human per-fiber myonuclei: UNKNOWN** (absence of evidence, not evidence of absence).
4. The precursor side is best served: DS02 measures precursor cells per fiber over 0–72 h and reveals that the frozen `S0` is 20×–80× below the measured resting pool.
5. Because the failed result is **insensitive to `k_loss`** and `k_loss` lacks a measurement anchor, no public dataset could have caught or explained the 52,35 % error. The data that would diagnose it does not exist publicly.
6. **No claim about functional retraining.** Nothing in the inventory measures retraining performance; nuclear retention is not motor memory.

### Next step

1. Longitudinal fiber-identified PCM1/DAPI with registered biopsy position — the only measurement that can identify a fiber over time and thus give the per-fiber time series `k_loss` needs. Does not exist publicly; must be done.
2. Ask the Cumming group (R1) whether per-fiber numbers exist as supplement tables. A request moves `M0` and `A0` from literature to measured.
3. Report the `Myonuclear domain` column in DS01 to its authors — with the specified fiber diameter, the values are three orders of magnitude wrong for µm² per nucleus.
4. Do not use DS02 to fit `S0`/`a_base`/`k_cycle` without stating the transfer assumption mouse→human and injury→training; it is assumption, not derivation.

### Korrigering mot prereg

`PREREG.md` wrote "26 parameters"; the frozen table in `inputs/Q127_PREREG.md` has **32** rows. C5 is evaluated over all 32. No criteria, thresholds or target lists were changed — the six extra rows are numerical settings (`dt`, `seed`, `stochastic_fibres`) and width/kinematics constants (`w_D`, `q_atrophy`, `p_protect`, `damage_sigma`, `repulsion`), none of them candidates for measured data.

### Filer

`PREREG.md`, `PREREG.sha256`, `DATA_SOURCES.json`, `results.json`, `RESULTS.md`, `samples/` (5 files + `SHA256SUMS.txt`). `harness/` contains the scripts that fetched metadata, downloaded the samples and generated `results.json`/`DATA_SOURCES.json` — runnable for checking. Released: nothing reused from `inputs/`, no code from other agents, no model run, no internal data.
