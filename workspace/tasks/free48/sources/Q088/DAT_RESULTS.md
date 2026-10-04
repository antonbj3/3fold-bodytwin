BT-DAT-Q088

# Public measured dataset constraining the Q088 model

## What was built on

`inputs/Q088_QUESTION.md` (Q088: ATP recovery, ion balance, pre-registered functional
outcome), `inputs/Q088_model.py` `PARAMETER_TABLE` (lines 94-233) and `REFERENCE_ANCHORS`
(lines 236-292), as well as `inputs/Q088_PREREG.md` "Required next resolution
step" (line 70-72), pinpointing exactly the missing class: matched longitudinal
measurements of O2/perfusion, PCr/ATP, an ion balance proxy, and function
during crossed hypoxia and work, with **measured** radius/CSA.

Previously aborted session (BT-HX-Q088) had left only `inputs/*` and an empty `agent.log`
(0 bytes). No dataset work existed. `PREREG.md` + `PREREG.sha256` were frozen before
first search (`e454e19017515f4de7f2da56e398238e3fb575476943b1790982f53e3460d585`).

Searched in Dryad, Zenodo, and Metabolomics Workbench via the repositories' APIs, so that
license, file size, and checksum come from the API response and not from memory. Full
list of sources, derivations and hypotheses separated is in `DATA_SOURCES.json` .

## Not redone

own sub-agent. **No model parameters changed, no fit** — `inputs/` is untouched
(sha256 unchanged, mtimes left from 25/9). No own measurement data. EnDA_siffror in this answer
is in `DATA_SOURCES.json` or in `sample_load_report.json`.

## Utvalda dataset (5 antagna, 6 uteslutna)

| # | Dataset | License | Size | Format | N | Limits |
|---|---|---|---:|---|---:|---|
| 1 | Dryad `10.5061/dryad.9p8cz8wwr` — fibrgeometri + capillary density, 23 women, rest vs trained leg | CC0-1.0 | 37 741 B | XLSX, 2 ark | 23 | `radius_m`, `oxygen_boundary_transfer_m_s` |
| 2 | Zenodo `7920956` DESANLIS — NIRS muskeloxygenering, ocklusion | CC-BY-4.0 | 168 264 B | XLSX | UNKNOWN | `oxygen_boundary_kpa`, reoxideringstid bakom `T90_F` |
| 3 | Zenodo `5005079` / Dryad `fh4vb2c` — ATP/O-efficiency, 12 arter | CC0-1.0 | 587 594 B | PDF | 0 (isolerade mitokondrier) | `oxygen_cost_mol_o2_mol_atp` (direct measurement of an assumed constant) |
| 4 | Zenodo `5047412` / Dryad `sqv9s4n49` — extra cellular K⁺, vilopotential, AP + force | CC0-1.0 | 259 851 B | 9 XLSX | 0 (mouse) | `ion_leak_rate_s`, `ion_pump_rate_s` (direction, not rate) |
| 5 | Zenodo `4936109` / Dryad `pp56ck7` — q-space fiber character, TA vs SOL | CC0-1.0 | 154 635 B | XLSX | UNKNOWN (mouse) | `radius_m`, non-invasive route |

Nearest misses: `nzv32s3v35` (sEMG + NIRS hemodynamics, CC-BY, **407 MB** → over 50 MB),
`22329660` MCMITO (software only), `19472107` "31P-MRS reproducibility" (MATLAB code only,
**no data**), `4061844` Sarcolab (2 astronauts, miRNA/cytokines → not parameter addressable),
`4726374` (crosslink kinetics → would require reinterpretation, not restriction).

## Test: loaded, form, units, coordinate frame

`GET https://datadryad.org/api/v2/versions/417431/download` → HTTP 200, 38 103 B →
`samples/dryad_9p8cz8wwr_v4.zip` (sha256 `6dbf8692…95a2f85` ). Loaded with
`python3 load_sample.py` → `sample_load_report.json`.

- **Form**: ark `CD+CF` (46 rader × 8 columns), ark `CSA+Per+SFI+CAF+CFi` (138 × 12),
  columns `Sample, Age, Leg, Fiber type, CSA, Perimeter, SFI, CAF, CFi, CFCE, CFPE, Shape adj CFPE`.
- **Units** (from the headers of the deposit): `Area (μm²)` , `CD (caps/mm²)` , `C:F (caps/fibre)`
  ; `CSA` and `Perimeter` have **no** unit in header → UNKNOWN, README specifies µm² and µm.
- **Koordinatram**: 2D transvers snittyta of kryosnittad biopsi; origin **UNKNOWN**,
  axis order not applicable (scalar per fiber/section), micrometer scale **UNKNOWN**,
  fiber orientation not recorded (only the dimensionless SFI). The time axis is cross-sectional:
  the only time-like the axis is `Leg` (Resta/Exercised), sampled 5 days after one
  heavy resistance training.
- **Geometry that was loaded**: CSA median 3 177 µm² (1 518–5 510), perimeter median 233 µm,
  SFI median 1.409, ekvivalentradie median **31.80 µm** (21.98–41.88), C:F median 1.276,
  capillary density median 369 caps/mm².

**N2 placebo control**: measured median radius 31.80 µm against the model's
assumed 25 µm (`radius_m`) = ratio 1.272, inside the frozen sweep 20–40 µm.
The assumption remains in terms of magnitude. Measured C:F 1.276 vs Piiper &
Scheid 1986 anchor "omkring 2" = ~1.6× lower. Reported, **not** rescaled.

**N1 provenance check**: both zip members' sha256 matches Dryad's published
digests exactly (XLSX `28cc4cc1…3b55dcb` , README `19669374…8f4d93` ).

## Privacy findings on the sample (reported, not fixed)

1. README specifies **mm²** for total tissue area, but the column header specifies **µm²**; the order of magnitude
   corresponds to µm². The value was left rescaled.
2. SFI recalculated from the deposited `CSA`/`Perimeter` gives median 1.374 against reported 1.409
   (median relative diff −2.03 %, range −14.9 % to −0.22 %) — columns are not adjacent
   exakta.
3. **The sheets share almost no sample-ID:n**: 17 ID:n is only in the geometry sheet, 16 only in
   `CD+CF`, despite both claiming to be participants × bone. A naive join on `Sample` quietly
   losing data; linking must be via `Age`/`Leg`.

## What wasn't fixed — documented gaps

- **Matched longitudinal time curve: not found.** No public dataset provides in the same
  cohort time-resolved PCr/ATP + O2/perfusion + ion proxy + measured radius + repeated functional
  outcome, crossed over hypoxia and work. Zenodo (type=dataset), Dryad and Metabolomics
  Workbench were screened; hittest title only deposits code. Corollary: `T90_F` , `V_PCr0`
  , `PCr_tau_s` and `t90_ion` remain **model predictions**, not measured contrasts.
- **Ion balance: only partial.** Only qualitative public restrictions exist; no
  public Na/K-ATPase rate or time-resolved human gradient recovery was found.
- **Geometry at the model's own scale: only partial.** Whole-fiber morphology (human) and
  q-space (mouse) are available; no per-territory diffusion radius with specified coordinate
  frame. `radius_m` is therefore still an assumption with order of magnitude support only.
- Finding no fully matching design is explicitly **not**
  an error according to PREREG but a documented gap.

## Reproduction

```
python3 load_sample.py     # -> sample_load_report.json
sha256sum -c PREREG.sha256 # PREREG.md: OK
```

1 thread, no cloud queuing, ~5 min wall clock. `matplotlib` is broken in the environment
(numpy 2 vs numpy 1) and was not used; `numpy/scipy/pandas/openpyxl` was used.
