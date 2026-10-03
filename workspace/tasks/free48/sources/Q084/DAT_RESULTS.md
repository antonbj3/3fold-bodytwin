BT-DAT-Q084

# Results — public datasets that constrain the Q084 model

## What it built on

`inputs/NIGHT_PREAMBLE.md`, `inputs/Q084_QUESTION.md` (the Q084 row) and `inputs/Q084_model.py`.
The model was read but **not** rerun and no frozen parameter changed. The interrupted session
had left `inputs/` with question, model and the preceding BT-HX-Q084 session's
PREREG/RESULTS/results.json; `samples/` and `scratch/` were empty and there was no
`DATA_SOURCES.json`, no `PREREG.md` and no `RESULTS.md` here. The work builds on
`inputs/Q084_model.py:27-44` (frozen parameters), `:114-173` (`local_coefficients`) and
`:280-378` (outputs), and on the source they were parameterised against: Destrian et al., PNAS 123(4)
e2519599123 (2026), DOI `10.1073/pnas.2519599123`.

`PREREG.md` was written and hashed (`c0f4fd75…f909f`, `PREREG.sha256`) before any download.
The hash is verified unchanged in `results.json`.

## Vad som hittades

Six datasets are included, of which **four measured**, one simulated and one incompletely interpreted.
All files ≤ 50 MB; the largest individual file is 8.73 MB, total 32.79 MB. Complete
URL/licence/size/format/subjects/parameter allocation is in `DATA_SOURCES.json`.

| ID | Dataset | Licence | Measured | Size | Load |
|---|---|---|---|---|---|
| DS1 | FCS/RICS D and α for eGFP and NCp7-eGFP, PLOS Figshare 1320326, DOI `10.1371/journal.pone.0116921.t002` | CC BY 4.0 | yes | 5.6 kB | LOADED |
| DS2 | Live-cell confocal, osmotisk stress, *D. discoideum* Ax2 GFP-Ran, Zenodo 14142946 | CC BY 4.0 | ja | 6.55 + 8.07 MB | LADDAD |
| DS3 | cryo-ET-segmenteringar *Chlamydomonas*, EMPIAR-11830 via Zenodo 15875786 | CC BY 4.0 | ja | 3.08 MB | LADDAD |
| DS4 | PolyD-NPC Brownian dynamics, Zenodo 14885286 | CC BY 4.0 | **nej (simulerad)** | 6.12 MB | LADDAD |
| DS5 | PDB 1EMA, GFPstruktur 1.9 O | CC0 (wwPDB-policy) | ja | 244 kB | LADDAD |
| DS6 | fluotracify FCS-TCSPC PEX5, Zenodo 8109282 | CC BY 4.0 | ja | 8.73 MB | DELVIS |

13 model parameters and 4 model outputs are affected. All four frozen criteria are met
(`results.json → acceptance.overall_passed = true`).

## Nyckeltal med fil

- **DS1** cytoplasmic eGFP `D = 34 ± 3 µm²/s` (FCS, `α = 0.92 ± 0.08`), nucleus `31 ± 1`
  (FCS) and `28 ± 3` (RICS). The binder construct NCp7-eGFP gives `4.5 ± 1` in cytoplasm against
  `34 ± 3` for free eGFP — that is the model's free/bound separation in measured numbers.
  → `samples_loaded.json` → `DS1_figshare1320326_fcs_rics_eGFP.table`.
  The model freezes `D_alpha0 = 20`; the cohort measures 34 ± 3.
- **DS2** each film is read as `(74, 2, 60, 60)` uint16, 148 TIFF pages = 2 channels × 74
  frames, `finterval = 5.0138 s`, length `366.01 s`, raw range 0–2730 counts. Channel 0 is
  background (≈3–18 counts), channel 1 GFP (≈916–1076 counts). 15 ROIs per archive. Measured A/B on
  the same cells and optics: the GFP channel's mean intensity `973.71` (distilled water) against
  `1041.20` (0.4 M sorbitol), ratio `1.0693`; time-variability CV `1.290 %` against `0.836 %`,
  ratio `0.6480`. Pixel size can **not** be recovered (ResolutionUnit = 1).
  → `samples_loaded.json` → `DS2b_hoffmann2024_osmotic_pair`.
- **DS3** fyra medlemmar (`649_Membrane/npc/Nucleosome/Ribosome.tiff`), vardera
  `shape_zyx = [511, 1024, 1024]` uint8 i en Deflate64-zip (`samples_loaded.json` →
  `DS3_empiar11830_chlamydomonas_cryoET_segmentation.members`). Voxel pitch 0.784 nm from
  the EMPIAR-11830 image set's metadata, thus field of view `802.8 × 802.8 × 400.6 nm`. The central plane
  z = 255 har icke-noll-andel 1.21 % / 0.31 % / 0.97 % / 1.17 %. 16×16-blockkarta (64 px =
  50.2 nm) gives spatial CV **1.93–4.93** with blocks from 0.0 to 66.9 % — the obstacle field is
  patchy, not sinusoidal.
- **DS4** 20 062 files in 2 068 directories. The sweep covers pore diameters 36–80 nm, cargo radius
  0.5–9.0 nm, `k_on` 1e-1/1e-2/1e-3 per ns, `k_of` 1e-2 per ns, attraktionsyta 0–4. Exempel:
  `pore_d40/leng_209/kon_1e-1/kof_1e-2/cargo_r01_5` ger `flux = 1.30469e-05`, 128 seeds,
  167 molecules across the first interface, `rate_ave.txt` with 9 rows (0.016 at interface
  0 to 0.822 at interface 8). The MSD file has 24 rows but **no documented
  time unit**, so no diffusivity is derived.
- **DS5** 1717 ATOM-siter, 149 HETATM (CRO/HOH/MSE), 236 aminosyror, cell 51.77 × 62.85 ×
  70.67 Å. Gyration radius around the centre of gravity `Rg = 16.18 Å`, sphere equivalent
  `R = Rg·√(5/3) = 2.09 nm` mot modellens frysta `tracer_radius_R = 2.3 nm`.

## What failed

- **The decisive source is not publicly available and licensed.** Destrian data is in an
  unlicensed Google Drive folder, the full material is offered "upon reasonable request" (>1 TB).
  These facts fail PREREG S1 and S2. No number from it is used here.
- **`phi_n` cannot be determined absolutely from DS3.** The record lacks a label legend, so it is
  unknown which uint8 value is which molecular class. The non-zero fractions above are therefore
  **not** an obstacle or packing volume fraction. `phi_m`, `phi_n` and
  `interaction_strength` remain `UNCONSTRAINED` in absolute numbers; only the **spatial shape**
  is measured.
- **DS6 does not count as numerical loading.** `.pqres` is a proprietary PicoQuant format; only
  the file header `PQRESLT.1.0.` and the filenames' time indications were read. 750 files, 250 per mode
  (FCS/Trace/TCSPC), duration 0–5447 s, median 2725 s. No photon or correlation arrays
  were decoded, so no `k_on`/`k_off` is acquired from here.
- **Dryad FRAP could not be acquired.** `doi:10.25349/D94C8M` (Sqh-GFP, 4.19 GB, CC0) has the right
  licence and the right subject, but Dryad's file download now responds `401 Unauthorized` without a bearer
  token. Zero bytes were acquired; no number from it is cited.
- **No dataset measures cytoplasmic viscosity**, so `k_visc` and
  `tortuosity_coefficient` lack a measurable counterpart. DS2's intensity ratio 1.0693 is raw counts,
  not calibrated concentration or volume, and is therefore not converted to `V_r`.
- **EMPIAR-10988** (S. pombe, promising cytosol label and coordinate lists) would be the best
  candidate for `phi_m`/`phi_n`; EMPIAR REST responded empty and EBI's FTP path gave 404, and
  no licence could be read. Registered as the next lead.
- **OpenCell** (CC BY-SA 4.0) would give `L` directly from the HEK293T cell volume 1 pL, but
  processed metadata CSVs could not be found in time. The number is not used.
- **BindingDB** has measured `kon`/`koff` only in TSV of 566 MB; REST gives only IC50/Ki/Kd.
  Outside the 50 MB budget.

## Assessment

The task is to find and show data, not to draw biological conclusions, so no such
verdict is stated. What can be said is structural: the model freezes `D_alpha0 = 20
µm²/s` while the only measured eGFP cohort here gives `34 ± 3`; `tracer_radius_R = 2.3 nm`
lies about 10 % above the structural estimate `2.09 nm` from 1EMA and is not contradicted by it; and
the model's obstacle field `phi_n = 0.18 + 0.08·cos(2πx/L)` with `L = 10 µm` does not describe the
field distribution DS3 measures, where non-zero blocks vary from 0 to 67 % with spatial CV up to
4.93 over a field of view of 0.8 µm. The model's geometric assumptions are thus not supplied by
any public licensed measurement found.

## Next step

1. EMPIAR-10988 through PDBj/EMDB API for the explicit cytosol label — the only candidate that
   can give `phi_m` and `phi_n` in absolute numbers.
2. ChlamyAnnotations (Zenodo 15615650, CC BY 4.0, 344 kB) for the label legend missing in
   DS3; it allows `phi_n` to be read in absolute numbers rather than only spatial shape.
3. An FCS dataset for free GFP in cytoplasm with raw curves, not a mean table, to obtain `D`
   and `α` per cell and a `D_eff/D_alpha0` that can be compared with the model outputs.
4. A FRAP or FCS dataset with a licence that can be downloaded without a token, for `k_on0`, `k_off0`
   and `bound_occupancy_*`.
5. The Destrian source: contact the authors about a licensed deposition. Without it,
   the only measurement the model is actually calibrated against is missing.

## Reproduktion

```
sha256sum -c PREREG.sha256
python3 -W ignore scratch/inspect_samples.py     # skriver samples_loaded.json
```

Requires `7z` on PATH: the cryo-ET zip uses Deflate64, which Python's `zipfile` cannot
decode. `scratch/read_chlamy_segmentation.py` is a standalone test of the same slice reading.
