BT-DAT-Q080

Assignment type: dataset identification. No new simulation, no retraining, no judgment words.
Interpretation note: "No judgment words" in `BRIEF.md` is read as *no value judgment* — no judgment about
whether the Q080 model is right, validated or biologically meaningful. Upstream's
`primary_pass` is reproduced only as a quoted field.

## What it built on

`inputs/Q080_QUESTION.md` (Q080 row, K01–K03), `inputs/Q080_model.py` (frozen equations,
`Parameters`, unit check), `inputs/Q080_results.json` (17 parameters with
`source_or_assumption` + `source_catalog`), `inputs/Q080_PREREG.md` (source list +
`next_resolution_step`), `inputs/Q080_RESULTS.md`. Read after `inputs/NIGHT_PREAMBLE.md`.
Upstream is a read copy: `inputs/` is unchanged, the model was not rerun, no
Q080 results were recomputed.

**Previous interrupted session:** the `BT-DAT-Q080` directory itself contained only `BRIEF.md`,
`inputs/`, `ALLOW_WEB`, `agent.log`, `AGENT_EXIT=1` — no runnable predecessor, only
`inputs/`. What existed was therefore the upstream model, not interrupted dataset work.
It was built on, not rebuilt.

## The catalog — 6 records, each with URL, license, size, format, subjects, constrained parameter

Full machine-readable catalog: `DATA_SOURCES.json`. Below per record; all numbers are
`measured` from `samples/_fetch_log.json` or recomputed from the file on disk.

| id | Dataset | License | Format | Unit | Coordinate frame | Subjects/count | Constrains | Bytes | Live |
|---|---|---|---|---|---|---|---|---:|---|
| D1 | ChEMBL bioactivity, human Adenosine A1 (CHEMBL226) | CC BY-SA 3.0 | REST-JSON | nM (per `standard_units`) | none | 1000 of 9789 matching records fetched (API page limit `limit=1000`), 74 source documents, 837 ligands, 236 assays; Ki 781, IC50 86, EC50 85, Kd 46, Koff 2 | **T4** | 1 473 789 | yes (200) |
| D2 | ChEMBL bioactivity, human β2-adrenergic (CHEMBL210) | CC BY-SA 3.0 | REST-JSON | nM | none | 1000 of 7286 matching records fetched, 127 assays; Ki 227, EC50 406, IC50 365, Kd 2 | T4 (weakly) | 1 504 305 | yes (200) |
| D3 | ChEMBL-kinetiksubset ADORA1 (Kon/Koff) | CC BY-SA 3.0 | REST-JSON | `1/nM.1/min` | ingen | 2 av 2 Koff-poster, **0 Kon-poster**, 1 assay | T4 (negativt resultat) | 2 909 | ja (200) |
| D4 | RCSB PDB **3SN6** β2AR–Gs–VHH, X-ray 3.2 O, P 2₁2₁1 | CC0 1.0 | mmCIF | O (`_atom_site.length_unit` missing; mmCIF standard applied) | kristallfasens ASU, right-handed Cartesius, origo in deposited cell origin — **not the membrane plane** | 1 struktur, 5 polymera kedjor, 10 274 atomsiter | **T1, T2** | 2 016 591 | ja (200) |
| D5 | Human Protein Atlas `proteinatlas.tsv` | CC0 1.0 | zipped TSV | per column in the column name (nTPM/pTPM/nCPM/NX) | none | 20 162 gene rows × 119 columns; ADORA1 and ADRB2 present | T1 — **nothing quantitative** | 7 460 135 | yes (200) |
| D6 | BioImage Archive/Europe PMC receptor topography in nm (S-EPMC9499342 AChR nanoclusters STED+STORM, S-EPMC6460894 CD4, S-EPMC7753144 MET) | CC0 1.0 (archive terms) | record metadata + PDF/AVI | claimed nm, not read from file | **UNKNOWN** | cell count not specified in the records | T2/T3/T6 in principle, **unusable** | 0 data files | **no** |

Largest download 7 460 135 bytes ≤ limit 52 428 800 bytes. Unit conversions and
field names per record are in `DATA_SOURCES.json → catalogue[].constrains[].unit_transformation`.

## The sample downloaded and displayed: `samples/3SN6_b2ar_Gs.cif`

- **Form:** mmCIF, 2 016 591 byte, sha256 `f17b5faa3097784a5d05604cb4946ac7e333a4e8fc731a1f909c265b044a6737`.
- **Enheter:** Oh. `_atom_site.length_unit` is *absent* in this entry, so mmCIF default
  used; the coordinates are consistent with: `_cell` 119.339/64.555/131.240 Oh.
- **Koordinatram:** deponerad kristallografisk asymmetrisk enhet i `P 1 21 1`,
  right-handed Cartesian, origin at the cell origin. **Not a membrane frame** — the xy projection
  below is a projection in the crystal frame. Every claim about membrane packing is therefore a
  *ceiling*, not a measurement of a membrane.
- **Content:** 10 274 atom sites, `auth_asym_id` A = Gsα, B = Gβ, G = Gγ, N = VHH,
  R = Endolysin+β2-adrenergi receptor (443 Cα).

### Derived quantities from the coordinates alone (`derived`, formulas in `results.json`)

| Quantity | Value | Formula |
|---|---:|---|
| Record chain's Rg | 3.167 nm | `sqrt(mean‖Cα−centroid‖²)·0.1` |
| Projected convex hull, xy | 23.7416 nm² | `ConvexHull(Cα_xy).area·0.01` |
| Implied max packing | 42 120.2 receptor/µm² | `1/(23.7416 nm²·10⁻⁶ µm²/nm²)` |
| Implied nearest neighbor | 5.236 nm | `sqrt(A·2/√3)` (hexagonal lattice) |
| Frozen `sigma_patch` / neighbor | **57.3×** | `0.30 µm·1000 nm/µm ÷ 5.236 nm` |
| Frozen `mean_receptor_density` / max packing | **0.00131** (0.13 %) | `55 ÷ 42120.2` |

The numbers are recomputed independently by `verify.py` from the same file and agree with `results.json`.
The page limit of 1000 records is a sample limit, not a population figure; the entire population is
9 789 and 7 286 records respectively and has not been fetched.

Two things follow from this and nothing more. (1) A 300 nm Gaussian patch corresponds to 57
receptor diameters — the frozen "patch" is an agglomerate of many receptors, not a
receptor. (2) 55 receptor/µm² is 0.13 % below the physical ceiling set by receptor size;
it is consistent with a sparse membrane and is therefore not a deviation in itself.

## Kontroller (frysta i PREREG.md)

- **Accession placebo, `pass`:** deliberately wrong accession `CHEMBL279` resolves to
  *Vascular endothelial growth factor receptor 2*, not the A1 receptor. The controls
  `CHEMBL226 → Adenosine receptor A1` and `CHEMBL210 → Beta-2 adrenergic receptor`
  confirmed. The placebo therefore triggered, which is the point.
- **The placebo triggered again:** `6G9J` was requested as β2AR–Gs–GRP101 and returned HTTP 200 with
  658 018 bytes, but the entity description is *Mitogen-activated protein kinase 1*. The file
  remains in `samples/` as a recorded negative and is not counted in anything.
- **Empty-result placebo, `pass`:** overrestrictive ChEMBL query with nonfinite standard_type
  gave 0 records (`api_total_count` 0), so filtering works and no numbers from that route
  are used.
- **Permutation placebo (seed 20260926, 200 draws):** the receptor chain's own Cα cloud,
  coordinate permutation. Observed hull area 23.7416 nm² against placebo median 34.138 nm²
  (range 32.178–35.765) — the cloud is more compact than the random order, i.e. the chain is a
  real compact body. **Known weakness, stated:** the same placebo is *invariant* for
  Rg, because columnwise permutation preserves each axis margin and thus the variance sum.
  The placebo therefore only has power for joint statistics (2D hull area). This is not hidden.

## Kriterierna C1–C8

| Kriterium | Utfall | Kort |
|---|---|---|
| C1 live access | **PASS** | 5/5 `verified_live` records: HTTP 200 and bytes online = bytes on disk = bytes in catalog |
| C2 license specified | **PASS** | 6/6 records have a license string |
| C3 T1–T6 coverage ≥5/6 | **FAIL** | strongly covered: T1, T2, T4 = **3/6**. T3, T5, T6 lack file-level datasets |
| C4 Kd anchoring | **PASS** | see below |
| C5 sample ≤50 MB downloaded and interpreted | **PASS** | 3SN6, 2 016 591 bytes, form/unit/frame reported |
| C6 ≥1 derived physical quantity | **PASS** | 6 of them, recomputed from the file by `verify.py` |
| C7 honest negatives | **PASS** | 7 recorded negatives + 1 `verified_live: false` record |
| C8 no invented measured data | **PASS** | all measurement counts recomputed from the files; upstream numbers quoted, never recomputed |

**C3 fails and is reported as a failure, not reinterpreted.** T3
(`D_R = 0.04 µm²/s`), T5 (`D_L = 25 µm²/s`, `k_clear`, pulse) and T6 (the spatial outputs)
have no file-level dataset in the delivery. The reason is D6: the only public records giving
in-membrane nanometer rats and receptor diffusion are Europe PMC mirrors without a
coordinate table (only PDF/AVI), so they must not count. Counting Briddon 2004, Suzuki 2005
and Calebiro 2013 as datasets would give 6/6, but that is exactly what PREREG's
*anchoring placebo* forbids: they are figure readings, not machine-readable releases.

## C4: where the frozen Kd lies in the measured distribution

`derived` ur `inputs/Q080_results.json`: `Kd = koff_s1/kon_s1_nM1 = 0.66/0.02 = 33.0 nM`.
`measured` ur D1, `standard_type = "Kd"`, `standard_relation = "="`, `standard_units = "nM"`:
n = 45, kvartiler 0.36 / 0.71 / 1.10 nM, min 0.24 nM, max 45 000 nM.

The frozen 33 nM is **outside the IQR, above the 75th percentile and far below max**.
It is an anchoring figure, not an approval. The max value 45 000 nM is an
aggregate affinity record, not a resolved small-molecule Kd, so the IQR is the meaningful
extract. All three frozen sensitivity points (S-a strict, S-b as fetched, S-c censored
out) give identical n = 45, because all fetched Kd records already have relation `=`.

## What failed

1. **C3: 3/6 instead of ≥5/6.** T3, T5, T6 are not delivered.
2. **T4 is only half covered.** 0 Kon records for ADORA1; the 2 Koff records carry the unit
   `1/nM.1/min`, not `s⁻¹`. `kon = 0.02 s⁻¹nM⁻¹` and `koff = 0.66 s⁻¹` remain
   **assumptions** and are labeled UNKNOWN against public data.
3. **D5 gave nothing.** `proteinatlas.tsv` is the per-gene summary, not
   the tissue matrix: 0 per-tissue value columns for ADORA1/ADRB2. A cohort spread of
   receptor density was therefore not computed, and no such spread is claimed.
4. **D6 is unusable.** No coordinate table fetched; `coordinate_frame` is UNKNOWN.
5. **No second independent structure.** 4SNI returned HTTP 404. D4 rests on a single
   crystal structure without cross-replication.
6. **BindingDB:** all three URL forms gave HTTP 404. No BindingDB numbers are reported.
7. **T5 completely without data.** No public measurement of effective 2D ligand spread
   (µm²/s) or ligand dose field (nM·s) was fetched. `D_L`, `k_clear`, `source_rate`,
   `pulse_duration` are and remain assumptions.

## What the next step is

1. **Create a real deposition for the D6 records** (coordinate table in nm, with
   `x, y[, z, frame, localization precision, cell id, time]` and a machine-readable license) so
   T2/T3/T6 go from `verified_live: false` to measurable. This is the single
   highest-value action and requires no new model run.
2. **Look for `kon` in `s⁻¹nM⁻¹` by a route other than ChEMBL** (stopped-flow/SPR original data per
   receptor, e.g. Supporting Information tables) — otherwise the pair `kon`/`koff` is
   permanently untested against public data.
3. **Fetch HPA's tissue matrix** (`rna_tissue_consensus`) if cohort variance in
   receptor expression is to be used; remember that nTPM is not receptor/µm².
4. **For the model**, unchanged: same-cell measured receptor coordinates plus calibrated
   ligand dose field, and a preregistered fit of patch width, effective ligand diffusion and
   kinetics on holdout cells.

## Reproduktion

```
python3 fetch.py                 # 11 requests, logs HTTP status and bytes in samples/_fetch_log.json
python3 analyze.py               # interprets the samples -> results.json
python3 build_data_sources.py    # -> DATA_SOURCES.json
python3 verify.py                # C1..C8 + recomputation of all derived numbers
```

`PREREG.sha256` = `843314f5fe1847d65542aaa10f9cd11b2a073e9c7a8c9317d76300b00fdaf91f`,
matched against `PREREG.md` in `results.json → integrity`. Upstream's `prereg_sha256`
(`6da63a31…`, quoted from `inputs/Q080_results.json`) is another document and is
not used as a requirement on this catalog.
