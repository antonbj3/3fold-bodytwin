BT-DAT-Q154

# Results — public measurement data that can constrain the Q154 model (D-glucose at the BBB)

## Kort svar

I found and downloaded **10 public sources**, all ≤ 50 MB, all HTTP 200, all
verified with `verify_samples.py` (form, units, reference frame) and SHA-256 in
`samples/_download_log.tsv`. **None of them can constrain the model's glucose parameters.**
Seven of ten sources are `constrains: not_applicable`, two are `order_only`, one is
negative evidence. Frozen criteria met for all files; the null model for H1–H4
could not run with data, so the hypotheses are reported as `NOT_TESTED_BY_DATA`
according to the countertest in `PREREG.md`. This is a negative but verified finding, not an
investigation that failed.

## What existed from the interrupted session

`PREREG.md` + `PREREG.sha256` existed, and the hash matches:
`0935ab6f17c928a861c819de74db5febe2818ae0e30b7d87bb9d89438254f5a8` = actual sha256 of `PREREG.md`.
`fetch_chembl.py` existed with the error that `target_synonym__icontains=SLC2A3`Resolvable to
**GLUT14 (CHEMBL4295906)**, not GLUT3. I replaced the lookup with UniProt→ChEMBL
at accession level (S10) and rewrote the script. `samples/` already contained BBBP, B3DB
(regression + LICENSE + README) and a broken ChEMBL lookup — these are verified
and reused, not downloaded again.

## Facit per hypotes

| Hypotes | Utfall | Nyckeltal | Fil |
|---|---|---|---|
| H1 logBB ≈ 0 for passive small polar substances | `NOT_TESTED_BY_DATA` | 0 of 1058 B3DB records contain glucose; median logBB −0,02, 33,0 % have \|logBB\| < 0,3 | `results.json` → `theochem_B3DB_regression.tsv.form` |
| H2 measured Km/Ki/IC50 for GLUT1/3, MRP1, P-gp | `NOT_TESTED_BY_DATA` | **0** glucose-substrate records against 4 targets; ChEMBL contains 13 glucose molecules (e.g. CHEMBL1614854) but no activity against these targets | `samples/chembl_glucose_vs_transporters.json` |
| H3 SLC2A1/3 > ABCC1/ABCB1 in endothelium | `PARTIALLY_MET`, `order_only` | vascular endothelial cells, nCPM: SLC2A3 769,0 > ABCB1 507,4 > SLC2A1 78,8 > ABCC1 16,3 | `results.json` → `hpa_rna_single_cell_type.tsv.zip.form` |
| H4 measured concentrations for `c_lumen_source_mM` / `c_met_lumen_mM` | `NOT_TESTED` | Metabolomics Workbench REST returns 0 assays for D-glucose | `DATA_SOURCES.json` → S11 |

## Source catalog (short)

| id | Source | License | Size | Form | Unit | Reference frame |
|---|---|---|---|---|---|---|
| S1 | MoleculeNet BBBP | `LICENSE_UNVERIFIED` (file lacks a license) | 149 kB | CSV, 2050 records | **missing** — binary class | not applicable |
| S2 | B3DB regression | **CC0 1.0** (LICENSE downloaded) | 378 kB | TSV, 1058 records | log10(C_hjisna/C_blod), dimensionless | brain vs blood, same species |
| S3 | ChEMBL SLC2A1 (CHEMBL2535) | CC BY-SA 3.0 | 1,41 MB | JSON, 829 records | nM (0,3–100 000) | in vitro, substrate identity |
| S4 | ChEMBL SLC2A3 (CHEMBL5215) | CC BY-SA 3.0 | 298 kB | JSON, 158 records | nM (8–100 000) | in vitro |
| S5 | ChEMBL ABCC1 (CHEMBL3004) | CC BY-SA 3.0 | 1,20 MB | JSON, 723 records | nM (1,5–98 000) | in vitro |
| S6 | ChEMBL ABCB1 (CHEMBL4302) | CC BY-SA 3.0 | 1,53 MB | JSON, 1000 of 3729 | nM (6–100 000) | in vitro |
| S7 | ChEMBL glucose × 4 targets | CC BY-SA 3.0 | 111 B | JSON, **0 records** | — | in vitro |
| S8 | HPA RNA single cell type | CC BY-SA 4.0 | 16,3 MB | ZIP→TSV, 3 087 080 rows | **nCPM** (not nTPM) | pan-tissue cell type, **not brain endothelium** |
| S9 | HPA per-gene JSON/TSV ×4 | CC BY-SA 4.0 | 178–243 B | JSON + TSV | qualitative + sparse nTPM | bulk tissue |
| S10 | UniProtKB accessions ×4 | CC BY 4.0 | 611 B | JSON | — | protein record |
| S11 | Metabolomics Workbench | `LICENSE_UNVERIFIED` | 2 B | `[]` | — | `FRAME_UNKNOWN` |

Full URL, `constrains`/`order_only`/`not_applicable` and justification per source:
`DATA_SOURCES.json`.

## Form, unit and frame check

- **Row count against published description:** B3DB 1058 = README's 1058 → `COUNT_VERIFIED`.
  BBBP 2050 = Martins et al. 2012 → `COUNT_VERIFIED`. ChEMBL SLC2A1/SLC2A3/ABCC1
  `complete: true`; **ABCB1 is 1000 of 3729 (26,8 %), `complete: false`** — this is an
  honest partial sample and is marked as such.
- **Units:** logBB dimensionless and declared in the column name; ChEMBL only `nM`
  (`mol/L = nM × 1e-9`); HPA `nCPM`, which is **not** the same as `nTPM` and not
  converted without a stated conversion.
- **Reference frame:** no spatial datasets were obtained, so the reference frame was tested: B3DB =
  brain/blood, ChEMBL = in vitro substrate, HPA = pan-tissue cell type. **HPA S8 is
  explicitly not brain capillary endothelium** and must not be read as a BBB measurement.
- **Size:** 31 files checked with `stat`; none over 50 MB; largest is S8 with
  16 346 850 B.

## What failed, and why

1. **H1 failed on substrate identity, not statistics.** B3DB is a good logBB dataset
   (median −0,02, 33 % within ±0,3) but lacks D-glucose. According to the frozen countertest it must
   not be used as a glucose proxy. H1 is thus not tested.
2. **H2 failed on a pure null finding.** ChEMBL has 1 710 records across the four targets and
   no glucose-substrate records. `km_luminal_in_mM = 1,5`, `vmax_luminal_efflux = 0,06`
   and `ki_comp_mM = 10` are **not measurable from public databases for D-glucose**.
   A first search with `canonical_smiles__icontains=glucose` gave 0 molecules; that
   search was invalid (glucose's SMILES does not contain the string "glucose") and
   was discarded. The negative instead rests on a `pref_name` search that found 13
   actual glucose molecules, so the zero is a real zero value, not a search artifact.
3. **H3 was met only partly and only as ordering.** SLC2A3 beats both pumps, but
   SLC2A1 (78,8) lies **below** ABCB1 (507,4). Since the unit is nCPM and the frame is
   pan-tissue this cannot set `vmax`.
4. **H4 was not tested.** Metabolomics Workbench REST gave 0 assays. An intermediate probe
   guessed `HMDB0000901` as D-glucose; the API answered "Valproic acid glucuronide".
   The guess was discarded and is not used in any number in this report.

## What the next step is

- **Substrate-specific logBB for D-glucose** is needed from a source that actually measures
  glucose at the blood–brain boundary. Open logBB sets do not contain it.
- **GLUT1/GLUT3 Km for D-glucose** must come from the literature's primary sources (measurement in
  erythrocytes or an expressed system), not from ChEMBL, which turns out to be empty for
  this particular substrate.
- **Brain-endothelium-specific expression** requires snRNA-seq with endothelium as its own cell type from
  brain tissue. HPA's pan-tissue panel is insufficient; this is the only expression answer
  the model needs and it is not openly available in that form.
- **Plasma glucose and lactate in mol/L** from an open human study, to replace the
  two explicitly assumed concentrations.

## Reproduktion

```
python3 verify_samples.py     # skriver results.json, skriver ut storlek/http/sha256
sha256sum samples/*
```

`PREREG.md` was not changed after downloading began; criteria, order
and countertest are unchanged since the hash was set.

## Disclaimer

The task frame (`≤ 45 min`, 1 thread) was delayed because HPA's cell-type archive first had to
be found at the right path. `~/projects/bodytwin` does not exist in this container, so no
node IDs from `MECHANISM_ANCHOR_GRAPH.json` can be cited and no internal datasets
(restricted model data, LHDL, OrthoLoad, OpenCap) have been used. No measurement values are manufactured,
interpolated or typical; every number above is traced to `results.json` or
`DATA_SOURCES.json`.
