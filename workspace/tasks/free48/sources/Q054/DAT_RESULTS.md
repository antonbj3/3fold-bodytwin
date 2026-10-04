BT-DAT-Q054

# Data sources that constrain the HX model in Q054

## What the build actually is

This is a **data-source task**, not an HX run. No new simulation, no change to
`inputs/Q054_model.py`. Question Q054 ("how do actual manufacturing deviations affect contact and
load distribution?") is answered by finding **published measured** datasets that can put a
measured value or measurement interval on the quantities HX's `parameter_table()` freezes in.

## Builds on

- `PREREG.md` (frozen, `sha256 78b8acc8…200a7`, verified OK against `PREREG.sha256` before any
  download in this task) — selection criteria U1–U6 and acceptance steps A1–A5.
- `inputs/Q054_model.py::parameter_table()` lines 326–354 — the 18Hx rows that freeze P1–P10 + O1–O4.
- **Previous interrupted session in the same folder:** `PREREG.md` + `samples/` (4 samples) +
  `verify_samples.py` were already written. `agent.log` (2026-09-25 22:03) ends with
  `UnknownError / Failed to execute statement` during a Zenodo search. Three things from that
  session needed correction, see below.

### Three errors from the interrupted session that I corrected

1. **`verify_samples.py` had never been run.** It built paths from `HERE / "ts_hipbone"` when
   the files are in `HERE / "samples" / …`; `rglob` found 0 files and the script died with
   `IndexError` on line 127. I added `SAMPLES = HERE / "samples"` and reran it. It now goes
   through all four samples.
2. **Incorrect source allocation for `fhn_taper`.** The interrupted session traced the files to
   Zenodo 4983030 *"hip joints of saurischian dinosaurs"* because the title matched the search for
   "taper". The record's 8 files are phylogenetic trees (`*.phy`, `*_characters.txt`) — **no taper data
   at all**, and dinosaur geometry is not measured human data (U4). The real source is
   **Zenodo 14981124**. This is reported as `R3_ZENODO_4983030` under `rejected`.
3. **Unit claims reviewed.** Unit and coordinate frame must not be guessed from filenames. The NIfTI unit
   is now read from `xyzt_units` (space code 2 = mm) and `sform_code`, and where the file itself is not
   enough, `UNKNOWN` stands (see GaitNDD below).

## Search paths (A4 null test)

12 registered searches in `api_cache/search_paths.json`, across **three independent catalogues**
(zenodo, physionet, websearch) — the acceptance requirement was at least two. S3 was found by two *different* phrases
("femoral+stem+taper+friction" and "hip+implant+taper+simulator+disassembly") to the same record,
which is a weaker null test than two separate repositories; it is stated as such.

## Accepted sources — 4

Complete fields, couplings and extracts from API responses: `DATA_SOURCES.json`.

| id | source | licence | size | n | format | constrains |
|---|---|---|---|---|---|---|
| `S1_TS_HIPBONE_CT` | [zenodo/18853791](https://zenodo.org/records/18853791) "Corrected hip bone segmentation for the TotalSegmentator small subset" | `cc-by-4.0` | 1 751 873 B | 104 hips | 104 NIfTI-1 masks | **P2** |
| `S2_MEDPELVIS3D` | [zenodo/21757686](https://zenodo.org/records/21757686) "MedPelvis3D" | `cc-by-4.0` | 202 240 B | 99 cases | 99 CSV, 57 landmarks | **P2** |
| `S3_FHN_TAPER` | [zenodo/14981124](https://zenodo.org/records/14981124) "Wear of femoral head taper connections…" | `cc-by-4.0` | 18 806 + 11 257 B | 24 stems (2×2×2, n=3) | 3 XLSX | **P10** |
| `S4_GAITNDD` | [physionet gaitndd 1.0.0](https://physionet.org/content/gaitndd/1.0.0/) (DOI 10.13026/C27G6C) | `Open Data Commons Attribution License v1.0` | 125 + 2×135 000 B | 64 records | WFDB | **P10** |

**S1** — human in-vivo CT, manually corrected segmentation; acetabulum and SI joint are the
most changed regions. Gives the **bone** surface, thus the outer edge of the cartilage layer, not the cartilage surface.

**S2** — 6 of 57 landmarks are acetabular (`Acetabulum_L/R`, `AcetAntWall_L/R`, `AcetPostWall_L/R`).
Gives shape, not a continuous surface.

**S3** — *explanted commercial implants*, measured in vitro: Taperloc Complete (Type 1 taper,
Zimmer Biomet) and Summit (12/14 taper, DePuy), CoCr and ZTA heads, load frame (LF) and
hip simulator (HJS). Disassembly force **1.76–3.44 kN**, mean 2.606, median 2.760 kN (n=24;
LF n=12 mean 2.633, HJS n=12 mean 2.578).

**S4** — 16 control / 20 hunt / 15 park / 13 ALS, 300 s @ 300 Hz per person, 2 shoe channels.
Body mass in kg is in `subject-description`.

## Rejected — 3 (`rejected`)

- `R1_ZENODO_4767469` — only `GD Poster.pdf` (3 549 710 B). No data → U1 fails.
- `R2_ZENODO_20807109` — `access_right='restricted'`, `files=[]`. Gated → may be listed but is
  never counted as verified (U6). Would have been relevant to P9/O1–O3.
- `R3_ZENODO_4983030` — wrong lead, see above.

## Samples loaded (A2)

`samples/` = 11 089 736 B total, largest individual file **1 751 873 B (1.67 MB)** — far below
50 MB. Complete demonstration report: `samples/load_report.json` (shape, size, units, frame).

| sample | shape | size | **unit** | **coordinate frame / orientation** |
|---|---|---|---|---|
| `S1` NIfTI | NIfTI-1 `.nii.gz`, 3D, int16, labels {0,1} | 311×311×431 voxels, 1.5 mm isotropic | **mm** — from `xyzt_units` space code 2 | **RAS+**, from `sform_code=1` (NIFTI_XFORM_SCANNER_ANAT); diagonal affine 1.5 mm, bbox 465×465×645 mm |
| `S2` CSV | CSV, UTF-8 BOM, 10 columns | 99 files × 57 landmarks, XYZ shape 57×3 | **mm** — *declared in the record, not in the file* | **patient-level LPS** — *declared in the record, not verified from the file* |
| `S3` XLSX | XLSX | Fig 10: 25×4; Fig 11: 45×22; dictionary 10×2 | **kN** in the column heading `Force (kN)`; Fig 11 without a declared unit | **not applicable** — tabular data without a spatial reference |
| `S4` WFDB | WFDB `.hea` + 12-bit `.rit`/`.let` | 2 channels × 90 000 samples @ 300 Hz = 300 s | **UNKNOWN in absolute terms** — FSR output, "roughly proportional to the force under the foot"; gain 3000 /mV without a stated unit, **not convertible to newtons** | **not applicable** — per-shoe force channel |

### Incident during the work: `samples/` was emptied externally

After everything was verified and before `RESULTS.md` was written, `samples/` was emptied outside my
commands: `ts_hipbone/` (the extracted NIfTI files) survived, but all downloaded zip/xlsx,
`gaitndd_*` and the two extracted directories disappeared. I ran no deletion.

Since `DATA_SOURCES.json` already carried URL, size, md5 and sha256 for every file,
recovery could be verified: I downloaded again from the same API URLs and got
**byte-identical** files. The check is now stronger than before:

- 4 Zenodo files: `md5sum -c` against `api files[].checksum` → all OK.
- 3 GaitNDD files: `sha256` against PhysioNet's own `SHA256SUMS.txt` → all OK.
- `subject-description.txt` (2362 B) was wrong the first time (`SUBJECT-README` gives 404);
  the right path is `/files/gaitndd/1.0.0/subject-description.txt`.
- The entire chain was rerun from the beginning; the acceptance steps pass unchanged.

## Coverage — the most important number

Calculated programmatically from the coupling table (`api_cache/binding_counts.json`):

- **2 of 14** frozen quantities have a binding source: **P2** (S1 + S2) and **P10** (S3 + S4).
- **0** parameters are *uniquely* bound to a single source.
- **12 of 14 remain UNBOUND:** O1, O2, O3, O4, P1, P3, P4, P5, P6, P7, P8, P9.

Of the four HX outputs that Q054 actually concerns — contact area, mean pressure, peak pressure,
indentation (**O1–O4**) — there are **zero** independent measurements. That is the most important result:
the very quantities that determine the answer to Q054 have no public standalone measurement on this
route, and none of the candidates could deliver them.

## What the sources must not be used for (symmetric scepticism)

- **P2 is not set to a number.** S1 gives bone surface, S2 gives 6 landmarks. Together, they give an
  *interval* for the surface that P2 and P9 jointly delimit. HX's 50.0 mm is a synthetic
  starting point and is nowhere written as measured (A3).
- **P10 is not set to 1800 N.** S3's 1.76–3.44 kN are *disassembly forces in a
  taper/trunnion joint* (neck–head), not the articulation surface and not in-vivo joint load. HX's
  1800 N lies directly in the interval, which is an order-of-magnitude check and nothing else.
  S4 cannot give newtons at all — the unit is unknown — and constrains only the load's *time shape*.
- **P4 remains unbound despite S3.** The taper geometry (d1, d2, contact length, cone angle) is in
  the paper's table, not in the downloaded data files. No diametral clearance is measured in the sample.
- **P1 (femoral head radius) remains unbound.** Searching for "femoral head radius" gave only
  clinical case series, no geometry dataset. None could be derived from S1/S2 either — they are
  acetabular, not humeral.

## Accepted gates

`api_cache/acceptance.json`, run by `check_acceptance.py` after everything was in place:

```
PASS  A1_fields          PASS  A2_sample_loaded     PASS  A3_no_synthetic_as_measured
PASS  A4_null_control    PASS  A5_byte_check        PASS  U5_coupling
==> ALLA FRYSTA STEG PASS
```

## Next steps

1. **P1 must be solved separately** — no public femur geometry was found. The Imperial-35-femur mesh
   (zenodo 167808, 52 112 646 B) is already on disk according to `NIGHT_PREAMBLE.md` §18 but is
   *surface shape without a measured value*; it gives no nominal diameter.
2. **O1–O4 require pressure plates or load–indentation curves.** Nothing found. This is a
   real void, not a search error — the search phrases hit the right field (MRI cartilage, smart garments)
   but both candidates were closed or empty.
3. **P7/P8 (cartilage modulus) and P9 (cartilage thickness)** should be bindable through published
   MRI cartilage segmentation with measured thickness; R2 is exactly that type of source, just gated.
4. `R2` (zenodo 20807109) and a request for access to gated catalogues are routes, not solutions.

## Reproduce

```
sha256sum -c PREREG.sha256
python3 verify_samples.py          # -> samples/load_report.json
python3 build_data_sources.py      # -> DATA_SOURCES.json + api_cache/binding_counts.json
python3 check_acceptance.py        # -> api_cache/acceptance.json
```

## Files in the folder

| file | contents |
|---|---|
| `PREREG.md` + `PREREG.sha256` | frozen before the first download, hash verified |
| `DATA_SOURCES.json` | 4 accepted + 3 rejected, with API extracts and couplings |
| `results.json` | machine-readable summary, cross-validated against `DATA_SOURCES.json` |
| `RESULTS.md` | this document |
| `samples/` | 4 downloaded samples + `load_report.json` |
| `verify_samples.py` | loads the samples and reports shape/size/unit/frame |
| `build_data_sources.py` | builds `DATA_SOURCES.json` from API cache + sha256 |
| `check_acceptance.py` | runs A1–A5 + U5 as a gate |
| `api_cache/` | raw API responses, search log, checksum, acceptance report |
