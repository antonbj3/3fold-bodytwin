BT-DAT-Q100

## Slutsats

Six public **measured** datasets constrain Q100’s geometry, operating point and kinematics, but **no**
public measured dataset constrains Q100’s layered material law. That is the answer to the question, not just
a list. `PREREG.md` + `PREREG.sha256` were created before the first retrieval and verify
(`32b81eee…75b2`). `results.json` carries all numbers with sha256, `DATA_SOURCES.json` all sources.

**Nyckeltal med fil**

| What | Value | File |
|---|---:|---|
| Sources meeting C1–C5 (of 6) | 6, of which 5 direct | `DATA_SOURCES.json` |
| Sources with a loaded sample that was read (C6) | 4 | `samples/sample_report.json` |
| Largest single sample | 28 311 680 B (27,0 MiB) | `samples/trajectories.npy` |
| Total network transfer (logged + estimated) | 39 083 865 B of 52 428 800 B | `results.json:byte_accounting` |
| Measured `psub` in the sample | 459,1–1494,1 Pa (n = 288, mean 997,08) | `samples/measures.csv` |
| Q100 protocol’s `p_sub` | 800 Pa — **inside** the measured interval | `inputs/Q100_model.py` |
| BAGLS sample, glottal area | 620 px² (mask 256×256, values {0, 255}) | `samples/bagls_0_seg.png` |
| BAGLS sample, frame rate | 4000 Hz → 250 µs/frame (declared in `.meta`) | `samples/bagls_0.meta` |
| VF-3D-MRI, NRRD-header | 224³ `double`, LPS, 0,383929 voxel, FOV 86,0001 mm, origin (−39,095; 39,526; 42,801) | `samples/vf3d_thick_nrrd_header.txt` |
| VF-3D-MRI, STL | 28 112 trianglar, bbox 15,04 × 22,41 × 29,56 | `samples/vf3d_thick_frame01.stl` |
| Ikuma Case2, MP4-container | 256×120 px, 904 frames, 158,2 s, tidsbas ger 5,714 fps mot deklarerade 2000 fps | `results.json:sample_loads.DS-04_ikuma_hsv` |

**What failed**

1. `H2` (public data for the layer law) is **disproved for public data**: the candidates
   `k_body`, `k_surface`, `c_body`, `c_surface`, `meniscus_force`, `meniscus_length`,
   `sheet_tension`, `mucus_height` got no hit in the Zenodo API (queries about mucosal wave,
   vocal fold stiffness, hemilarynx, surface tension gave only journal articles). Therefore
   `H3 = NO`: Q100’s layered-to-two-mass contrast **cannot** be decided on public data.
2. `VocalFold11b`/Zheng M5 was not found as a public database record and is also an idealized
   model, not measurement data (C4) — the `F_air` constraint via 3D flow is therefore not covered.
3. `GIRAFE` (Zenodo 13773163): 0 files, no license → failed C2/C3.
4. Aichinger 2016 (375 HSV + audio, 120 people): paper open, database not downloadable.
5. No video frame was decoded (no ffmpeg/cv2/imageio) — the MP4 sources are container-verified
   only, and are stated as such.
6. `DS-02` is ex vivo pig, not human; `DS-03` is n = 1; `DS-06` (DVTD) has only two speakers and
   no vocal folds → weakest entry, marked `indirect`.

**What the next step is**

1. Set `p_sub` per subject from DS-02/DS-03 instead of the protocol’s 800 Pa, with
   training/validation split per person.
2. Scale the BAGLS masks’ px² to mm² per recording and measure open/contact fraction against Q100’s
   `contact_percent` 30,792 % (baseline) / 14,221 % (layered).
3. Retrieve DS-03’s NRRD volume (one is 89 916 629 B) for thickness and volume per phase and identify
   `k_body`/`k_surface` with a published modulus prior — the only route keeping the layer hypothesis testable.
4. Mucosal wave speed remains the declared proxy `c_m = sqrt(T_sheet/(rho_mucus*h_m))`; measurable
   wave speed was not publicly available and was not found.

## What it built on

`inputs/Q100_QUESTION.md` (mekanism + fyra utdata), `inputs/Q100_model.py` (41 `Parameters` med
enheter, `F_air`/`F_contact`/`F_meniscus`, proxyformeln), `inputs/Q100_RESULTS.md` +
`inputs/Q100_results.json` (`measured_inputs_used=false`, `geometry_measured=false`,
`mucosal_wave_empirical_validation=UNKNOWN`; f0 302,148/138,910 Hz), `inputs/NIGHT_PREAMBLE.md`
(source vs derivation, explicit UNKNOWN, write only here). `DATA_SOURCES.json` +
`results.json` + `samples/` + `scripts/` is all new in this directory; no previous
BT-DAT-Q100 file existed, so no interrupted work was resumed — however, this was built on
the BT-HX-Q100 run’s own gaps (`data_status`).

## Sources (license · size · format · subjects → what they constrain)

| # | Source | License | Size | Format | People | Constrains |
|---|---|---|---:|---|---|---|
| DS-01 | BAGLS, Zenodo 3762320, [10.5281/zenodo.3762320](https://doi.org/10.5281/zenodo.3762320) | CC-BY-NC-SA-4.0 | 17 571 433 690 B | PNG + mask + JSON + mp4 | 640 recordings, 7 clinics, 380 healthy/262 ill/50 unknown, 1000–10 000 Hz, 256×120–512×512 px | `area0`, `area_contact`, `m1`, `m2`; outputs f0 and contact fraction |
| DS-02 | Ex vivo pig, subglottal pressure, Zenodo 10640031, [10.5281/zenodo.10640031](https://doi.org/10.5281/zenodo.10640031) | CC-BY-4.0 | 4 340 255 187 B | mp4 + CSV + npy | 6 larynges, 288 recordings, gap 0/1/2 mm | `p_sub`, `area0`, `area_slope`, `area_contact`, `flow_resistance`; contact fraction |
| DS-03 | Dynamic 3D-MRI of vocal folds, Zenodo 19629778, [10.5281/zenodo.19629778](https://doi.org/10.5281/zenodo.19629778) | CC-BY-4.0 | 5 228 064 118 B (159 zip members: 60 nrrd, 60 stl, 12 wav, 6 png, 2 mp4, 2 gif) | NRRD + STL + WAV | 1 person, 6 phonation types × 10 phases | `k_body`, `k_surface`, `m1`, `m2`, `area0`, `area_slope`, `area_contact`; f0, contact fraction |
| DS-04 | Ikuma et al., bifurcations, Zenodo 4928535, [10.5281/zenodo.4928535](https://doi.org/10.5281/zenodo.4928535) | CC-BY-4.0 | 6 212 491 B (4 mp4) | mp4 2000/4000 fps | 4 cases: healthy, unilateral paresis, polyp ×2 | f0, contact fraction; falsifies `k_contact_surface` |
| DS-05 | Ultrafast 3D-MRI phase II, pathologies, Zenodo 19915187 | CC-BY-4.0 | 72 124 458 B (5 mp4) | mp4, audio removed | 5 cases: polyp, MTD, scar, nodules, paresis | f0, contact fraction per pathology |
| DS-06 | DVTD, figshare 11897187, [10.6084/m9.figshare.11897187](https://doi.org/10.6084/m9.figshare.11897187.v1) | CC0 | 1 129 467 424 B | MRI-STL + FEM + aeroacoustics | 2 speakers, 22 German sounds | only f0 indirectly (`indirect`) |

## Proven laddar (form, enheter, koordinatram)

- **DS-02** `trajectories.npy`: ndarray `float64`, shape **(288, 2048, 6)**, axes = recording ×
  frame × glottal wall point; the 6 columns are left posterior/medial/anterior and right
  posterior/medial/anterior (the record’s own column order, dictionary). Values −36,157…34,044.
  **Units and coordinate frame are not declared** by the record → `UNKNOWN` (scale O(10) makes mm
  plausible but unverified). `measures.csv`, however: `psub` in **Pa**, `gap` in **mm**,
  `glottis_length` in **px**, adduction in **mNm**. md5 matches the record’s metadata.
- **DS-01** `bagls_0.png`: PNG, PIL mode RGB, **(256, 256, 3)** `uint8`, intensity 0–253,
  image grid in px, resolution 256×256 declared in `.meta` (`Color=false` but the file has three
  identical channels — noted). `bagls_0_seg.png`: **(256, 256)** `uint8`, values {0, 255},
  glottal area **620 px²**; px² cannot be converted to mm² without the endoscope working distance that
  the record lacks. `.meta` gives `Sampling rate (Hz) = 4000` → 250 µs/frame.
- **DS-03** `vf3d_thick_frame01.stl`: binary STL, 28 112 triangles, file size = 84 + 50·n
  (consistent), vertex array (28112, 3, 3), bbox 15,04 × 22,41 × 29,56; **STL declares no
  unit**. `vf3d_thick_nrrd_header.txt`: NRRD0004, `type: double`, `sizes: 224 224 224`,
  `space: left-posterior-superior`, `space directions` 0,383929 per axis, `space origin`
  (−39,095104; 39,526318; 42,801071) → FOV 86,0001 mm per axis. The unit mm is **derived** from
  spacing and the filename token `08mm`; the header lacks `space units` → marked as unverified.
  The entire volume (89 916 629 B) was not loaded, only the header (1 237 B) via byte-range.
- **DS-04** `ikuma_case2_paralysis_2000fps.mp4`: MP4 with `ftyp/free/mdat/moov`, 256×120 px,
  904 frames (`stsz`), `media timescale 10240`, `stts`-delta 1792 → **5,714 fps** and 158,2 s,
  while the record declares 2000 fps. The conflict is reported, not smoothed over: **f0 must not
  be calculated from this file** without an external time base. No pixels were decoded.

## Method and reproducibility

1. `scripts/fetch_bagls_sample.py` — reads BAGLS `test.zip` (207 236 572 B) with HTTP Range:
   4 000 000 B tail for the central directory (10 501 entries), then `test/0.png`, `test/0_seg.png`,
   `test/0.meta`; 4 028 979 B transferred.
2. `scripts/fetch_vf3d_members.py` — reads ZIP64 tail (200 000 B) to find
   `PK\x06\x06`, extracts an STL member (1 405 684 B) and the header (1 237 B) from an
   89 916 629 B NRRD member; 840 639 B transferred.
3. `scripts/check_samples.py` — reads everything in `samples/`, writes `samples/sample_report.json`
   with shape/dtype/value range/units/coordinate frame + md5/sha256.
4. `scripts/make_results.py` — builds `results.json` from `PREREG.sha256`,
   `DATA_SOURCES.json` and `samples/sample_report.json`; no numbers are written by hand.

Run: `nice -n 19 python3 scripts/fetch_bagls_sample.py && python3
scripts/fetch_vf3d_members.py && python3 scripts/check_samples.py && python3
scripts/make_results.py`. `samples/network_usage.log` logs every range retrieval; two previous
iterative retrievals (4 028 979 B + 828 351 B) occurred before the log existed and are explicitly counted
in `results.json:byte_accounting` (39 083 865 B total, under 50 MiB).

## Source status per claim

- **Source** (verbatim from retrieved record/paper page): all URLs, licenses, file sizes,
  subject counts, sampling rates, resolutions, the NRRD header, md5/sha256.
- **Derivation** (computed by script from downloaded bytes): glottal area 620 px², FOV 86,0001 mm,
  bbox, 250 µs/frame, 5,714 fps, 158,2 s, psub interval.
- **Hypothesis**: H1 supported, H2 refuted for public data, H3 = NO.
- **UNKNOWN** (explicit): mucosal wave speed as a measurement; DS-02 unit/frame; DS-01 px→mm;
  DS-03 space unit; DS-04 actual frame rate; whether any subject matches Q100’s
  operating point (none does). No measurement data have been copied into any Q100 parameter, and
  `p_sub = 800 Pa` is untouched.
