BT-DAT-Q090

# Results — public measured datasets that can constrain BT-HX-Q090

## Status: negative dat|result, reported as a defeat

No license-clear public measured data could be downloaded and opened within the sample budget of 50 MB. Falsification conditions **F1** and **F2** triggered, **F3** partly. I report this rather than fill the gaps with values from articles.

## What was built on

`inputs/Q090_QUESTION.md`, `inputs/Q090_model.py`, `inputs/Q090_RESULTS.md`. The model has 39 parameter rows, of which **32 are tagged `Antagen`**, 3 `Fryst`, 3 numerical and 1 `OVERIFIERAD, UR MINNET` (`reference_mobility_m = 8.0e-5 m`). Zero rows are measurement-backed. The model's own sensitivity table guides prioritization: `t_PDL` 94.856507 %, `E_PDL` 48.698176 %, `E_TMJ_disc` 0.958975 %.

## Key numbers with file

- Sample: `samples/zmk-tooth-cohort-0.8.zip`, 22 005 589 B, sha256 `387ad734…26f06e`, md5 `b1e6d338…` **matches** the repo's published checksum. `zipfile.testzip()` OK, 11 entries.
- Form: two JPEG renderings 2803×2176 and 2404×3756 px, two notebooks (67 cells), `README.md`, `requirements.txt`.
- **Units: `UNKNOWN`. Coordinate frame: `UNKNOWN`.** No scale, voxel size, affine matrix or origin exists in the archive. Not assumed.
- **Measured data in the archive: none.** The public product is a loader package, not the cohort. The only external file the notebooks reference is an unrelated breast MRI example. Therefore **no** derived geometry was counted.
- Candidates: 3 license-clear and human-measured but too large for a ≤50 MB sample (Poseidon3D 1 029 192 000 B, STS-Tooth 31 823 540 000 B, AlphaDent 4,9 GB); 1 license-restricted (Teeth3DS+); 2 access-excluded (PhysioNet requires DUA+CITI, NLM returned HTTP 403); 3 that are only article PDFs without data upload.

All numbers: `results.json` and `DATA_SOURCES.json`.

## What failed

`pdl_thickness_m` (0.0002 m) remains `UNCONSTRAINED` — F1. No time-series or cohort dataset was found. The anchor 0.08 mm @ 100 N cannot be repaired with data but with the original publication; it remains `UNKNOWN` externally. A prereg typo ("37 of 39") was discovered in review and is reported without changing the frozen file.

## Next step

Request specific files from Poseidon3D/STS-Tooth per patient instead of the entire zip; look for microscopic PDL fiber angles in direction/measured PDL thicknesses (no hits in this search); and verify Bien & Topp (1970) in the original.
