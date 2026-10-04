BT-DAT-Q021

## Results (short answer)

Seven public candidates found, of which **one is admissible and downloaded**: `doi:10.11588/DATA/ZPEZPX`
(heiDATA, CC BY 4.0, 19 COPD patients, 76 measurements, 12 regions). It loads, has units and
coordinate frame in the file's own data dictionary, and gives a measured number compared against the model's frozen
values. **No model parameter was changed** (`results.json:model_parameters_changed = 0`); the four
failed criteria from `inputs/Q021_RESULTS.md` remain failed.

### What was built on

- `PREREG.md` from the interrupted session (SHA-256 `8aa013b8…`, `sha256sum -c` → `PREREG.md: OK`
  both before and after the work). Freeze underway: no dataset was opened before the freeze.
- Helper `_dc2.py` (DataCite query) from the interrupted session, plus new `_fetch.py` that saves
  every publisher response under `provenance/` (no. 1–5 in PREREG F5).
- `inputs/Q021_model.py:57` parameter table as mapping vocabulary, `inputs/Q021_RESULTS.md` as
  **frozen** comparators (0.321475 / 0.540052 / 6.667433 / 2.441848 mmHg), Hall 2014 0.690/0.784.
- `samples/perfusion_results.tab` + `samples/data_dictionary.txt` (554 870 + 1 275 B).

### Key numbers with file

| Quantity | Value | File |
|---|---:|---|
| Samples total | 556 145 B (ceiling 52 428 800 B) | `results.json:sample.files`, `results_sample.json:C1_size` |
| Records / columns / region columns | 1 976 / 46 / 42 | `results_sample.json:parse` |
| Patients / measurements | 19 / 4 per patient (76) | `results_sample.json:parse` |
| Units from the file | MTT `[s]`, PBF `[ml/100ml/s]`, PBV `[ml/100ml]`, volume `[mm^3]` | `results_sample.json:units_from_file` |
| Coordinate frame | 12-region anatomical lung index (left/right × upper/middle/lower third × front/back), **read from the file's own dictionary**; body position is **not** stated | `results_sample.json:coordinate_frame` |
| Measured upper/lower flow ratio, room air (n=38) | median **1.0528** (min 0.4461, max 1.9694, sd 0.3876) | `results.json:c5` |
| Same, in the Hall unit ml·min⁻¹·ml⁻¹ (×0.6) | median **0.6317** | `results.json:c5` |
| Model's non-dependent/dependent, rest / heavy load | 0.321475 / 0.540052 | `inputs/Q021_RESULTS.md` |
| Hall 2014, rest / exercise | 0.690 / 0.784 | `inputs/Q021_RESULTS.md` |
| Flow ratio O₂/air (hypoxic response, n=38) | median 0.5753 | `results.json:c5` |

C1–C4 pass (`results.json:frozen_checks`). C5 is reported **regardless of outcome**:
the model's 0.321475 lies below the measurement's minimum 0.4461, and Hall 2014's 0.690 lies
*within* the interval — i.e. the measured ratio is consistent with Hall and **not** with the model's rest output.
No parameter was touched (PREREG F3/K3).

### What failed / counterchecks

- **K2 unit trap triggered.** Dictionary states PBF in `[ml/100ml/s]`, but the file's numbers give ~879 L/min
  for the six merged thirds — ~2 orders of magnitude wrong. Absolute flows are therefore
  not used at all; only dimensionless ratios are compared. Declared format `text/tab-separated-values` is
  also wrong: the file has 0 tab characters (comma). Publisher-declared size 562 063 B ≠ 554 870 B
  transferred, but MD5 matches the publisher exactly — reported, not hidden.
- **G4 lacks an admissible dataset.** Best candidate, PhysioNet BOLD (`10.13026/phvt-3277`), is
  *Credentialed Access* → E1 failed, listed with the access condition written out.
- **G5 lacks an admissible dataset.** Zenodo `14049935` (26 COPD, exercise protocol) has an empty
  file list in the API and `/files` gives HTTP 403 → E4 failed.
- **G7 best candidate unavailable.** Ottawa CT/SPECT/lobar NIfTI (`10.5281/zenodo.15650888`,
  100 studies) has an empty file list; size could not be verified → E1/E4 failed, marked
  `single_source` according to K1.
- **E2 held strictly.** Only one animal dataset was found (deer, `10.6084/m9.figshare.13348439.v2`):
  explicitly rejected as a human source for `hpv_gain`.
- **F2 no whole-lung averages as "regional".** DS2/DS3 are over-covering as *lobar*/EIT-regional.
- S1: one dataset was loaded (2 files from the same record: measurement table + publisher dictionary that C3/C4
  require). If S1 is read as "one file", the dictionary is metadata, not a second sample.

### Next step

1. Retrieve DS4's file list via the Borealis/Ottawa repo (or DataCite-v2) — the only candidate with
   voxel frame that can set `n_regions`/`lung_height_m` against measurements.
2. G4: either apply for PhysioNet credentialing or find an open Hb/P50 cohort; until then
   `p50_mmHg`/`hill_n`/`hb_g_dl` are entirely unconstrained by data.
3. G3: the measured O₂/air ratio 0.5753 is a *hyperoxia* response, not hypoxia — it cannot
   calibrate `hpv_gain` without a hypoxia arm and a directed model response.
4. Body position is missing from the file: before any number is compared with Hall, "upper third" = non-dependent
   must be confirmed for these scans.

Reproduce: `python3 check_sample.py` (1 thread, no network) → `results_sample.json`.
