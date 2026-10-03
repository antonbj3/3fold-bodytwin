BT-DAT-Q044

**Builds on.** Interrupted session in this folder (`AGENT_EXIT`=1) that only had `BRIEF.md`,
`ALLOW_WEB`, `agent.log` and `inputs/`. The model being constrained is `inputs/Q044_model.py` (five
concentric layers, 21 electrodes on a sphere, average/mastoid reference, RDM), and its run
BT-HX-Q044 in `inputs/Q044_RESULTS.md` (topography norm 1,1449×10⁻⁹ V; skull conductivity ±50 %
→ RDM 0,0939/0,0604; skull thickness ±50 % → 0,0254/0,1103; reference choice → RDM 0,9387). The night's
requirement §1 (grep in `~/projects/3fold-workspaces/bodytwin`) could not be run: the harness denies
`external_directory`, so `notes/`, `docs/DATASET_MAP.md` and other `BT-DAT-*` are unknown here.
Everything below is our own work in this folder. `PREREG.md` was frozen before download, digest
`54366bed…` is in `PREREG.sha256` and is identical in `results.json`.

**Public sources that constrain the model** (complete table with URL, license, size,
format, subjects, verification: `DATA_SOURCES.json`)

| id | license | size | subjects | constrains |
|---|---|---|---|---|
| CHB-MIT v1.0.0 (PhysioNet) | ODC-Attribution v1.0 | 42,6 GB, 664 EDF | 22 patients, 1,5–22 years | electrode vector/reference, potential scale, time series |
| OpenNeuro ds007358 | CC0 | 26 005 files, 16,08 GB | 2 000 participants, 2 sites | electrode montage, reference, cohort spread |
| TotalSegmentator v2 (Zenodo) | CC BY 4.0 | 23,58 GB, 1 file | 1 228 CT | skull thickness, skull radius, nonsphericity |
| IXI | CC BY-SA 3.0 | 27,37 GB (archive) | ~600 healthy, 3 hospitals | r_wm, r_gm, CSF/skull/head boundary |
| BrainWeb (McGill) | UNKNOWN (no license text on the page) | 7 109 137 B sample | 0 (phantom) | all five layer boundaries, r_scalp |
| CamCAN | CC BY 4.0 | UNKNOWN | ~N=700, 18–87 years | anatomy + sensor position per age |

**Downloaded samples** (`samples/`, all ≤ 50 MB; `sha256` and extraction in `results.json:samples`)
`chb01_03.edf` 42 399 744 B (= remote `Content-Length`, placebo 2), `brainweb_normal_crisp.raw`
7 109 137 B (= 181·217·1, exact), `ds007358_participants.tsv` 89 290 B, plus two small
CHB-MIT files. Form: EDF 23 channels × 256 samples per 1 s record, int16 → physical **µV**,
256 Hz, 3 600 records = 1 h; **no coordinate frame exists in the file** (no x/y/z, no fiducials).
The phantom: uint8, (181, 217, 181) in z/y/x, 1 mm isotropic, **mm**, origin x −90, y −126, z −72 mm.
Placebo 1 (two independent parses) = identical.

**What the samples say about the model** (all numbers in `results.json`)
- Geometry is not a question of a better sphere: in the phantom, the white matter surface is at median 67 mm
  (IQR 60–75) against the model's `r_wm` 50 mm, and the head surface at ≥87 mm (censored by the volume,
  |x|≤90) against `r_scalp` 95 mm. The model's scalpel therefore does not work as the "head's" sphere.
- `skull_thickness` 7 mm against 3 mm median (IQR 2,75–4,25, max 7) in the phantom: the right order of magnitude
  is half a factor, and a single value hides an interval.
- The channel set is 23 bipolar pairs (FT9/FT10, **no M1/M2**). The model's
  `reference="mastoid"` therefore cannot be built from this montage, and a referential u
  requires an explicit montage→reference transform first; the label `T8-P8` is also repeated.
- Scale: measured AC-RMS 29,9 µV (median per channel) and AC peak 191 µV against the model's 1,1449×10⁻⁹ V
  → 1,67×10⁵ times. The file is a seizure file (listed in `RECORDS-WITH-SEIZURES`), so this is
  no normal interictal reference. RDM is scale-free, so this does not affect RDM but all
  absolute outputs.
- `ds007358`: the `age` column exists but is "n/a" in all 2 000 rows → age-dependent geometry cannot
  be constrained from that file.

**What failed.** BrainWeb has no license text found → the field is UNKNOWN, not assumed. IXI
returned HTTP 403 from this machine, so no IXI sample was downloaded. TotalSegmentator is a single
23,6 GB file, so no CT sample. `age` in ds007358 was unusable. My own first RMS computation was
DC-contaminated (the channels are at ≈ −800 µV) and reported 803 µV; it is corrected to
AC statistics and the error is mentioned here rather than hidden.

**Gaps.** (1) No public measured conductivity data was found — especially not
skull conductivity, the model's own dominant axis; MR-EPT/EPG publishes papers, not verified
open endpoints. `sigma_skull`=0,010 S/m remains a literature anchor. (2) No electrode impedance
time series. (3) Skull thickness rests on a simulated sample plus an uncited CT cohort.

**Next step.** Fetch per-subject CT from a cohort with an open file server (TotalSegmentator
v2 or an openly shared head CT) and measure the skull band's thickness distribution per scan; this gives
the single axis the model itself points out. In parallel: build a 10-10 realization cohort from
CHB-MIT/ds007358 with explicit montage→average-reference transform and test whether the sphere gives
RDM > 0,05 against it. No measured data has been used here to validate the model.
