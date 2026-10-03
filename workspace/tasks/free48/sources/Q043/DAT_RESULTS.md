BT-DAT-Q043

## Result

Public data closes the **geometry, force-capacity and rate-coding** side of the Q043 model; it
does **not** close the electrode-geometry side. 9 records admitted (8 human, 1 animal flagged),
6 candidates rejected or unreachable. `DATA_SOURCES.json` holds the full fields; `results.json`
holds every number below. No Q043 parameter was re-estimated here.

## Sample downloaded and verified (cap 50 MB)

`samples/` — primary: **10.5281/zenodo.10069266** (Caillet et al., PLOS Comput Biol 2023, CC BY
4.0). The file is 80,250,770 B, so it was **not** fetched whole: HTTP Range (206) on the zip
central directory, then 5 members, **914,465 B transferred**, all CRC32-verified, sha256 per
file in `samples/manifest_10069266.json`. Secondary: `10.5281/zenodo.18922801`, two CSVs,
176,453 B.

| file | shape | units | frame |
|---|---|---|---|
| `2_ankle_tib_ant_(400fibres_plus_tendon).osim` | OpenSim 4.00 XML: 401 `Thelen2003Muscle` (`TibialisAnterior_0..399` + `tib_ant_r_tendon`), 6400 `MovingPathPoint`, 4 `Mesh` | SI: m, N, s, rad | subject-specific right ankle, OpenSim global frame (`tibia_r`, `talus_r`, `calcn_r`) |
| `tib_ant_tendon.obj` | 4952×3 vertices, 9904 triangles | mm in file; `scale_factors 0.001 0.001 0.001` → m | right limb (`_r`), global frame after scaling |
| `2_S15_Segmentation.nii.gz` | `dim = [3,512,464,528]`, `pixdim = 0.44922, 0.44922, 1.00002` | mm (`xyzt_units=2`) | scanner anatomical **RAS** (`qform_code=1`, quaternion `(0.001081, 0.002629, 0.999996)`, `qoffset = (187.368, 57.149, -362.484)` mm), uint16, 16 labels |
| `motorUnitResultsUpdate_Raiteri.csv` | 3840×6 (`pid, decomp, intensity, condition, muid, outcome`) | Hz (median MU discharge rate) | no spatial frame; index = pid × decomp × intensity × condition × muid |

## What the data says about the model

Read from the downloaded files: `n_mus` **400** (not 120); `total_peak_force_n` **905 N** (not
100); activation/deactivation time constants **0.01 / 0.04 s** (model: 0.030–0.090 s);
tibialis anterior + tendon bbox **67.18 × 54.90 × 180.81 mm** (model box 60 × 70 mm);
optimal fibre length 0.098 m, tendon slack 0.223 m, pennation 0.0873 rad = 5.0°. In the 17
participants of DS2, 9–37 motor units resolved per participant × intensity (median 21.5); median
discharge **12.9 Hz**, p05 **9.2 Hz**, max **23.5 Hz** (11.9 Hz at 20% MVT → 14.4 Hz at 40% MVT),
sd 2.631 Hz. So `min_rate_hz = 8` is plausible; `first_peak_rate_hz = 45` is **not** tested —
the cohort stops at 40% MVT.

## Falls / open

Still assumptions, unconstrained by any admitted dataset: `geometry_length_mm = 12`,
`action_potential_duration_ms = 5/20`, `first_peak_rate_hz = 45`, `peak_rate_step_hz = 10`,
`isi_cv = 0.2`, `force_noise_fraction = 0.02`, `emg_event_amplitude_mv = 0.035` (DS2 states no
EMG amplitude unit). Participant counts for DS3/DS4 are `UNKNOWN` from the records. Rejected:
Zenodo 10936952 (simulated HDEMG, not measured), 5700435 "Data of AdHere" (advertising dataset,
title collision), 6953924 (network weights). Unreachable from this sandbox: NinaPro DB1–5 (DNS),
OpenCap LabValidation. Next: fetch DS8 (FHS, 250 subject-specific OpenSim models, 39.7 MB) for
the geometry cohort, and DS5/DS6 (CC0, 21 and 20 participants, EMG + torque + PCSA) as the
held-out test of the rank-1 EMG/force assumption. No answers invented; no measured value derived
from the synthetic run.
