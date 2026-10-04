# Video analysis: inventory and innovation hunt for BodyTwin (V0)

V0, 2026-09-22. Everything has been read without writing anything in the sources. Nothing has been mounted and no old waves have been run. The numbers are in `results/V0/`:

- `inventory.json`: clusters, paths, counts, GB, code, reports, status, gaps and audit notes.
- `candidates.json`: ranked candidates.
- `hjc_video_vs_mocap.json` and `smoke_test_KA.json`: the smoke test.

This document builds on `MECHANISM_TRANSFER_INVENTORY.md` (O5) and does not repeat its analysis of tracker, retarget, IK bridge, splats and camera twins.

Counts and sizes come from `find -maxdepth 6` per directory and `du` on selected directories. **Not accessible:**

- `external_media Volume` (sdb2): FUSE interruption, *Transport endpoint is not connected*. `~/datasets` points there.
- `external_mount`.
- `external_mount`: linux_data.img is not mounted and contains football_pose_v1.
- INTENSO sdb1: not connected.

We make no statements about the contents there.

## 1. Inventory per cluster

| Cluster | What | Source / license | Scope | Code | Reports (key numbers) | Status | BodyTwin relevance |
|---|---|---|---|---|---|---|---|
| **C1 OpenCap LabValidation** `external_media` | 10 people, 5 synced iPhones with intrinsics/extrinsics, mocap-TRC, video-TRC (HRNet/OpenPose ×2/3/5 cameras, with R/L_HJC and LSTM-augmented `*_study` markers in the same file), OpenSim IK/ID/JR/SO, force plates, EMG, height/mass | Public (SimTK 2385); code Apache-2.0; data license file missing locally | 1628 videos, 16,97 GB (23 GB total); 16 IK trials per person | `bodytwin/scripts/msk/markerless_vs_marker_validation.py`, `pose_to_opensim_ik.py`, ~40 `msk/*` on the mocap arm | Joint-angle RMSE median 4,65° (HRNet 2 cameras, n = 540 pairs), shuffle-null 16,3°; no comparison of force, events or HJC before V0 | Active | **Highest.** Only local source with synced video, markers and force plate |
| C1b video2kin `…/mechanism_data/video2kin` | Own Pose2Sim run (subject2 walking1) | Derived from C1 | 18 files | Pose2Sim-venv, `labVal_script.py` | "4,0–5,5° against mocap" appears only in a document; V0 found no JSON | Stale | Replaced by C1 as reference |
| C2 IG training clips `local_path` → `bodytwin/data/msk_pose` | MediaPipe → TRC → OpenSim IK, monocular | **Personal data** (public accounts, individuals); contents not quoted | 893 videos, 2,09 GB; 267 pose directories, 0,425 GB | `msk/pose_extract.py`, `run_pose_batch.py` | 228/264 pass the gates; marker RMS 0,085 m; ankle FAIL | Stale-active | Low for mm/rad; stress test for QC |
| C3 World Cup 2026 broadcasts `external_media`, `external_media`, `segmentation_dev` | 4K matches, SR wide shots, DB with 66 videos, 16556 shots, 1823 segments, 633 players | Public professional sport via Internet Archive, broadcaster's copyright | 310 + 1513 + 52 videos; 267 + 20,3 + 2,8 GB | `twin_capture_lane/pipeline_play_pull`, CS-I `video_index` (89 files), `football_tracking`, `pitch_awareness/PnLCalib` | Longest tracks 36/12/5 frames, 9–13 % ID switches (proxy), TPR 0,09 for the same team, height bias 1,2455, corner view 57 % error | Corpus exists; biomechanically negative | Low. Methods for identity, time base and pitch scale |
| C4 body_twin_collection `external_media` | YouTube: weightlifting with known load, force-plate tests, plyometrics, slow motion; 18 World Cup excerpts | Public videos, elite/professionals; license not checked | 282 videos, 14,14 GB, 70 transcripts | CS-I `iwf_*`, `powerlifting_ledger_cell.py`, `video_metric/disc_diameter.py` | No kinetics JSON found | Collected, unanalyzed | Medium. Known external load and scale anchor |
| C5 weigh-ins `…/mechanism_data/athlete_video_corpus` | Weigh-ins with known mass, one DEXA comparison | Public people in a professional role | 6 videos, 0,31 GB | `msk/body_composition_from_video.py` | RMSE 4,0 %BF and 2,7 kg on the regression; 5/541 frames usable | Stale | n = 6 is not enough |
| C6 CS video corpus `external_media` | 201 admitted segments: football 80, CameraBench 70, Kansas bridge 51; humancap_scratch, football pitches | Mixed public/open | 87 videos, 6,78 GB (22 GB total) | CS-I `human_capture_chain_v0/v1`, `hum_retarget` | Leg-length CV median 0,14 (19 clips); retarget RMS 0,179 m | Stale | Low-medium (admission, countertest) |
| C7 Mocap/force without video `mechanism_data/{knee_gc_measured, fukuchi_check, k2muse_full, openbiomechanics_*, weightlift_grf, …}`, `external_media[1-6].Data*Competition-latest.zip` | In-vivo knee force (Grand Challenge), running, K2MUSE, pitching/hitting etc. | Public, licenses vary | GC zip ~8,5 GB (no top-level videos; nested zips not opened), K2MUSE 7,6 GB | `msk/grand_challenge_multisubject_verify.py`, `validate_*_force.py`, `opensim_jam_build` | `gait_mechanics_results.json` exists (not read) | Active for the force chain | High for \|R\| reference, but reaches video only via C1 |
| C8 Body composition `mechanism_data/{DXX*, BMX*, sarcopenia, vat_test, totalseg, totalseg2}` | NHANES DXA/anthropometry, TotalSegmentator CT | Public | totalseg 76 MB + 299 MB | `body_composition_from_video.py`, `vat_test/cv_models.py` | See C5 | Partly used | Medium-high for mass/inertia prior |
| C9 Splats `~/projects/splat_capture`, `src/cad_to_sim/splat_*` | Phone→3DGS, ESDF, certificates | Own object capture + ETH3D | See O5 | See O5 | σ_min AUC 0,88 < scale 0,91; prospective 0,549 | Stale/negative | Low now (no metric scale, rigid scene) |
| C10 Camera models/metric video CS-I `src/video_metric`, `feem_*`, `p20_camera_*` | Scale anchors, plane pose, disc diameter, PSF/sensor | Own code; feem = own object filming (112 videos, 5,6 GB in `anton_flytt…/Downloads/Feem`) | – | See What column | p20 self-test 4/4; graf3 c9 OPEN_NEGATIVE | Stale | Medium. Building blocks for scale and σ |
| C11 Humanoid gait simulation CS-I `physics_exp/hum_gait_*.py` | Newton humanoid; momentum-GRF F = m(g + a_y) replaced a defective contact sensor | Own simulation | 10 scripts | `hum_gait_momentum_remeasure_v0.py` | Root velocity 0,147 m/s against the literature's 0,15–0,30 | Stale | The method can transfer to video-COM → GRF |
| C12 Own filming / personal | `sdc1-tmp/PERSONAL`, `Osorterat/Videos`, `8838…/Users`, `Foton-Xiaomi`, `lite random videos`, `~/captures` | **Personal data**; only path and volume | 109 / 51 / 130 / – / 21 / 3 videos | – | – | – | Not used. No own lab filming with a reference |
| C13 Not body-relevant | Courses (1848 videos, 134 GB), cad_course_corpus CAD tutorials (687, 19,4 GB), cs_video physics/optics (274, 21,7 GB), film, sharp_football | – | – | – | – | – | None |

## 2. Building blocks missing in BodyTwin

In `references/current_bodytwin`, there is no code for camera, video, pose or movement. In `~/projects/bodytwin/scripts/msk`, the following are missing:

1. **Gait events from video.** `gait_kinematics_fidelity.py` gets heel strike from GRF.
2. **Ground reaction from video-COM.** The momentum estimator exists only in the humanoid simulation.
3. **Measured distribution for joint-center error.** K1 tests assumed ±10 mm. Measured: see §4.
4. **Segment length with σ from video** to I ∝ m·L².
5. **Camera covariance propagated to 3D points.**
6. **Population σ for regional mass** (`region_mass_v1` has density → mass but no σ).

## 3. Candidates, ranked by what they unlock per cost

The full decomposition (idea → mechanism → equation → operation → representation → assumptions, route, adaptation, test) is in `results/V0/candidates.json`. A summary follows here.

| # | Candidate | Domain coupling | Equation / operation | Unlocks | First experiment: hypothesis · baseline · countertest · falsified if | Cost |
|---|---|---|---|---|---|---|
| 1 | **K-A: joint-center uncertainty → the collaborator's \|R\|** | Camera twin and markers → the collaborator chain (K1) | e = R_pelvis·(p_video − p_ref); ΔR/R = f_K1(e) | the collaborator's most sensitive leaf gets a measured error distribution instead of an assumption | Median \|ΔR/R\| > 10 % · marker regression through the same K1 · permute axes/people · < 5 % with an extended grid | Low. **Smoke-tested** |
| 2 | K-B: video-COM → vertical GRF | Humanoid simulation (momentum-GRF) → actual video → force plate | F_z = M(g + c̈_COM), c = Σ m_i c_i / M with de Leva fractions or `region_mass_v1` | Ground reaction without a force plate for jumping, squats and rising | MAE < 15 %BW in DJ/squat · F = Mg and the same estimator on mocap · time shift 0,2 s + wrong person's mass · does not beat Mg in DJ | Low |
| 3 | K-C: gait events from video | Sports tracking → gait events | Zeni: HS = argmax (heel − pelvis)·v̂; reference F_z > 20 N | %GC clock without a force plate; makes C2/C4 video analyzable per cycle | \|error\| < 30 ms HS, < 50 ms TO · Zeni on mocap · half-cycle shift · HS error > 50 ms | Low |
| 4 | K-D: segment length → mass/inertia with σ | Video QC (leg-length CV) → scaling | σ_I/I ≈ √((σ_m/m)² + (2σ_L/L)²) | Individual scaling in `uncertainty_v1` | \|bias\| < 10 mm, CV < 3 % · L = 0,245·H · permute people · does not beat height regression | Low |
| 5 | K-E: camera twin → measurement uncertainty | Triangulation Fisher (splat experiment) → marker σ | F = Σ(I − ddᵀ)/(σ_px² r²/f²); coverage P(\|e\| < 2σ) | σ per marker can be labeled `measured` | ρ(σ_pred, \|e\|) > 0,3 · constant σ · swap camera geometry between sessions · ρ < 0,1 or coverage < 50 % | Low-medium |
| 6 | K-F: regional DXA → mass prior | Body shape → mass/inertia | m_seg/M = β·x + ε | `region_mass_v1` with population σ | R² > 0,6 for leg mass/M · de Leva constant · permute covariates · R² < 0,2 | Low |
| 7 | K-G: disc diameter 450 mm → bar force | Metric video (`disc_diameter.py`) → inverse dynamics with known load | F_bar = m_bar(g + s·ÿ_px), s = 450/D_px | Monocular scale and a plausibility test with known external load | Impulse balance within 10 % in ≥ 5/8 clips · scale from the athlete's height · wrong diameter 400 mm · no difference between correct and wrong scale | Medium (decoding) |
| 8 | K-H: jersey/re-ID → identity over time | Football re-ID → gate in multi-person video | Hungarian on IoU, pose, velocity, number | Prevents mixing two people's data | ID switches < 2 % · current 9–13 % · mix numbers between tracks · no improvement | Medium-high |

The ranking is justified as follows:

- **K-A–K-D** use only text files in C1, which are already synced and have a force plate. Each candidate fills a missing node in the collaborator's chain.
- **K-E** requires interpreting OpenCap's pickle format.
- **K-F** has no video dependence. It remains because mass/inertia is its own chain step.
- **K-G and K-H** require video decoding and give less for individual biomechanics.

## 4. Smoke test K-A (`results/V0/hjc_video_vs_mocap.py` → `hjc_video_vs_mocap.json`, `smoke_test_KA.json`)

The run used CPU with one thread, at most 130 MB RSS and about 6 minutes. The source is untouched. The scope was 1530 TRC files, 3060 trial sides and 0 read errors. The code runs.

| Measure (median across trial sides) | All configurations | HRNet, 2 cameras |
|---|---|---|
| 3D-RMS video-HJC against mocap-R/L_HJC | 39,8 mm | 41,1 mm |
| Systematic \|error\| ant/sup/lat | 19,2 / 11,8 / 13,1 mm | 19,9 / 13,7 / 14,0 mm |
| Signed (HRNet 2 cameras) | – | −19,9 ant / +9,1 sup / −11,5 lat mm (video lies posteriorly and medially) |
| Jitter-SD within trial | 12,7 / 13,8 / 7,0 mm | – |
| Offset removed (OpenCap's `_offsetRemoved`) | 34,7 mm | 36,5 mm |
| **Floor: marker regression R_HJC_reg against R_HJC** | **27,8 mm** (lat −17,7 mm) | same |
| Video against R_HJC_reg | 46,7 mm | 45,9 mm |
| Fraction outside K1's ±10 mm on at least one axis | 96 % | 96 % |
| Best configuration | OpenPose_highAccuracy with 5 cameras, 35,6 mm | – |

Interpretation:

- **Video-HJC is not sufficient as input to the collaborator's chain.** The systematic error is about twice K1's largest tested shift. According to K1, −10 mm anteriorly already gives +22 % \|R\|.
- **The K1 translation is not interpretable yet.** The clipped sum, about 0,116, consists of +0,220 and −0,126 that partly cancel. K1 L1 must be rerun with ±40 mm.
- **The more important discovery concerns the reference itself.** Two marker definitions of HJC in the same data differ by 27,8 mm RMS. The choice of HJC method is therefore itself larger than K1's entire test interval. A claim of a joint center accurate to mm must therefore state which definition it concerns.
- **One calibration per person is not enough.** Individual bias varies less between people (SD 8,9/11,4/5,9 mm) than within the same person across trials and configurations (12,2/14,7/14,7 mm).
- **The time shift of 0,5 s is a weak countertest.** It gave 634,6 mm, but global translation dominates. It does not carry the conclusion.
- **Next steps:**
  1. Extend the K1 grid to ±40 mm.
  2. Run K1 with the HJC distributions (video, regression and reference) as Monte Carlo.
  3. Add a functional HJC (sphere fitting on thigh markers across movement) as a third marker method.

## 5. Audit notes

- `bodytwin/reports/probes/markerless_vs_marker_validation.json`, G4, passes the gate only thanks to 5 % tolerance. Median RMSE **increases** from 4,651 to 4,807 and 4,971° with 2, 3 and 5 cameras.
- The number "4,0–5,5° against mocap" for video2kin in `MECHANISM_MOVEMENT_BASELINE_CORPUS.md` has no JSON source that we have found.
- The Grand Challenge zip files are not unpacked. Whether the nested *Raw Motion Data* contain video has not been checked. If they do, it is the only source with video and in-vivo joint force.
