# Video analysis, addendum: what V0 missed (V0b)

V0b, 2026-09-23. Addendum to `VIDEO_ANALYSIS_OVERVIEW.md` (V0) and `MECHANISM_TRANSFER_INVENTORY.md` (O5). Everything was read without writing anything to the sources. Nothing is mounted and no pipelines have been run.

The material is in `results/V0b/`:

- `census.json`: word census, 94 search terms (Swedish, English and Anton’s own) against 98 directories, with hits per directory and subdirectory and counts in fold ledgers.
- `clusters.json`: 22 clusters with paths, numbers and relevance, plus the five innovation candidates in machine-readable form.

## Conclusion first

V0 searched BodyTwin, mechanism and the data volumes, but almost none of the cad-to-simulation trees. Most of the video work is there. Three findings change the situation for the ongoing agents:

1. **V3 has a ready counterpart.** `cad-to-simulation-J` has a physics plausibility certificate for COM: the flight phase must be ballistic and the contact phase must have a realisable GRF within the friction cone. The certificate is validated on CMU MoCap and has an adapter from SMPL-X joints. It is not used in BodyTwin. (N2)
2. **Gravity gives both metric scale and time base.** The L ledger has a ballistic depth anchor, Z = f·g/|a_bild|. A QC shows that monocular SMPL-X scale is a nullspace that a body-height prior cannot break. No one has yet used gravity on the human’s own flight phase. (N3)
3. **V0’s claim that kinetics JSON is missing for the C4 corpus is incorrect.** There are measurements of force-plate jumps, g from barbell drops, 74 IWF lift attempts with weight in kg, radar and audio contact. (N7)

## State of the volumes

| Volume | State 2026-09-23 | Consequence |
|---|---|---|
| `external_media Volume` (sdb2, where `~/datasets` points) | FUSE interruption (*Transport endpoint is not connected*). It is still **not inventoried**. | EuRoC (`~/datasets/euroc/V1_01_easy`) is here and, according to memory, the 85 GB render_match_benchmarks, patent corpus and z24. Memory’s "/dev/sda2" is the same disk; the device name has changed. |
| `external_mount` | FUSE interruption | – |
| `external_mount` (`linux_data.img` on nvme1n1p3) | Not mounted. Contents reconstructed from references in docs and memory. | References include `corpus/football_pose_v1`, `v2cad_synth/{gt,handheld}`, `v2cad_labels`, `twin_baseline/venv_probe`, `scratch/kitchen_splat_v0`, `datasets/{camera_sensor_ptc, humanoid_hand_grasp_assets, lens_design_references, illumination_spectra, textile_weave_reflectance, display_projection_optics, panel_teardown_evidence}` and `feem_sfm`. |
| `external_mount`, `external_mount` | Mounted | Contain only dental videos and entertainment films. |
| `sde1` (vfat) | Mounted | Contains 3D-printer data and no video. |

## Kluster

"Missed by V0" means the cluster is absent from `results/V0/inventory.json`. "Partial" means V0 mentioned the path but not the contents, or judged the cluster irrelevant. Full paths are in `clusters.json`. CS-I stands for `local_path`.

| # | Cluster | Path (main) | What | Scope | Status | Key numbers | Relevance (BodyTwin · V3 · V4 · FV1) | Missed by V0 |
|---|---|---|---|---|---|---|---|---|
| N1 | Video modality inventory and corpus index | CS-I `docs/VIDEO_MODALITY_INVENTORY_V1.md`, `data/video_index_tables/` | 33 modalities with existing cells, a corpus index in three layers and `corpus_query.py` | 476 cells; 557 videos in the manifest; semantics for 261 body_twin clips in 47 groups | Active | Actual gaps: multi-video co-registration and event camera. Container fps is at most 60 in the entire body_twin corpus. | The checklist can be used as a gate for the V4 report. "Slow motion" in the corpus is an export effect, so V3 needs fps truth. | Yes |
| N2 | Physics plausibility certificate for COM | `cad-to-simulation-J/src/cad_to_sim/{floating_base_plausibility, pose_to_biomechanics}.py`, `docs/BIOMECHANICS_CERT_V1.md` | Requires ballistic flight phase and GRF within the friction cone in the contact phase. Has an adapter for SMPL-X/CMU → segment COM (Winter). | 6 CMU jumps; running; jitter test | v1 ready since 2026-07-16, not imported | Residual 0,011–0,076. 0 % false enclosure at 8 mm jitter. Abstain gate at 3,5 g. | **V3 directly.** V4 gets a new measure: fraction of IMPLAUSIBLE windows video versus mocap. FV1 can use the abstain pattern. | Yes |
| N3 | Metric scale from gravity and the scale nullspace | `cad-to-simulation-L/data/FOLD_LEDGER_L.jsonl` (L387–L396), CS-I `reports/probes/barbell_drop_g.json`, git 6328484a37 | Z = f·g/\|a\| without baseline. Monocular SMPL-X scale is a first-moment nullspace that a body-height prior cannot break. | Synthetic, barbell drop and IWF | Partially verified | Z error 0,004 with fixed camera. A following camera diverges (error 0,43 at g_t = 0,3). Barbell g: 12,7 % error in 1 of 3 clips, IWF 0 of 8 windows. SMPL-X scale 0,806. | **V3 directly:** gives scale and clock. V4: an anchor must be decorrelated from the body prior. | Yes |
| N4 | SMPL-X certificate fit and observability σ_min | CS-I `reports/probes/smplx_certified_fit_v0.json`, `cad-to-simulation-C/scripts/physics_exp/spine_*motion_from_video*_C.py` | Certificate fit per joint. Observability = σ_min(∂2D/∂pose). Monocular projection has rank 2, and an occluded point kills observability for distal joints. | 55 frames of an upper-body clip; 6 theoretical theorems | MEASURED v0. QC: observability is not corroboration. | Visible joints 18,5 mm, poorly observed 212 mm, hip 0,20–0,30 m. Fisher trace remains 608,70 under 0,3 m bias. | **FV1 directly:** a formal model for occlusion. V4 gets a predicted error per joint. | Yes |
| N5 | Humanoid hand and grasp (HUM-*) | CS-I `docs/HUM_HAND_V0.md`, `GRASP_*`, `scripts/hand_leg/`; `3fold-motion-engine/.../grasp_wrench_cert.py` | Hand with 21 DoF and SMPL-X bone lengths, grasp synthesis v1–v5, 8 grasp classes from the factory’s manual operations, force-closure certificate and contact grounding through object motion | 30 scripts, around twenty reports | GREEN/PARTIAL per node | Grasp v4/v5 fail the double gate (penetration 4,2 mm). Contact grounding AUC is 1,0, but only synthetic. | BodyTwin: grip force as input to the MECHANISM hand. FV1: contact is a nullspace in appearance. | Yes |
| N6 | Disassembly, maker and wave2 video | `~/corpus_incoming/{disassembly,maker,wave2}`, CS-I `reports/probes/{nh35_direction_pair, eta_single_take, render_match_v03, hum_hand_video_retarget*}.json` | 14 YouTube clips: 4 watch servicing (NH35 disassembly and assembly as a pair, ETA 2824-2), 4 maker and 6 repair/cooking. Metric render match and hand retarget. | 1,77 GB | Part-unit-level sequence exists, but only from text/OCR. Sequence from pixels, grasp and force was not found. | ETA: order reverses between disassembly and assembly, ρ = −1,0 (p = 0,042). NH35 between two videos: ρ = −0,03. NH35 scale 0,107 mm. Hand retarget 69 → 50,2 mm RMS. | Broad. MediaPipe hand scale error is 1,21–1,56 per finger. | Yes |
| N7 | Kinetics from the body_twin corpus | CS-I `reports/probes/{forceplate_jump, barbell_drop_g, iwf_attempt_ledger, radar_speed_validation, audio_contact_kinematics}.json` | Flight and contact times, g closure, lift attempts with kg via OCR, radar against OCR, audio against contact | GT suite: 9 positive, 5 partial, 9 negative | 2026-07-20 | 11 flight phases with median 0,50 s. Jump height 0,31 m, contact time 250 ms. 74 IWF attempts. Audio contact residual 697 ms (negative). | V3: flight time is an independent COM check. V4: ready negative controls. | Yes (V0 wrote that it was missing) |
| N8 | Player identity: face and shirt number | `~/projects/index_extract_proto/`, `twin_capture_lane/{FACE_INDEX_HANDOFF, MENTAL_MODELING_PLAN}.md` | ArcFace and shirt-number OCR per track in World Cup segments (public professional sport), plan for belief state | 69 095 face detections, 9 507 track embeddings, 361 segments | Prototype run | 5185 shots with 26 labels | Identity over time; track seams at occlusion | Partial |
| N9 | Upscaling: football FSRCNNX, browser and IPTV | `twin_capture_lane/UPSCALER_HANDOFF.md`, `~/.var/app/dev.fredol.open-tv/config/mpv/`, `browser_video_enhancer/` | Fine-tuned real-time SR in open-tv and a paused browser extension | 251 clips, 13 754 image pairs | v1 and v2 delivered | +0,485 dB against the original shader | Low. SR before pose requires a hallucination certificate (D12). | Partial (V0: irrelevant) |
| N10 | Video to exact CAD / bodywork | CS-I `docs/{VIDEO_TO_EXACT_CAD_PLAN_V1, VAGEN_FRAN_VIDEO_TILL_GEOMETRI_V1, KAROSSGRIND_MOT_VIDEON_V1}.md` | SfM, cue fusion and CAD fit through analysis by synthesis, with a preregistered gate on the worst frame | Volvo 240 from diffusion video | Active 2026-09 | The silhouette identifies 4 of 14 parameters (20–40 mm). Shape parameters move < 0,5 mm. Articulation from video: NO-BUILD. | The same limitation applies to body shape: the silhouette gives proportions, not mass distribution. | Yes |
| N11 | Optics, camera twin, reflection and image formation | CS-I `data/CAMERA_MODEL_MOTO_EDGE50.json`, `scripts/physics_exp/{camera_pipeline_sim, inverse_optics_cell, lens_raytrace_cell, rs_mains_flicker_cell, reflection_tracks, polarization_leg}.py`, `docs/{INVERSE_OPTICS_CELL_V0, CAMERA_IN_LOOP_CALIB_VERDICT, RENDER_MATCH_HARNESS, DIFFUSION_HANDOFF}.md`, FEEM_* | Phone camera model, ray tracing and inverse optics for a 365 mm triplet, PSF in render match, camera-in-loop projector calibration, rolling shutter and mains flicker, MTF/bokeh, specular reflection/BRDF, polarisation and sensor noise | Cells: 21 for RS/flicker, 64 specular, 29 polarisation, 24 bokeh/PSF | Partially stale | Camera-in-loop 0,036 px (synthetic). Inverse optics: EFL 0,67 %, but R1 25 % (failed). DepthPro 6,5 % median error: may RANK but not MEASURE. | V4: σ_px from PSF and noise model. V3: rolling shutter in fast limbs. FV1: separate reflections through motion fields. | Partial (V0 C10 had only p20/graph3) |
| N12 | Video authenticity, time base and motion floor | CS-I `src/cad_to_sim/temporal_flow_coherence.py`, `scripts/video_index/fps_truth_*`; coordinator2 memory | Authenticity cascade with 10 legs on EuRoC. Interpolated frames passed until a residual leg was added. Motion floor: amplitude ε/π, frequency ceiling fps/2. | EuRoC; stream_eyes on 16 streams | MEASURED | Seam detection 88 % | **V3:** interpolated "slow motion" gives false accelerations. | Yes |
| N13 | Sensor-rich video (GPMF-IMU, GPS, shaker, DIC) | `8838…/multimodal_corpus/cs_video/metric_telemetry/`, CS-I `gpmf_vision_crossvalidation_cell.py` | GoPro with embedded IMU/GPS, checkerboard/ChArUco calibration and load cell | 10 GPMF, 7 checkerboard, 4 shaker, 4 DIC | Partially used | GPMF 10/10 verified | V3: the camera’s own acceleration can leak into a_COM. V4: calibration reference answers. | Partial (V0: not body-relevant) |
| N14 | YouTube provenance, moved and deleted | `twin_capture_lane/video_provenance.tsv`, `8838…/_neutralize_quarantine/`, `NEUTRALIZE_HANDOFF.md`, `~/research/{video_corpus_grounding_2026-07-15.md, storage_cleanup_20260922/}` | URL lists, quarantined metadata, deletions and moves | See table below | – | 556 URLs in 25 categories. 376 .info.json and 490 .vtt. 42 partial downloads (7 GB) deleted. | The entire body_twin corpus can be downloaded again. | Yes |
| N15 | MECHANISM hand and foot close-up (OpenSim) | `~/projects/bodytwin/docs/MECHANISM_HAND_*.md`, `MECHANISM_FOOT_CLOSEUP_POSE.md`, `data/msk_models/subject2_hand_complete.osim` | Hand with five fingers and annular pulleys. MediaPipe’s hand model is used on a bare foot. | 44 bodies, 257 muscles | HYPOTHESIS, awaiting QC | Annular pulleys retain the moment-arm sign in 13 of 13 cases. The moment arm agrees with finite differences in 12 of 12. | Receives grip force (N5). Explains why the ankle failed in C2: it has only 3 landmarks. | Partial |
| N16 | Own filming (personal data) | Feem (V0), planned `~/corpus_incoming/room_components/`, git d8ed99aba0 (operator_hand_v1), CS-I `data/HAND_*_V0.json`, apartment splat | Paths only. Contents not read. | – | – | – | Not used | Partial |
| N17 | Arena splat from broadcast and SfM/3DGS framework | CS-I `reports/probes/arena_flyaround_v2.json`, `docs/{GSPLAT_READY, GPU_RECON_READY}.md` | 3DGS arena built from a selected pan in the full match | 82 of 82 frames registered | MEASURED | 22,25 dB, 25 532 gaussians, reprojection 0,33 px | Ground plane and scale around a person; static occluding geometry for FV1 | Yes (V0 C9 only splat_capture) |
| N18 | Models and teachers on disk | CS-I `reports/probes/gpu_topdown_pose_bench.json`; SAM2.1 in `omniverse_extension/vision_models`; `content_platform/models/{smplx, dino…}`; YOLO11 | ViTPose-B, SAM2.1, DINO, CLIP, SMPL-X, UniRig and YOLO11. Licence rule PROPOSE/RANK/MEASURE. | – | On disk | ViTPose-B: 1,21 s per frame with 46 figures, yield 0,76 | V4: fourth 2D detector. FV1: SAM2 masks. | Yes |
| N19 | Industrial YouTube subtitle corpus | `robot_lab-spine/data/youtube_subtitles/` | .vtt and analysis JSON from 11 CNC/manufacturing channels | 2745 files, 152 MB | Stale | – | No direct relevance | Yes |
| N20 | Isaac Assist: scene_eyes, Isaac camera, teleop, GR00T and grasp poses | `omniverse_extension/{workspace, docs, scripts}` | Measurement eye, replicator camera, teleop and humanoid | Files with hits: scene_eyes 165, isaac camera 173, teleop 156 | Stale since June | GR00T directories are almost empty | **FV1:** a sim camera can render known occlusions | Yes |
| N21 | content_platform: SMPL-X → retopology → UniRig → retarget | `content_platform/{core, scripts, research, vendor/UniRig, models/smplx}` | Character pipeline for games | 132 files with SMPL-X, 54 with retarget | Abandoned | – | Source of the SMPL-X weights | Yes |
| N22 | Monocular 3D human reconstruction (SOTA review) | `~/projects/bodytwin/docs/MECHANISM_HUMAN_3D_RECONSTRUCTION.md`; `MECHANISM_DAVID_TWIN_BUILD.md` (personal data, path only) | Primary sources for TokenHMR, HMR2.0, WHAM, GART and others | – | Review complete. The recommended next step is not performed. | Pose PA-MPJPE 43,8/49,8 mm. Shape PVE 84,6/104 mm, so metric shape from monocular video is rejected. | V4: external baseline. BodyTwin: ceiling for body shape. | Yes |

## What existed and whether it can be recreated

| Corpus | Subject | Now | Code that used it | Reproducible |
|---|---|---|---|---|
| body_twin_collection (V0 C4) | Weightlifting, jumps, sprint, slow motion | Exists on 8838…/multimodal_corpus (282 videos). Neutralisation 2026-07-19 renamed it `movement_mechanics_corpus`, but today it is named `body_twin_collection`. | CS-I `video_index/*`, N7 cells | Yes: `video_provenance.tsv` has URL and title. Quarantine has .info.json and .vtt. |
| cs_video (V0 C13) | Optics, fluid, free fall, telemetry | Exists (274 videos) | N11, N13 | Yes: `video_provenance.tsv` |
| Partial downloads (.part/.ytdl/.aria2) | Mixed | 42 files (7 GB) deleted 2026-07-19 | – | Yes via URL, if the source remains |
| Risk quarantine | Violence and weapons | 11 in `_neutralize_quarantine/risky_videos/`. Handoff states 30 total. | – | – |
| Disassembly/maker/wave2 | Watch servicing, crafts, repair | Exists in `~/corpus_incoming` | N6 | Yes: MANIFEST.jsonl |
| World Cup 2026 (Internet Archive) | Football | Exists (V0 C3). `ia_wc2026_480p/download.list` has 70 entries. | twin_capture_lane | Yes |
| Industrial subtitles | CNC and manufacturing | Only text remains (N19) | robot_lab | Yes via video ID |
| Ingest queue (CameraBench, Kansas, DroneZaic) | Camera motion, bridge inspection | The stream-process-discard pipeline has never run on real videos. 2 entries queued, 1 skipped for licence reasons. | `scripts/video_ingest.py` | – |
| render_match_benchmarks and EuRoC | Rendering, VIO | Moved 2026-09-05 to New Volume, now disconnected | CS render match, C/L lanes EuRoC | Requires remounting the disk. Anton will do it. |
| Isaac Sim assets 5.0 | Sim assets | Deleted 2026-09-22, manifest exists | Isaac Assist | Yes: NVIDIA URL in README |
| MECHANISM foot close-ups | Foot | Never downloaded (YouTube bot-blocked) | `foot_closeup_pose.py` | – |

## Ordcensus i korthet

Each row shows the number of files with hits and the five directories with the most hits. Noise from common words (frame, depth, match, hand, sharp, reflect, gaussian) is described in `census.json`. The CS trees are worktrees with partly the same contents, so the numbers must not be summed across them.

| Word | Files | Clusters V0 did not search |
|---|---|---|
| humanoid | 1483 | 3fold-motion-engine 408, omniverse_extension 235, CS-I 190 |
| teleop | 427 | omniverse_extension 156, kimate_backup 73 |
| scene_eyes | 3594 | CS-I 1070, omniverse_extension 165 |
| disassembly/disassembly | 1473 | CS-I 272, robot_lab-spine 159, CS-J 77 |
| smpl | 731 | CS-I 135, content_platform 132 |
| ansikte/face | 204 | content_platform 82, CS-I 52 |
| euroc | 2218 | CS-I 354, CS-C 242, CS-G 120 |
| kameratvilling | 267 | CS-I 85 |
| rolling shutter | 655 | CS-I 152, CS-C 88 |
| brdf / polarisation / psf | 1689 / 3058 / 1799 | CS-I 456 / 544 / 322 |
| isaac camera | 235 | omniverse_extension 173 |
| youtube/yt-dlp | 7747 | cad_course_corpus 4457 (CAD courses), robot_lab-spine 327 |
| egocentric / HOI / EPIC / Ego4D | 48 / 13 / 0 / 0 | Thin. No egocentric corpus exists locally. |
| wham / gvhmr / 4dhumans / tram / vitpose / rtmpose | 69 / 50 / 32 / 34 / 55 / 38 | Mostly literature review (N22) and a benchmark (ViTPose, N18) |

Two counts in fold ledgers: CS-I has 976 rows, of which 104 concern hand or grasp and 54 humanoid or gait. The L ledger has 704 rows, of which 239 concern optics and 67 EuRoC.

## Innovation candidates made possible only by the new findings

The same candidates are in `clusters.json` → `innovation_candidates`.

**I1. Gravity as ruler: metric body scale from the flight phase** (N3, N2 and C1; for V3 and V4)
During flight the COM follows a parabola with known curvature g. For a fixed camera it gives scale Z = f·g/|a_bild| without a body-height prior, and thus segment lengths in mm.
- **First test:** project mocap COM in C1’s drop-jump trials through one camera’s `cameraIntrinsicsExtrinsics.pickle` and add 2–4 px noise. The flight phase is taken where the force plate shows < 20 N.
- **Baseline:** scale from SMPL height prior.
- **Countertest:** contact phase instead of flight phase must give large errors.
- **Falsified if** median scale error is > 5 % at 3 px noise, or if the method does not beat the height prior.
- **Cost:** CPU, no decoding.

**I2. Physics plausibility as a measure of the video chain** (N2 and C1; for V4 and V3)
Run `plausibility_residual`/`certify_motion` on video TRC and mocap TRC with the force plate’s contact mask. The fraction of windows where the ground would need to pull, or where GRF falls outside the friction cone, becomes a physics measure that needs no reference answer for joint angles.
- **Hypothesis:** video gives at least 15 percentage points more IMPLAUSIBLE windows than mocap. The residual correlates with V0’s HJC error with ρ > 0,3.
- **Countertest:** contact mask shifted 0,2 s.
- **Falsified if** the difference is < 5 percentage points or ρ < 0,1.

**I3. Observability σ_min predicts where the video errs** (N4 and C1; for V4 and FV1)
Compute σ_min(∂2D/∂pose) per joint from OpenCap’s actual setups with 2, 3 and 5 cameras. Compare against |video − mocap| across V0’s 3060 trial pages.
- **Hypothesis:** Spearman ρ > 0,3.
- **Baseline:** constant σ.
- **Countertest:** permuted setups.
- **Falsified if** ρ < 0,1.
- **Interpretation:** if σ_min increases with more cameras while measured error also increases (V0’s G4 anomaly), the error is dominated by detector bias rather than geometry. Then V4 should put the resource into the detector.

**I4. Gravity as clock: expose retimed "slow motion" before kinetics** (N3, N7, N12 and N1; for V3)
With known pixel scale, for example a 450 mm barbell plate or body height, estimated g ∝ (fps_sann/fps_deklarerad)².
- **Hypothesis:** real-time clips with verified fps give g within 15 %. The slow-motion groups deviate by a retiming factor (4, 8 or 16), and the deviation covaries with the interpolation residual leg.
- **Data:** forceplate_jump (11 flight phases), known_load_lifts and slow_motion.
- **Falsified if** spread in real-time clips is > 25 %, or if the deviation does not covary with the interpolation leg.
- **Cost:** decoding, so the job must go through the shared queue.

**I5. Assembly order from disassembly video, then grip force to the hand model** (N6, N5 and N15; broad, according to Anton’s addition)
The ETA clip shows at text level that order reverses exactly between disassembly and assembly (ρ = −1,0). The same test must now be done at pixel level.
- **Hypothesis:** the "part removed" order from image segmentation in the NH35 disassembly, reversed, agrees with "part inserted" in assembly with Kendall τ > 0,6.
- **Baseline:** text’s ρ = −0,03.
- **Countertest:** random segment order must give τ ≈ 0.
- **Falsified if** τ < 0,3.
- **Then:** render-match scale (0,107 mm) and contact grounding through object motion can give the grasps forces that the MECHANISM hand receives.

## Limitations

- New Volume, `external_mount` and `external_mount` are still unread. We make no absence claims about their contents.
- The census reads text files up to depth 5 and smaller than 1,5 MB. Videos and binaries are counted only at volume level.
- Personal data is marked with a path and not cited. This concerns `ig_downloads` and BodyTwin’s `scripts/msk/inherited_hand`, which runs against the IG corpus, as well as Feem, operator_hand, the apartment splat and `MECHANISM_DAVID_TWIN_BUILD.md`.
