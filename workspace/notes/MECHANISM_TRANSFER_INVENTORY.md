# Mechanism → BodyTwin: transfer inventory (motion, video, splats, camera twin)

O5, 2026-09-22. Read-only in the source projects. No source graphs written, no old waves run. Two existing CPU self-tests were run (see §5). The hash comparison is in `results/O5/sha256_mechanism_vs_bodytwin.tsv`.

## 1. Hash status: mechanism compared with bodytwin

Urval: filer under `scripts/ src/ docs/ reports/probes/` vars namn matchar splat, camera, pose, football, motion, track, video, human_capture, retarget, optic, pitch_calib, smpl eller anthropo.

| Utfall | Antal | Kommentar |
|---|---|---|
| Byte-identical | 715 | All inputs in the assignment, including `football_tracking`, `video_index/*`, `p20_camera_twin.py`, `splat_*`, `l_splat_*`/`l_certified_splat*` with JSON, `human_capture_chain_v0/v1`, `hum_retarget`, `pitch_calibration.py` and `temporal_split_cert.py`. |
| Different | 3 | `scripts/mechanism_video_analyze.py`: mechanism has a newer version (1268 lines) that adds the pitch anchor; bodytwin has 854 lines. `scripts/msk/pose_to_opensim_ik.py`: bodytwin has the newer version (811 lines) with flags for occlusion asymmetry, low frame coverage and the model's joint limits; mechanism has 577 lines. `scripts/msk/body_composition_from_video.py`: bodytwin adds a `reference_body` block that names a real private individual. That block must not be carried forward. |
| Only in bodytwin | 77 | Mainly `scripts/msk/bonemet_*`, `scripts/cad/*` and two reports from `mechanism_video_analyze_*`. |
| Only in mechanism | 0 in the selection | `cognition/probe_pitch_anchor.py` (sha `e62d63ef…`) is outside the selection and exists only in mechanism. |
| `data/MECHANISM_ANCHOR_GRAPH.json` | different | mechanism `01c5d753…`, bodytwin `2c0c4b17…`. The difference has not been analysed and no graph is to be overwritten. |

Conclusion: almost all code for motion, video and splats exists in a single version. The evidence may therefore only be imported once, with the mechanism path as provenance. The three files with different content require a choice per file (see §4). A file's presence in cad-to-simulation-I/-L says nothing about its presence in bodytwin: `camera_twin_v1.py` and `FOOTBALL_TRACKS_ADJUDICATION_L381.md` are missing in both.

## 2. What the code actually does (read in the code)

### 2.1 Motion and football

- **`scripts/football_tracking/build_football_tracking_v0.py`** (`02b082e8…`)
  - **In:** `/mnt/data_root/corpus/football_pose_v1/*.json`, med COCO-17-key points in pixels and bbox per frame. The installation is missing today, so the script cannot be rerun.
  - **Method:** Hungarian matching (`scipy.linear_sum_assignment`) on three terms: bbox-IoU, pose distance normalised by the bbox diagonal and centroid velocity per frame. The terms are scaled so that the median for positive pairs becomes 0 and the median for negative pairs becomes 1. The scaling is calibrated against a track with seven transitions that a greedy tracker itself produced. The calibration is thus **not independent ground truth**. The gate is 0,5. `MAX_GAP=20` frames.
  - **Ut:** `reports/probes/football_tracking_v0.json`.
  - **Units and frames:** only pixels in the image. No camera motion or metric scale. `dt` is the difference in frame index, not seconds.
  - **Result:** longest track 36/12/5 (fg_open/fg_corner/lg) versus 8 for retarget v0. Matching rate 0,84–0,92.
  - **Counter-test:** shuffled order collapses the track length in two clips out of three. In `lg`, 1 of 12 shuffles beats the normal longest track, so the gate fails there.
  - **Proxy for ID switches:** 9–13 % flagged transitions.
  - **L's `temporal_split_cert`:** ABSTAIN on all tracks, because the tracks are shorter than 60 frames.
- **`scripts/hum_retarget/build_hum_retarget_v0.py`** and `docs/HUM_RETARGET_V0.md`
  - **Method:** COCO-17 in 2D → simplified capsule rig with SMPL-X topology → Newton/MuJoCo with PD servo. Ten SMPL-X joints are unobserved, including pelvis, back, neck and feet.
  - **Result:** 8 frames (267 ms at 30 fps). Kinematic RMS 0,179 m. Shuffled order gives **identical** RMS, so the measure lacks temporal information. Root tracking RMS in physics is 0,362 m versus 0,367 m with shuffled order, which is marginal. Height proxy 1,495 m, implausibly low for a professional player. The player is 28×69 px.
- **`scripts/human_capture_chain_v0.py` / `v1.py`**, med rapporterna `human_capture_v0/v1.json`
  - **Method:** MediaPipe BlazePose (Apache-2.0, CPU) for 3D landmarks in hip-centred metres.
  - **v0:** bone-length CV has median **0,14** and p90 0,31 in 19 clips. Visibility and relative error have negative rank correlation (ρ −0,22) in 89 % of the clips.
  - **v1:** association with IoU reduces CV from 0,106 to 0,059 in one multiperson clip, but gives no improvement in another (0,074 versus 0,0745). Only three clips are truly multiperson.
  - **Global scale** via foot point on the pitch and pitch homography: synthetically exact with true homography and 0,07–0,09 % error with B2 line fitting. **The corner camera gives 57 % error at confidence 0,776**, and the confidence did not warn. Without a known metric pitch in the image, scale is missing.
- **`cognition/probe_pitch_anchor.py`** (only mechanism) and `scripts/mechanism_video_analyze.py` (newer in mechanism)
  - **Method:** the pitch (105×68 m) with PnLCalib gives the player's height from ankle to nose.
  - **Result:** after two bug fixes, the raw data has standard deviation 0,10 m but a systematic bias. The bias is corrected at population level against Wikipedia heights for 51 players, with the factor **1,2455** (CI 1,221–1,269). This is an ensemble estimator, not an individual measure.
  - **Bodytwin's runs of `mechanism_video_analyze`:** both clips selected `learned_body_size_prior` and have `chosen_is_true_measurement=false`.
- **`scripts/msk/pose_to_opensim_ik.py`** (nyare i bodytwin)
  - **Method:** MediaPipe-`.trc` → OpenSim IK on the scaled LaiArnold2017 model for subject2, thus **another person's** anatomy. No per-person scaling and no global translation, because the landmarks are recentred on the hip in each frame.
  - **Result:** marker-RMS **0,085 m**, 6,5× the floor for the same person with real markers (0,013 m).
  - **Gates:** knee and hip PASS, ankle FAIL. For hip flexion, the PASS ceiling coincides with the model's mechanical range, so a PASS can mean the solver was locked against the limit. 228/264 corpus clips PASS.
- **`scripts/msk/body_composition_from_video.py`**
  - **Method:** silhouette width → ellipse with assumed depth/width ratio *k* → circumference → regression for %BF and mass.
  - **Validation:** against 252 men with hydrostatic weighing. Test-RMSE 4,0 %BF and 5,97 lb (≈2,7 kg). Mass is 2,7× more sensitive to the assumed height value than %BF.
  - **Application to real video:** 5 of 541 frames usable.
- **`scripts/video_index/football_scene_verifier.py`** (`4bca1d09…`)
  - **Method:** a second signal for clip boundaries, SSIM on greyscale at 6 Hz, confirms, moves or rejects candidates from the Tier-0 histogram. Crossfades are sought via accumulated drift. Decoding via ffmpeg pipe; `cv2` gave 0 frames on AV1.
  - **Evidence:** `I-FOOTBALL-SCENES`. 16/16 samples agreed on visual inspection, in 22 clips and 48 scenes. n is small and the reviewer is the same agent.
  - **Use in BodyTwin:** preprocessing (segmentation into scenes), not biomechanics.
- **`scripts/video_index/fps_truth_probe.py`** with the tables `data/video_index_tables/fps_truth`
  - **Method:** determines the real recording rate from at least two independent declarations: title, ASR, OCR and a manually read sign. The container's fps explicitly does not vote, because slow-motion material is retimed.
  - **Use in BodyTwin:** this is the only existing check of the **time base**. Velocities and angular velocities depend directly on it.
- **`cad-to-simulation-L/docs/FOOTBALL_TRACKS_ADJUDICATION_L381.md`** (`927d0a63…`)
  - **Content:** 29/31 segments flagged for an ID switch within the same team are admitted anyway. The appearance leg has TPR 0,09 at 5 % FPR for the same shirt.
  - **Decision:** ADMIT, with the gap disclosed. Source manifest L317 is not in the tree (lost scratchpad material).
  - **Consequence for BodyTwin:** identity over time is not certified for people who look alike.

### 2.2 Splats

- **`src/cad_to_sim/splat_world.py`** (`277e9e7c…`)
  - **Method:** INRIA-3DGS `.ply` → voxel-ESDF for collision tests in robot planning. Floaters are removed with an opacity gate. Occupancy is either isotropic with dilation or as anisotropic Mahalanobis ellipsoids.
  - **σ measure:** the eikonal residual `|‖∇sd‖−1|` has only rank correlation 0,37 against real ESDF error and is not a metre-calibrated σ. Coverage-σ is `1/(1+support)`.
  - **Task:** free or occupied space. Neither dimensions nor material.
- **`src/cad_to_sim/splat_quality_cert.py`** (`3b215fbe…`)
  - **Method:** GT-free certificate per Gaussian with the classes needle, ghost and giant. The tail is computed as median + k·MAD, plus absolute limits. Input contract: scale after `exp()` and opacity after `sigmoid()`, otherwise ABSTAIN.
  - **What it measures:** the health of the representation, not geometric accuracy.
- **`selftest_splat_world.py`**
  - **What it measures:** conversion errors (false free or false occupied) against the raw Gaussians. Capture errors are not measured. The reference capture from splat_capture has `gate2=False` (failed run).
  - **Not run:** the script writes to `mechanism/reports/splat_world.json`.
- **`scripts/local/l_splat_*` and `l_certified_splat*`** (ETH3D pipes, 14 views, 2473 SfM points, gsplat on GPU)
  - **Idea:** triangulation Fisher `F=Σ(I−ddᵀ)` and σ_min as a certificate for floaters.
  - **Results that hold:**
    - The direction of σ_min coincides with the line of sight (0,999).
    - Reprojection error and σ_min are uncorrelated (−0,035).
    - AUC against distance to surface is 0,88. The scale baseline is **better**, 0,91, and the combination 0,93.
  - **Results, weak or negative:**
    - Prospective AUC 0,549, labelled as "breakthrough".
    - Cropping based on the fence **worsens** held-out L1 for every fraction; best is not to crop at all.
  - **Metric scale:** missing. The COLMAP scale is arbitrary.
- **`~/projects/splat_capture`** (8,6 GB)
  - **Content:** phone capture → glomap SfM → 3DGS. The gates are photometric: PSNR ≥ 27, SSIM ≥ 0,85, LPIPS ≤ 0,25.
  - **Outcome:** the latest runs in `pipeline_runs.jsonl` have status failed or error. `scene_normalized.ply` is scale-normalised, thus without metric scale.

### 2.3 Kameratvillingar

The camera twins are two different things with the same name. Neither is a measurement-camera model.

- **`scripts/physics_exp/p20_camera_twin.py`** (`0435aed9…`)
  - **Content:** photometric and optical forward model. Airy-PSF via FFT of the pupil, Zernike defocus, sensor chain (QE, dark current, Poisson, full well, read noise, gain) and photon transfer curve.
  - **Missing:** projection, intrinsics/extrinsics, distortion, rolling shutter and timestamps.
  - **The calibration is circular:** it recovers parameters in its own simulated sensor. The "inversion" is the Maréchal formula applied backwards to the same simulated PSF.
  - **The PSF gate** compares 5 px against 4,88 px with a resolution of 1 px.
  - **Label error:** the line `gain … DN/e⁻` should be e⁻/DN, as the summary states.
- **`cad-to-simulation-I/reports/probes/graf3_kameratvilling_v1.json`** with `scripts/graf/camera_twin_v1.py` (152 lines)
  - **Content:** rendering camera in WebGPU. FOV from sensor and focal length (0,67 px error), CoC for thin lens (GPU against analytical < 0,05 px), exposure adaptation and cost of DoF.
  - **Purpose:** image synthesis, not inverting real images.
  - **Missing:** distortion, timestamps and uncertainty.
  - **c9 is OPEN_NEGATIVE:** the focus distance is less than the focal length.

### 2.4 Other directories

- **`~/projects/sharp_football`:** odds model for betting, Dixon–Coles/Poisson. **No motion or image data. Not relevant.**
- **`~/projects/twin_capture_lane`:** pipeline that fetches and segments World Cup 2026 matches in 4K from Internet Archive, with SQLite DB, face index and an FSRCNNX upscaler trained on football.
  - **Value:** corpus and temporal levels.
  - **Warning:** frames after upscaling are not measurement data. Licences and personal data require checking before use in BodyTwin.

## 3. The chain raw data → quantity → computation, and the gaps

| Link | What exists | Quantified error or status | What is missing for mm or radian level |
|---|---|---|---|
| Time base | `fps_truth_probe` (declaration voting), ffprobe duration, `dt` = frame index in the tracker | Timestamps per frame are missing. VFR is flagged but not modelled | Hardware timestamps or sync LED/audio for multiple cameras, rolling-shutter time per row and known exposure time (motion blur) |
| Identity over time | IoU/pose/velocity with Hungarian, IoU tracker in v1 | 9–13 % proxy flags, TPR 0,09 for the same shirt, tracks < 60 frames | Independent identity ground truth, ReID, one person per protocol in the lab |
| Occlusion | Visibility per landmark, L/R flags in IK (bodytwin version) | The ankle gate FAIL is partly interpreted as occlusion | Multiple views and an explicit occlusion model per joint |
| Camera motion | `camera_motion_certified` tables (not read in depth), homography per frame for the pitch | Corner view: 57 % height error at conf. 0,776 | Extrinsics per frame with covariance and a fixed calibration object |
| Scale | Pitch homography (only football), population bias 1,2455, body-size prior | Synthetic 0,1 %, real 0,10 m std before bias correction | Known object or scale bar in the image, stereo with known baseline, or IMU |
| Joints and segments | MediaPipe world-landmarks → OpenSim IK on scaled model of another person | Marker-RMS 85 mm, bone-length CV median 14 %; radian error against ground truth not measured | Motion-capture reference (same trial), static calibration per person, calibration of multiple cameras with distortion |
| Body shape → mass/inertia | Silhouette → ellipse (*k* assumed) → circumference → %BF/mass | 4,0 %BF and ≈2,7 kg on the regression link. *k* is an unobserved nullspace | Depth from multiple views or splats with metric scale and **individual validation** (DXA/weighing), segment densities |
| Splat → geometry | ESDF and certificates for floaters/needles | Photometric gates; σ_min AUC 0,88 < scale 0,91; no metric scale | Validation of distances in mm against a known body/phantom, a body-specific task, and awareness that a moving non-rigid body violates SfM's rigidity assumption (`mechanism_video_analyze`: "WRONG tool" for static camera and moving subject) |
| Camera model → uncertainty | p20 (photometry), graf3 (rendering), triangulation Fisher in the splat experiments | None of them propagates an intrinsics/extrinsics covariance to 3D joints | Covariance for a pinhole and distortion model → propagation to 3D points, segment length and joint angle |

**Where it can feed BodyTwin:** BodyTwin uses `mm` and a provenance vocabulary: `region_mass_v1`, `uncertainty_v1` with first-order quadrature, `frame_v1` and `mesh_ingest_v1`. The published reference (`references/current_bodytwin`) has **no** code for camera, video, pose, motion or splats. The hits concern optics for tissue photons (NIRS) and video citations in the sprint cell.

Closest to connecting:

1. **Segment length and body height** as `assumption` or `calibrated` with σ, for scaling segment mass and inertia (I ∝ m·L²). The gap: σ_L/L ≈ 0,06–0,14 gives σ_I/I ≈ 0,12–0,3 for the L term alone.
2. **Joint-angle series from IK** as `derived` and never `measured`, with the flags for coverage and joint limits.
3. **Mass from body shape** only as a prior with regression-RMSE, never as an individual measurement.
4. **Triangulation Fisher** (σ_min and redundancy) as a general certificate for reconstruction from multiple views, for example markers or landmarks from multiple cameras.

## 4. Transfer table

| Source (sha256 prefix) | BodyTwin counterpart | Unique change | Experimental support | Scope | Recommended use |
|---|---|---|---|---|---|
| `football_tracking/build_football_tracking_v0.py` (02b082e8, identical) | None | Hungarian tracker with cost scaled pos/neg and counter-test with shuffled order | Shuffling collapses 2/3 clips; lg fails | The calibration track is self-generated, only pixels, the corpus mount is missing | Method pattern (counter-test with shuffling, honest dt) for tracking in lab video. No data import |
| `hum_retarget/build_hum_retarget_v0.py` (3e6713f0, identical) | None | COCO→SMPL-X mapping with list of unobserved joints | RMS 0,179 m; shuffling = identical (negative); physics-RMS 0,36 m | 8 frames, 28×69 px, implausible height 1,495 m | Only the joint mapping and the list of unobserved joints. The results must not be used |
| `human_capture_chain_v0/v1.py` (4a178948 / 896df2c2, identical) | None | Bone-length CV as a self-consistency measure, foot point on pitch → metric height | CV median 0,14; pitch height synthetically ≤ 0,1 %, corner 57 % | No motion-capture ground truth, scale only with known pitch | CV per segment as **QC gate** before IK and scaling; the pitch method as synthetic reference in E1 |
| `cognition/probe_pitch_anchor.py` (e62d63ef, **only mechanism**) and `mechanism_video_analyze.py` (8707fb59 mechanism, newer than bodytwin 36a01c84) | None | Height from pitch and player extent, bias 1,2455 calibrated on the population | KS p = 0,28 against roster; 0/200 null controls | Ensemble, not individual; 51 public heights as ground truth | Cautionary example: a geometrically exact measure got 25 % systematic bias on real data. The method must not be used for individual scaling |
| `msk/pose_to_opensim_ik.py` (2ce9fb8a mechanism / **eb43d678 bodytwin, newer**) | None in the public reference | Bridge MediaPipe → OpenSim IK; the bodytwin version has flags for occlusion, coverage and joint limits | 228/264 PASS; marker-RMS 85 mm; ankle FAIL | Model for another person, no translation, no per-person scaling | **Start from the bodytwin version**. The angles are labelled `derived` with RMS as σ floor. Prohibited as a basis for mm or radian claims until a motion-capture reference exists |
| `msk/body_composition_from_video.py` (e5dca233 / e0417085) | `region_mass_v1` (density → mass) | Silhouette → ellipse circumference → %BF/mass | Test-RMSE 4,0 %BF, ≈2,7 kg (n = 252 men) | The nullspace *k* is unobserved; one population (men); the regression says nothing about one's own body | Prior with σ. **Do not include** the `reference_body` block in the bodytwin version, because it names a private individual |
| `video_index/football_scene_verifier.py` (4bca1d09, identical) | None | SSIM as a second signal for clip boundaries | 16/16 in visual sample | Small n, same reviewer | Preprocessing: no computations across clip boundaries |
| `video_index/fps_truth_probe.py` (a6182031, identical) | None | Voting on recording rate, where the container does not vote | Tables in `data/video_index_tables/fps_truth` | Declarations, not physical measurement | Mandatory time-base gate for all web video before velocity or angular velocity |
| `L381` decision (927d0a63) | None | Quantified gap in ID switches within the same team | 29/31 admitted; TPR 0,09 | Manifest lost | Registered as a known identity limitation in multiperson video |
| `physics_exp/p20_camera_twin.py` (0435aed9, identical) | `geometry/optics` concerns tissue optics, not camera | PSF, sensor and PTC chain | Self-test 4/4 rerun (§5) | No geometry or time; circular calibration | Only photometric noise link, e.g. noise for landmarks as a function of light. Not as a geometric camera twin |
| `graf3_kameratvilling_v1` and `camera_twin_v1.py` (e73eb32a / 76c6dfff, only CS-I) | None | FOV and CoC for rendering, verified on GPU | c1–c8 PASS, c9 OPEN_NEGATIVE | Rendering camera without inversion or uncertainty | Synthetic forward model to generate test images with known parameters in E1, if rendering is needed |
| `splat_world.py` (277e9e7c) and `splat_quality_cert.py` (3b215fbe, identical) | `mesh_ingest_v1` (surface in mm) | ESDF and health certificates for Gaussians | Self-test for the certificate rerun (§5); ESDF-σ ρ = 0,37 | Not metric; the task is robot collision | Only the certificate as input gate if splats ever provide body surface. No biomechanical use now |
| `l_splat_*` and `l_certified_splat*` (identical) | None | Triangulation Fisher σ_min and redundancy per point | AUC 0,88 (scale 0,91); prospectively 0,549; cropping worsens | ETH3D (rigid scene), GPU, arbitrary scale | The Fisher and redundancy certificate is transferable to landmarks from multiple cameras; the "breakthrough" headings must not come along |
| `splat_capture` | None | Capture pipeline and photometric gates | Latest runs failed | Rigid scene, no scale | Not now |
| `twin_capture_lane` | None | Corpus with World Cup video and SR upscaler | – | Licences and personal data; SR is not measurement data | Only as a source of candidate clips after licence checking. Never SR frames as measurement data |
| `sharp_football` | – | Odds model | – | – | Irrelevant |
| `data/MECHANISM_ANCHOR_GRAPH.json` (01c5d753 versus 2c0c4b17) | BodyTwin's own graphs | Not analysed | – | – | No merging; targeted comparison of individual nodes as needed |

## 5. Rerun self-tests (CPU, ≤ 4 threads, no GPU; write targets checked in advance)

- `p20_camera_twin.py`: writes nothing, only stdout. rc = 0, 4/4 ✓, 0,18 s.
  - Result: Airy 5 px versus 4,88; Strehl 0,77 versus 0,78; gain 2,48 (true 2,5), read noise 5,9 (6,0), full well 19346 (20000); σ_rec 0,113 (0,110).
  - Log: `results/O5/p20_camera_twin_selftest.txt`.
- `splat_quality_cert._selftest()`: no file IO. Returned `True`, 0,45 s.
  - Log: `results/O5/splat_quality_cert_selftest.txt`.
- **Not run:**
  - `selftest_splat_world`: writes to `mechanism/reports/`.
  - Tracker, retarget and human_capture: corpus mounts and venv are missing, and they require MediaPipe.
  - `l_splat_*`: GPU.

The runs only confirm that the code runs and reproduces its own synthetic gates. They say nothing about validity against real cameras.

## 6. Proposed transfer experiments (not run)

**E1: Camera uncertainty → segment length → scaled inertia and moment arm (CPU, synthetic and then real).**

- **Setup:** put a known skeletal model (segment lengths in mm) in a pinhole camera with distortion.
- **Perturbations:** focal length ±1–5 %, camera height and tilt according to the B2 covariance, 2D landmark noise σ = 1–5 px at 70–700 px body height, and errors in the foot point.
- **Chain:** reconstruct with pitch homography (`pitch_calibration.py`, unchanged) → segment length → scaling according to de Leva: m_seg ∝ M, I ∝ m·L², moment arm ∝ L. Propagate with `uncertainty_v1` (quadrature) and with Monte Carlo.
- **Baselines:**
  - (a) Fixed model for another person, subject2, as in today's IK bridge.
  - (b) Body-size prior.
- **Expected:**
  - σ_L/L ≈ 1–3 % with nearby camera and well-anchored pitch, and ≫ 10 % in the corner view.
  - I uncertainty roughly twice as large as L uncertainty.
  - First-order quadrature underestimates in the corner view, because the error is nonlinear and correlated.
- **Counter-test:** coverage. The fraction of true values within propagated 1σ and 2σ must be close to 68/95 %. Measure coverage especially in the case with 57 % error; the current confidence measure did not warn there.
- **Decision rule:** if coverage in the corner view is below 50 %, the camera anchor must not propagate a σ without an extrapolation gate.
- **Real step:** requires video with known individual height and segment lengths, thus a reference with motion capture and anthropometry.

**E2: QC for bone-length consistency before IK, validated against motion-capture reference.**

- **Setup:** run `pose_to_opensim_ik` (bodytwin version) on a public dataset with synchronised video and motion capture and with a licence that permits the use. Licence and access are checked first.
  - Hypothesis: segment-CV > τ and L/R asymmetry predict joint-angle error against the reference.
- **Baseline:** all PASS clips are used.
- **Expected:** RMSE in joint angle at degree level (10°+) for ankle and hip rotation, lower for knee flexion. Rank correlation between CV and angle error > 0.
- **Counter-test:**
  - Switch identity between tracks and permute the flags: the relationship should disappear.
  - Check that a PASS at the joint limit is not counted as correct.
- **The outcome determines** whether video-based angles can be labelled `calibrated` with σ, or only `derived`.

## 7. Overall assessment

- No existing Mechanism chain gives mm or radian quality for an individual.
- Best measured levels:
  - 85 mm marker-RMS on a model for another person.
  - 14 % median-CV for bone length.
  - Scale only with known pitch: synthetic < 0,1 %, real 25 % bias at population level.
  - %BF ±4.
- The camera twins model photometry or rendering, not a measurement camera with covariance and time base.
- What is transferable is the methods: counter-tests with shuffling and null tests, time-base gate, triangulation Fisher and redundancy, provenance flags, QC patterns for IK and a list of unobserved joints. Results and corpus data must not be transferred.
- The next valuable step is E1: cheap, on CPU and synthetic, and it connects camera assumptions to BodyTwin's `uncertainty_v1`.
