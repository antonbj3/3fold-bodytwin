# Researcher scenarios, tested workflows and modality map

Agent O4, 2026-09-22. Code: published BodyTwin `references/current_bodytwin` (commit `b95a8dc`). All runs were made in the session's scratchpad with `BODYTWIN_OUT` outside the source tree. The source tree was not changed (`git status` clean). The contract for combined interventions is in [INTERVENTION_CONTRACT.md](INTERVENTION_CONTRACT.md).

Labeling:

- **DOCUMENTED [Sx]**: the need or workflow is stated in a published source, tool documentation or correspondence.
- **ASSUMED**: our hypothesis about what the role needs. No real user has confirmed it.
- **RUN**: tested against the code in this session.

Sources labeled "retrieved" were checked on the web 2026-09-22. Sources labeled "not checked" are standard references that I cite from memory. They need to be checked before being used externally.

---

## 1. Forskarroller

### 1.1 the collaborator / biomekanisk modellering (dokumenterat fall)

**Dokumenterat.** the collaborator'The group works with digital human twins with a biomechanical focus. "geometry manager"/"geometry engine". He wants to see how we model physiological geometry and parametric shapes. He could not distinguish implemented from planned in our repo. Source: mail 99, 100 and 107 (`references/john_context.md`).

His tools already have three steps for scaling, **DOKUMENTERAT [S1, retrieved]**: affinity landmark scaling (12 DOF), RBF-morphing with thin plate, and surface based STL-morphing. They run in succession, where each step inherits the previous transform. Muscle attachments accompany the bone morph. He has published landmarks/RBF-morphing (2014) and adaptation of torso skeleton to biplanes low dose X-rays with biomechanical conditions (2022)Source: `JOHN_BODYTWIN_CONTEXT_2026-09-22.md`.

**Questions**

1. How do attachments, joint centers and joint axes follow a parameterized shape change with stable identity, and what is the uncertainty in the moment arm? (ASSUMED: his "geometry manager" may mean identity and interface management rather than the morph itself.)
2. Which shape parameters change a calculated joint load more than the measurement uncertainty, and which are negligible? (ASSUMED)
3. Can a geometry with several tissue regions (bone/cartilage/muscle/fat) and interfaces be fed to his musculoskeletal model and to FE analysis without losing identity? (ASSUMED, supported indirectly by his CAD, FE and optimization focus.)


**Data.** Grand Challenge knee loads, SimTK: markers, force plates, EMG from 14 muscles, measured tibial contact force and bone and implant geometry registered from CT. License "MIT Use Agreement" [S3, retrieved]. CAMS-Knee: 6 people with instrumented knee implants, fluoroscopy, markers, GRF and EMG. Data on request, license unknown [S4, retrieved]. OrthoLoad: hip contact forces from 10 people, AVER75/HIGH100 for download [S5, retrieved].

**Missing connections in BodyTwin (RUN, see flow 2):**

- No representation of attachment points, joint centers or joint axes.
- No function from mesh to moment arm.
- No inverse dynamics cell.
- The inertia tensor is calculated (`region_mass_v1.mass_moments_surface`) but is not used by any consumer.


### 1.2 Rehabilitation and movement research

**Questions**

1. Is a changed gait pattern after injury consistent with reduced muscle force, pain-avoiding redistribution of joint torques or changed joint stiffness, and which measurement distinguishes them? (ASSUMED)
2. Is markerless video (2+ smartphones) sufficient to track joint load over rehabilitation weeks, or must a force plate or IMU be added? (Partly DOCUMENTED: OpenCap is described for screening, intervention effects and rehabilitation decisions [S6, retrieved].)
3. How much of the difference between two visits is measurement noise and how much is real change? (DOCUMENTED as a research question: repeatability and variation between sessions for Theia3D [S7, retrieved].)

**Today's workflow.** Marker-based mocap (Vicon) with force plates and OpenSim. Markerless through OpenCap (video → pose → OpenSim IK → dynamics [S6]) or Theia3D (commercial). According to [S7], Theia3D has poorer agreement for pelvic tilt (CMD ≤ 0,48) and for rotation angles in hip and knee.

**Data.** AddBiomechanics 1.0: 273 people, over 70 h of movement and force plates, CC BY 4.0, approximately 4 h download at 50 Mbit/s [S8, retrieved]. Grand Challenge [S3]. GRABMyo: forearm EMG, 43 people × 3 days, CC BY 4.0, 9,4 GB [S9, retrieved].

**Missing connections (RUN, flow 3):**

- BodyTwin has no cell reading video, pose, IMU or force plates.
- `sprint_spring_mass`, `cartilage_poroelastic_relaxation` and `achilles_elastic_energy_return_running` read only embedded literature values.
- `vestibular_vor` refers to a sister cell `postural_sway_pendulum.py` with PhysioNet force plates. That cell does not exist in the published code.
- Mechanism's football tracks are 2D detections per frame without ReID between frames. The longest continuous track is 8 frames. There is no metric scale.

**A useful answer must contain:**

- Competing explanations as parameterized models.
- Which observation distinguishes them.
- An error budget per modality: angle error per joint and plane, timebase error.
- A test–retest boundary for when a change is real.

### 1.3 Tissue and material research

**Questions**

1. Which combination of aggregate modulus, permeability and thickness explains a measured creep or relaxation trajectory, and which test (confined vs unconfined, strain rate) distinguishes them? (ASSUMED, based on standard biphasic theory.)
2. Can absolute material parameters be identified from existing tendon tests, or only normalized shape? (DOCUMENTED in our own code: `tendon_relaxation_v1` states that the UniPD series lacks displacement and specimen geometry, so only normalized relaxation is identifiable.)
3. How do mineral fraction and collagen organization change the modulus, and which measurement channels are independent? (The cell `bone_composite_modulus_bounds` exists.)

**Today's workflow.** Material testing, then fitting in FEBio (MIT license, biphasic and hyperelastic materials [S10, retrieved]). Specimen-specific joint models in Open Knee(s): MR plus joint mechanical testing for 8 knees, CC BY [S11, retrieved].

**Saknade kopplingar:**

- `thickness_v1` propagates thickness uncertainty to τ (dτ/τ = 2·dL/L, RUN).
- No cell connects a measured force–displacement trajectory to FE or biphasic inversion.
- `tendon_paired_analysis_v1` is the analysis pattern, but paired data are missing.
- Materials are global literature constants and not spatial fields.

**A useful answer:** an identifiability table (which parameters a given protocol determines) plus a proposal for the most informative next test. The OED engine that exists for optics (`profiled_information_v1`) may be transferable. This is ASSUMED and untested.

### 1.4 Physiology and systems biology

**Questions**

1. Can glucose load, oxygen demand and cardiac output be explained by a coherent chain, and which step dominates uncertainty? (Partly answered: `metabolic_chain_v1` shows that the blood oxygen step dominates, with half-width 2,32 of 3,05 L/min. The reference body is synthetic.)
2. Does a published CellML model reproduce its own figures, and where does it break against independent data? (DOCUMENTED problem: approximately half the models in the CellML repository were curated to reproduce their papers, and known cases deviated [S12, retrieved].)
3. Which observation (blood test, ECG, breathing) would distinguish two mechanisms for the same trajectory? (ASSUMED)

**Today's workflow.** CellML/SBML models in PMR and BioModels, simulated in OpenCOR or COPASI. Data from PhysioNet with open or credentialed licenses. This is standard practice, not checked in the session.

**BodyTwin today.** Many cells run. Six producers (`bone_stress`, `bone_stress_consequence`, `femoral_neck_stress`, `fluid_compartments`, `metabolic_cost`, `muscle_perfusion`) are unpublished and replaced by synthetic JSON. According to grep, 20 cell files mention their result names, including `cardiac_output`, `thermoregulation` and the bone family. The edge list `tests/cell_dependencies.json` (55 records) does not include these producers as edges, so the dependency is not visible there.

**A useful answer:** the chain with uncertainty per step, the boundaries of the validity domain and a list of unbound inputs.

### 1.5 Imaging and measurement technology

**Questions**

1. How much of a measured cartilage thickness change over 2 years can be due to segmentation, sequence or reconstruction? (DOCUMENTED: method comparison of cartilage thickness in OAI [S13, retrieved search result].)
2. How does segmentation uncertainty from TotalSegmentator or MONAI affect a derived biomechanical quantity such as moment arm, inertia or cortical thickness? (ASSUMED)
3. Which part of a markerless joint angle is a calibration or occlusion artifact? (DOCUMENTED problem [S7].)

**Today's workflow.** 3D Slicer (BSD) for segmentation and registration. TotalSegmentator: 117 CT classes, dataset CC BY 4.0 with 1 228 CT, weights Apache-2.0 [S14, retrieved]. MONAI (Apache-2.0) for own training. Slicer and MONAI were not checked in the session.

**Data.**

- TotalSegmentator CT [S14].
- SKM-TEA: knee MR with raw data, 155 patients, 1,6 TB. The license has not been established [S15, retrieved].
- OAI: 4 796 people, access through NIMH Data Archive with an agreement [S13, retrieved].
- BodyParts3D: CC BY-SA 2.1 JP, already used.

**BodyTwin idag.**

- `mesh_ingest_v1` requires mm and closed surfaces.
- `geometry_validity_v1` rejects unevaluated self-intersection.
- `uncertainty_v1` exists as vocabulary.
- No segmentation or image error model is connected to the geometry's uncertainty.
- The only public anatomy is a kidney mesh, with volume error 0,0975 % at 2 mm. According to the field handoff, volume is an insufficient control metric.

**A useful answer:** an error chain from image to quantity (voxel, segmentation, mesh, parameter) with separate contributions.

### 1.6 Intervention and treatment research

**Questions**

1. How does a combination of two interventions differ from the sum of the individual ones, and does the answer depend on the order? (ASSUMED, the assignment from Anton.)
2. What model structure is required for a simulated combination effect to be credible according to ASME V&V40 (context of use, verification, validation, UQ)? (DOCUMENTED framework [S16, retrieved].)
3. Which measurable intermediate step (biomarker, strain, image metric) tests the prediction earliest? (ASSUMED)

**Today's workflow.** In silico trials with credibility assessment [S16]. Pharmacokinetics/pharmacodynamics and mechanistic ODE models in SBML or CellML.

**BodyTwin.** The contract and an executable example are in [INTERVENTION_CONTRACT.md](INTERVENTION_CONTRACT.md) (RUN).

### 1.7 Surgical and prosthetic research

**Questions**

1. How does an osteotomy or implant position (e.g. femoral offset, neck angle) change joint contact force and bone strain, and what is measurable afterward? (Partly DOCUMENTED: instrumented implants measure contact force in vivo [S3, S4, S5].)
2. How are interfaces between bone and implant and material partitions preserved through deformation and meshing? (ASSUMED, shared with dental.)
3. Which part of the load change after the intervention is due to geometry and which to a changed movement pattern? (ASSUMED)

**Today's workflow.** CT, then Slicer segmentation, CAD implants, FE (FEBio/Abaqus) and a musculoskeletal model (the reference model/OpenSim). Validation against OrthoLoad and Grand Challenge.

**BodyTwin.**

- The bone chain (flow 1) has the shape parameters neck angle, offset and cortical thickness in a beam model.
- The producers are unpublished and fed with synthetic substitutes.
- No connection exists from mesh to cross-section.

---

## 2. Three simulated workflows through today's code

### Flow 1 — "How do femoral neck angle and offset change the bone's remodeling regime?" (surgery, prosthetics, tissue, the collaborator)

Command (RUN):

```
cp -r examples/synthetic_inputs/* $OUT/
BODYTWIN_OUT=$OUT python src/bodytwin/cells/musculoskeletal/{bone_remodeling,bmu_turnover_kinetics,bone_strain_bounded_classifier_proposal}.py
```

The chain is shape parameter (NSA, offset, outer diameter, cortical thickness) → beam stress → strain σ/E (E = 17 GPa) → mechanostat regime → BMU turnover.

**Resultat.**

- `bone_remodeling` gives 2 000 µε medially and 1 647 µε laterally, with the sweep 1 529–2 412 µε. Everything is classified as "formation, not maintenance".
- The cell itself warns: the beam is dominated by bending (σ_axial 3 MPa against σ_bend 31 MPa), because muscle forces are missing. The verdict is therefore probably too high.
- `bmu_turnover_kinetics` gives 7/8 neck configurations in overload. Three of ten gates are False, but the process exits with 0.

**Where the flow breaks:**

| Brott | Typ | Detalj |
|---|---|---|
| The producers `bone_stress` and `femoral_neck_stress` are unpublished | missing cell | The consumers read handwritten synthetic JSON |
| The synthetic flag is lost | missing contract | Inputs have `"synthetic": true`. Outputs from `bone_remodeling` and `bmu_turnover_kinetics` contain 0 occurrences of "synthetic" and are labeled as results. `result_envelope_v1` is only connected to `metabolic_cost`. |
| Muscle forces are missing | missing mechanism | The beam is loaded only by joint force with offset. Bending is exaggerated, which the cell itself reports. |
| Mesh to cross-section is missing | missing link | Outer diameter and thickness are scalars. No connection exists from `mesh_ingest`/`thickness_v1` to the beam section. |
| Exit code does not reflect gates | contract | Gates that are False still produce exit 0 in `bmu_turnover_kinetics` |
| Conflicting measurements | incomplete anchors | `bone_remodeling` notes that "habitual gait in maintenance window" is not supported by Frost 1987. Cortical turnover is outside a factor of 2 against the stated anchor. |

**Shared building block.** Shape parameter → stress/strain is used by surgery (implant position), tissue (mechanostat) and the collaborator (geometry to load).

### Flow 2 — "How does a shape parameter change muscle force, moment arm and joint reaction?" (the collaborator, dental, surgery)

Composition in the scratchpad (RUN, `flow2_flow3.py`). The chain uses only exported functions:

1. `region_mass_v1.mass_moments_surface` (mesh → volume, mass, inertia)
2. `region_mass_v1.physiological_cross_section_mm2` (PCSA = V·cosθ/L_f)
3. `muscle_pcsa_crossbridge_specific_tension_rebuild.sigma_n_per_cm2` (specific tension 25,35 N/cm² from the crossbridge budget)
4. `tmj_lever_model.sagittal_equilibrium` (sagittal torque balance about the condyle)

The muscle is a **synthetic ellipsoid** (12×20×30 mm, as a masseter stand-in) scaled isotropically with s = 0,9/1,0/1,1. Fiber length (25·s mm) and pennation (15°) are ASSUMED. The moment arm is held at 22,5 mm and the bite distance at 50 mm (center of the cell's sweep).

| s | V (mm³) | massa (g) | I_zz (kg·mm²) | PCSA (mm²) | F per sida (N) | F_bett molar (N) | R_kondyl (N) |
|---|---|---|---|---|---|---|---|
| 0,9 | 21 928 | 23,2 | 2,04 | 941 | 239 | 215 | 262 |
| 1,0 | 30 079 | 31,9 | 3,46 | 1 162 | 295 | 265 | 324 |
| 1,1 | 40 036 | 42,4 | 5,58 | 1 406 | 357 | 321 | 392 |

**Check.** PCSA scales as s², and inertia as s⁵. Both follow analytically.

**Where the flow breaks:**

1. **Moment arm and bite distance do not follow the shape parameter.** With isotropic scaling of the entire head, both r_muskel and d_bett scale with s. The ratio F_bett/F_muskel is then invariant: condylar reaction per bite force does not change. The geometry's specific contribution to joint load thus lies in *anisotropic* changes of attachments and joint centers, and that representation is missing. This is precisely the gap that the collaborator's "geometry manager" reasonably refers to (ASSUMED).
2. No muscle mesh exists in the repo, and no fiber length or pennation is measured (ultrasound/DTI missing).
3. The inertia tensor has no consumer, because there is no inverse dynamics.
4. The contract works: `ThicknessDistribution` rejected provenance `'synthetic, not MR'` and required one of `measured/calibrated/literal_cited/assumption/UNKNOWN`. This is good behavior.

**Shared building block.** Mesh → PCSA/inertia → torque is needed by the collaborator, dental (TMJ, bite force), rehab (strength) and surgery.

### Flow 3 — "Is video sufficient to say anything about cartilage or tendon load in walking or running?" (rehab, tissue, imaging)

RUN with `thickness_v1.relaxation_time_s` and assumed thicknesses 1,8–3,0 mm (H_A = 0,73 MPa, k = 4,3e-15). τ became 1 032–2 867 s. A relative thickness uncertainty of 10 % gives 20 % in τ.

**Where the flow breaks:**

1. **Video → movement is missing in BodyTwin.** No cell reads pose, IMU or force plates. Mechanism's tracks are monocular 2D detections in broadcast TV without metric scale and with short identity (8 frames). The camera twin `p20_camera_twin` models lens and sensor but not calibration against a moving person.
2. **The question is insensitive to the video.** Stance phase is 0,7–1,0 s while τ is ≥ 400 s. The margin is over 100× for all reasonable thicknesses (`cartilage_poroelastic_relaxation` G3). Video thus adds no information to *this* cartilage question. Cartilage thickness (MR) and permeability (test) dominate. This is a concrete example of a modality needing to be chosen against a quantity and not added generally.
3. Where video **can** add information: `sprint_spring_mass` (leg stiffness from contact and flight time) and Achilles energy (force from ID). Both require a timebase with known error and a metric scale. Today's cells read embedded literature constants.
4. Thickness lacks a measurement source. MR data (OAI with agreement, SKM-TEA with unknown license) are not loaded. Open Knee(s) with CC BY is the closest usable option.

---

## 3. Modalitetskarta

"Frame" means the coordinate frame in which the measurement is given. "Derived" means that the quantity requires a model or reconstruction beyond the instrument. The error models are typical values from the literature and must not be read as checked numbers, except where a source is given.

| Modality | What is actually measured | Unit | Frame | Error model (main terms) | Calibration | Measured / derived | BodyTwin quantity it can constrain |
|---|---|---|---|---|---|---|---|
| CT | X-ray attenuation | HU | scanner/patient (DICOM) | noise, partial volume, beam hardening, metal artifact; voxel 0,5–1 mm | phantom (water/air), density phantom for BMD | HU measured; density, modulus and segmentation derived | bone geometry (`mesh_ingest`), cortical thickness, bone density → modulus (`bone_composite_modulus_bounds`) |
| MR | Proton signal (T1/T2/PD weighting), quantitative T2 maps | arbitrary intensity; T2 in ms | scanner | B0/B1 inhomogeneity, geometric distortion, movement, segmentation bias | phantom, shimming | signal measured; thickness, volume and T2 derived through segmentation/fitting | cartilage thickness (`thickness_v1` → τ), muscle volume → PCSA, soft tissue regions |
| Ultrasound (B-mode) | Echo time and amplitude | mm (via assumed sound speed 1 540 m/s) | probe (moving) | assumed sound speed, probe pressure, 2D plane, operator-dependent | phantom; probe positioning through optical tracking | fascicle length and pennation derived | fiber length and pennation (flow 2), tendon elongation → stiffness |
| Biplanar X-ray (EOS / dual fluoroscopy) | 2 projections | px → mm | lab/scanner | model–image registration error, limited soft tissue | calibration object | 3D bone pose derived through 2D–3D registration | joint kinematics without skin artifact, bone shape through SSM (the collaborator's 2022 work) |
| Video, multicamera, markerless | Pixel intensity → 2D keypoints → 3D | px, then m, grader | camera → world frame | calibration (intrinsics/extrinsics), sync, occlusion, pose network bias; rotation angles and pelvic tilt weakest [S7] | checkerboard/ChArUco; sync | everything except pixels is derived | segment kinematics, timebase (contact phases), leg stiffness (SLIP), input to ID |
| 3D reconstruction/splats | Images → Gaussian scene | arbitrary scale until it is metric | reconstruction frame | scale ambiguity, "floaters", needles and ghosts (`splat_quality_cert`), no tissue information | known distances/markers | entirely derived | external body shape and volume (with scale) → segment mass and inertia through a density assumption |
| IMU | Acceleration, angular velocity (magnetic field) | m/s², rad/s | sensor frame | bias, drift, soft tissue artifact, sensor-to-segment calibration | static pose, functional movements | orientation and position derived (integration) | segment angular velocity, events (heel strike), inertial forces |
| Force plate | Force and torque | N, N·m | plate frame → lab | offset, drift, CoP error at low load, cross-coupling | static load, CoP calibration | GRF measured; CoP derived | external load for ID, stance phase, `sprint_spring_mass`, postural sway |
| EMG (surface/HD-) | Electrical potential | µV | electrode (skin location) | crosstalk, electrode displacement, impedance, normalization choice | MVC normalization | activation and timing derived; force requires a model | `motor_unit_recruitment` (recruitment/rate coding), muscle activation limit in ID |
| Instrumented implant | Contact force in prosthesis | N | implant frame | strain gauge calibration, telemetry loss | factory calibration | measured (in vivo) | **independent truth** for joint load (flow 2 and the bone chain) |
| Lab tests (blood, urine) | Concentration | mmol/L, pM | time (sample) | preanalytics, assay variation, circadian rhythm | calibrators | measured | endocrine and metabolic cells (`glucose_insulin_minimal_model`, PTH → `bone_rankl_opg_lemaire2004`) |
| Tissue tests (mechanics) | Force, displacement, time | N, mm, s | specimen/machine | grip slip, unknown specimen geometry, hydration, temperature | load cell, extensometer | force measured; modulus and permeability derived | `tendon_paired_analysis_v1`, `cartilage_poroelastic_relaxation` (H_A, k) |
| Histomorphometry | Labeled surfaces and distances | µm, µm/dag | biopsy | sampling site, 2D section | — | MAR measured; turnover derived | `bmu_turnover_kinetics`, `bone_remodeling_transient_arithmetic` |

### Combinations that add information

- **Video + force plate.** The force plate gives stance phase and GRF independently of the image. It tests the video's timebase and makes ID possible. The errors are independent: optics against load cell.
- **Markerless/markers + instrumented implant (Grand Challenge, CAMS-Knee).** The only independent truth for joint contact force. It tests both geometry (moment arm) and muscle distribution.
- **MR (volume) + ultrasound (fiber length/pennation).** PCSA requires both, and they have different error sources.
- **CT (geometry) + biplanar X-ray (kinematics).** Separates shape error from movement error.
- **EMG + ID torques.** EMG constrains the muscle distribution left undetermined by ID. It is independent of moment arm geometry, except for crosstalk.

### Combinations sharing the same error or source (do not count as independent)

- **Video kinematics + splat body shape from the same cameras.** Calibration, sync and scale are shared. A scale error affects both segment lengths and masses.
- **IMU kinematics + IMU events.** The same sensor and the same drift.
- **Moment arm from a morphed generic model + muscle force from the same model's static optimization.** Validation against "reasonable muscle forces" is circular. Only implant force or EMG is independent.
- **Segmented volume + segmented thickness from the same MR and the same network.** Segmentation bias is shared.
- **Literature anchors from the same cohort.** Examples are Frost thresholds cited through secondary summaries (reported in `bone_remodeling`) and "Parfitt/Eriksen" turnover. The same underlying source must be counted once.
- **Synthetic substitutes + cells already based on the same textbook constants.** E = 17 GPa is used both in the stress substitute and in `bone_remodeling`. The consistency gate is then self-fulfilling ("bit-for-bit consistency" in `bmu_turnover_kinetics` is true by construction).

---

## 4. Combined and sequential interventions

See [INTERVENTION_CONTRACT.md](INTERVENTION_CONTRACT.md). Smallest executable example: `bone_rankl_opg_lemaire2004.rhs` with inputs I_P (PTH), I_O (OPG) and I_L (RANKL).

RUN result, osteoclast AUC excess against the baseline trajectory:

- PTH alone: +0,166 pM·dag.
- OPG alone: −0,0028.
- Simultaneous: +0,150, which is a non-additivity of −0,0127 (approximately −8 % of the PTH effect).
- Order effect (PTH→OPG against OPG→PTH): +0,0008.

The PTH peak reproduces the cell's own +304 % (4,04×). **How strong the OPG arm is was our own choice and is not calibrated.**

---

## 5. Three broad experiments for the day, ranked

Each experiment requires at least two independent data sources or mechanisms and has a counter-test set before the run.

### 1. Shape parameter → moment arm → knee contact force against measured implant force (Grand Challenge, MIT Use Agreement)

- **Roles:** the collaborator, surgery/prosthetics, rehab.
- **Data:** CT-registered bone and implant geometry, markers, GRF, EMG and measured medial/lateral tibial force [S3].
- **Execution:** Build the missing link, that is, attachment points and joint centers with stable IDs on the bone mesh, following a parameterized deformation (affine → RBF, the same class as the reference model [S1]). Calculate moment arms and run ID/static optimization (OpenSim is Apache-2.0 and executable without a license). Compare contact force with generic and geometry-specific moment arms.
- **Preregistered counter-test:** If the geometry-specific moment arm does not reduce RMS error against measured contact force more than the difference between two generic scalings, the hypothesis that geometry is the limiting factor fails. Muscle distribution then limits, which the EMG arm tests.
- **Negative control:** Isotropic scaling must produce an unchanged ratio of contact force/external load, according to flow 2.
- **Benefit:** Directly closes the gap in flow 2 and gives the collaborator demo an external truth.

### 2. What markerless movement can constrain: timebase and leg stiffness from kinematics against force plates (AddBiomechanics, CC BY 4.0)

- **Roles:** rehab, imaging/measurement technology, sports.
- **Execution:** Take running and walking trials with both kinematics and GRF. Derive contact and flight time and SLIP leg stiffness (`sprint_spring_mass`) *only* from kinematics, with contact estimated from the foot/toe trajectory. Compare with GRF-derived values. Propagate the measured timebase error to Achilles energy (`achilles_elastic_energy_return_running`) and to the cartilage margin (`cartilage_poroelastic_relaxation`).
- **Counter-test:** If the kinematic timebase error produces a leg stiffness error greater than the variation between individuals in the same dataset, video without force plates cannot distinguish individuals in this quantity.
- **Negative control:** The cartilage margin must be insensitive (over 100×). If it proves sensitive, the chain is wrong.
- **Shared source to avoid:** AddBiomechanics kinematics are themselves fitted against GRF (dynamic consistency). Use raw marker IK, not the dynamics-adjusted trajectories, otherwise the source is shared.

### 3. Mechanical load + biochemical intervention in bone turnover: combination, order and contract gate

- **Roles:** intervention research, tissue, physiology.
- **Execution:** Connect the mechanostat's strain (`bone_wolff_law_mechanostat.mechanostat_rate` or strain from the beam cell) to the RANKL/OPG input in `bone_rankl_opg_lemaire2004`, following a published mechanoregulation form (e.g. Scheiner/Pivonka lineage; not checked in the session). Add a bone mass state through `bone_remodeling_transient_arithmetic`. Run the contract's six arms (baseline, A, B, A+B, A→B, B→A) from [INTERVENTION_CONTRACT.md](INTERVENTION_CONTRACT.md). As a second measurement source: the histomorphometry anchors (MAR, BMU period) and the PTH response from Lemaire 2004.
- **Counter-test:**
  1. Non-additivity and order effect must exceed baseline drift (C-AUC 1,5e-5) and solver tolerance by at least 10×, otherwise they are reported as zero.
  2. Without mechanical coupling, the load arm must produce exactly zero osteoclast effect (zero case).
  3. Every synthetic input must carry through as a flag in the output. This remedies the break in flow 1.
- **Benefit:** The first intervention chain across two mechanisms. Tests the contract and flag propagation simultaneously.

**Order and rationale.** Experiment 1 has the highest documented need (the collaborator, the Grand Challenge tradition) and an independent truth. Experiment 2 is the cheapest and determines what video can actually constrain, which directs the Mechanism transfer. Experiment 3 is the newest for BodyTwin but has the weakest external validation. Calibration against clinical intervention data is missing.

---

## Sources

- [S2] Seth et al. 2018, OpenSim, PLOS Comput Biol 14(7):e1006223; Apache-2.0 — https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006223 (retrieved)
- [S3] Grand Challenge Competition to Predict In Vivo Knee Loads — https://simtk.org/projects/kneeloads ; Fregly et al. 2012 J Orthop Res (retrieved)
- [S4] Taylor et al. 2017, CAMS-Knee, J Biomech — https://pubmed.ncbi.nlm.nih.gov/29037443/ ; https://movement.ethz.ch/data-repository/cams-knee-project.html (retrieved; license not established)
- [S5] Bergmann et al. 2016, Standardized Loads Acting in Hip Implants, PLOS ONE 11(5):e0155612; orthoload.com (retrieved)
- [S6] Uhlrich et al. 2023, OpenCap, PLOS Comput Biol 19(10):e1011462 (retrieved)
- [S7] Theia3D validations: Kanko et al. 2021 (cited in search results); systematic review 2025 — https://www.sciencedirect.com/science/article/pii/S0933365725002672 ; healthy/clinical comparison — https://www.nature.com/articles/s41598-024-80499-8 (search results, not full text)
- [S8] AddBiomechanics Dataset 1.0, CC BY 4.0 — https://addbiomechanics.org/download_data.html ; Werling et al. ECCV 2024 (retrieved)
- [S9] GRABMyo v1.1.0, PhysioNet, CC BY 4.0, 9,4 GB — https://physionet.org/content/grabmyo/1.1.0/ (retrieved)
- [S10] Maas et al. 2012, FEBio, J Biomech Eng 134(1):011005; MIT license (search result)
- [S11] Open Knee(s), SimTK; data CC BY — https://pmc.ncbi.nlm.nih.gov/articles/PMC7890148/ (search result)
- [S12] Lloyd et al. 2008, CellML Model Repository, Bioinformatics 24(18):2122; reproducibility cases PLOS ONE e0254749 (search results)
- [S13] OAI through NIMH Data Archive — https://nda.nih.gov/oai ; cartilage thickness method comparison — https://pmc.ncbi.nlm.nih.gov/articles/PMC10076900/ (search results)
- [S14] TotalSegmentator, Wasserthal et al. 2023 Radiol AI; dataset CC BY 4.0, weights Apache-2.0 — https://github.com/wasserth/TotalSegmentator (search result)
- [S15] SKM-TEA, Desai et al. NeurIPS 2021 — https://github.com/StanfordMIMI/skm-tea (license not established)
- [S16] Viceconti et al. 2021, In silico trials: VV&UQ, Methods 185:120–127 — https://pubmed.ncbi.nlm.nih.gov/31991193/ (search result)
- Not checked in the session: Fedorov et al. 2012 (3D Slicer), Cardoso et al. 2022 (MONAI), Updegrove et al. 2017 and Vascular Model Repository (SimVascular), Goldberger et al. 2000 (PhysioNet), Lemaire et al. 2004 J Theor Biol, Scheiner/Pivonka mechanoregulation.
