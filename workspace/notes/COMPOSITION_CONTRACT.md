# Composition Contract for BodyTwin (C1, 2026-09-22)

Contract and inventory were written before the run (PREREG sha256 in `results/C1/PREREG.sha256`); §3–§5 after. All numbers in §3–§5 are in `results/C1/c1_results.json` or `results/C1/out/*/stage*.json`.

## 1. Inventory: how the parts represent the same thing

| Part | Geometry | Identity | Device | Frame | Material | Uncertainty | Provenance/status |
|---|---|---|---|---|---|---|---|
| BodyTwin framework `result_envelope_v1` | none (scalars) | quantity key (dotted string), `region_id` | free string `units`, checked against consumer expectation (`validate_for_consumer`, row 492) | `frame` = `"scalar"` or JSON-descriptor (row 167, 418), canonically compared; ** none transform** | any `material_region_id` | `uncertainty` any scalar; required only with `require_uncertainty` | 7 term vocabulary `PROVENANCE_VOCABULARY` (line 90); `value=None` = abstention, refused (row 598) |
| BodyTwin `region_envelope_v1` | region-payload | `region_id` against frozen registry | mm, mm², kg, J/K | copied descriptor, never invented | roster `I2_MATERIAL_REGION_IDS` | scalar | `derive_provenance` = minimum trust (row 116, order row 97) |
| BodyTwin geometry `mesh_ingest_v1` , `region_mass_v1` | triangle mesh | mesh-sha256 | **only mm** (`mesh_ingest_v1.py:30` throws otherwise), density kg/mm³ | implicit (dataset's) | `DENSITY_REGISTRY` (`region_mass_v1.py:126`); `lookup_density` returns `None` , never default value (line 155) | none | per material property (`literal_cited`, `assumption` ) |
| D1 `geomgr` | shares + landmarks as linear functional of vertices | `Landmark.id` stable, **no version** | **only mm** (`landmarks.py:34`) | `frame='dataset'` string (`landmarks.py:30`), no transform | no | TPS-interpolation error per landmark (`model.py` `interp_error_mm`) | own 4 term vocabulary `measured/atlas/derived/synthetic` (`landmarks.py:17`), **Not the same as the framework** |
| H5 `h5_moments` | closed surface → (V,S,T) | none | mm, kg/mm³, kg·mm² (`h5_moments.py:3`) | CT-coordinates; `affine_moments` takes F,t but no frame identity | density per shell, caller returns | none | none (README says inner shell is SYNTHETIC, code doesn't carry it) |
| M1 motion engine / pinocchio | URDF links | link name (`femur_l`) | m, kg, kg·m² | **pinocchio assumes z up**, OpenSim y up (`m1_engine.py:24`, corrected by hand after 946 % error) | none | none | none |
| The field motor `faltkarna_v1_mesh_to_sdf` | dense/sparse SDF, sample point = `lo + i·pitch` (not cell center) | none; `kropps_id_falt_v1` has owner index per block | mm | local window `origin_l` | owner channel, not material | waterproof gate, `ytkorrektion` ±pitch/2 | none |
| Graph engine `graph_store` (compose example in motion-engine) | – | node ID `n.<domain>.<namn>`, without version | free numbers in `numbers` | – | – | – | `status` string (`VERIFIED-FRESH`/`FAILED`) |
| K1 components (50 leaf) | – | see `representation.identity`: e.g. "vertex index (unstable between runs)", "marker name", "joint ID" | mixed mm/m/rad | scanner frame LPS/RAS, bone frame, global, joint coordinate | – | sensitivity base, 14 of 50 UNKNOWN | node status |
| Mechanism/splats | Gaussians, camera images | track-ID | pixels; MediaPipe hip-centered meters; COLMAP arbitrary scale | camera, hip centered | – | eikonal residual (rank 0,37 against true error) | ABSTAIN at wrong input contract |
| VSD-data (this run) | Femur_L-mesh | subject `z001`, landmark name + **corner index from v2-mesh** which does not fit in v3 | mm | CT LPS | – | 5 rater × 4 attempt | measured (CT) |

### Incompatibilities that actually exist

1. **mm against m.** The geometry side (BodyTwin geometry, geomgr, H5, the field engine) is mm-only and refuses otherwise. The dynamics side (pinocchio, OpenSim, the reference model) is m. Inertia separates 10⁻⁶ (kg·mm² → kg·m²). No part carries a conversion rule; M1 recalculated by hand (`M1/README.md`, "What did not work", point 7).
2. **y up against z up.** OpenSim and gait2392 have y up, pinocchio's default gravity has z up. M1 got 946 % error in `pelvis_ty` before gravity was set by hand.
3. **Frame as string without transform.** The framework compares frame descriptors but cannot transform; geomgr writes `'dataset'`; H5 has no frame identity. A CT-vector and a segment vector have the same type.
5. **two status vocabularies.** geomgr `measured/atlas/derived/synthetic` against framework `measured/calibrated/calibrated_public_model/literal_cited/tuned/synthetic/unknown`. `derived` is an operation, not an origin.
6. **Flags that are dropped.** `bone_remodeling`/`bmu_turnover_kinetics` never reads `synthetic` (X3: 0 occurrences in output). H5: s synthetic inner shells are only worn in README. Field engine's A188: UNKNOWN was never written, thin region got 0 cells with `transfer_ok=true`.
7. **Subject identity.** Geometry from one subject (VSD z001, 1,80 m, 87 kg) and movement from another (gait2392 subject01, 72,6 kg) can be combined without saying anything.
8. **Osecurity as a scalar.** No part can say "the known part of the error budget is X, but term Y is unknown".

## 2. The contract (extension of `result_envelope_v1`, not a parallel one)

Each link in a chain takes in and out **quantities** and puts them into a regular `ResultEnvelope` (`make_envelope`), sealed with `seal_chain` against the upstream link's sealed artifact. The novelty lies in what a quantity must carry and in six rules that a link maintains before it counts.

### 2.1 Entities

`Body` (subject), `Segment`, `Region` (material region in a segment), `Landmark`, `Joint`, `Material`, `Field` (SDF/ownership field), `Trajectory` (q, q̇, q̈ over time), `Load` (external force, GRF), `Observation` (raw measurement, eg an assessor's landmark). Everyone has `id@version` and `body` (which subject they belong to).

### 2.2 Mandatory fields per quantity

| Field | Rule |
|---|---|
| `id` | `<entitet>:<body>/<namn>@<version>`, e.g. `Landmark:vsd.z001/Femur_L.FHC@v3.0+c1` |
| `value` | number/array or `None`. `None` = UNKNOWN (abstention). |
| `units` | from a closed unit table (mm, m, kg, kg/mm³, kg/m³, kg·mm², kg·m², N, N·m, row, s, 1). Conversion only via `convert()`, which enters the factor in the derivation. |
| `frame` | frame-ID out of a frame registry. Frames are only exchanged via a registered `Transform(src→dst)` which itself is a quantity with status and uncertainty. |
| `status` | framework vocabulary. geomgr's `atlas` → `calibrated_public_model`; `derived` is no status: derived value gets **least confidence among its inputs** (`derive_provenance`). |
| `uncertainty` | dict `{term: value}` in the unit of quantity; a term may be `"UNKNOWN"` (the total becomes a lower bound with the list of unknown terms) or `"ENSEMBLE"` (the term has a different unit, such as the angle of the frame, and is evaluated by rerunning the chain per observation). |
| `assumptions` | list, the union of all upstream assumptions. |
| `derivation` | upstream quantity id + link source sha256. |

### 2.3 Rules, maintained by the link (not the producer)

- **R1 Status is propagated downward.** Output status = `derive_provenance(indata)`. A producer cannot raise it. Synthetic in ⇒ synthetic out.
- **R2 UNKNOWN never becomes a number.** A link that receives `value=None` in an input it needs throws `ChainStop` and writes no output number. (Uncertainty `"UNKNOWN"` does not stop; it makes the total a lower bound.)
- **R3 Device and frame are checked on each input.** The link declares expected unit and frame per input. Wrong unit ⇒ `UnitError`, wrong frame ⇒ `FrameError`. Raw arrays without quantity wrappers are refused.
- **R4 Invariant on the link boundary.** Each link that changes representation (mesh → moment, mm → m, CT → segment frame) checks at least one physical invariant against independent data (length ratio, mass sign, positive definiteness, right-handedness, anatomical direction). A mislabeled number (mm numbers labeled m) is caught here, not by the label.
- **R5 Identity.** Quantities from different `body` may only be combined via a declared `Mapping` link (registration/scaling) that is in `assumptions` and `uncertainty`.
- **R6 Conservation per region.** A representation switching link (mesh → field) counts volume per declared region. A region with 0 cells becomes UNKNOWN (A188/N7), never a neighbor label.

### 2.4 Graph node

The last link writes a node in graph-engine's `{nodes, edges}` format (same as `3fold-motion-engine/examples/compose`), with `status` = contract status, `numbers` with entities in the key name, `uncertainty` and `assumptions`, and a `produces` edge from each link.

## 3. Thin implementation (`results/C1/compose/`)

| File | Contents |
|---|---|
| `contract.py` | `Q` (quantity), closed `UNITS`-table + `convert`, `Frame`-registry, `Transform` (rejects non-orthonormal, left-handed and between-subjects), `Mapping` (R5, scales per dimension: L¹, M/L³ s³, ML² s⁵), `Link` (`check` = R2/R3/R5 in order form → UNKNOWN → unit → frame → identity; `out` = R1 via `derive_provenance`; `envelope` = `make_envelope` + `seal`/`seal_chain`), `quantities_from_sealed` (verifies receipt and expected upstream artifact). Numpy + `bodytwin.framework`, works in both venvs. |
| `stage1_geometry.py` (`.venv-bodytwin`) | L0–L5: VSD `Femur_L` z001 → 20 rater observation → geomgr barycentric anchors → anatomical frame (R4: right-handed, anterior/right/up to LPS) → mapping to subject01 (R4: aspect ratio) → geomgr TPS + certificate → H5 exact torque + field motor `surface_raster_and_flood` 2 mm (R6 regions). |
| `stage2_dynamics.py` (`.venv-motion`) | L6–L8: reads L5 via sealed envelope only; thigh = gait2392 − leg(0) + leg(p) (R4: the soft part rest should be a physical body); the knee joint is moved; M1's `build_state`/`wrench_world`/`gen_forces` + pinocchio `rnea`, gravity as quantity; hip intersegmental force out of `data.f`; graph node in graph-engine format, read back with `GraphStore`. |
| `check_knee_shift_mechanism.py` | Independent Newton–Euler sum (without rnea) for the knee displacement contribution. |
| `summarize.py`, `run_all.sh` | The whole run (about 20 s per variant of step 1, a few seconds per variant of step 2; maximum 283 MB RSS measured) → `results/C1/c1_results.json`, log `results/C1/run_all.log`. |

The graph node: `results/C1/graph_node.json` (nominal) and `graph_node_synthetic.json`. Nothing written in any graph.

Additions to the contract during construction, both forced by data:
- Osecurity terms that have units other than the quantity (frame angle in degrees, origin in mm) may not be summed. They are labeled `ENSEMBLE` and are evaluated by rerunning the entire chain for the 20 rater observation. First version RSS-summed degrees and meters; it was detected in the output and corrected before final execution.
- The envelope key is the entire quantity id. The first version used the last part of the id, and mass, center of mass, and inertia overwrote each other in the envelope. The link now refuses two outputs with the same id.
- K The source code hash of the receipt caught that step 1 had changed after the seal (`stale result: producer code changed since seal`), so step 2 refused to read. It is intended behavior; the final run was done with `run_all.sh` from blank.

## 4. Outcome (one run, `run_all.sh`)

The chain: VSD z001 Femur_L (21 380 corner, certificate approved) → HJC–KJC 414,2 mm (20 observations: SD 3,05 mm, frame rotation versus mean 1,9° on average, at most 6,3°) → scale to subject01 s = 1,096 (gait2392 hip–knee 0,4541 m) → bone mass 1,175 kg (1,9 g/cm³) → mapped 1,548 kg → p = +10 mm length via TPS (certificate approved, knee center moved 9,997 mm) → thigh 8,984 → 9,010 kg, soft tissue remainder 7,44 kg (positive definite) → ID above left HS→HS (1,247–2,46 s).

| Output (left hip, p = +10 mm) | Value |
|---|---|
| peak hip flexion moment at p = 0 | 48,03 N·m (same as M1/OpenSim) |
| max \|Δτ\| hip flexion | 2,43 N·m, in support (t = 1,95 s); in swing 0,64 |
| of which leg moment / knee displacement | 0,30 / 2,13 N·m (additive to 1,4·10⁻¹⁴) |
| max \|Δτ\| adduction / rotation | 0,56 / 0,16 N·m |
| hip intersegmental force, peak at p = 0 | 654 N (0,92 BW; without muscle forces) |
| relative change in force peak | −0,03 % |
| status | `literal_cited` (the density); strict consumer refuses, demo consumer accepts |

Error budget for Δτ hip flexion: rater spread throughout the chain (20 runs) SD 0,026 N·m; 2 mm voxel against exact 0,0024 N·m; engine numerics against OpenSim 5·10⁻⁹ N·m; known RSS 0,026 N·m. **Lower bound**, five terms UNKNOWN: density dispersion (registry indicates none; ±10 % gives 0,030 N·m), A_map, A_soft, A_kin, A_frame.

Mechanism check: the knee contribution calculated as pure Newton–Euler sum over the distal bodies (5,09 kg) matches the chain rnea result to 1,4·10⁻¹³ N·m. My first check left out that the acceleration of distal bodies changes as the knee moves and just explained 0,65 by 2,13 N·m; that version is reported in JSON.

### Counter sample

| Test | Outcome |
|---|---|
| (a) synthetic flag into L0 | **passed.** Status `synthetic` in all seven sealed envelopes (L0–L5, L7) and in the graph node; strict consumer refuses. The numbers are identical to nominal run, only the status differs. |
| (b) density UNKNOWN (`MSK-BONE-CORTICAL` is not in the registry) | **approved.** `ChainStop` in L5a, no envelope L5 and no minutes written. |
| (b2) thin region 0,8 mm at 2 mm pitch | **approved.** 0 cells (whole femur 77 273), `ChainStop` in L5b with "UNKNOWN; ingen grannetikett" (A188 not repeated). |
| (c1) inertia in kg·mm² to L6 | **approved**, `UnitError`. |
| (c2) mm numbers labeled m | **passed.** Label passes unit check, R4 catches: length ratio 0,0011 outside [0,714, 1,4], `ContractError`. |
| (c3) knee vector in CT-frame | **approved**, `FrameError`. |
| (c4) gravity in z-up frame | **approved**, `FrameError`. Past the contract, the error had given 946,7 % RMS-wrong in pelvis_ty (same as M1 saw) and 19,9 N·m in hip flexion. |
| (c5) z001 bone mass without mapping | **approved**, `IdentityError`. |
| (d) p = 0 versus M1's `run_id` without composition layers | **approved.** Max \|Δ\| 3,6·10⁻¹⁵ over all 23 generalized forces; 17 of 23 bit identical. Thigh at p = 0 differs 0 kg in bulk and 2,8·10⁻¹⁷ kg·m² in inertia from gait2392. |

### Predictions

Held: P1 length (414 mm), P2 bone mass (1,17 kg), P3 hip force < 0,5 % (−0,03 %), P4 rater spread (0,026), voxel (0,0024) and density (0,030 N·m), P5 status, P6 field volume (2,7·10⁻⁵ and 1,1·10⁻³).

Fell:
- P1-the dispersion: 3,05 mm, over predicted 3.
- **P3 the moment: 2,43 N·m against predicted 0,05–1,0.** I predicted through the path of inertia (H5/M1 reasoning), which gave 0,30 N·m. The location of the knee center, i.e. a landmark that becomes a joint location, accounts for 88 % of the effect. It is correct with the K1's ranking (joint center most sensitive) but applies here to ID-moment, not just the muscle recruitment.

## 5. New questions that the connection opens

1. **Joint placement or inertia: which geometry link governs ID-error?** When landmark → joint → ID is linked, it can be seen that 10 mm femur length acts 7 times more via knee placement (2,13 N·m) than via leg inertia (0,30). Experiment: run the chain for all 19 VSD-subject with manual femur landmarks, with ∂τ/∂(knee center) and ∂τ/∂(bone inertia) exactly via the TPS linearity. Prediction: the ratio is > 5 for all subjects in support. Forged by a subject with quotient < 2.
2. **Donor geometry: how big is the error of A_map?** The contract forced that the z001 leg at subject01 time is a declared mapping. Experiment: within VSD, where each subject has own legs, replace each subject's femur with another subject's isotropically scaled femur (leave-one-out, same gait data) and measure Δτ and Δ knee placement against the subject's own leg. Prediction: the donor error in hip flexion is below 1 N·m, thus below M1's knee coupling error (1,8 N·m RMS, M1 A1k). Falsified if median exceeds 1 N·m. Then the subject's own geometry (CT, SSM) is worth its cost of ID.
3. **Which UNKNOWN-term should be measured first?** The chain has five unknown terms and one known RSS at 0,026 N·m. Experiment: give each unknown term an open, documented interval (A_soft: soft tissue mass ±20 % with center of mass ±10 mm along femur; A_frame: 0,65° plus assessors' 6,3°; A_kin: IK if with extended femur in OpenSim) and count the Δτ band per term through the same chain. Prediction: A_kin dominates (IK with longer femur changes the angles, and the angles carry the hip moment in support), followed by A_soft. Falsified if A_soft or A_frame is greater. The outcome determines whether the next measurement should be the body surface (splat with metric scale → soft tissue) or the movement (markers → IK).

A fourth question opens but is not testable yet: shape parameter → joint reaction → tissue damage with error budget all the way. The intersegmental force (654 N) lacks muscle forces. The next link is K1's local recruitment, and the S-N exponent m ≈ 10–15 (K1, not controlled) amplifies a loading error of 1 % to 10–16 % in damage.
