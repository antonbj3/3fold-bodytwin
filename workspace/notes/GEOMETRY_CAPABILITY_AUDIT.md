# BodyTwin: geometric capability today (O2, 2026-09-22)

The question concerns what exists for the chain **real anatomical mesh → physiologically meaningful shape parameter → traceable landmarks/regions/interfaces/attachment points → measurable consequence in a computation**. The answer is based on code run today, not on documentation.

Sources, all read without modification:

| Source | Path | Version |
|---|---|---|
| Published | `references/current_bodytwin` → `external_mount` | `b95a8dc` "Add geometry contracts and refuse unapproved chained results" |
| Staging | `../3fold-bodytwin` | `2ec6cb2` |
| Private | `source_repository/` | `2c42ec52f9` |
| Field engine | `local_path` | `3ffdb0e` |
| Dental session's cell inventory | `../dental/notes/COMPUTE_CELL_INVENTORY_bodytwin.md` | Not duplicated here. It covers 232 cells, `tmj_lever_model` and the geometry modules, §2b. |

The published working tree already had three modified `reports/verification_*.json` with mtime 2026-09-17. I have not modified them.

---

## 1. What was run

Everything ran on CPU with `CUDA_VISIBLE_DEVICES=''`, `OMP_NUM_THREADS=4` and Python in `3fold_staging/.venv-bodytwin` (numpy 2.5.3, scipy 1.18.1, trimesh 5.1.0, warp 1.13.0). Runs were made in a copy of the published tree, `$S/pub` (`S=/tmp/coordinator-1000/-home-anton/b240ebb1-…/scratchpad`), because `compose_mesh_to_field.py` writes to `examples/anatomy/` in the source tree (`out=ROOT/'examples/anatomy'`). Logs are in `results/O2/`.

| Run | Exact command | Outcome |
|---|---|---|
| Atlas → seam → SDF field | `results/O2/compose_cmd.txt`: `THREEFOLD_ROOT=local_path WARP_CACHE_PATH=$S/warpcache … python examples/anatomy/compose_mesh_to_field.py` | exit 0 in 5,7 s. All 6 gates PASS. Voxel volume is 103 640 mm³ against the mesh's 103 741,155 mm³, a relative error of 9,75e-4. The field is 34×26×58 with 2 mm pitch. The field's SHA is `b1f26106…d7d9`, the same as in RUNNING.md. Result JSON and npz are **byte-identical** to the published files (`compose_before_sha.txt`). |
| Geometry tests (20 files: frame, geometry_validity, i2_*, material_roster, region_*, thermal, thickness, uncertainty, cross_module) | `results/O2/geometry_tests_cmd.txt` (requires `PYTHONPATH=src`) | **306 passed, 2 skipped** in 14,8 s. The two skipped lack external reference files (`BODYTWIN_REF_PROPERTY_REGISTRY_DIR` and `BODYTWIN_REF_AREA_DESCRIPTOR`). |
| Same tests without `PYTHONPATH=src` | `geometry_tests_no_pythonpath.log` | 5 collection errors (`No module named 'bodytwin'`). RUNNING.md specifies `pytest -q tests` without PYTHONPATH and without package installation. Neither `pyproject` nor `setup` exists. **The documented test line therefore does not work as written** for the geometry tests. |
| O2 probe on the kidney atlas and OpenMandible | `results/O2/o2_probe_cmd.txt` → `o2_geometry_probe.json` | exit 0 in 1,3 s. See §4. |
| O2 landmark probe | `python results/O2/o2_landmark_feasibility.py` → `o2_landmark_feasibility.json` | exit 0. See §4. |

Warp wrote "CUDA error 100: no CUDA-capable device" to stderr. This is expected with the GPU hidden and is not a GPU failure.

---

## 2. Capability table, published `b95a8dc` (`src/bodytwin/geometry/`)

Conventions throughout the package: mm is explicitly required (`units='mm'`, otherwise ValueError). There is no unit guessing, no repair and no resampling. Missing data gives ABSTAIN/UNKNOWN and is never replaced by a default value.

Status abbreviations from RUNNING.md: VF = VERIFIED-FRESH, SO = SYNTHETIC-ONLY, OGF = OWN-GATE-FAIL.

| Capability | File / function | In → out (units, frame) | What the tests check | Evidence | Limit |
|---|---|---|---|---|---|
| Strict surface seam | `mesh_ingest_v1.surface_mesh(v,f,units='mm',weld_exact=)` | (N,3) float in mm and (M,3) int → `MeshArrays` (read-only, int32). `motion_arrays()` gives m (×0,001). | watertight, consistent orientation, positive volume, no duplicates and exact byte welding; 8 invalid inputs rejected | VF; kidney today: 4 vertices welded, all vertex coordinates unchanged | no certification of self-intersection; no anatomical validity |
| Tet seam | `mesh_ingest_v1.tetrahedral_mesh` | tet (K,4) → oriented boundary surface | positive orientation, non-manifold surface, parity for internal surfaces; 1/3-volume test | VF (contract) | tet generation missing, the mesh must be supplied |
| Validity certificate | `geometry_validity_v1.geometry_certificate(..., evaluate=True)` | surface → dict with accepted, reason and three-valued `self_intersections` | open, nonconforming and self-intersecting surface rejected; `None` = abstains | tests passed; **kidney today: accepted=True, self_intersections=False, 0,17 s** | O(n) BVH plus pairwise test, not tested on mesh with over 40 000 triangles |
| Mass, center of gravity and inertia | `region_mass_v1.mass_moments_surface`, `mass_moments_tetra` | closed surface or tet (mm) plus `MaterialProperty` → kg, mm and kg·mm² about the center of gravity in the global frame | analytical box, surface = tet, rotational covariance, L³ scaling, missing density → UNKNOWN | test passed; kidney today: see §4 | homogeneous density per region; the density registry only has MSK-MUSCLE, NEU-TISSUE, CARD-VESSEL, MSK-BONE and GENERIC (assumption) |
| Inertia in frame | `frame_v1.Frame`, `frame_region_bridge_v1.moments_in_frame` | origin plus orthonormal right-handed basis (`axes`: **columns = local basis vectors**) → centroid and inertia in local frame | identity, translation, rotation = axesᵀ I axes, missing density | test passed; today: principal-axis frame gave off-diagonal 7,8e-15 | no fiber directions or anisotropy fields; the frame is rigid, without deformation |
| Volume, area and interfaces per region | `region_mass_v1`: `area_mm2`, `interface_area_by_pair_mm2`, `region_boundary_volume_mm3`, `voxel_region_volume_mm3`, `tetra_volume_mm3`, `functional_exchange_area_mm2` (must be declared) | front/back-labeled triangles → mm² and mm³ | interfaces counted once, partition preserved, pitch contract | tests passed | requires conforming interfaces; exchange area ≠ wall area |
| Cross section (muscle) | `anatomical_cross_section_mm2(V,L)`, `physiological_cross_section_mm2(V,L_f,θ)` | mm³, mm and rad/deg → mm² | ACSA and PCSA kept distinct, angle unit must be specified | test passed | fiber length and pennation must be specified; no connection to mesh yet |
| Interfaces from tet | `tetra_interfaces_v1.tetra_interfaces(v,cells,labels)` | conforming tet with label ≥ 1 → oriented triangles with front/back | 4 valid and 2 invalid checks | SO 6/6 | conformity assumed; nonadjacent intersections not detected |
| Interfaces from voxel | `label_interfaces_v1.label_interfaces(uint8)` | label grid → quad triangles | – | **OGF** 4/5 (a non-manifold edge at contact between multiple labels) | origin 0 and pitch 1, no frame |
| Thickness | `thickness_v1.ThicknessDistribution`, `relaxation_time_s` | samples in mm → distribution; τ = L²/(H_A·k) | statistics, units, exponent 2 | test passed | distribution supplied by caller; no measurement on mesh |
| Ray traversal in tet | `tetra_ray_walk_v1`, `_batch`, `_cuda` | tet partition plus ray → cell sequence and parameter intervals | 512 rays exact, CPU = CUDA | SO | first exit only; singular starts rejected |
| Synthetic fixtures | `layered_box_v1`, `layered_cylinder_v1` | – | – | SO | not anatomy |
| Payload, registry and uncertainty | `i2_payload_builder/validate/region_registration`, `material_roster_v1`, `uncertainty_v1`, `thermal_capacity_v1` | region record with units, provenance and KNOWN/ABSENT/MISSING | contract tests (included in the 306) | passed | bookkeeping; no geometric operation |
| Mesh → field | `examples/anatomy/compose_mesh_to_field.py` plus the field engine's `faltkarna_v1_mesh_to_sdf.mesh_to_sdf_del` | kidney in mm → SDF and solid with 2 mm pitch | 6 gates | VF; **byte-identical today** | rasterized distance transform, not exact distance; no contact integration |
| Optics | `geometry/optics/*` | photon transport in tet/voxel | – | VF/SO/OGF per variant | see dental inventory |

Motion and muscle cells in published code:
- `cells/musculoskeletal/tmj_lever_model.py` is the only cell with a moment arm. The moment arms are **scalar sweep parameters** (`r_muscle_mm` 15–30, `d_molar_mm` 40–60, `d_incisor_mm` 95–120) and are not read from any geometry.
- `sagittal_equilibrium(F, r, d)` solves bite force and condylar reaction.
- According to the cell's own derivation, the molar/incisor bite force ratio depends only on `d_incisor/d_molar`.
- No published code does any of the following: landmarks, shape parameters, deformation, attachment points on a surface, muscle paths, scaling of any type other than affine, or tet generation.

## 3. Staging and private

**Staging** (`2ec6cb2`) is an older subset of the published code. It has `mesh_ingest`, `label/tetra_interfaces`, `tetra_ray_*`, `layered_*` and optics, and shared files are identical. The following are missing in staging: `frame_v1`, `frame_region_bridge_v1`, `geometry_validity_v1`, `region_mass_v1`, `material_roster_v1`, `thickness_v1`, `thermal_capacity_v1`, `uncertainty_v1`, `i2_*`, `examples/geometry` and their tests. Staging only has `examples/optics` in addition to published code. Geometrically, there is therefore nothing in staging that is not already published.

**Private** (`source_repository/`). Most of the actual musculoskeletal geometry is here, via **OpenSim 4.6** (`.venv-msk`). Everything is private and nothing is published.

| Capability | Where | Has it been run? | Comment |
|---|---|---|---|
| Moment arms via `Muscle.computeMomentArm` | `scripts/msk/moment_arm_validation.py`, `moment_arm_lit_compare.py` | yes: `data/msk_smoketest/moment_arm_validation/` (2 810 raw rows, 2026-07-21) plus literature comparison **11 PASS / 6 FLAG / 9 PLAUSIBLE** of 26 rows (FLAG = adductor magnus) | The model is LaiArnoldModified2017, scaled to an OpenCap LabValidation subject. The sign is measured, not assumed. Graph nodes `HOLE-MOMENT-ARM-SENSITIVITY-UNWIRED-AT-KNEE-CELL` and `AUDIT-MOMENT-ARM-TAUTOLOGY-…` are OPEN. |
| Segmentwise scaling (OpenSim ScaleTool: marker pairs → scale factors, the muscle's path points follow the body) | `scripts/msk/subject_specific_scaling.py`, `knee_jw_transplant/make_scale.py`, `data/msk_models/*_scaled.osim` (35 models) | yes | Scaling is linear per segment. Fmax is not scaled (verified in the C++ source). With `preserve_mass_distribution=true`, mass becomes a global factor. |
| Attachment points and muscle paths | `add_arm_muscles.py`, `add_trunk_flexors.py`, `add_erector_spinae.py`, `anatomical_hand.py`, `attach_band.py` (PathPoint/PathSpring) | yes (evidence JSON exists) | The attachment points are points in the segment frame, not barycentric points on a bone geometry. |
| Contact and cartilage thickness on mesh | `cartilage_contact.py`, `cartilage_thickness_material.py` (JAM `Smith2018ContactMesh`, lenhart2015 Apache-2.0), `hip_cartilage_contact.py` | yes (`data/msk_smoketest/cartilage_*`, `hip_cartilage_contact`) | Own ray-based thickness computation exists. |
| Shape from video (bone length and ratios, not SMPL) | `scripts/david_shape_batchfit.py`, `…_identifiability_audit.py` | yes (`data/david_shape_*`) | `smplx` is not installed; the script says so explicitly. |
| Parametric body mesh (LBS shape plus residual) | `scripts/physics_exp/d_human_mesh_observability.py` | synthetic | Identifiability certificate, no real anatomy. |
| Remesh with projection to nearest point | `scripts/physics_exp/spine_isotropic_remesh.py` | yes, on robot mesh | numpy and scipy |
| RBF/thin-plate | 19 hits, all in `scripts/cad/*` and `prosth_hand_v1.py` (CAD and prosthesis) | – | no anatomical morphing |
| Statistical shape model/PCA-SSM, TotalSegmentator processing, BodyParts3D, Z-Anatomy, tetgen, the reference model | 0 code hits (the reference model only mentioned in a docstring) | – | `MSKEL-GEOM-OVERDET` is **designed** (ASSUMED, HIGH risk) but not built. `AUTO-TMJ-LEVER-CERT-V2-PIN-GEOMETRY…` (OPEN) is exactly the chain in the assignment: "pin geometry from a digitized CT/CBCT mandible". |

The dental inventory, row 30, states that no OpenSim model contains a mandible. Thus there is no private TMJ geometry either.

**Field engine, kernel and motion** (`3fold_public/*`): 0 hits for RBF, thin-plate, landmarks or morph.
- The field engine has `closest_point` (trimesh) in the watertightness gate and its own barycentric hit during ray winding. These are internal tools, not an API.
- `kropps_id_falt_v1` concerns owner indexes in contact fields (robot cells), not anatomy.
- The motion engine has inertia for analytical primitives and URDF inertia, not for mesh.
- Mass inertia for mesh already exists in BodyTwin (`region_mass_v1`).

## 4. Today's probe on real meshes, with published functions without modification

**Kidney** (BodyParts3D FJ3145/FMA7205, 1 493 vertices after welding):
- The certificate is accepted and self-intersection is evaluated as False.
- Kidney density is missing in the registry, so `mass_moments_surface(density=None)` gives `known=False, missing=['density']`. This is correct ABSTAIN.
- With GENERIC density (1,05 kg/L, *assumption*), mass is 0,1089 kg and principal inertias 22,8 / 81,2 / 89,9 kg·mm².
- An affine elongation by 10 % along the longitudinal axis gives a volume and mass factor of exactly 1,100 and inertia factors of 1,100 / 1,311 / 1,291. Vertex identity is trivially preserved. This shows that a shape parameter can give a measurable consequence with existing code, but only for affine mappings.

**OpenMandible** (`external_mount`, git `e1f8cef`):
- `misc/mandible_model_used_for_prescribing_muscles_insertion.stl` has 128 976 raw vertices and 42 992 triangles. After exact welding, **it passes the strict seam** (`surface_mesh` PASS) with volume 23 387 mm³ and extent 107×101×54 mm, which is reasonable in mm.
- All **14 muscle attachments on the mandible** (`bc_directed/*_mandible_insertion.stl`: superficial and deep masseter, temporalis anterior/middle/posterior, medial and lateral pterygoid, left and right) are **exact vertex subsets** of the mandible mesh (fraction 1,0 for all 14). The attachments can therefore be carried as vertex-id sets: no projection, and barycentric coordinates become trivial.
- The origin surfaces on the skull exist as 14 separate STL files, and `misc/skull_model…stl` also exists.
- Landmarks can be read from labeled parts without manual pointing:
  - `Cartilage.stl` is the condylar cartilage. Both sides are in one file, so the centroid ends up in the midline, (138,4; 185,4; 75,7) mm. A split per side requires a split into components.
  - The enamel for incisor L1/R1 and molar L6/R6 gives bite points.
- **Preliminary, with unchecked axis orientation and the crown's centroid instead of the occlusal contact point:**
  - The condyle–incisor distance in y–z projection is 86,1–86,6 mm; the cell's sweep is 95–120.
  - The condyle–L6 distance is 63,9–65,4 mm; the cell's sweep is 40–60.
  - The ratio d_inc/d_mol ≈ 1,33, against 1,6–3,0 in the cell's sweep.
  - This points to `tmj_lever_model` needing to be anchored in geometry. It is precisely the unbound node `AUTO-TMJ-LEVER-CERT-V2`. The numbers should not be interpreted until axes and contact points are defined.

## 5. Missing building blocks for the demonstration chain

| Link | Exists | Missing |
|---|---|---|
| Real mesh input | `surface_mesh`, the certificate and field seam (VF). Kidney and OpenMandible pass the seam. | reading multiple parts/labels as a region set (STL per part → label per vertex/triangle) |
| Shape parameter | affine only (today in the probe); OpenSim segment scaling private | **landmark-driven non-affine deformation**: thin-plate-RBF (`scipy.interpolate.RBFInterpolator(kernel='thin_plate_spline')` exists in venv) that moves vertices and retains topology and indexes; a **physiologically named parameter** (e.g. ramus height, mandibular body length or gonial angle) defined as landmark displacement; a gate for inverted triangles and self-intersection after deformation (reuse `geometry_certificate(evaluate=True)`) |
| Landmarks, regions and attachment points | labeled interfaces (front/back); attachments as vertex subsets in OpenMandible | a **landmark type**: vertex-id or (triangle-id, barycentric coordinates) plus derivation rule (centroid of label set, extremum along axis) plus unit and frame; a split into components for bilateral labels; an **anatomical frame** from landmarks (e.g. intercondylar axis, occlusal plane) built on `frame_v1.Frame` |
| Consequence in computation | mass and inertia (`region_mass_v1`), PCSA formula, `tmj_lever_model.sagittal_equilibrium` | **moment arm from geometry**: the muscle's line of action from the attachment's centroid on the mandible to the origin's centroid on the skull, moment arm = ‖(p − c) × û‖ about the condylar axis, plus d_bite from condyle to contact point; a thin adapter that feeds `sagittal_equilibrium` and condylar reaction per shape parameter value |

**Minimum new code**, about 200–300 lines of numpy and scipy, all additive in `geometry/`:
1. `landmarks_v1`: type and rules, where a label set gives centroid/extremum and barycentric evaluation on deformed mesh.
2. `tps_morph_v1`: source landmarks and target landmarks → vertex displacement. Gates: identity at zero displacement, exact interpolation at landmarks, rejection at inversion, same vertex and label indexes before and after.
3. `moment_arm_v1`: line of action plus joint axis → moment arm in mm, with analytical tests (right angle, parallel line → 0).
4. An example that chains 1–3 with `region_mass_v1` and `tmj_lever_model.sagittal_equilibrium`.

The field engine and motion engine have none of this. The barycentric hit in the field engine is internal and should not be bound in.

## 6. Recommended minimum demo chain

**Dataset: OpenMandible base model**, with mandible, 14 attachments, 14 origins on the skull, condylar cartilage and teeth per material. It already passes the strict seam, and the attachments are exact vertex subsets.
- License: the code is GPL-3.0, but **the license for model data is not specified** (UNKNOWN, according to `../dental/notes/orient_C_code_chain.md`). This is enough for an internal demonstration. Before publication, clarification from the author is required, or a switch to BodyParts3D mandible (CC BY-SA 2.1 JP / CC BY 4.0, but attachments are missing there) or Open-Full-Jaw (CC BY-NC-SA 4.0, 17 patients, conforming tet and tooth/PDL/bone, but without muscle attachments).

**Chain:** mandible STL → `surface_mesh`/certificate → landmarks from labels (condyle per side, incisor and molar contact, attachment centroids) → TPS change of a named parameter (e.g. ramus height ±10 %) → certificate after deformation → moment arms for masseter, temporalis and medial pterygoid about the intercondylar axis, d_incisor, d_molar and mass/inertia for the mandible → `sagittal_equilibrium` → **measurable consequence**: bite force per N muscle force, condylar reaction and the molar/incisor ratio as a function of the shape parameter. Compare with the cell's sweep intervals and the Waltimo anchor (909/382 N for men).

Decorrelated anchor: OpenMandible's own Ansys/FE boundary conditions (`02 Use cases`) and published moment arms for jaw muscles. The private `moment_arm_lit_compare.py` is a ready pattern for the literature comparison (±30 % and ±20° with external justification).

Limb alternatives (not recommended first):
- TotalSegmentator cases on disk (`external_media{,2}`, 3+8 CT with femur, hip and gluteus segmentation; CC BY 4.0 according to the private source registry) require marching cubes and surface cleaning before the seam.
- The JAM/Grand Challenge knee (`data/external/opensim_jam_grand_challenge_DM/Geometry`) only has bone surfaces. The attachments are in `.osim` as points.

## 7. Other observations

- `compose_mesh_to_field.py` writes to the source tree. The documented `pytest -q tests` requires `PYTHONPATH=src`.
- The frame convention in `frame_v1` is easy to misuse: `axes` has the basis vectors as **columns**. My first probe gave rows and got off-diagonal 17,7 kg·mm² instead of ~1e-14. This was a caller error, not the library's fault, but a demo chain should test it explicitly.
- The disk inventory is in `results/O2/mesh_inventory_raw.txt`. Other anatomical meshes on disk:
  - `mandibular_defect_147` (338 STL, `external_mount`)
  - Teeth3DS (CC BY-SA 4.0)
  - STS-3D-Tooth
  - MMDental (license discrepancy noted in `KATALOG_TILLAGG_2026-09-22.md`)
  - OpenSim's bundled `Geometry/*.vtp` (gait models)
  - multisegment_foot STL
  - JAM cartilage
  - lenhart2015 (Apache-2.0)
- No Z-Anatomy or complete BodyParts3D archive was found. The archive's SHA is in `examples/anatomy/THIRD_PARTY.md`.
