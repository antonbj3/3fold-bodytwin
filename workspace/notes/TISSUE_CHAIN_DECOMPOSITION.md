# The tissue chain in the smallest computable components

Agent K2, 2026-09-22. The graph is in `results/K2/components.json` (72 nodes, 71 edges, 11 links), built by `results/K2/build_components.py`. The preregistration is `results/K2/PREREG.md` (sha256 in `PREREG.sha256`). All numbers below are in `results/K2/out/*.json`.

## 1. The chain

Geometry → region partition → material → CT density → density–modulus → load → solver → strain/SED → mechanobiology → observable quantity → measurement noise.

Each link is broken down into idea, mechanism, equation, operation, representation and assumption, until the leaves are individual operations. Each leaf has the fields `have` (file:function on our side or null), `status`, `neg` (negative results N1–N14 in `FIELD_RESULTS_FOR_BODYTWIN.md` §3), `sens` (sensitivity to the final question "local strain → bone adaptation → BMD/BV/TV") and `dental` (shared with the dental session).

| link | leaves (selection) | exists on our side | status | sensitivity |
|---|---|---|---|---|
| L1 geometry | load surface, frame, validity, segmentation threshold | `geometry/mesh_ingest_v1`, `frame_v1`, `geometry_validity_v1` | tested; 5 of the first 10 Keast tibias are not watertight in `*-cortical.stl` | M |
| L2 region | cortical shell = two nested surfaces; thickness along normal; surface → voxel; surface → tet; lost region → UNKNOWN | `tetra_ray_walk_v1`, `thickness_v1`, `label_interfaces_v1`; K2 `k2_common.label_volume` | **tested (T2)**; UNKNOWN output missing (A188) | H |
| L3 material | roster KNOWN/ABSENT/MISSING; E_kortex; anisotropy; ν; viscoelasticity; poroelasticity | `material_roster_v1`, `bone_composite_modulus_bounds`, `cartilage_poroelastic_relaxation`, `tendon_relaxation_v1` | roster only covers density and heat capacity; anisotropy ABSENT | H |
| L4 CT density | HU → ρ_QCT (phantom); ρ_ash → ρ_app; partial volume | nothing | **ABSENT: no CT with HU in our datasets** | H |
| L5 density → modulus | choice of relation; density type; extrapolation to cortex; X3's n = 2 | K2 `t1_density_modulus.REL`; `results/X3/model.py` | **tested (T1)** | H |
| L7 solver | voxel-hex8; tet-FE; homogenization of thin layer; h convergence | K2 `k2_voxfe`; private `warpfem_*`; dental F1 (ccx) | voxel-FE patch tested; **homogenization tested (T3)** | M/H |
| L8 response | strain, principal strain, SED, averaging volume | K2 `k2_voxfe.elem_strain` | the averaging volume for stimulus is an open assumption | H |
| L9 mechanobiology | Frost thresholds; RANKL/OPG + Scheiner; BMU | `bone_wolff_law_mechanostat`, `bone_remodeling`, `bone_rankl_opg_lemaire2004`, `results/X3` | X3 calibrated against Lang 2004 and Cummings 2009 | H |
| L10 observation | BV/TV; aBMD/vBMD; synthetic image | private `bmd_femur_partial_recovery.py` | no forward model from ρ field to DXA | H |
| L11 noise | DXA precision; QCT artifacts; uncertainty propagation | `uncertainty_v1` | DXA precision does not exist in code | H |

## 2. Choice of leaves to test

Three leaves have high sensitivity and were untested:

1. **L5.choice**: which density–modulus relation is used. E differs 10–23× between published relations at the same density.
2. **L2.vox**: whether cortical thickness is preserved in the surface → voxel conversion. Here the field engine already has two negative results, A188 (thin region disappeared silently) and A194 (the sum check was blind).
3. **L7.homog**: which rule is used to mix a thin layer in a coarse voxel.

Data: Keast et al. 2023, tibia SSM (CC BY 4.0). The dataset contains CT-segmented surfaces for both outer cortex and endosteal boundary. The cortical shell is therefore measured and not synthetic. There is no CT density in any of our datasets, so the density in T1 and T3 is synthetic and declared as such.

## 3. T2: cortical thickness through the surface → voxel conversion (real geometry)

Method: 5 tibias and up to 3000 rays per tibia along the inward normal. The reference is the hit on the endosteal surface. The same line is then sampled in the voxel label. The control with an analytical cylindrical shell reproduced the thickness with max error 2.3e-4 mm.

| h | median \|error\| (mm) | lost cortex, t_ref < 1 mm (n = 503) | lost, 1–2 mm | \|ΔV\|/V cortex |
|---|---|---|---|---|
| 2 mm | 0.58–0.64 | **41 %** (grid phase 0) / **48 %** (phase h/2) | 13 % | ≤ 0.74 % |
| 1 mm | 0.29–0.31 | 9–11 % | ≤ 0.3 % | ≤ 0.13 % |
| 0.5 mm | 0.14–0.15 | 0–0.8 % | 0 | ≤ 0.03 % |

Outcome against PREREG:
- **P2a holds.** At 0.5 mm, median \|error\| is ≤ 0.153 mm, below the limit 0.25 mm.
- **P2b confirmed.** At 2 mm, 41 % of the thin cortex is lost, against the threshold 30 %.
- **P2c holds** (volume error ≤ 0.74 %). Together with P2b, this means the volume check does not see the local loss. This is the same blindness A194 showed, now on real anatomy.

The lost cortex receives the neighboring label TRAB or 0 without warning, which is the same error as A188. The grid phase moves the loss fraction by 6 percentage points.

## 4. T1: choice of density–modulus relation → subchondral strain (real geometry, synthetic density)

Method: proximal 60 mm of tibia 102480 in voxel-FE with h = 1 mm (83 792 elements). Cortex has E = 17 GPa. Trabecular bone has homogeneous ρ_app. The load is 2000 N on the plateau's joint surface and the cut is fixed-clamped. Five relations are compared: Morgan 2003 tibia, pooled and femoral neck, Keller 1994 and Carter–Hayes 1977.

| ρ_app | E spread S0 | spread in subchondral strain S | damping D = ln S / ln S0 | spread in cortical strain | spread in stiffness |
|---|---|---|---|---|---|
| 0.20 | 22.9 | 5.00 | 0.51 | 1.56 | 2.51 |
| 0.30 | 14.9 | **5.92** | **0.66** | **1.78** | 2.56 |
| 0.45 | 9.6 | 5.41 | 0.75 | 1.98 | 2.25 |

Outcome against PREREG (ρ = 0.3):
- **P1a holds**, S = 5.92 ≥ 2.
- **P1b holds**, D = 0.66 < 0.8. Cortex carries part of the load, but only part.
- **P1c failed.** The cortical strain spread is 1.78, not below 1.5. The choice of relation therefore also reaches cortex.
- Countertest (i): a homogeneous body gives S = S0 within 1.3e-12. Countertest (ii): h = 1.5 mm changes the spread by −2.1 %.

Beyond PREREG (reported but does not gate):
- Within the Morgan family, the spread is only 1.33–1.45×. Keller and Carter–Hayes give 2.7–4.3× more strain than Morgan pooled.
- At ρ = 0.3, the Morgan relations give 443–616 µε and Keller/Carter–Hayes 1806–2621 µε. The repo's mechanostat threshold is at 1500 µε (`bone_remodeling.THRESHOLD_VARIANTS`). **The choice of relation therefore determines the regime**: maintenance or formation. This leaf precedes X3's mechanobiology and every BMD prediction.

## 5. T3: homogenization of a thin cortical layer at h = 2 mm (real geometry, synthetic density)

Method: same crop and load as T1, with relation R2 at ρ = 0.3. The reference is binary voxel-FE with h = 0.5 mm (663 793 elements), solved matrix-free with Jacobi-PCG: 6249 iterations, residual 1.0e-8, 1.08 GB RSS. The coarse model has h = 2 mm (10 387 elements). The volume fractions per coarse voxel are computed from 64 subsamples.

| rule | ΔQ3 stiffness | ΔQ1 subchondral strain |
|---|---|---|
| majority label | −4.1 % | −11.9 % |
| Voigt | **+0.01 %** | −30.9 % |
| Reuss | −20.1 % | **+2.5 %** |
| HS upper | −4.5 % | −24.2 % |
| countertest: homogeneous fine against coarse | +3.1 % | +0.2 % |

Outcome against PREREG:
- **P3a holds.** Voigt gives stiffness within 10 %.
- **P3b holds** with a small margin: Reuss −20.1 % against the threshold −20 %.
- **P3c failed** in the preregistered form. Not all rules miss Q1 by more than 25 %: Reuss is at +2.5 % and majority at −11.9 %.

The important part lies beside the criteria. **No rule hits both stiffness and local strain.** The rule that is best for stiffness (Voigt) is second worst for subchondral strain, and vice versa for Reuss. The choice of rule moves Q1 by 33 percentage points, from −30.9 % to +2.5 %. The discretization error for a homogeneous body is only 0.2–3.1 %, so the spread comes from the mixing rule and not from the mesh.

The coarse model also loses 11.5 % of the cortical volume in the element mask (φ_ben ≥ 0.5), while total bone volume only changes +0.6 %. The sum check is therefore blind here too (N5). Caveat: the Q1 set is 456 coarse voxels against 24 912 fine ones, a comparable volume but not identical points.

## 6. Shared with dental (coordination, no duplication)

| leaf | BodyTwin | dental |
|---|---|---|
| thickness along normal (`tetra_ray_walk_v1` + `thickness_v1`) | cortex, cartilage | PDL width (GEOM-ANCHOR-PDL-BONE) |
| thin layer disappears at coarse h, lost region → UNKNOWN | cortex < 1 mm (T2) | PDL 0.1–0.3 mm is more serious. The cylinder control shows that the loss stops only at h ≤ t: t = 0.5 mm is lost to 100 % at h = 2, 38.5 % at h = 1 and 1.5 % at h = 0.5. PDL therefore requires h ≲ 0.1–0.2 mm, or a tet route |
| two-phase bounds Voigt/Reuss/HS (`bone_composite_modulus_bounds`) | bone, T3 | enamel, dentin, jawbone (MAT-TOOTH-TISSUES, MAT-BONE-HU-MODULUS) |
| density–modulus, HU → ρ | T1; CT missing | jawbone from CBCT, where gray values are not HU calibrated |
| PDL constitutive model, ν ≈ 0.5, relaxation | — | dental F1 owns the leaf |
| cartilage/disc biphasic τ = L²/(H_A k) | knee cartilage | TMJ disc |
| RANKL/OPG under compression | X3 | orthodontic movement (hypothesis) |
| region vocabulary | 7 id | needs DENT regions |

Dental uses Open-Full-Jaw tet and CalculiX for PDL (F1). K2 has not touched those leaves.

## 7. Deviations from PREREG

- T2: five of the first ten cases (103559, 103862, 107813, 108421, 113033) have a `*-tibia-cortical.stl` that is not watertight (4–25 bodies). Parity is then undefined, so the cases were replaced with the next watertight ones. Cases used: 102480, 102924, 107215, 112802, 116070.
- T1: the run was interrupted by systemd-oomd at 22:06 and resumed. ρ = 0.2 was already done and was reused from JSON. Same code and same input.
- T3: the fine reference (0.5 mm) is solved matrix-free with Jacobi-PCG instead of AMG-CG. The discretization is the same and the residual requirement is the same (1e-8). The switch was made after the oomd interruption to keep memory below 1 GB.
- No tet route was tested (declared in PREREG).
