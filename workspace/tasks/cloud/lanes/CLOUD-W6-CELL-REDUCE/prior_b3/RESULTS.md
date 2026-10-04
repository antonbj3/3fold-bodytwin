# CLOUD-W4-B3-CELL: How cell shape affects whether diffusion-reaction parameters can be identified (partial run)

**Status:** the numerical criterion passes for the cases that finished (λ = 25, 50 and 100 nm). The λ = 200 and 400 nm cases are **UNKNOWN**. Biological calibration is **UNKNOWN**.
The sweep was stopped at the user's request after about 30 minutes. The results below keep only the completed cases.
This is a **synthetic method check** run on a public segmentation geometry. It is **not** an empirical validation: no transport measurement of any kind was used.

## Data (public, verified in this session)
- Bucket (anonymous HTTPS): `https://janelia-cosem-datasets.s3.amazonaws.com`
- Array: `jrc_hela-2/jrc_hela-2.zarr/recon-1/labels/groundtruth/crop1/mito/{s1,s2}`
  - Format: zarr v2, zstd compression, uint8.
  - Encoding: 0 = absent, 255 = unknown. Only 0 and 1 occur in the data.
  - CellMap annotation version: `0.1.1`, read from `crop1/.zattrs`.
- Voxel scale (z, y, x) in nm, read from the OME-NGFF multiscales metadata:
  - s1: 5.24 × 4 × 4, shape 100 × 500 × 500, zarray ETag `22d08d651d119275292b5695eb7696c6`, array sha256 `6436d1ad…b28f`
  - s2: 10.48 × 8 × 8, shape 50 × 250 × 250, zarray ETag `90c64ea077797609eb6b168f6c0106ca`, array sha256 `2dd499ca…76fd`
- Organelle: one connected mitochondrion.
  - Volume: 0.1306 µm³ (194,693 voxels at s2).
  - Radius of the equal-volume sphere: R = 314.7 nm.
  - **The organelle touches the crop boundary on three faces: z-min, z-max and x-min.** The crop is only 524 nm deep in z. This matches the brief's warning that the previous crop touched its boundary.
- Background references, which I did not open in this session because the proxy blocked doi.org, so treat them as unverified:
  - Heinrich et al. 2021, Nature, doi:10.1038/s41586-021-03977-3 (COSEM segmentation)
  - Xu et al. 2021, Nature, doi:10.1038/s41586-021-03992-4 (OpenOrganelle)
  - None of the numbers in this report depend on these papers.

## Model (synthetic)
- Equation: −D∇²c + k·c = 0 inside the mitochondrion, with c = c₀ on the membrane (a Dirichlet boundary, standing for uptake from the cytosol).
- Scaling and units: D = 1 and c₀ = 1. Lengths are in nm. k = 1/λ², where λ is the penetration length.
- Discretisation: cell-centred finite volumes on the anisotropic voxel grid. Membrane faces use a half-cell conductance.
- Observables:
  - total uptake J, in units of D·c₀·nm
  - effectiveness η = J/(k·V·c₀), which is dimensionless
- Crop-edge closure. Voxels of the organelle that lie on the edge of the crop are treated in two ways:
  - **exposed**: c = c₀, as if the crop edge were membrane
  - **reflecting**: zero flux, i.e. the organelle continues as a mirror image
  - The true value lies outside the crop and is not observed.
- Refinement: each s2 voxel is subdivided f = 1, 2, 3 times per axis. This gives h = 8, 4 and 2.67 nm in-plane, and 10.48, 5.24 and 3.49 nm in z.
- Solver: pyamg smoothed-aggregation preconditioned CG, tolerance 1e-12.
- Baseline comparison: an equal-volume sphere, voxelised on the same grids, compared with the analytic result η = 3/φ²·(φ·coth φ − 1), where φ = R/λ.

## Results (finest grid, f = 3)
| λ (nm) | η exposed | η reflecting | η sphere (analytic) | crop closure: exposed/reflecting − 1 | crop-edge share of J (exposed, h = 8 nm) | refinement change f2→f3 (exp / refl) | k̂/k if an equal-volume sphere is assumed (exp / refl) |
|---|---|---|---|---|---|---|---|
| 25  | 0.2768 | 0.2304 | 0.2194 | +20.1 % | 20.4 % | 0.87 % / 1.03 % | 0.60 / 0.90 |
| 50  | 0.4798 | 0.4158 | 0.4009 | +15.4 % | 20.0 % | 0.38 % / 0.49 % | 0.63 / 0.91 |
| 100 | 0.7287 | 0.6730 | 0.6539 | +8.3 %  | 19.3 % | 0.15 % / 0.21 % | 0.66 / 0.90 |
| 200 | UNKNOWN | UNKNOWN | 0.834 (analytic only) | UNKNOWN | — | UNKNOWN | UNKNOWN |
| 400 | UNKNOWN | UNKNOWN | 0.949 (analytic only) | UNKNOWN | — | UNKNOWN | UNKNOWN |

Flux balance, |J_boundary − ∫k·c| / J, was at most **9.2e-13** across all completed runs.

## Independent checks
1. **Sphere against the analytic solution** (a check of the finite-volume discretisation):
   - The relative error is +1.0 %, +0.8 % and +0.6 % at λ = 25 nm for f = 1, 2, 3, and +0.7 %, +0.4 % and +0.3 % at λ = 50 nm.
   - At λ = 100 nm it is +0.35 % and +0.19 % for f = 1 and 2. The f = 3 case is UNKNOWN because the sweep was stopped.
   - The error shrinks monotonically as the grid is refined.
2. **Second solver:** Jacobi-preconditioned CG from scipy, with no AMG, agrees with AMG to about 1e-15 relative on the h = 8 nm grid for λ = 25, 50 and 100 nm.
   - I tried a SuperLU direct solve, but it was too slow and I abandoned it. It is not a result.
3. **Segmentation resolution:** at λ = 50 nm (exposed closure), the native s1 mask gives η 1.15 % lower than the s2 mask subdivided to the same grid spacing.
   - This measures geometry uncertainty from segmentation resolution, which is separate from discretisation error.

## Preregistered criterion (frozen, not changed)
- **Flux balance < 1e-5: PASS.** The maximum was 9.2e-13.
- **Refinement change < 2%: PASS for completed cases.** The maximum was 1.03 % (λ = 25 nm, reflecting). The cases at λ = 200 and 400 nm are UNKNOWN.
- **Crop-boundary contribution reported: yes.**
  - The crop edge supplies about 19–20 % of uptake under the exposed closure.
  - Switching closure changes η by 8–20 %. This is far larger than the discretisation error.
- **Biological calibration: UNKNOWN.** There are no same-cell transport measurements.

## Interpretation, counterexamples and failure cases
- **The crop boundary is the dominant uncertainty.** Its effect (8–20 % on η, and a factor of 0.60–0.66 against 0.90 on the inferred k) is 10–100 times the numerical error.
  - A shape-resolved model built from this crop therefore does not identify k better than the closure assumption allows.
- **Counterexample to "shape matters most".** With the reflecting closure, the segmented mitochondrion differs from the equal-volume sphere by only 2.9–5.0 % in η. Fitting the sphere model would bias k by about 10 %.
  - Under the exposed closure, the same fit biases k by 34–40 %.
  - So on this crop, most of the apparent effect of shape is actually an effect of the crop artifact.
- **Structural non-identifiability** (standard limits, not computed here for λ ≥ 200 nm):
  - When λ ≫ organelle size, η → 1 for any shape, so J ≈ kVc₀ and shape cannot be identified.
  - When λ ≪ organelle size, J ≈ c₀·S·√(Dk), so only S·√(Dk) can be identified: surface area and k are confounded.
  - J alone never separates D from k without a second observable.
- **Failure cases and limitations:**
  - Voxel staircasing overestimates the surface area. This is benign here only because λ ≥ 25 nm ≥ 3h.
  - Only one crop and one organelle were used.
  - The reaction is first order, with a Dirichlet membrane: no permeability and no cristae compartments.
  - The mito_mem and mito_lum labels were not used.
  - The λ = 200 and 400 nm cases and the f = 3 sphere case at λ = 100 nm were not run.
  - The crop-edge share was measured only on the h = 8 nm grid.

## Reproduce
```
pip install numpy scipy pyamg zstandard
python3 fetch_cellmap.py s1 s2   # ~10 s, writes data/*.npy + meta (sha256 above)
python3 fv_shape.py              # full sweep, prints run_log lines (~30+ min; stopped early here, log kept as run_log.txt)
python3 supplement.py            # ~1 min: crop-edge share, independent CG, s1-vs-s2 geometry
python3 analyze.py               # builds results.json from run_log.txt + supplement.json
```
Uncertainty:
- The discretisation error is below 1.1 %, judged from the f2→f3 change and the sphere-versus-analytic error.
- The geometry error from segmentation resolution is about 1.2 %.
- The error from the crop closure is 8–20 %, and it dominates.
