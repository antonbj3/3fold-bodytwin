BT-DAT-Q052

## What was built on

`inputs/Q052_QUESTION.md` (coupling porosity → mechanics → transport for **the same sample**),
`inputs/Q052_model.py` (28 frozen parameters, unchanged — sha256 before == after, `MODEL_INPUT.sha256`)
and `inputs/NIGHT_PREAMBLE.md` §3 (precision, source/derivation/hypothesis, UNKNOWN). The register is
`DATA_SOURCES.json`; `PREREG.md` + `PREREG.sha256` were frozen before the first download. No
internal BodyTwin material was used.

## Public measured datasets that constrain the model

Seven datasets, of which four downloaded (16.47 MB total, budget 50 MB), all CC-BY-4.0, metadata
verified through the Zenodo API in this session:

| # | Dataset | What it constrains |
|---|---|---|
| DS-01 | 11061947 µCT trabecular bone (NIfTI, 9 volumes) | φ, d, d_ref, b_tau, S_v |
| DS-02 | 7516228 3D pore geometry → measured D_eff (M-factor, 4608 rows) | D_0 (the form φ/τ²) |
| DS-03 | 4641712 τ–φ compilation, 2204 points | b_tau (parsing PARTIAL) |
| DS-04 | 14617243 BAM Ti-6Al-4V, E/G vs T (6 samples) | E_s |
| DS-05 | 20856007 AM lattice, compression (NOT_FETCHED) | χ_E, m_E, m_d |
| DS-06 | 167808 35 volunteers' femur/tibia (NOT_FETCHED, 153 MB) | L (scale context) |
| DS-07 | 10635546 apparent E cranial cancellous bone | E_eff (bone side) |

## The sample loads: shape, units, coordinate frame

`samples/BMLPL_001_REF_17_SEG_SUB.nii`: NIfTI-1, `sizeof_hdr=348`, magic `n+1`, dim 100³, int16,
`vox_offset=352`, labels {0,1}. `xyzt_units=0` → **the unit is UNKNOWN in the header**; `pixdim=0.01741`
is interpreted as mm (17.41 µm; the register text states 17.59 µm) — documented assumption. `qform_code=1`
(quaternion 0,0,1), `sform_code=0` (srow zeros) → scanner coordinates, no anatomical axes.
Extent 1.74 mm.

## Key numbers (all in `results.json`)

- φ = 0.518–0.902, mean **0.692**; 3 of 8 volumes lie in the frozen [0.70, 0.90] → the model's
  φ interval is not covered by measurements at its lower end.
- Pore chord length 0.124–0.453 mm, mean **0.223 mm** = 0.28 × `d_ref` → `d` in the model is 2.7–7.2 × larger
  than measured pore geometry. The identity S_v=6φ/d gives 0.62 mm: **`d` is not a chord length** (ratio 2.78).
- S_v = 2.86–9.27 mm⁻¹, mean 6.69 mm⁻¹ vs the model's 6.0 mm⁻¹ (1.12×).
- τ (6-connected BFS, by definition an upper bound) 1.31–1.69 vs the model's 1.154 at measured φ → b_τ
  too small; DS-02 gives τ=1.50 from the M-factor and b_τ≈0.16.
- Frozen K_g, b_τ with **measured** geometry: K = 8.58e-10 m² = 0.046 × nominal 1.87e-8 and 5.7 ×
  the literature bone 1.5e-10 → K is geometry-sensitive, not validated.
- E_s: BAM 115–117 GPa at 24 °C vs frozen 110 GPa (+4–6 %).

## What failed

Coverage: **7 of 28** parameters have an absolute measured constraint (requirement ≥ 6 met), 1 only
relatively, 1 scale context, **19 lack public support**; 6 of them are the response law (hypothesis). Placebo
(64 permutations, seed 52) gives median 1, max 3 → coverage is not due to chance. All eight frozen
criteria pass.

## Gaps and next steps

No open measured dataset was found for K (requires pressure/flow), for λ, g_p, μ, σ_app, or for
the coupling process setting s → (φ, d) on the same sample. DS-05/DS-06 have complete metadata but are not acquired.
Next steps: acquire DS-05 (8.6 MB) for χ_E/m_E/m_d, and seek measured flow data; until then, Q052's
`d`, `b_tau` and response law should not be interpreted as validated.
