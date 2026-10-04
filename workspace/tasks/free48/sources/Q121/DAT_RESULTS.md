BT-DAT-Q121

**What was built on.** Upstream `BT-HX-Q121` (`inputs/Q121_model.py`) met its frozen
acceptance test (222.465 mL vs 227.0 ± 60.7 mL/24 h) with **36 of 46** parameters still
marked `assumption` in `PARAMETER_META` — every residence time, gas-transfer coefficient,
compliance and residual gas volume. This run maps public measured data onto those 36 keys.
PREREG.md frozen first, sha256 `2adf08d8…` (pre-erratum `c6a340e4…`, both in `results.json`;
the erratum corrected one estimated count, no criterion).

**8 public datasets found, 7 anonymously retrievable** (`DATA_SOURCES.json`, full fields):

| # | dataset | lic. | size | n | constrains |
|---|---|---|---|---|---|
| DS1 | MRE bowel segments, Zenodo 13839321 | CC BY 4.0 | 441.5 MB | 114 | `v0_*`, `initial_liquid_*` |
| DS2 | TotalSegmentator v2, Zenodo 10047292 | CC BY 4.0 | 23.6 GB | 1228 | `v0_*`, `compliance_*` |
| DS3 | MSD Task10_Colon | CC BY-SA 4.0 | 6.24 GB (HEAD) | 190 | indirect (RAS frame/scale) |
| DS4 | Colonic high-res manometry, SPARC 33v2 | CC BY 4.0 | 400.83 MB | 95 | `k_forward/backward_pressure_*`, `rectal_pressure_threshold_mmhg`, `motility_gain` |
| DS5 | H2/CH4 breath test CIPO, Zenodo 16910802 | CC BY 4.0 | 19.0 kB | UNKNOWN | `k_fermentation_h`, `k_michaelis_g`, `yield_h2_mmol_g`, `k_microbial_h` |
| DS6 | Lactulose H2–CH4, Zenodo 20578413 | CC BY 4.0 | 26.1 kB | 90 | same, on the model's own substrate |
| DS7 | WMC transit table, figshare 5453194 | CC BY 4.0 | 5.6 kB | UNKNOWN | `tau_proximal_h`, `tau_distal_h`, `tau_rectal_h` |
| DS8 | UW-Madison serial abdominal MRI, Kaggle ds/3577354 | CC BY-NC 4.0 | UNKNOWN | 107 | between-day spread of `v0_*` (needs account) |

**Sample downloaded and proven to load.** Case 23 of DS1, pulled from the 441.5 MB archive by
HTTP Range (central directory at byte 441,443,907), CRC32-matched per member, **2,923,353 B
total in `samples/`** vs a 52,428,800 B limit. Form: NIfTI, shape (300, 320, 27), uint16,
1.40625 × 1.40625 × 4.8 mm voxels, **LPS, mm (xyzt code 2), sform_code 1**, affine in
`samples/23_data.nii.gz`. Intensity units are **UNKNOWN** (MR signal, not HU). Measured label
volumes: colon 600.5 mL vs small intestine 647.4 mL; mapped onto the model's three
compartments 297.3 / 231.8 / 71.4 mL against nominal `v0` 40 / 50 / 25 mL (ratios 7.43 /
4.64 / 2.86). **Scale bracket only** — these are mannitol-distended luminal *content* volumes
in one IBD patient, not gas volumes, and the mapping is a stated segmentation choice.

**Second sample, the most direct constraint.** DS7 loads as a 7×6 xls: healthy-control CTT
median 21.5 h (3.1–64.6), WGTT 28.5 h (9.3–73.6). The model's Στ = 42 h (10+14+18) falls
inside the healthy WGTT range, 1.47× its median. Statistic types differ (median transit vs mean
residence of a serial chain), and no per-segment taus are in the deposit, so the three taus
cannot be separated by this file.

**Checks.** C4: R recomputes to 62.363598 mL·mmHg/(mmol·K) vs 62.3637 coded; V_molar at
273.15 K/760 mmHg = 22.41397 vs 22.414 — both constants are dimensionally sound. `h_h2_mmol_ml_mmhg`
remains **UNKNOWN**: temperature-dependent, unsourced.

**What fell / gaps.** 4 of 6 families covered. **No admissible dataset found for
phase-transfer/dissolution** (`k_phase_transfer_h`, `h_h2_mmol_ml_mmhg`,
`k_dissolved_diffusion_h`) — the very terms the model's separation claim rests on — **nor for
water/secretion** (`base_water_rate_ml_h`, `water_loss_h`). MSK label tar returned 404, so
DS3's label path is UNKNOWN; DS8's byte size is UNKNOWN. Two intestinal-gas cohorts (Azpiroz
n=37; McWilliams n=207) are listed separately as publication tables with no deposit, and are
**not** counted as coverage.

**Next step.** Pull DS4's pressure traces and DS5/DS6 breath time courses to bound the
pressure-conductance and fermentation constants, and treat the phase-transfer gap as the
priority — it is the one that would let retention be re-attributed to dissolution.
