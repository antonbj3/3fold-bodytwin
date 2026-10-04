BT-HX-Q012

## Conclusion

A first executable first-principles model of the chain **dietary precursor → colonic microbial TMA →
hepatic FMO3 → TMAO → renal excretion** is in `model.py` with 11/11 PASS in
`dimensional_check()`, 29/29 PASS in `test_model.py`, and all numbers in `results.json`.
Frozen at `PREREG.sha256` = `1607b3c9885ecea9d966abd4302239c5488613a505427eb0f111d6912e245647`.

Prediction at the frozen protocol (70 kg, 430 mg L-carnitine = 250 mg d₃ + 227 g steak at 180 mg,
cf. Koeth 2013 Fig. 1): **C_ss = 0.657008 µM**, **Φ = 0.988652** (degree of microbial
dependence), capacity ceiling **108.2326 µM** with margin **164.74×**, exposure contrast
**3.3184×** across the microbiome capacity's 25th–75th percentile.

## Review measures (frozen, PREREG §4)

| Measure | Value | Threshold | Outcome |
|---|---:|---:|---|
| G1 degree of microbial dependence Φ | 0.988652 | ≥ 0.95 | **PASS** |
| G2 exposure contrast | 3.3184× | ≥ 2× | **PASS** |
| G3 identifiability of `p = f_esc·y` | 34.11 | < 0.10 | **FAIL** |
| G4 capacity margin | 164.74× | ≥ 3× | **PASS** |

G1 is conditional on `r_host_umol_h = 0.062 µmol/h`, a free assumption that has not been
looked up. Sensitivity is weak: at `r_host = 0.20 µmol/h`, Φ becomes 0.9643, still
above 0.95; G1 fails only at `r_host = 0.30 µmol/h` (Φ = 0.9474).

## What failed

- **G3 FAIL.** Relative 90 % width for the product is 34.11 against frozen threshold 0.10; individual
  factors have 6.198 (`f_colonic_escape`) and 18.312 (`y_microbial_conversion`). With the
  declared priors, the product cannot be determined to 10 %. The product width is
  also *not* smaller than the individual factors', so H2's strong form fails: from a
  single steady-state plasma concentration without labelling, `f_esc` and `y` are separately
  unidentifiable — only the product is identifiable.
- **A prereg defect found and reported, not silently fixed.** PREREG §6 placebo required
  Φ to remain unchanged when `y` varies. That is false: Φ = 1 − R_host/(R_host+R_micro)
  depends on `y` by definition (Φ(0.001) = 0.3034 → Φ(1.0) = 0.9977). The true
  invariance — unchanged Φ under *joint* scaling of `r_host` and `y` — is tested instead
  and holds (spread 0.000e+00 over 20× rescaling). The PREREG file is unchanged.
- **Numerical antibiotic suppression ratio: UNVERIFIED.** Koeth 2013 Fig. 1b–e
  verbally states "near complete suppression" and "virtually no detectable formation"; no
  numerical ratio exists in the acquired text (PMC3086762 and PMC3701945 gave HTTP
  500/404 from Europe PMC REST). The threshold 0.95 in G1 is a declared interpretation of the
  verbal expression.
- **Absolute plasma TMAO levels: UNVERIFIED.** No comparison against published
  concentrations is made. `C_ss = 0.657 µM` is a model property, not a measured level.
- **K03/K11/K14 have no labels in the workspace** (grep in `../../notes`, `../../data`,
  `../../docs`, `../../scripts`, `../../START*.md`, `../../tasks` → zero hits; no
  `results/K03|K11|K14` among 3105 directories). The mapping K03=diet, K11=substance/metabolite,
  K14=microbiome function is positional, from `inputs/QUESTION.md`'s own column.

## Uncertainty that remains

Monte Carlo, 20000 draws, seed 20260925, 7 lognormal priors: median 0.654965 µM,
5th percentile 0.109349 µM, 95th percentile 4.199799 µM — **90 % interval width 38.407×**.
OAT log-variance shares (`variance_contributions`, interactions neglected): microbiome
0.6664, diet 0.2836, host 0.0501. `vd_l` has share 0.0 because Vd does not enter the
steady-state expression — it appears only in the half-life 4.0842 h.

## Which variables are needed (Q012)

Greek selection (`needed_variables`) until the 90 % width ≤ 2×:

| Step | Variable | Block | 90 % width afterwards |
|---:|---|---|---:|
| — | start (everything at median) | — | 38.407× |
| 1 | `y_microbial_conversion` | microbiome | 9.052× |
| 2 | `f_colonic_escape` | diet | 2.864× |
| 3 | `f_tma_to_hepatic` | microbiome | 2.265× |
| 4 | `f_fmo3` | host | 1.645× |

Four variables, of which two microbiome, one diet, one host. `gfr_ml_min_per_kg`, `vd_l` and
`r_host_umol_h` are not needed for the 2× target. Four more add nothing.

## Sensitivity ±50 % (frozen selection)

| Parameter | −50 % | +50 % |
|---|---|---|
| `y_microbial_conversion` | 0.3322 µM (−0.494) | 0.9818 µM (+0.494) |
| `f_colonic_escape` | 0.3322 µM (−0.494) | 0.9818 µM (+0.494) |
| `gfr_ml_min_per_kg` | 1.3140 µM (+1.000) | 0.4380 µM (−0.333) |

`y` and `f_esc` are exactly equivalent in this respect (+0.494 each) — the product governs.
`gfr` is inverted; its absolute effect is largest despite the smallest log-variance share, because
C_ss ∝ 1/CL.

## Null models and placebo

N0-A (diet-proportional, k calibrated at the frozen point) gives 0.6570 µM at y = 0 where
the mechanism gives 0.007456 µM — 88.1× error, so the null model cannot express
microbial knockout. N0-B (microbiome is everything) gives 2.4132 µM, 3.67× above. Placebo
(Km → 0, same code path) reproduces the closed-form expression with relative difference 0.000e+00.

## What the next resolution step requires

1. **Separate `f_esc` from `y`.** Requires labelled passage through the colon, i.e. d₃-L-carnitine →
   d₃-TMA in plasma and 24-h urine, as in Koeth 2013 Fig. 1. Without this, the product
   `f_esc·y` is all that is identifiable (G3 FAIL).
2. **Bind `r_host` with a numerical value.** The antibiotic arm's remaining TMAO as a
   fraction of the baseline value, measured on the same individuals before/after. This determines whether G1 is
   0.9887 or 0.845.
3. **Acquire an absolute plasma TMAO value with table/figure**, to anchor `CL` and
   `V_d` to data rather than the assumption GFR = 1.8 mL/min/kg.
4. **Measure FMO3 `Km`/`Vmax`.** `f_fmo3` is rank 4; the capacity limit is at
   68.335 g/d precursor per day, which is far above the protocol but not if
   the dietary peak spikes (e.g. carnitine supplements).
5. **Acquire the actual K03/K11/K14 inputs** so that the mapping can be checked.

## Verification

`python3 test_model.py` → 29/29 PASS, 1.4 s, 1 thread, < 1 GB. `python3 model.py` →
`results.json`, 1.2 s. Analytical limiting cases: stoichiometry 1:1 (TMA/TMAO = 1.000000000000,
no 3:1 — FMO3 adds one O per TMA), exact dose linearity (doubling gives
2.000000000000×), capacity plateau (103.233 µM at Km = 5 µM, 108.233 µM at Km → 0,
analytical ceiling 108.233 µM), microbial zero (y = 0 leaves 0.007456 µM host source),
mass balance `C·CL = R` closed to 1e-12, ODE convergence towards both closed-form expressions
(relative difference 2.374e-13).
