# RESULTS — CLOUD-W3-OED-CORREL: correlated measurement error in a minimal protocol


## Sources / versions
- External data: **none**. No URLs or DOIs were used because every number was generated inside `oed_correl.py`. The error magnitudes (landmark total SD 7 mm, imaging total SD 9 mm, mean offsets +8 / −3 mm, prior SDs 18/24/22 mm) are **assumed, illustrative values**. They were not taken from any checked publication.
- Code: `oed_correl.py` uses only the Python 3.11 standard library, with seed 20260924. Preregistration `PREREG.md` sha256 `113206f9…f571` matches `PREREG.sha256`.
- Exact command: `python3 oed_correl.py` (about 25 s). It writes `results.json`.

## Design
- Parameters θ = [pelvis width, thigh length, shank length] (mm), with a Gaussian prior (θ2–θ3 correlation 0.6).
- Candidate measurements (8), each measuring one component:
  - landmark palpation: L1 pelvis, L2 thigh, L3 shank, L2b and L3b as repeat palpations;
  - imaging: I1 pelvis, I2 thigh, I3 shank.
- Error model: `e_j = mu_g + b_{g,subject} + eps_j`. The per-subject **systematic bias** `b_g` is shared by every measurement in modality g, which gives a within-modality correlation ρ.
  - S5 changes the imaging bias to a multiplicative 2 % scale error, so it is not exactly the Gaussian model that is fitted.
- Minimal protocol: budget M = 3 of the 8 candidates. The design is chosen by exhaustive c-optimal search, minimising the posterior variance of the target c'θ. The target is leg length θ2+θ3, except in S4, which uses the contrast θ2−θ3.
- Posterior (linear-Gaussian): `P = (Σ0⁻¹ + HᵀR⁻¹H)⁻¹`, `m = PΣ0⁻¹μ0 + PHᵀR⁻¹(y − μ̂)`. The 90 % interval is `c'm ± 1.645·sqrt(c'Pc)`.
- Calibration vs held-out:
  - R and μ̂ are estimated from a **calibration cohort** of 100 synthetic subjects with known truth.
  - Covariance-aware R: shared + independent variance per modality, estimated by moments from the residual covariance.
  - Naive R: the diagonal of the same sample covariance. The marginal variances are identical; only the correlation is dropped.
  - Coverage is measured on a **separate held-out cohort** of 5000 subjects drawn from an independent RNG stream.
- Uncertainty is reported three ways:
  - Wilson 95 % CI on the held-out coverage;
  - exact analytic coverage given the calibrated R and the true R;
  - 200 replicate calibrations, each tested on 1000 held-out subjects.

## Baseline and main results (target θ2+θ3, 90 % nominal)

| Scenario | Design (aware) | Design (naive) | Held-out coverage aware [95 % CI] | Held-out coverage naive [95 % CI] | Actual SD aware / naive (mm) | Naive claimed SD (mm) |
|---|---|---|---|---|---|---|
| S1 ρ=0.8 (prereg case) | L2, L3, I3 | L2, L3, L3b | **0.889** [0.880, 0.898] | **0.717** [0.704, 0.729] | 10.46 / 12.69 | 8.41 |
| S2 ρ=0.5 | L2b, L3b, I2 | L2, L2b, L3b | 0.873 [0.863, 0.882] | 0.739 [0.727, 0.751] | 9.94 / 11.25 | 7.83 |
| S3 ρ=0 (baseline, no correlation) | L3, L2b, L3b | same | 0.885 [0.876, 0.894] | 0.877 [0.867, 0.885] | 8.35 / 8.35 | 7.95 |
| S5 ρ=0.8, multiplicative imaging bias | L2, L3b, I3 | L2, L3, L3b | 0.873 [0.863, 0.882] | 0.680 [0.667, 0.693] | 10.50 / 12.73 | 7.92 |

Replicate calibrations (200 replicates × 1000 held-out subjects each), for S1:
- aware coverage median 0.886 (2.5–97.5 %: 0.846–0.919); 95.5 % of replicates fall within 85–95 %;
- naive coverage median 0.706 (0.641–0.776); 0 % of replicates fall within 85–95 %;
- the two methods selected the same design in 0 % of replicates.

**How naive independence changes the selected measurement.** With ρ = 0.8 the naive model assumes that palpating again (L3b) averages the landmark error away. It therefore picks three landmark measurements, and the shared soft-tissue bias counts twice in θ2+θ3.
- The covariance-aware model replaces the repeat with an imaging measurement (I3). The two modalities' biases are independent, so this adds genuinely new information.
- Result: the actual error SD falls from 12.69 mm (naive design) to 10.46 mm (aware design). The oracle design under the true R has SD 10.43 mm.
- The naive model claims SD 8.41 mm, so it is overconfident by a factor of about 1.5.
- Separating the two effects: applying covariance-aware inference to the *naive* design restores coverage (0.907) but not precision (SD 12.49 mm). So correlation affects **both** the choice of measurements and whether the intervals are calibrated.

## Preregistered criterion
"PASS if covariance-aware 90 % intervals cover 85–95 % and naive intervals fail in high-correlation case across ≥1000 trials."
- High-correlation case S1, 5000 held-out trials: aware coverage 0.889, which is inside 85–95 %; naive coverage 0.717, which is outside.
- **PASS (synthetic).** It also holds in 95.5 % of replicate calibrations. In the other 4.5 %, aware coverage fell below 0.85 because of the small calibration cohort.

## Counterexamples and limitations
1. **S4, contrast target θ2−θ3 with ρ = 0.8.**
   - The shared bias cancels, so the naive intervals are *over*-conservative: coverage 0.998 (claimed SD 7.97 mm vs actual 4.48 mm). "Naive fails" does not always mean undercoverage.
   - The **covariance-aware method also fails here**: coverage 0.827 [0.816, 0.837], analytic 0.833, and only 39 % of replicate calibrations fall within 85–95 %.
   - Cause: the independent-error variance is estimated as a small difference of noisy moments from 100 subjects. The design search then exploits the underestimate (a winner's-curse effect). Using plug-in covariance estimates is therefore not reliably calibrated for bias-cancelling targets. Fixes such as a larger calibration cohort, propagating uncertainty in the estimated R, or cross-validated design selection were not tested.
2. **S3 (ρ = 0):** naive and aware agree and are equally calibrated (0.877 vs 0.885). Naive independence is harmless when the errors are in fact independent.
3. **Coverage sits slightly below nominal everywhere** (about 0.87–0.89). This comes from plug-in estimation with N_cal = 100 plus selecting the design on the same estimates. When R is known exactly, the intervals are calibrated by construction. With S1's calibrated R, the exact analytic coverage is 0.901, so the empirical 0.889 there mainly reflects Monte Carlo error.
4. **Assumptions that were not tested on data:**
   - linear measurement operators;
   - Gaussian prior and errors;
   - bias shared within a modality but independent across modalities;
   - no bias correlation within a subject across sessions;
   - one bias per modality, not per landmark.

   Real palpation and imaging errors (for example landmark-specific soft-tissue artefact, or scale errors correlated with body size) may be structured differently. S5 probes one such deviation only.
5. **Everything is synthetic.** The magnitudes are illustrative, not measured. Empirical validity on real musculoskeletal protocols is **UNKNOWN**.

## Files
- `oed_correl.py`: simulation, calibration, design search, held-out evaluation.
- `results.json`: all numbers above, including per-scenario analytic checks and replicate summaries.
