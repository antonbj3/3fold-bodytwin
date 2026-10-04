# CLOUD-W3-VASCULAR-SAMPLE — Results

**Preregistered verdict: FAIL (synthetic).** The Fisher-rank part passed. The coverage part failed: PS coverage was 0.972, above 0.95. **The real-world optimal schedule is UNKNOWN.**

## What was done
- `code/model.py` is an original closed two-compartment model (plasma and interstitium). Its states are volumes Vp and Vi and solute masses Mp and Mi.
  - Fluid filtration is J = kf·(ΔVp/Vp_eq − ΔVi/Vi_eq), a linearised surrogate for Starling filtration.
  - Solute flux is F = PS·(Cp−Ci) + (1−σ)·J·C_upwind.
  - Total volume and total solute mass are conserved exactly. Over 0–240 min the measured drift was 7e-15 L and 2e-13 g.
- Scenario: a 1 L solute-free bolus into plasma at t=0.
  - Known values: Vp_eq=3 L, Cp_pre=40 g/L, σ=0.9.
  - Unknown θ, estimated on a log scale: kf, PS, Vi_eq, Ci0. Nominal values are 0.03 L/min, 0.003 L/min, 12 L and 20 g/L.
  - Known noise is additive Gaussian with SD 0.4 g/L for both plasma and interstitial samples.
- Candidate sampling times ran from 5 to 240 min in 5-min steps.
- Designs were chosen by greedy D-optimality at the nominal θ.
- Fisher rank used an SVD threshold of σ_min/σ_max > 1e-6. This threshold was fixed in the code before the full run.
- Coverage used 1000 simulations with held-out truths drawn as θ_true = θ_nom·exp(N(0, log 1.25)). These truths were never used for design.
  - Each fit was nonlinear least squares.
  - Intervals were 90% Wald intervals on log θ, using the FIM at the estimate.
- Commands: `pip install numpy scipy` (numpy 2.4.6, scipy 1.17.1, Python 3.11.15), then `cd code && python3 run.py 1000`. The run takes about 4 min with seed 20260924.

## Results
| Design | Plasma times (min) | Interstitial times (min) | FIM rank | σmin/σmax |
|---|---|---|---|---|
| P6 plasma-only | 30,85,115,165,200,240 | – | 3 | 1.4e-7 |
| P8 plasma-only (equal budget) | +35, +240 | – | 3 | 1.3e-7 |
| All 48 plasma candidates | 5…240 | – | 3 | 9.7e-8 |
| P6 + I1 | P6 | 5 | 4 | 2.6e-5 |
| **P6 + I2** | P6 | 5, 120 | **4** | 1.3e-4 |
| Minimal full-rank design found | 85,165,240 | 5 | 4 | 1.2e-5 |

- With plasma samples only, the weakest parameter direction was about 0.82·logPS + 0.58·logCi0. The data cannot separate the exchange rate from the unknown interstitial starting concentration.

Held-out coverage for P6+I2 (1000 of 1000 fits converged; Wilson 95% intervals in brackets):

| Parameter | Coverage | Wilson 95% | Median 90% CI width (log) | In 0.85–0.95? |
|---|---|---|---|---|
| kf | 0.915 | [0.896, 0.931] | 0.50 | yes |
| PS | **0.972** | [0.960, 0.981] | 1.94 | **no** |
| Vi_eq | 0.862 | [0.839, 0.882] | 3.40 | yes (the interval crosses 0.85) |
| Ci0 | 0.922 | [0.904, 0.937] | 0.07 | yes |

## Preregistered criteria
1. **Added interstitial points raise Fisher rank: PASS.** Rank went from 3 to 4. Even all 48 plasma samples stayed at rank 3.
2. **90% coverage within 85–95% in ≥1000 simulations: FAIL.** This rule was fixed in the code as applying to every parameter. PS over-covers at 0.972, and Vi_eq is borderline at 0.862.
3. **Real-world optimal schedule: UNKNOWN.** No empirical data were used.

## Counterexamples and caveats
- **Rank depends on the threshold.** The plasma-only σ ratio is about 1e-7, so it is near-singular, not exactly singular. A threshold below 1e-7 would give rank 4, and the rank criterion would then fail. The practical non-identifiability is clear anyway: P8 plasma-only gave median log-CI widths of 13 to 13,800, from 300 simulations.
- **Plasma-only intervals are uninformative, not accurate.** Their coverage of 0.97–1.00 is meaningless because the intervals are so wide.
- **Nearly no fluid shift breaks the design.** At kf=1e-6, P6+I2 drops back to rank 3 (σ ratio 1.2e-13). Identifiability here depends on the fluid transient.
- **Early-only interstitial samples** at 5 and 10 min still reach rank 4, but with lower information: log-det 14.8 versus 16.8 for P6+I2.
- **Five times the noise** gives over-coverage of 0.96–0.99 (300 simulations). The Wald intervals are miscalibrated there because the model is nonlinear.
- **Why PS over-covers:** the Wald interval ignores the curvature of the model and the prior spread of the truth around the design point.
- **Pilot run disclosure:** a 20-simulation pilot was run before the full run to check the pipeline. No settings were changed after it.
- **The interstitial measurement is idealised.** A noise-free-bias Ci sample, from something like microdialysis or lymph, is an assumption. So are the linear filtration law, the fixed σ, the closed system with no lymph return or renal loss, and all parameter magnitudes.

## Sources
No external data, papers or datasets were retrieved or used. No URLs or DOIs were checked, so none are cited. Parameter magnitudes are the author's assumptions and have not been checked against literature. This is purely a synthetic demonstration, and model output does not establish real-world validity.
